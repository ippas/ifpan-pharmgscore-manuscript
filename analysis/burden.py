import os
from itertools import count

import pandas as pd
import hail as hl
from hail.plot import show
from bokeh.layouts import gridplot

from analysis.utils.load_spark import hl_init
from analysis.utils.load_spark import wd, localfs_path, scratch_path
from analysis.utils.annotations import find_annot_desc, load_annotations
from analysis.sample_control import annotate_genes

if int(os.getenv('SLURM_NNODES')) > 1:
    tmp_path = os.path.join(scratch_path, 'tmp/')
else:
    tmp_path = os.path.join(localfs_path, 'tmp/')

def str_count():
    counter = count()
    while True:
        yield f"{os.getenv('SLURM_JOBID')}-{str(next(counter))}"
dir_counter = str_count()


hl_init()
hl.plot.output_notebook()

# # get random control: analysis/sample_control.py
# # it creates vcf matrix table ('data/pvcf/vcf-adr-ctrl-better.mt')
# # saves patients' eid and group at adr-control-patients.csv
# 
vcf = hl.read_matrix_table(wd + 'data/pvcf/vcf-adr-ctrl-better.mt')
cadd = hl.read_table(os.path.join(
    '/net/archive/groups/plggneuromol/imdik-zekanowski-gts',
    'data/external-data/cadd-full.ht'
))
# 
# # %% pca
# vcf = vcf.repartition(5000)
# vcf = vcf.checkpoint(scratch_path + 'tmp/f9859c')

# eigenvalues, pcs, _ = hl.hwe_normalized_pca(vcf.GT, k=20)
# vcf = vcf.annotate_cols(eig_scores=pcs[vcf.s].scores)

# df = pd.read_csv(wd + 'data/adr-control-patients.csv', dtype={'f_eid': str})
# adr_eid = set(df[df.adr].f_eid)
# vcf = vcf.annotate_cols(
#     adr_group=hl.if_else(hl.literal(adr_eid).contains(vcf.s), 'adr', 'control')
# )
# vcf.write(wd + 'data/pvcf/vcf-adr-ctrl-better-eig.mt')
# vcf = hl.read_matrix_table(wd + 'data/pvcf/vcf-adr-ctrl-better-eig.mt')

p = hl.plot.scatter(
    vcf.eig_scores[0], vcf.eig_scores[1], label=vcf.adr_group,
    title='PCA', xlabel='PC1', ylabel='PC2', size=2,
)
show(p)
p = hl.plot.scatter(
    vcf.eig_scores[2], vcf.eig_scores[3], label=vcf.adr_group,
    title='PCA', xlabel='PC3', ylabel='PC4', size=2,
)
show(p)

vcf = vcf.filter_cols(vcf.eig_scores[0] > -0.02)
p = hl.plot.scatter(
    vcf.eig_scores[0], vcf.eig_scores[1], label=vcf.adr_group,
    title='PCA', xlabel='PC1', ylabel='PC2', size=2,
)
show(p)

p = hl.plot.scatter(
    vcf.eig_scores[2], vcf.eig_scores[3], label=vcf.adr_group,
    title='PCA', xlabel='PC3', ylabel='PC4', size=2,
)
show(p)

vcf = hl.split_multi_hts(vcf)
vcf = vcf.explode_rows(vcf.within_gene)
vcf = vcf.annotate_rows(
    maf=vcf.info.AF[vcf.a_index - 1],
    cadd_score=cadd[vcf.row_key].cadd_score,
)

genes_list = hl.import_table(wd + 'raw/genes-list.csv', delimiter=',')
dg = hl.dict(list(
    zip(genes_list['Gene'].collect(), genes_list['Group'].collect())
))
vcf = vcf.annotate_rows(
    gene_group=dg.get(vcf.within_gene),
    one_group='one_group'
)

control_list = pd.read_table(wd + 'raw/control-genes-lists.csv', sep='\t')
# duplicated - remove
set(genes_list.Gene.collect()) & set(control_list.skeletal_muscle.dropna())
set(genes_list.Gene.collect()) & set(control_list.long_genes.dropna())
set(control_list.skeletal_muscle.dropna()) & set(control_list.long_genes.dropna())
control_list['skeletal_muscle'] \
    [control_list['skeletal_muscle'] == 'ADCY2'] = None
control_list['long_genes'] \
    [control_list['long_genes'].isin(['NRXN3', 'AGBL1', 'ANKS1B'])] = None
all_group_genes = (
    genes_list.Gene.collect()
    + list(control_list.short_genes.dropna())
    + list(control_list.long_genes.dropna())
    + list(control_list.skeletal_muscle.dropna())
    + list(control_list.intestine.dropna())
)
assert len(all_group_genes) == len(set(all_group_genes))
gene_list_dict = {}
for gene in genes_list.Gene.collect():
    gene_list_dict[gene] = 'list'
for col in control_list.columns:
    for gene in control_list[col].dropna():
        gene_list_dict[gene] = col
gene_list_dict = hl.literal(gene_list_dict)
vcf = vcf.annotate_rows(
    gene_list=gene_list_dict.get(vcf.within_gene)
)

vcf.write(wd + 'data/pvcf/vcf-adr-ctrl-better-eig-long.mt')
##########

vcf_g = hl.read_matrix_table(wd + 'data/pvcf/vcf-adr-ctrl-better-eig-long.mt')

# burden
group_col = ['within_gene', 'gene_group', 'gene_list', 'one_group'][0]
"""
within_gene - for all genes
gene_group - CYP, CNS, other
gene_list - list, intestine, short, long, ...
one_group - genes from the list as one group
"""
genes_on_list = False  # within_gene but only for genes on the list
AF_TH, CADD_TH = 0.05, 1
vcf_fr = vcf_g.filter_rows(
    (vcf_g.maf < AF_TH)
    & (vcf_g.cadd_score > CADD_TH)
    & ~hl.is_missing(vcf_g[group_col])
)
if group_col == 'within_gene' and genes_on_list:
    vcf_fr = vcf_g.filter_rows(
        vcf_g['gene_list'] == 'list'
    )

if group_col != 'within_gene':
    vcf_b = (
        vcf_fr
        .key_rows_by(*vcf_fr.row_key.keys(), group_col)
        .distinct_by_row()
        .key_rows_by(*vcf_fr.row_key.keys())
    )
vcf_b = (
    vcf_fr
    .group_rows_by(vcf_fr[group_col])
    .aggregate(
        burden_score=hl.agg.sum(vcf_fr.GT.n_alt_alleles() * vcf_fr.cadd_score)
        # burden_score=hl.agg.sum(vcf_fr.GT.n_alt_alleles() * -hl.log10(vcf_fr.maf))
    )
)

burden_mt_path = f'{scratch_path}/tmp/cadd-burden-{next(dir_counter)}'
vcf_b = vcf_b.checkpoint(burden_mt_path)

vcf_b_entries = vcf_b.entries()
v_stats = (
    vcf_b_entries
    .group_by(vcf_b_entries[group_col], vcf_b_entries.adr_group)
    .aggregate(stats=hl.agg.stats(vcf_b_entries.burden_score))
)
v_stats = v_stats.repartition(100)
v_stats = v_stats.checkpoint(f'{tmp_path}/fc0a33-{next(dir_counter)}')
v_stats.show()

# %% logisit regression
vcf_b = hl.read_matrix_table(burden_mt_path)

mt = vcf_b.annotate_cols(y=hl.if_else(vcf_b.adr_group == "adr", 1, 0))

LTEST, NCOV = 'lrt', 3
result_ht = hl.logistic_regression_rows(
    test=LTEST,
    y=mt.y,
    x=mt.burden_score,
    covariates=[1, *[mt.eig_scores[i] for i in range(NCOV)]],
)

# annotate with gene group name
genes_list = hl.import_table(wd + 'raw/genes-list.csv', delimiter=',')
result_ht = result_ht.annotate(
    group=genes_list.key_by('Gene')[result_ht.key].Group
)
result_ht = result_ht.checkpoint(f'{tmp_path}/044894-{next(dir_counter)}')

# annotate with burden sum
v_adr = v_stats.filter(v_stats.adr_group == 'adr').key_by(group_col)
v_control = v_stats.filter(v_stats.adr_group == 'control').key_by(group_col)
result_ht = result_ht.annotate(
    adr_burden_sum=v_adr[result_ht.key].stats.sum,
    control_burden_sum=v_control[result_ht.key].stats.sum,
)
result_ht = result_ht.annotate(
    burden_comp=hl.if_else(
        result_ht.adr_burden_sum > result_ht.control_burden_sum,
        'adr > control',
        'adr ≤ control',
    )
)

result_ht = result_ht.checkpoint(f'{tmp_path}/a7a3a5-{next(dir_counter)}')

p = hl.plot.qq(
    pvals=result_ht.p_value,
    title=f"Burden CADD >{CADD_TH}; genes list; {NCOV} cov; {LTEST}",
    label=result_ht.group,
    size=4,
)
p2 = hl.plot.qq(
    pvals=result_ht.p_value,
    title=f"Burden CADD >{CADD_TH}; genes list; {NCOV} cov; {LTEST}",
    label=result_ht.burden_comp,
    legend=True,
    size=4,
)
show(gridplot([p, p2], ncols=2))

result_ht.order_by(result_ht.p_value).show()


# result_ht.to_pandas().to_csv(
#     wd + 'results/adr-burden/burden-2/burden-log-maf-list.csv',
#     index=False
# )

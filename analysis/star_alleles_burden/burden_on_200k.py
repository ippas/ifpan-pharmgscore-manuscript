import os
import re
from collections import defaultdict, namedtuple
from typing import Counter
from functools import reduce

import pandas as pd
import numpy as np
import hail as hl
from hail.plot import show
from bokeh.io import show, output_notebook
from bokeh.layouts import gridplot
from IPython.display import display
output_notebook()

from analysis.utils.load_spark import hl_init
from analysis.utils.load_spark import wd
from analysis.utils.load_spark import tmpdir_path_iter
from analysis.utils.pathogenicity_scores import PharmGScore
from analysis.utils.annotations import query_biomart
from analysis.star_alleles_burden.pharmvar import PharmVar
from analysis.utils.annotations import query_biomart, reorder_cols
from analysis.utils.annotations import find_annot_desc, load_annotations
from analysis.utils.annotations import translate_icd10_codes
from analysis.utils.variant_filtering import VCFFilter


def codes_to_set(codes_array, first_n=4):
    codes_set = hl.set(hl.map(lambda x: x[:first_n], codes_array))
    codes_set = codes_set.remove(hl.missing(hl.tstr))
    return codes_set


hl_init()
hl.plot.output_notebook()
tmpdir_path = tmpdir_path_iter()
hl.spark_context().setLogLevel('WARN')

############
## READ DATA
############
CYPS = [
    'CYP1A2',
    'CYP2B6',
    'CYP2C19',
    'CYP2C8',
    'CYP2C9',
    'CYP2D6',
    'CYP2E1',
    'CYP3A4',
    'CYP3A5',
]

list_genes = (
    hl.import_table(
        wd + 'raw/genes-list.csv',
        key='Gene',
        delimiter=','
    )
    .rename({'Gene': 'gene_name', 'Group': 'group'})
    .key_by('gene_name')
)

fafb = hl.read_table(wd + 'data/full-annots-for-burden.ht')
fafb.count()

burden = hl.read_matrix_table(wd + 'data/pharmacogenetic-score/ukb-pgs.mt')
burden_norm = hl.read_matrix_table(wd + 'data/pharmacogenetic-score/ukb-pgs-norm.mt')

annot = load_annotations(wd + 'data/ukb-annotations.ht')
a = annot.select(
    'f_41202', 'f_41204', 'f_41270', 'f_21022', 'f_189', 'f_22009', 'f_22006'
)
a = a.annotate(
    f_21022=hl.int(a.f_21022),
    f_189=hl.float(a.f_189),
    f_22009=hl.flatmap(lambda x: hl.map(hl.float, x), a.f_22009),
)
a = a.checkpoint(next(tmpdir_path))


#############
## REGRESSION
#############

fafb2 = fafb.select(
    fafb.f_31,  # sex
    fafb.f_189,  # townsend
    fafb.f_21022,  # age at the recruitment
    fafb.f_22006,
    fafb.f_22009,
    fafb.group,
    fafb.intention,
    fafb.self_harm,
    fafb.drug_abuse,
    fafb.mental_health_inpatient,
    fafb.is_match_6,
    fafb.match_6,
    fafb.is_match_17,
    fafb.match_17
)
fafb2 = fafb2.transmute(
    sex=hl.if_else(fafb2.f_31 == 'Male', 1, 0),
    townsend=hl.float(fafb2.f_189),
    age=hl.int(fafb2.f_21022),
    genetic_pc=hl.flatmap(lambda x: hl.map(hl.float, x), fafb2.f_22009),
    genetic_ethnic_grouping=fafb2.f_22006,
)
fafb2 = fafb2.checkpoint(next(tmpdir_path))

burden2 = burden.semi_join_cols(fafb2)
mt = burden2.annotate_cols(
    **fafb2[burden2.col_key]
)
mt = mt.checkpoint(next(tmpdir_path))

## Comparisons

comps = {}
# Comparison 1:
mt0 = mt.annotate_cols(
    group=hl.if_else(mt.group == 'adr', 1, 0),
)
mt0 = mt0.checkpoint(next(tmpdir_path))
comps[('adr', 'control')] = mt0

# Comparison 2:
mt1 = mt.filter_cols(mt.is_match_6 | (mt.group == 'adr'))
mt1 = mt1.annotate_cols(
    group=hl.if_else(mt1.group == 'adr', 1, 0),
)
mt1 = mt1.checkpoint(next(tmpdir_path))
comps[('adr', 'is_match_6')] = mt1

# Comparison 3:
mt2 = mt.filter_cols(mt.is_match_17 | (mt.group == 'adr'))
mt2 = mt2.annotate_cols(
    group=hl.if_else(mt2.group == 'adr', 1, 0),
)
mt2 = mt2.checkpoint(next(tmpdir_path))
comps[('adr', 'is_match_17')] = mt2

# Comparison 4:
anx_dep = ['anxiety', 'both', 'depression']
mt3 = mt.filter_cols(
    hl.literal(anx_dep).contains(mt.mental_health_inpatient)
    | (mt.group == 'adr')
)
mt3 = mt3.annotate_cols(
    group=hl.if_else(mt3.group == 'adr', 1, 0),
)
mt3 = mt3.checkpoint(next(tmpdir_path))
comps[('adr', 'mental health inpatient (anx/dep/both)')] = mt3

# Comparison 5:
mt4 = mt.filter_cols((mt.self_harm == 'with_drugs') | (mt.group == 'adr'))
mt4 = mt4.annotate_cols(
    group=hl.if_else(mt4.group == 'adr', 1, 0),
)
mt4 = mt4.checkpoint(next(tmpdir_path))
comps[('adr', 'self harm with drugs')] = mt4

# Comparison 6:
mt5 = mt.filter_cols(mt.group == 'adr')
mt5 = mt5.annotate_cols(
    group=hl.if_else(mt5.intention == 'intentional', 1, 0),
)
mt5 = mt5.checkpoint(next(tmpdir_path))
comps[('adr intentional', 'adr accid/therap/unspec')] = mt5

# Comparison 10:
adr_ctr = hl.import_table(
    wd + 'data/adr-control-patients.csv',
    delimiter=',',
    types={'f_eid': hl.tstr, 'adr': hl.tbool},
).key_by('f_eid')
mt9 = mt.annotate_cols(
    group=hl.int(adr_ctr[mt.col_key].adr)
)
mt9 = mt9.filter_cols(hl.is_defined(mt9.group))
mt9 = mt9.checkpoint(next(tmpdir_path))
comps[('first adr', 'first control')] = mt9

# Comparison 11:
mt10 = mt.annotate_cols(
    group=hl.if_else(mt.group == 'adr', 1, 0),
)
mt10 = mt10.filter_cols(
    mt10.genetic_ethnic_grouping == 'Caucasian'
)
mt10 = mt10.checkpoint(next(tmpdir_path))
comps[('adr cauc', 'control cauc')] = mt10

# Comparison 12:
genes_list = pd.read_csv('raw/genes-list.csv')
cyps = genes_list[genes_list.Group == 'CYP'].Gene
cyps_without_2d6 = cyps[cyps != 'CYP2D6']
mt2 = mt.annotate_entries(
    pgs=hl.if_else(mt.pgs < 50, 0, mt.pgs)
)

mt_genes = mt2.filter_rows(
    hl.literal(cyps.to_list()).contains(mt2.gene_name)
)
mt_genes = mt_genes.annotate_cols(
    pgss=hl.agg.sum(mt_genes.pgs),
    burden_group='cyps'
)
mt_genes = mt_genes.checkpoint(next(tmpdir_path))
m1 = mt_genes.cols()
mt_genes = mt2.filter_rows(
    hl.literal(cyps_without_2d6.to_list()).contains(mt2.gene_name)
)
mt_genes = mt_genes.annotate_cols(
    pgss=hl.agg.sum(mt_genes.pgs),
    burden_group='cyps_without_2d6'
)
mt_genes = mt_genes.checkpoint(next(tmpdir_path))
m2 = mt_genes.cols()

mt12 = m1.union(m2)
mt12 = mt12.drop('group')
mt12 = mt12.to_matrix_table(
    ['burden_group'],
    ['s', 'intention', 'age', 'sex', 'townsend', 'genetic_pc']
)
adr_accid = ['accidental', 'therapeutic', 'unspecified']
mt12 = mt12.annotate_cols(
    group=hl.if_else(hl.literal(adr_accid).contains(mt12.intention), 1, 0),
)
mt12 = mt12.checkpoint(next(tmpdir_path))
comps[('adr accidental', 'all other')] = mt12

results, plots = {}, {}
Comparison = namedtuple('Comparison', ['name', 'score', 'covs', 'ltest'], defaults=['wald'])

runs = [
    Comparison(
        name=('adr accidental', 'all other'),
        score='pgss',
        covs=(17, 'age', 'sex', 'townsend'),
        ltest='wald'
    ),
    Comparison(
        name=('adr accidental', 'all other'),
        score='pgss',
        covs=(0, ),
        ltest='wald'
    ),
    Comparison(
        name=('adr', 'control'),
        score='pgs',
        covs=(17, 'age', 'sex', 'townsend'),
        ltest='wald'
    ),
    Comparison(
        name=('adr intentional', 'adr accid/therap/unspec'),
        score='pgs',
        covs=(17, 'age', 'sex', 'townsend'),
        ltest='wald'
    ),
    Comparison(
        name=('adr', 'self harm with drugs'),
        score='pgs',
        covs=(17, 'age', 'sex', 'townsend'),
        ltest='wald'
    ),
    Comparison(
        name=('adr', 'mental health inpatient (anx/dep/both)'),
        score='pgs',
        covs=(17, 'age', 'sex', 'townsend'),
        ltest='wald'
    ),
]
for run in runs:
    if run in results:
        continue
    mtr = comps[run.name]
    n_cov, *covs_str = run.covs
    for k in mtr.__dict__:
        if not k.startswith('_') and mtr[k] is mtr[run.score]:
            score_name = k[6:] if k.startswith('score') else k
    title = (
        f'{score_name.upper()};'
        f' {" vs ".join(run.name)} (NCOV={n_cov} + {"/".join(covs_str)});'
        f' {run.ltest}'
    )
    print('-----------------------------------------------------')
    print(title)
    print('-----------------------------------------------------')
    result = hl.logistic_regression_rows(
        test=run.ltest,
        y=mtr.group,
        x=mtr[run.score],
        covariates=[
            1,
            *[mtr[k] for k in covs_str],
            *[mtr.genetic_pc[i] for i in range(n_cov)]
        ],
    )
    result = result.annotate(
        group=list_genes[result.key].group
    )
    file_name = '-'.join(run.name)
    file_name = re.sub('[/() ]+', '_', file_name)
    file_name = re.sub('_+$', '', file_name)
    file_name += f'-{run.covs[0] + len(run.covs[1:])}covs'
    result = result.annotate_globals(title=title)
    result = result.repartition(1000)
    result.write(wd + f'data/gwas-pgs/{file_name}.ht')
    result = hl.read_table(wd + f'data/gwas-pgs/{file_name}.ht')
    results[run] = {
        'title': title,
        'result': result,
    }

pruns = [
    *runs
]

for run in pruns:
    print(run.name)
    result = results[run]['result']
    p = hl.plot.qq(
        pvals=result.p_value,
        title=results[run]['title'],
        size=4,
        collect_all=True,
    )
    show(p)
    df = result.order_by('p_value').drop('fit').to_pandas()
    if 'gene_name' in df:
        idx_column = 'gene_name'
    else:
        idx_column = 'burden_group'
    dff = df.set_index(idx_column)
    dff.index.name = ''
    dff.columns.name = idx_column
    display(dff.drop('group', axis=1).head(5))
    display(dff[dff.group == 'CYP'].drop('group', axis=1))


###################################
## Contingency table test for codes
###################################
code_overrep_data_path = 'data/diagnostic-code-overrepresentation/'
drug_related_codes_3 = (
    [f'T{i}' for i in range(36, 66)]
    + [f'X{i}' for i in range(40, 45)]
    + [f'Y{i}' for i in range(40, 60)]
    + ['Z88']
)
drug_related_codes_4 = (
    [f'T{i}{j}' for i in range(36, 66) for j in range(10)]
    + [f'X{i}{j}' for i in range(40, 45) for j in range(10)]
    + [f'Y{i}{j}' for i in range(40, 60) for j in range(10)]
    + [f'Z88{j}' for j in range(10)]
)
psychiatric_codes_3 = (
    [f'F{i:02}' for i in range(100)]
    + ['Z91']
)
psychiatric_codes_4 = (
    [f'F{i:02}{j}' for i in range(100) for j in range(10)]
    + [f'Z91{j}' for j in range(10)]
)

genes_list = pd.read_csv('raw/genes-list.csv')

n_codes = 4
b = a.select(
    all_codes=codes_to_set(a.f_41270[0], n_codes),
    genetic_ethnic_grouping=a.f_22006
)
b = b.checkpoint(next(tmpdir_path))

# filter out and annotate required table for analysis
cyps = genes_list[genes_list.Group == 'CYP'].Gene
cyps_without_2d6 = cyps[cyps != 'CYP2D6']
genes_to_analyse = (
    genes_list.Gene.to_list()
    + [cyps.to_list()]
    + [cyps_without_2d6.to_list()]
)
for i, gene in enumerate(genes_to_analyse):
    print(i, gene)
    if isinstance(gene, str):
        if gene not in mt.gene_name.collect():
            print("Gene not present in dataset")
            continue
        th = 50
        genes = [gene]
    else:
        th = 100
        genes = gene
    genes_id = ','.join(genes)
    result_path = f'{code_overrep_data_path}/{genes_id.replace(",", "-")}-{th}.ht'
    if os.path.exists(result_path):
        continue
    mt_genes = mt.select_cols()
    mt_genes = mt_genes.filter_rows(
        hl.literal(genes).contains(mt_genes.gene_name)
    )
    mt_genes = mt_genes.annotate_cols(
        **b[mt_genes.col_key]
    )
    mt_genes = mt_genes.filter_cols(
        mt_genes.genetic_ethnic_grouping == 'Caucasian'
    )
    mt_genes = mt_genes.entries()
    if not isinstance(gene, str):
        mt_genes = mt_genes.annotate(
            pgs=hl.if_else(
                mt_genes.pgs < 50, 0, mt_genes.pgs
            )
        )
    mt_genes = (
        mt_genes
        .group_by(mt_genes.s, mt_genes.all_codes)
        .aggregate(pgs=hl.agg.sum(mt_genes.pgs))
        .key_by('s')
    )
    mt_genes = mt_genes.annotate(genes_id=genes_id)
    mt_genes = mt_genes.annotate(
        high_pgs=mt_genes.pgs > th,
    )
    mt_genes = mt_genes.repartition(100)
    mt_genes = mt_genes.checkpoint(next(tmpdir_path))
    patients_counts = (
        mt_genes
        .group_by(
            mt_genes.genes_id,
            pgs=hl.if_else(mt_genes.high_pgs, 'non-functional', 'functional')
        )
        .aggregate(n=hl.agg.count())
    )
    patients_counts = (
        patients_counts
        .to_pandas()
        .set_index('genes_id')
        .pivot(columns='pgs')
        .droplevel(0, axis=1)
        .to_dict('index')
    )
    if 'non-functional' not in patients_counts[genes_id]:
        continue
    print(pd.DataFrame(patients_counts).T.iloc[:, [1, 0]])

    # explode codes and compute statistics
    mt_ct = mt_genes.explode('all_codes')
    mt_ct_grouped = mt_ct.group_by(group_code=mt_ct.all_codes)
    mt_ct = mt_ct_grouped.aggregate(
        nonfun_case=hl.agg.sum(mt_ct.high_pgs),
        fun_case=hl.agg.sum(~mt_ct.high_pgs)
    )
    mt_ct = mt_ct.annotate(
        nonfun_noncase=
            patients_counts[genes_id]['non-functional'] - mt_ct.nonfun_case,
        fun_noncase=
            patients_counts[genes_id]['functional'] - mt_ct.fun_case,
    )
    mt_ct = mt_ct.transmute(
        nonfun_case=hl.int32(mt_ct.nonfun_case),
        fun_case=hl.int32(mt_ct.fun_case),
        nonfun_noncase=hl.int32(mt_ct.nonfun_noncase),
        fun_noncase=hl.int32(mt_ct.fun_noncase),
    )
    # filter out code with small counts
    mt_ct = mt_ct.filter(
        hl.min(hl.array([
            mt_ct.nonfun_case,
            mt_ct.nonfun_noncase,
            mt_ct.fun_case,
            mt_ct.fun_noncase
        ])) > -1
    )
    mtp = mt_ct.annotate(
        **hl.contingency_table_test(
            mt_ct.nonfun_case,
            mt_ct.nonfun_noncase,
            mt_ct.fun_case,
            mt_ct.fun_noncase,
            min_cell_count=100
        )
    )
    mtp = mtp.annotate(
        d=hl.if_else(mtp.odds_ratio > 1, '^', '_')
    )
    mtp2 = mtp.checkpoint(wd + result_path)

    code_meaning = translate_icd10_codes(mtp2.group_code.collect())
    mtp2 = mtp2.annotate(
        meaning=
            code_meaning[mtp2.group_code].meaning,
        parent_meaning=
            code_meaning[mtp2.group_code.first_match_in(r'(.{3})')[0]].meaning
    )
    mtp2 = mtp2.checkpoint(next(tmpdir_path))
    if n_codes == 3:
        drug_related_codes = drug_related_codes_3
    elif n_codes == 4:
        drug_related_codes = drug_related_codes_4
    mtp2_filtered = mtp2.filter(
        hl.literal(drug_related_codes).contains(mtp2.group_code)
    )
    for m in [mtp2, mtp2_filtered]:
        show(hl.plot.qq(
            pvals=m.p_value,
            title=f'{genes_id} > {th}',
            size=4,
            collect_all=True,
        ))
        m.order_by(hl.asc('p_value')).show(20)


the_5_codes = ['T424', 'T432', 'T402', 'Y450', 'T426']
dfs = []
for ht_path in os.listdir(code_overrep_data_path):
    *gene, th_raw = ht_path.split('-')
    gene = '-'.join(gene)  # e.g. 'HLA-A' or gene group
    th = int(th_raw.split('.')[0])
    is_group = gene == '-'.join(cyps) or gene == '-'.join(cyps_without_2d6)
    if th == 51 and not is_group or th == 100 and is_group:
        pass
    else:
        continue
    ht = hl.read_table(wd + f'{code_overrep_data_path}/{ht_path}')
    ht = ht.filter(
        hl.literal(the_5_codes).contains(ht.group_code)
    )
    ht = ht.annotate(
        gene_name=gene,
    )
    ht = ht.checkpoint(next(tmpdir_path))

    df = ht.to_pandas()
    df = df.pivot('gene_name', 'group_code', ['p_value', 'odds_ratio'])

    dfs.append(df)

print(len(dfs))
df = pd.concat(dfs)

dff = df.reorder_levels(['group_code', None], axis=1)
dff = dff.sort_index(level=[0, 1], axis=1, ascending=[True, False])
dff.columns = (
    dff.columns.get_level_values(0)
    + '-'
    +  dff.columns.get_level_values(1)
)
dff.insert(0, 'gene_name', dff.index)
dff['sort_genes'] = (
    2 * dff.gene_name.str.startswith('CYP')
    - dff.gene_name.str.match('(CYP(.*)-?){2,}')
)
dff = dff.sort_values(by='sort_genes', ascending=False).drop('sort_genes', axis=1)
dff.to_csv(f'{code_overrep_data_path}/pgs-burden-genes-50-100.csv', index=False)

translate_icd10_codes(the_5_codes).show()


##  Export PGS burden results for histogram plotting
n_codes = 4
b = a.select(
    all_codes=codes_to_set(a.f_41270[0], n_codes),
    genetic_ethnic_grouping=a.f_22006
)

genes = pd.read_csv('raw/genes-list.csv').Gene.to_list()
mt_genes = mt.select_cols()
mt_genes = mt_genes.filter_rows(
    hl.literal(genes).contains(mt_genes.gene_name)
)
mt_genes = mt_genes.annotate_cols(
    **b[mt_genes.col_key]
)
ht = mt_genes.entries()
ht.drop('all_codes').export(wd + 'data/genes-list-pgs-burden.tsv')


## Export PGS burden for group of CYP genes
# (Figure 4 - barplots with mean PGS burden of CYPs for different
# comparison groups, ROC for other groups of genes)

genes_list = pd.read_csv('raw/genes-list.csv')
n_codes = 4
b = a.select(
    all_codes=codes_to_set(a.f_41270[0], n_codes),
    genetic_ethnic_grouping=a.f_22006
)
b = b.checkpoint(next(tmpdir_path))

cyps = genes_list[genes_list.Group == 'CYP'].Gene
cyps_without_2d6 = cyps[cyps != 'CYP2D6']
burden_genes = [
    ('cyps', cyps.to_list()),
    ('cyps_without_2d6', cyps_without_2d6.to_list()),
    ('list', genes_list.Gene.to_list()),
    ('all-protein', None)
]
hts = []
for i, (name, gene_group) in enumerate(burden_genes):
    mt_genes = mt
    if gene_group:
        mt_genes = mt_genes.filter_rows(
            hl.literal(gene_group).contains(mt_genes.gene_name)
        )
    mt_genes = mt_genes.annotate_cols(
        **b[mt_genes.col_key]
    )
    mt_genes = mt_genes.annotate_entries(
        pgs=hl.if_else(mt_genes.pgs < 50, 0, mt_genes.pgs)
    )
    mt_genes = (
        mt_genes
        .annotate_cols(pgss=hl.agg.sum(mt_genes.pgs))
        .select_cols(
            'group', 'mental_health_inpatient', 'self_harm', 'intention',
            'pgss'
        )
        .cols()
        .rename({'pgss': 'pgs'})
        .annotate(burden_group=name)
        .checkpoint(next(tmpdir_path))
    )
    hts.append(mt_genes)

ht = reduce(lambda hta, htb: hta.union(htb), hts)
ht.export(
    wd + 'results/manuscript/pgs-burden-groups-comparison-groups.tsv'
)

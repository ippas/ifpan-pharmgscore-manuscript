import os

import pandas as pd
import hail as hl
from hail.plot import show
import statsmodels.stats.multitest as sms

from analysis.utils.load_spark import wd
from analysis.utils.load_spark import hl_init
from analysis.utils.load_spark import tmpdir_path_iter
from analysis.utils.pathogenicity_scores import PharmGScore
from analysis.utils.annotations import translate_icd10_codes


hl_init()
hl.plot.output_notebook()
tmpdir_path = tmpdir_path_iter()
hl.spark_context().setLogLevel('WARN')


## Supplementary Table A - PharmGScore

pgs = PharmGScore.read_table()
pgs = pgs.select('pgs')
pgs.export(wd + 'results/supplementary/pgs.tsv.bgz')


## Supplementary Table E and Figure 1 - GWAS PGS burden on genes, qqplots

for file_path in os.listdir('data/gwas-pgs'):
    result = hl.read_table(wd + 'data/gwas-pgs/' + file_path)

    tables_dict = {
        'all': result,
        'cyps': result.filter(result.group == 'CYP'),
        'adrel': result.filter(hl.is_defined(result.group))
    }

    for key, res in tables_dict.items():
        print(file_path, key)
        res = res.filter(hl.is_defined(res.p_value))
        p = res.p_value.collect()
        _, p_fdr = sms.fdrcorrection(p)
        df_fdr = pd.DataFrame({
            'gene_name': res.gene_name.collect(),
            'p_value_fdr': p_fdr
        })
        ht_fdr = hl.Table.from_pandas(df_fdr, key='gene_name')
        res = res.annotate(**ht_fdr[res.key]).order_by('p_value')
        res = res.order_by('p_value')
        res.export(
            wd + f'results/supplementary/gwas-pgs/{file_path[:-3]}-{key}.tsv'
        )

        if key != 'cyps':
            p = hl.plot.qq(
                pvals=res.p_value,
                title=res.title.collect()[0],
                size=4,
                collect_all=True,
            )
            show(p)


## Supplementary Table D - diagnostic code overrepresentation

genes_list = pd.read_csv('raw/genes-list.csv')
cyps = genes_list[genes_list.Group == 'CYP'].Gene
cyps_without_2d6 = cyps[cyps != 'CYP2D6']
the_5_codes = ['T424', 'T432', 'T402', 'Y450', 'T426']

code_overrep_data_path = 'data/diagnostic-code-overrepresentation/'

genes_to_analyse = (
    list(zip(genes_list.Gene, genes_list.Group))
    + [(cyps.to_list(), 'CYP')]
    + [(cyps_without_2d6.to_list(), 'CYP')]
)
feature_results = []
for i, (gene, group) in enumerate(genes_to_analyse):
    print(i, gene)
    if isinstance(gene, str):
        th = 50
        genes = [gene]
    else:
        th = 100
        genes = gene
    genes_id = '-'.join(genes)
    ht = hl.read_table(wd + f'{code_overrep_data_path}/{genes_id}-{th}.ht')
    ht = ht.drop('d')
    ht = ht.annotate(feature=genes_id, group=group)
    feature_results.append(ht)
diag_code_overrep_ht = hl.Table.union(*feature_results)

diag_code_overrep_ht = diag_code_overrep_ht.filter(
    hl.literal(the_5_codes).contains(diag_code_overrep_ht.group_code)
)
diag_code_overrep_ht = diag_code_overrep_ht.select(
    'feature', 'p_value', 'odds_ratio', 'nonfun_case',
    'fun_case', 'nonfun_noncase', 'fun_noncase', 'group'
)

code_meaning = translate_icd10_codes(ht.group_code.collect())
diag_code_overrep_ht = diag_code_overrep_ht.annotate(
    meaning=
        code_meaning[diag_code_overrep_ht.group_code].meaning,
    parent_meaning=
        code_meaning[diag_code_overrep_ht.group_code.first_match_in(r'(.{3})')[0]].meaning
)
diag_code_overrep_ht = diag_code_overrep_ht.key_by('group_code', 'feature')
diag_code_overrep_ht = diag_code_overrep_ht.checkpoint(next(tmpdir_path))

col_order = (
    'group_code', 'feature', 'p_value', 'p_value_fdr', 'odds_ratio',
    'nonfun_case', 'fun_case', 'nonfun_noncase', 'fun_noncase',
    'meaning', 'parent_meaning'
)
ht_cyp = diag_code_overrep_ht.filter(
    diag_code_overrep_ht.group == 'CYP'
)
df_code_fdrs = []
for g_code in set(ht_cyp.group_code.collect()):
    ht = ht_cyp.filter(ht_cyp.group_code == g_code)
    p = ht.p_value.collect()
    _, p_fdr = sms.fdrcorrection(p)
    df_code_fdrs.append(
        pd.DataFrame({
            'group_code': ht.group_code.collect(),
            'feature': ht.feature.collect(),
            'p_value_fdr': p_fdr
        })
    )
    df_fdr = pd.concat(df_code_fdrs)
ht_fdr = hl.Table.from_pandas(df_fdr, key=('group_code', 'feature'))
ht_cyp = ht_cyp.annotate(**ht_fdr[ht_cyp.key])
ht_cyp = ht_cyp.order_by('group_code', 'p_value')
ht_cyp = ht_cyp.select(*col_order)
ht_cyp.export(
wd + 'results/supplementary/diagnostic-code-overrepresentation-sheetA.tsv'
)

ht_no_cyp = diag_code_overrep_ht.filter(
diag_code_overrep_ht.group != 'CYP'
)
df_code_fdrs = []
for g_code in set(ht_no_cyp.group_code.collect()):
    ht = ht_no_cyp.filter(ht_no_cyp.group_code == g_code)
    p = ht.p_value.collect()
    _, p_fdr = sms.fdrcorrection(p)
    df_code_fdrs.append(
        pd.DataFrame({
            'group_code': ht.group_code.collect(),
            'feature': ht.feature.collect(),
            'p_value_fdr': p_fdr
        })
    )
    df_fdr = pd.concat(df_code_fdrs)
ht_fdr = hl.Table.from_pandas(df_fdr, key=('group_code', 'feature'))
ht_no_cyp = ht_no_cyp.annotate(**ht_fdr[ht_no_cyp.key])
ht_no_cyp = ht_no_cyp.order_by('group_code', 'p_value')
ht_no_cyp = ht_no_cyp.select(*col_order[:-2], 'group', *col_order[-2:])
ht_no_cyp.export(
    wd + 'results/supplementary/diagnostic-code-overrepresentation-sheetB.tsv'
)

## QQ plots data

result = hl.read_table(
    wd + 'data/gwas-pgs/adr-control.ht'
)
result = result.select('p_value', 'group')
result = result.annotate(
    adr_control=result.p_value,
    adr_control_ad_relevant=hl.if_else(
        hl.is_defined(result.group), result.p_value, hl.missing(hl.tfloat)
    )
)

result.export('data/gwas-pgs/pvals-for-qqplot.tsv')

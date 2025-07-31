"""Computes all scores and burden for PharmVar"""
# %%
import hail as hl
import pandas as pd
import numpy as np

from analysis.utils.pharmvar import PharmVar, is_sub
from analysis.pharm_ml_score.pmls import pmls_pl2
from analysis.utils.pathogenicity_scores import annotate_vcf_with_scores


# %%

pv_version = '6.0.10'
# PharmVar - skips star alleles without any variants
pv_db = PharmVar(pv_version)
pv = pv_db.read_table()

pv = annotate_vcf_with_scores(pv)  #  *with all scores
pv_data = pv.filter(
    pv.is_sub
    & (pv.dist_to_exon <= 60)  # filter out missing
    & (hl.len(pv.alleles[0]) == 1)
    & (hl.len(pv.alleles[1]) == 1)
)
pv_data = pv_data.repartition(48)
pv_data.write(pv_db.table_filtered_path, overwrite=True)


# %%
pv_data = PharmVar(pv_version).read_table_filtered()

dataset = pv_data.transmute(
    cadd_nn=pv_data.scores.cadd_norm,
    fathmm_xf_nn=pv_data.scores.fathmm_xf_norm,
    pair0_nn=pv_data.scores.sift_norm,
    pair1_nn=pv_data.scores.phylop100_norm,
)
dataset_df = dataset.to_pandas()
dataset_df = dataset_df.drop(['locus', 'alleles'], axis=1)
out = pmls_pl2(
    dataset_df.copy(), id_columns=['gene_name', 'star_allele', 'function'],
    models_dir_prefix='data/pmls-models/sift-phylop100-',
    genes=set(dataset_df.gene_name)
)
out.to_csv(f'tmp/pharmvar-{pv_version}-pmls_pred.ht', index=False)

# %% gather all data
pv_data = PharmVar(pv_version).read_table_filtered()

data_roc = pv_data.transmute(
    pgs=pv_data.scores.pgs,
    adme=pv_data.scores.adme,
    adme2=pv_data.scores.adme2,
    cadd=pv_data.scores.cadd,
    fathmm_xf=pv_data.scores.fathmm_xf,
    ma=pv_data.scores.ma,
    provean=pv_data.scores.provean,
    sift=pv_data.scores.sift,
    phylop100=pv_data.scores.phylop100,
)

roc_scores = ['pgs', 'adme', 'adme2', 'cadd', 'fathmm_xf', 'ma', 'provean', 'sift', 'phylop100']
data_roc = data_roc.group_by('gene_name', 'star_allele', 'function').aggregate(
    burden=hl.struct(**{
        score: hl.agg.sum(data_roc[score]) for score in roc_scores
    }),
    all_na=hl.struct(**{
        score: hl.agg.all(
            hl.is_missing(data_roc[score])
            | hl.is_nan(data_roc[score])
        ) for score in roc_scores
    }),
)
data_roc = data_roc.annotate(
    burden=hl.struct(**{
        score: hl.if_else(data_roc.all_na[score], -np.inf, data_roc.burden[score]) for score in roc_scores
    })
)
data_roc = data_roc.drop('all_na')
data_roc = data_roc.transmute(
    **data_roc.burden
)
roc_df = data_roc.to_pandas()


pv = PharmVar(pv_version).read_table()
pv_data_all_sub = pv.filter(is_sub(pv))
all_sub_df = (
    pv_data_all_sub
    .key_by('gene_name', 'star_allele', 'function')
    .select()
    .distinct()
    .to_pandas()
)

out = pd.read_csv(f'tmp/pharmvar-{pv_version}-pmls_pred.ht')
final_df = all_sub_df.merge(
    roc_df.merge(out[['star_allele', 'pmls_pred']], on='star_allele'),
    how='left'
)
final_df = final_df.fillna(-np.inf)
final_df.to_csv(f'data/pharmvar-{pv_version}-final.csv')


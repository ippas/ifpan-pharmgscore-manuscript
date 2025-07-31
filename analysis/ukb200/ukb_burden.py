"""Compute burden on UKB, annotate with polygenic"""
# %%
import hail as hl
import polars as pl

from analysis.pharm_ml_score.pmls import pmls_pl2
from analysis.utils.pharmvar import PV_GENES
from analysis.utils.ukb_const import PHENOTYPE_MAP


# %%
# minimal in a ukb200-pgx
NEUTRAL_SCORE = {
    'pgs': 0,
    'adme': 0,
    'adme2': 0,
    'cadd': -2.63,
    'fathmm_xf': 0,
    'provean': -5.78,
    'ma': -3.04,
    'sift': -1,
    'phylop100': -14,
}

# %%
mt = hl.read_matrix_table(vcf_data_path + 'ukb200-pgx.mt')
mt = mt.filter_rows(
    (mt.dist_to_exon <= 60)
    & (hl.len(mt.alleles[0]) == 1)
    & (hl.len(mt.alleles[1]) == 1)
    & hl.literal(PV_GENES).contains(mt.gene_name)
).persist()
mtr = mt.rows().persist()

compute_burden = lambda mt, score: hl.agg.sum(score * mt.GT.n_alt_alleles())
compute_burden = lambda mt, score, name: hl.if_else(
    hl.agg.all(hl.is_missing(score)),
    hl.literal(NEUTRAL_SCORE)[name],
    hl.agg.sum(score * mt.GT.n_alt_alleles())
)

mt_burden = mt.group_rows_by('gene_name').aggregate(
    n_alt_alleles=hl.agg.sum(mt.GT.n_alt_alleles()),
    pgs=compute_burden(mt, mt.scores.pgs, 'pgs'),
    burden_adme=compute_burden(mt, mt.scores.adme, 'adme'),
    burden_adme2=compute_burden(mt, mt.scores.adme2, 'adme2'),
    burden_cadd=compute_burden(mt, mt.scores.cadd, 'cadd'),
    burden_fathmm_xf=compute_burden(mt, mt.scores.fathmm_xf, 'fathmm_xf'),
    burden_provean=compute_burden(mt, mt.scores.provean, 'provean'),
    burden_ma=compute_burden(mt, mt.scores.ma, 'ma'),
    burden_sift=compute_burden(mt, mt.scores.sift, 'sift'),
    burden_phylop100=compute_burden(mt, mt.scores.phylop100, 'phylop100'),
)
ht_burden = mt_burden.entries()
ht_burden.write(vcf_data_path + 'burden-scores.ht')

compute_nn_score = lambda mt, score_expr: hl.if_else(
        hl.agg.all(hl.is_missing(score_expr)),
        hl.missing(hl.tfloat64),
        hl.agg.sum(score_expr * mt.GT.n_alt_alleles()),
    )
ht_non_ref = mt.entries()
ht_burden_nn = ht_non_ref.group_by('gene_name', 's', 'locus', 'alleles').aggregate(
    cadd_nn=compute_nn_score(ht_non_ref, ht_non_ref.scores.cadd_norm),
    fathmm_xf_nn=compute_nn_score(ht_non_ref, ht_non_ref.scores.fathmm_xf_norm),
    pair0_nn=compute_nn_score(ht_non_ref, ht_non_ref.scores.sift_norm),
    pair1_nn=compute_nn_score(ht_non_ref, ht_non_ref.scores.phylop100_norm),
)
ht_burden_nn.export(vcf_data_path + 'burden-nn.tsv')

# %%
dataset_df = pl.read_csv(
    vcf_data_path + 'burden-nn.tsv',
    separator='\t',
    schema_overrides={'s': pl.String},
    null_values='NA'
)
ukb_df = dataset_df.to_pandas()
out = pmls_pl2(
    ukb_df,
    id_columns=['gene_name', 's'],
    models_dir_prefix='data/pmls-models/sift-phylop100-',
    genes=set(ukb_df.gene_name)
)
out_ht = hl.Table.from_pandas(out, key=['gene_name', 's']).persist()
out_ht.aggregate(hl.agg.min(out_ht.pmls_pred))


# %%
ht_burden = hl.read_table(vcf_data_path + 'burden-scores.ht')
out_ht = hl.read_table(vcf_data_path + 'pmls_pred.ht')
poly_ht = hl.read_table('data/polygenic/polygenic-results.csv')

burden = ht_burden.annotate(
    pmls_pred=out_ht[ht_burden.key].pmls_pred,
    poly_phenotype=poly_ht[ht_burden.key].poly_phenotype,
    poly_label=poly_ht[ht_burden.key].poly_label,
)
burden = burden.annotate(
    pmls_pred=hl.if_else(hl.is_missing(burden.pmls_pred), -2.7, burden.pmls_pred)  # -2.72 min score for pmls
)
burden.show()


# %%
burden = hl.read_table(vcf_data_path + 'burden.ht')
mt = burden.to_matrix_table(
    row_key=['gene_name'],
    col_key=['s'],
)
mt = mt.annotate_entries(
    poly=hl.literal(PHENOTYPE_MAP).get(mt.poly_phenotype),
)
mt.write(vcf_data_path + 'burden.mt', overwrite=True)

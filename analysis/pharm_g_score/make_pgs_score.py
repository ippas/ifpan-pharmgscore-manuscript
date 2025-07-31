# %%
import os

import numpy as np
import polars as pl
import hail as hl
import pandas as pd
from tqdm import tqdm
from scipy.stats import rankdata

from analysis.utils.pathogenicity_scores import PharmGScore

# %%
make_score_path = os.path.join('tmp', 'make-score.csv')


# %% Compute PGS -----------------------------------------------------------------

BURDEN_SCORE_NAMES = [
  "cadd_raw", "fathmm_xf", "provean_converted", "mutation_assessor"
]
def normalize_quantiles(x, target=None):
    """
    Quantile normalization of a vector based upon a specified target distribution.
    """
    x = rankdata(x, method='dense', nan_policy='omit')
    x = rankdata(x, method='min', nan_policy='omit')  # true n_q
    y = (x - np.nanmin(x)) / (np.nanmax(x) - np.nanmin(x))

    return np.exp(y * 10 - 5) - np.exp(-5)


def make_score_pl(gene_names, gene_scores_path, dist_exon_th=np.inf, dist_gene_th=np.inf):
    gene_scores = []

    for g_name in tqdm(gene_names):
        g_path = os.path.join(gene_scores_path, f'{g_name}.tsv.bgz')

        if os.path.exists(g_path):
            df = pl.read_csv(g_path, separator='\t', null_values='NA')

            df = df.with_columns(
                pl.lit(g_name).alias('gene_name')
            )

            for column in BURDEN_SCORE_NAMES:
                df = df.with_columns(
                    pl.col(column).alias(f'{column}_oryg')
                )

            df = df.filter(
                (pl.col('dist_to_exon') <= dist_exon_th) &
                (pl.col('dist_to_gene') <= dist_gene_th)
            )

            for column in BURDEN_SCORE_NAMES:
                df = df.with_columns(
                    pl.Series(
                        f'{column}_pgs',
                        normalize_quantiles(df[column].cast(pl.Float32)),
                        nan_to_null=True
                    ).cast(pl.Float32)
                )

            df_selected = df.select(
                ['locus', 'alleles', 'gene_name'] \
                + [f'{c}_pgs' for c in BURDEN_SCORE_NAMES]
            )

            gene_scores.append(df_selected)
        else:
            raise FileNotFoundError(f"File: {g_path} doesn't exist")

    all_scores = pl.concat(gene_scores)

    return all_scores


if __name__ == '__main__':
    hl_init()

    genes_positions_path = 'data/pharmacogenetic-score/genes-position.tsv'
    genes_positions = pd.read_csv(
        genes_positions_path,
        sep='\t',
        low_memory=False
    )
    genes = set(genes_positions['Gene name'])

    dist_th = 60
    df_pl = make_score_pl(
        genes,
        gene_scores_path='data/gene-scores',
        dist_exon_th=dist_th
    )
    df_pl = df_pl.with_columns(
        pl.col('alleles').str.extract(r'([ACGT]+).*([ACGT]+)', 1).alias('allele_1'),
        pl.col('alleles').str.extract(r'([ACGT]+).*([ACGT]+)', 2).alias('allele_2'),
    )
    df_pl = df_pl.drop('alleles')
    df_pl.write_csv(make_score_path, separator='\t')

    pgs = hl.import_table(
        make_score_path,
        missing='',
        types={
            'locus': hl.tstr,
            'allele_1': hl.tstr,
            'allele_2': hl.tstr,
            'gene_name': hl.tstr,
            'cadd_raw_pgs': hl.tfloat32,
            'fathmm_xf_pgs': hl.tfloat32,
            'provean_converted_pgs': hl.tfloat32,
            'mutation_assessor_pgs': hl.tfloat32,
        },
    )
    pgs = pgs.annotate(
        scores=hl.array([
            pgs.cadd_raw_pgs,
            pgs.fathmm_xf_pgs,
            pgs.provean_converted_pgs,
            pgs.mutation_assessor_pgs,
        ])
    )
    pgs = pgs.transmute(
        pgs=hl.if_else(
            hl.all(hl.map(hl.is_missing, pgs.scores)),
            0,
            hl.mean(pgs.scores)
        )
    )
    pgs = pgs.annotate(
        **hl.parse_variant(
            pgs.locus + ':' + pgs.allele_1 + ':' + pgs.allele_2
        )
    )
    pgs = pgs.drop('allele_1', 'allele_2')
    pgs = pgs.key_by('locus', 'alleles')
    pgs = pgs.persist()

    # PharmGScore is per gene
    pgs_collect = pgs.group_by(pgs.locus, pgs.alleles).aggregate(
        **{
            col: hl.agg.collect(pgs[col]) for col in [
                'gene_name',
                'cadd_raw_pgs',
                'fathmm_xf_pgs',
                'provean_converted_pgs',
                'mutation_assessor_pgs',
                'pgs',
            ]
        }
    )
    pgs_collect = pgs_collect.persist()

    # However as stated in hail docs 'The element order of the resulting array is
    # not guaranteed, and in some cases is non-deterministic.' We need to check it.

    compare_pgs = pgs.annotate(
        collected=pgs_collect[pgs.key]
    )
    compare_pgs = compare_pgs.annotate(
        pgs_s=compare_pgs.pgs,
        pgs_c=compare_pgs.collected.pgs[compare_pgs.collected.gene_name.index(compare_pgs.gene_name)],
    )
    compare_pgs = compare_pgs.annotate(
        same_pgs=(compare_pgs.pgs_s == compare_pgs.pgs_c)
    )
    compare_pgs = compare_pgs.persist()
    compare_pgs.show()

    compare_pgs_count = compare_pgs.filter(~compare_pgs.same_pgs).count()
    if compare_pgs_count == 0:
        pgs_collect.write(PharmGScore.table_path, overwrite=True)

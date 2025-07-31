"""Save genes and their scores for a later quantile normalization"""

import pandas as pd
import hail as hl
from tqdm import tqdm

from analysis.utils.pathogenicity_scores import CADD, FathmmXF, dbNSFP
from analysis.utils.annotations import query_biomart


all_scores_path = 'data/all-pgs-raw-scores.ht'
genes_positions_path = 'data/pharmacogenetic-score/genes-position.tsv'


def annotate_all_possible_variants_with_scores(all_scores_path):
    """Annotate all possible variants with scores (based on CADD table)"""
    dbnsfp = dbNSFP.read_table(single=True)
    fathmm_xf = FathmmXF.read_table()
    cadd = CADD.read_table()

    gnomad = hl.read_table(
        RESOURCES_PATH + '/gnomad/gnomad.genomes.v3.1.1.sites.ht'
    )
    nfe = gnomad.freq_meta.collect()[0].index({'group': 'adj','pop': 'nfe'})

    cadd = cadd.select('score_raw')
    template = cadd.rename({
        'score_raw': 'cadd_raw',
    })
    all_scores = template.annotate(
        mutation_assessor=dbnsfp[template.key].single.MutationAssessor_score,
        provean_converted=-dbnsfp[template.key].single.PROVEAN_score,
        fathmm_xf=fathmm_xf[template.key].score,
        gnomad_nfe_AF=gnomad[template.key].freq[nfe].AF,
        in_gnomad=hl.is_defined(gnomad[template.key]),
    )

    all_scores.write(all_scores_path)


def download_genes_and_exons_positions(genes_positions_path):
    genes_positions = query_biomart(
        dataset='hsapiens_gene_ensembl',
        filters={
            'biotype': 'protein_coding',
        },
        attributes=[
            'chromosome_name',
            'external_gene_name',
            'start_position',
            'end_position',
            'exon_chrom_start',
            'exon_chrom_end',
        ],
        # version: 'http://mar2022.archive.ensembl.org/biomart/'
        url='http://ensembl.org/biomart/'
    )
    genes_positions = genes_positions[
        ~genes_positions['Chromosome/scaffold name'].str.match('CHR_|MT|GL|KI')
        & ~genes_positions['Gene name'].isna()
    ]
    genes_positions.to_csv(genes_positions_path, sep='\t', index=False)


if __name__ == '__main__':
    annotate_all_possible_variants_with_scores(all_scores_path)
    download_genes_and_exons_positions(genes_positions_path)

    chrom # set chrom
    chrom = {'23': 'X', '24': 'Y'}.get(chrom, chrom)
    print(f'Chromosome: {chrom}')

    genes_positions = pd.read_csv(
        genes_positions_path,
        sep='\t',
        low_memory=False
    )
    genes_positions = genes_positions[
        genes_positions['Chromosome/scaffold name'].isin([chrom])
    ]

    # save scores for each gene in a separate file (around exons)
    gene_offset = 2500
    gene_names = genes_positions['Gene name'].unique()
    all_scores = hl.read_table(all_scores_path)

    for gene_name in tqdm(gene_names):
        gene_positions = genes_positions[genes_positions['Gene name'] == gene_name]
        gene_intervals = []
        exon_intervals = []
        exon_breaks = []
        for i, row in gene_positions.iterrows():
            contig = row['Chromosome/scaffold name']
            gene_start = row['Gene start (bp)']
            gene_end = row['Gene end (bp)']
            exon_start = row['Exon region start (bp)']
            exon_end = row['Exon region end (bp)']
            exon_breaks.extend([exon_start, exon_end])
            gene_interval_offset = [
                hl.parse_locus_interval(
                    f'chr{contig}'
                    f':{gene_start - gene_offset}'
                    f'-{gene_end + gene_offset + 1}'
                )
            ]
            gene_interval = [
                hl.parse_locus_interval(
                    f'chr{contig}:{gene_start}-{gene_end + 1}'
                )
            ]
            exon_intervals.append(
                hl.parse_locus_interval(
                    f'chr{contig}:{exon_start}-{exon_end + 1}')
            )
        gene_scores = hl.filter_intervals(all_scores, gene_interval_offset)
        gene_interior = hl.filter_intervals(gene_scores, gene_interval)
        exon_interior = hl.filter_intervals(gene_scores, exon_intervals)
        gene_interior = gene_interior.annotate(
            is_in=True,
        )
        exon_interior = exon_interior.annotate(
            is_in=True
        )
        gene_scores = gene_scores.annotate(
            dist_to_gene=hl.min(
                hl.abs(
                    [gene_scores.locus.position - x for x in [gene_start, gene_end]]
                )
            ),
            dist_to_exon=hl.min(
                hl.abs(
                    [gene_scores.locus.position - x for x in exon_breaks]
                )
            ),
            is_gene_interior=gene_interior[gene_scores.key].is_in,
            is_exon_interior=exon_interior[gene_scores.key].is_in,
        )
        gene_scores = gene_scores.transmute(
            dist_to_gene=hl.if_else(
                hl.is_defined(gene_scores.is_gene_interior),
                0,
                gene_scores.dist_to_gene
            ),
            dist_to_exon=hl.if_else(
                hl.is_defined(gene_scores.is_exon_interior),
                0,
                gene_scores.dist_to_exon
            )
        )
        gene_scores.export(
            f'data/gene-scores/{gene_name}.tsv.bgz'
        )

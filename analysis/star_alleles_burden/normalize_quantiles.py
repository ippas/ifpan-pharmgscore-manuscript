"""Save genes for a later quantile normalization"""

import os
import re
from collections import defaultdict
from operator import itemgetter

import pandas as pd
import numpy as np
import hail as hl
from hail.plot import show
from tqdm import tqdm

from analysis.utils.load_spark import hl_init
from analysis.utils.load_spark import wd
from analysis.utils.load_spark import tmpdir_path_iter
from analysis.utils.pathogenicity_scores import CADD, FathmmXF, dbNSFP, DANN
from analysis.utils.annotations import query_biomart
from analysis.star_alleles_burden.pharmvar import PharmVar
from analysis.utils.annotations import query_biomart

hl_init()
hl.plot.output_notebook()
tmpdir_path = tmpdir_path_iter()

# annotate all possible variants with scores (based on CADD table)
all_scores_path = wd + 'data/pharmacogenetic-score/all-scores.ht'
if not os.path.exists(all_scores_path):
    dbnsfp = dbNSFP().convert_subset(
        dbnsfp_subset_path=wd + 'data/scores/dbNSFP4.3a-subset.ht'
    )
    fathmm_xf = FathmmXF.read_table()
    dann = DANN.read_table()
    cadd = CADD.read_table()

    cadd = cadd.select('score_raw')
    template = cadd.rename({
        'score_raw': 'cadd_raw',
    })
    all_scores = template.annotate(
        mutation_assessor=dbnsfp[template.key].mutation_assessor,
        provean=dbnsfp[template.key].provean,
        fathmm_xf=fathmm_xf[template.key].score,
        dann=dann[template.key].score,
    )

    # TODO: next time move below annotations to the above function
    gnomad = hl.read_table(
        '/net/archive/groups/plggneuromol/resources/gnomad/gnomad.genomes.v3.1.1.sites.ht'
    )

    nfe = gnomad.freq_meta.collect()[0].index({"group":"adj","pop":"nfe"})
    gnomad = gnomad.annotate(
        gnomad=True,
        AF=gnomad.freq[nfe].AF
    )

    all_scores = all_scores.annotate(
        gnomad_nfe_AF=gnomad[all_scores.key].AF,
        in_gnomad=gnomad[all_scores.key].gnomad,
    )
    all_scores = all_scores.annotate(
        in_gnomad=hl.if_else(
            hl.is_missing(all_scores.in_gnomad),
            False,
            all_scores.in_gnomad
        )
    )

    all_scores.write(all_scores_path)
else:
    all_scores = hl.read_table(all_scores_path)

# download genes and exons positions
genes_positions_path = wd + 'data/pharmacogenetic-score/genes-position.tsv'
if not os.path.exists(genes_positions_path):
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
else:
    genes_positions = pd.read_csv(
        genes_positions_path,
        sep='\t',
        low_memory=False
    )

# SLURM data
chrom = os.environ.get('SLURM_ARRAY_TASK_ID')
chrom = {'23': 'X', '24': 'Y'}.get(chrom, chrom)
print(f'Chromosome: {chrom}')
genes_positions = genes_positions[
    genes_positions['Chromosome/scaffold name'].isin([chrom])
]
#SBATCH --account plgdepresja
#SBATCH --partition plgrid
#SBATCH --array=1-24
#SBATCH --ntasks-per-node=2
#SBATCH --mem=61GB
#SBATCH --time 4:00:00

# save scores for each gene in a separate file (around exons)
gene_offset = 2500
gene_names = genes_positions['Gene name'].unique()
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
        wd + f'data/pharmacogenetic-score/gene-scores/{gene_name}.tsv.bgz'
    )

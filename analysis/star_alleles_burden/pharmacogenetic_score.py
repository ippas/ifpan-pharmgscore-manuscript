"""Creates new score optimized for ADME genes"""

import os
import re
from collections import defaultdict

import pandas as pd
import numpy as np
import hail as hl
from hail.plot import show

from analysis.utils.load_spark import hl_init
from analysis.utils.load_spark import wd
from analysis.utils.load_spark import tmpdir_path_iter
from analysis.utils.pathogenicity_scores import CADD, FathmmXF, dbNSFP, DANN
from analysis.utils.annotations import query_biomart
from analysis.star_alleles_burden.pharmvar import PharmVar
from analysis.utils.annotations import query_biomart


def phred_score(ht, score_col_name, order='desc', tmp_dir=None):
    ht = ht.filter(hl.is_defined(ht[score_col_name]))
    if order == 'desc':
        ht = ht.order_by(hl.desc(ht[score_col_name]))
    elif order == 'asc':
        ht = ht.order_by(hl.asc(ht[score_col_name]))
    else:
        raise ValueError('Order must be either "desc" or "asc"')
    ht = ht.add_index('score_index')
    if tmp_dir:
        ht = ht.checkpoint(tmp_dir + '/ht1a')
        ht = ht.repartition(ht.n_partitions())
        ht = ht.checkpoint(tmp_dir + '/ht1b')

    ht_ranks = (
        ht
        .group_by(ht[score_col_name])
        .aggregate(
            score_min_rank=hl.agg.min(ht.score_index) + 1,
        )
    )
    if tmp_dir:
        ht_ranks = ht_ranks.checkpoint(tmp_dir + '/ht2')
    ht = ht.drop('score_index')

    n = ht_ranks.aggregate(hl.agg.max(ht_ranks.score_min_rank))
    ht = ht.annotate(
        percent_rank=ht_ranks[ht[score_col_name]].score_min_rank / n,
    )
    ht = ht.transmute(
        phred=-10 * hl.log10(ht.percent_rank)
    )
    if tmp_dir:
        ht = ht.checkpoint(tmp_dir + '/ht3')

    return ht


def position_in_interval(locus, interval):
    """Return True if locus is in interval"""
    return (locus.position > interval.left) & (locus.position <= interval.right)


def case_builder(intervals, locus):
    """Build cases for gene annotation"""
    cases = hl.case()
    for (contig, gene), c_intervals in intervals.items():
        for interval in c_intervals:
            cases = cases.when(
                (locus.contig == contig) & position_in_interval(locus, interval),
                gene
            )
    return cases.or_missing()


hl_init()
hl.plot.output_notebook()
tmpdir_path = tmpdir_path_iter()


dbnsfp = dbNSFP().convert_subset(
    dbnsfp_subset_path=wd + 'data/scores/dbNSFP4.3a-subset.ht'
)
pv_vcf = PharmVar.read_table()
fathmm_xf = FathmmXF.read_table()
dann = DANN.read_table()
cadd = CADD.read_table()

# write pharmacogenetic gene start and end positions
intervals_path = 'data/pharmvar/intevals.list'

gene_positions = query_biomart(
    dataset='hsapiens_gene_ensembl',
    filters={
        'external_gene_name': list(set(pv_vcf.gene.collect()))
    },
    attributes=[
        'chromosome_name',
        'external_gene_name',
        'start_position',
        'end_position',
    ],
    # version: 'http://mar2022.archive.ensembl.org/biomart/'
    url='http://ensembl.org/biomart/'
)
gene_positions = gene_positions[
    ~gene_positions['Chromosome/scaffold name'].str.startswith('CHR_HSCHR')
]

offset = 2500  # consist all variants from PharmVar
intervals = defaultdict(list)
for _, gene_info in gene_positions.iterrows():
    contig = 'chr' + gene_info['Chromosome/scaffold name']
    gene_name = gene_info['Gene name']
    i = pd.Interval(
        gene_info['Gene start (bp)'] - offset,
        gene_info['Gene end (bp)'] + offset
    )
    intervals[(contig, gene_name)].append(i)

with open(intervals_path, 'w') as f:
    for (contig, gene), c_intervals in intervals.items():
        for interval in c_intervals:
            f.write(f'{contig}:{interval.left}-{interval.right}\n')

interval_table = hl.import_locus_intervals(intervals_path)

# Phred scores for pharmacogenetic variants
table_score_zip = [
    (dbnsfp, 'mutation_assessor', 'mutation-assessor'),
    (dbnsfp, 'provean', 'provean'),
    (fathmm_xf, 'score', 'fathmm-xf'),
    (dann, 'score', 'dann'),
    (cadd, 'score_raw', 'cadd'),
]
wd_ps = wd + 'data/pharmacogenetic-score/scores-phred/'
for score_ht, col_name, file_name in table_score_zip:
    for kind in ['all', 'interval']:  # order is important
        if kind == 'interval':
            score_ht = score_ht.filter(
                hl.is_defined(interval_table[score_ht.locus])
            )
        order = 'asc' if col_name == 'provean' else 'desc'
        score_phred_ht = phred_score(
            score_ht, col_name, order, tmp_dir=next(tmpdir_path)
        )
        (
            score_phred_ht
            .key_by('locus', 'alleles')
            .write(wd_ps + f'{file_name}-phred-{kind}.ht', overwrite=True)
        )

# annotate PharmVar variants
scores_dict = defaultdict(dict)
for _, _, file_name in table_score_zip:
    for kind in ['interval', 'all']:
        scores_dict[file_name][kind] = \
            hl.read_table(wd_ps + f'{file_name}-phred-{kind}.ht')

pv_vcf = pv_vcf.annotate(
    mutation_assessor= \
        dbnsfp[pv_vcf.key].mutation_assessor,
    mutation_assessor_phred= \
        -10 * hl.log10(1 - dbnsfp[pv_vcf.key].MutationAssessor_rankscore),
    provean= \
        dbnsfp[pv_vcf.key].provean,
    provean_phred= \
        -10 * hl.log10(1 - dbnsfp[pv_vcf.key].PROVEAN_converted_rankscore),

    mutation_assessor_phred_inter= \
        scores_dict['mutation-assessor']['interval'][pv_vcf.key].phred,
    mutation_assessor_phred_all= \
        scores_dict['mutation-assessor']['all'][pv_vcf.key].phred,
    provean_phred_inter= \
        scores_dict['provean']['interval'][pv_vcf.key].phred,
    provean_phred_all= \
        scores_dict['provean']['all'][pv_vcf.key].phred,

    fathmm_xf=fathmm_xf[pv_vcf.key].score,
    cadd_raw=cadd[pv_vcf.key].score_raw,
    cadd_phred=cadd[pv_vcf.key].score_phred,
    dann=dann[pv_vcf.key].score,

    fathmm_xf_phred_inter=scores_dict['fathmm-xf']['interval'][pv_vcf.key].phred,
    fathmm_xf_phred_all=scores_dict['fathmm-xf']['all'][pv_vcf.key].phred,
    cadd_phred_inter=scores_dict['cadd']['interval'][pv_vcf.key].phred,
    cadd_phred_all=scores_dict['cadd']['all'][pv_vcf.key].phred,
    dann_phred_inter=scores_dict['dann']['interval'][pv_vcf.key].phred,
    dann_phred_all=scores_dict['dann']['all'][pv_vcf.key].phred,
)

pv_vcf = pv_vcf.checkpoint(next(tmpdir_path))
pv_vcf.show()
pv_vcf.export('data/pharmacogenetic-score/star-alleles-scores-phred.tsv.bgz')


# annotate healthy patients
healthy_path = os.path.join(
    '/net/archive/groups/plggneuromol/imdik-zekanowski-sportwgs/',
    'data/joint/full-healthy.mt'
)
healthy = hl.read_matrix_table(healthy_path)

healthy = healthy.filter_rows(
    hl.is_defined(interval_table[healthy.locus])
)
healthy = healthy.filter_entries(healthy.GT.is_non_ref())
healthy = healthy.select_rows('rsid', 'a_index', 'was_split', 'within_gene')
healthy = healthy.select_cols('group')
healthy = healthy.select_entries('GT')
healthy = healthy.repartition(500)
healthy_s = healthy.checkpoint(next(tmpdir_path))

hr = healthy_s.rows()
hr = hr.checkpoint(next(tmpdir_path))

dbnsfp = dbNSFP().convert_subset(
    dbnsfp_subset_path=wd + 'data/scores/dbNSFP4.3a-subset.ht'
)
fathmm_xf = FathmmXF.read_table()
dann = DANN.read_table()
cadd = CADD.read_table()
# speed up annotation
dbnsfp = dbnsfp.join(hr.select())
cadd = cadd.join(hr.select())
fathmm_xf = fathmm_xf.join(hr.select())
dann = dann.join(hr.select())

healthy_pg = healthy_s.annotate_rows(
    dbnsfp=hl.struct(
        mutation_assessor=dbnsfp[healthy_s.row_key].mutation_assessor,
        provean=dbnsfp[healthy_s.row_key].provean,
    ),
    cadd_raw=cadd[healthy_s.row_key].score_raw,
    fathmm_xf=fathmm_xf[healthy_s.row_key].score,
    dann=dann[healthy_s.row_key].score,
    gene_name=case_builder(intervals, healthy_s.locus)
)
healthy_pg.write(wd + 'data/pharmacogenetic-score/healthy-pg.mt', overwrite=True)

healthy_pg = hl.read_matrix_table(wd + 'data/pharmacogenetic-score/healthy-pg.mt')
healthy_pg = healthy_pg.entries()
healthy_pg = healthy_pg.annotate(
    n_alt=healthy_pg.GT.n_alt_alleles()
)
healthy_pg.to_pandas().to_csv(
    wd + 'data/pharmacogenetic-score/healthy-pg.csv',
    index=False
)

# save healthy patients and pharmvar variants with allele frequency annotations
# to tsv.bgz file
gnomad = hl.read_table(
    '/net/archive/groups/plggneuromol/resources/gnomad/gnomad.genomes.v3.1.1.sites.ht'
)
pv_vcf = PharmVar.read_table()
nfe = gnomad.freq_meta.collect()[0].index({"group":"adj","pop":"nfe"})

pv_vcf = pv_vcf.annotate(
    AF=gnomad[pv_vcf.key].freq[nfe].AF,
)
pv_vcf.export('data/pharmacogenetic-score/pharmvar.tsv.bgz')

healthy_path = os.path.join(
    '/net/archive/groups/plggneuromol/imdik-zekanowski-sportwgs/',
    'data/joint/full-healthy.mt'
)
healthy = hl.read_matrix_table(healthy_path)
healthy = healthy.annotate_rows(
    gt_stats=hl.agg.call_stats(healthy.GT, healthy.alleles)
)
healthy = healthy.transmute_rows(
    pop_af=healthy.gt_stats.AF[1],
)
healthy_nref = healthy.filter_entries(healthy.GT.is_non_ref())

healthy_nref = healthy_nref.annotate_entries(
    n_alt=healthy_nref.GT.n_alt_alleles()
)
healthy_nref = healthy_nref.select_rows('within_gene', 'pop_af')
healthy_nref = healthy_nref.select_cols('group')
healthy_nref = healthy_nref.select_entries('n_alt')
# healthy_nref = healthy_nref.key_cols_by()
healthy_nref = healthy_nref.entries()
healthy_nref = healthy_nref.filter(healthy_nref.pop_af < 0.2)

healthy_nref.export(wd + f'data/pharmacogenetic-score/healthy-non-ref-filtered.tsv.bgz')

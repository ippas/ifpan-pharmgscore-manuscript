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
from analysis.utils.pathogenicity_scores import PharmGScore
from analysis.utils.annotations import query_biomart
from analysis.star_alleles_burden.pharmvar import PharmVar
from analysis.utils.annotations import query_biomart, reorder_cols
from analysis.utils.annotations import find_annot_desc, load_annotations
from analysis.utils.annotations import translate_icd10_codes
from analysis.utils.variant_filtering import VCFFilter


hl_init()
hl.plot.output_notebook()
tmpdir_path = tmpdir_path_iter()

# Convert VCF to hail matrix table
vcf_annotated_path = '/net/scratch/people/plgjacekh/pgs-tmp/part-x-0'
mt = hl.import_vcf(
    # [
    #     wd + 'raw/pvcf-200k/ukb23156_c[1-9]_b*_v1.vcf.gz',
    #     wd + 'raw/pvcf-200k/ukb23156_c1[0-8]_b*_v1.vcf.gz',
    # ],
    # [
    #     wd + 'raw/pvcf-200k/ukb23156_c19_b?_v1.vcf.gz',
    #     wd + 'raw/pvcf-200k/ukb23156_c19_b[1-2]?_v1.vcf.gz',
    #     wd + 'raw/pvcf-200k/ukb23156_c19_b3[0-1]_v1.vcf.gz',
    # ],
    # [
    #     wd + 'raw/pvcf-200k/ukb23156_c19_b3[2-9]_v1.vcf.gz',
    #     wd + 'raw/pvcf-200k/ukb23156_c19_b[4-9]?_v1.vcf.gz',
    #     wd + 'raw/pvcf-200k/ukb23156_c2[0-2]_b*_v1.vcf.gz',
    #     wd + 'raw/pvcf-200k/ukb23156_c{X,Y}_b*_v1.vcf.gz',
    # ],
    force_bgz=True,
    array_elements_required=False,
)

# Filter VCF
mt_filter = VCFFilter()
mt = mt_filter.min_mean_read_depth(mt)
mt = hl.variant_qc(mt)
mt = mt_filter.min_variant_missingness(mt)
mt = hl.split_multi_hts(mt, permit_shuffle=True)
mt_filtered = mt_filter.min_allele_balance(mt)

# Annotate VCF
pgs_genes_path = wd + 'data/pharmacogenetic-score/intermediate/pgs_genes.ht'
if not os.path.exists(pgs_genes_path):
    pgs = PharmGScore.read_table()
    pgs_genes = pgs.group_by('locus', 'alleles').aggregate(
        gene_names=hl.agg.collect_as_set(pgs.gene_name)
    )
    pgs_genes.write(pgs_genes_path)
pgs_genes = hl.read_table(pgs_genes_path)

mt_filtered = mt_filtered.annotate_rows(
    **pgs_genes[mt_filtered.row_key]
)
mt_filtered_3key = (
    mt_filtered
    .explode_rows('gene_names')
    .rename({'gene_names': 'gene_name'})
    .key_rows_by('locus', 'alleles', 'gene_name')
)

mt_scores = mt_filtered_3key.annotate_rows(
    **PharmGScore.read_table()[mt_filtered_3key.row_key]
)
mt_scores = mt_scores.annotate_rows(
    pgs_scores=hl.array([
        mt_scores.cadd_raw,
        mt_scores.mutation_assessor,
        mt_scores.provean,
        mt_scores.fathmm_xf
    ])
)
mt_scores = mt_scores.annotate_rows(
    pgs=hl.if_else(
        hl.all(hl.map(hl.is_missing, mt_scores.pgs_scores)),
        0,
        hl.mean(mt_scores.pgs_scores)
    )
)
mt_scores.write(vcf_annotated_path)  # chr1-chr4 - 4 nodes - 19h

paths = [
    '/net/scratch/people/plgjacekh/pgs-tmp/chr1-chr4-0',
    '/net/scratch/people/plgjacekh/pgs-tmp/chr5-chr9-0',
    '/net/scratch/people/plgjacekh/pgs-tmp/chr10-chr14-0',
    wd + 'data/pharmacogenetic-score/intermediate/chr15-chr19a-0',
    wd + 'data/pharmacogenetic-score/intermediate/chr19-to-end-0',
]
for vcf_annotated_path in paths:
    vcf_scores = hl.read_matrix_table(vcf_annotated_path)
    compute_burden = lambda mt, score: hl.agg.sum(score * mt.GT.n_alt_alleles())
    vcf_burden = vcf_scores.group_rows_by('gene_name').aggregate(
        n=hl.agg.sum(vcf_scores.GT.n_alt_alleles()),
        pgs=compute_burden(vcf_scores, vcf_scores.pgs),
        score_cadd=compute_burden(vcf_scores, vcf_scores.cadd_raw_oryg),
        score_fathmm_xf=compute_burden(vcf_scores, vcf_scores.fathmm_xf_oryg),
        score_provean=compute_burden(vcf_scores, vcf_scores.provean_oryg),
        score_mutation_assessor=compute_burden(
            vcf_scores, vcf_scores.mutation_assessor_oryg
        ),
    )
    vcf_burden.write(vcf_annotated_path[:-1] + '1')  # chr1-chr4 - 2 nodes - 54mins

# Merge PharmGScore
burdens_a = [
    hl.read_matrix_table(file)
    for file in [
        '/net/scratch/people/plgjacekh/pgs-tmp/chr1-chr4-1',
        '/net/scratch/people/plgjacekh/pgs-tmp/chr5-chr9-1',
        '/net/scratch/people/plgjacekh/pgs-tmp/chr10-chr14-1',
        wd + 'data/pharmacogenetic-score/intermediate/chr15-chr19a-1',
    ]
]
burden_a = hl.MatrixTable.union_rows(*burdens_a)

burden_b = hl.read_matrix_table(
    wd + 'data/pharmacogenetic-score/intermediate/chr19-to-end-1'
)

a_cols = burden_a.s.collect()
b_cols = burden_b.s.collect()
col_names = set(a_cols) & set(b_cols)
burden_a = burden_a.filter_cols(hl.literal(col_names).contains(burden_a.s))
burden_b = burden_b.filter_cols(hl.literal(col_names).contains(burden_b.s))
burden_b = reorder_cols(burden_b, burden_a)
burden = burden_a.union_rows(burden_b)
burden = burden.repartition(10000)

burden.write(wd + 'data/pharmacogenetic-score/ukb-pgs.mt')

annot = load_annotations(wd + 'data/ukb-annotations.ht')
coding19 = hl.import_table(wd + 'raw/dataset/data-codings/coding19.tsv', key='coding')
icd_10_fields = {'main': 'f_41202', 'secondary': 'f_41204', 'all': 'f_41270'}
a = annot.select('f_41202', 'f_41204', 'f_41270')

pgs = PharmGScore.read_table()
vcf = hl.read_matrix_table(wd + 'data/pvcf/vcf-adr-ctrl-better.mt')

test_vcf = hl.import_vcf(
    wd + 'data/pharmacogenetic-score/test-files/test.vcf',
    array_elements_required=False
)

# annotate all patients
PharmGScore.compute_burden(vcf, next(tmpdir_path)).write(
    wd + 'data/pharmacogenetic-score/vcf-adr-ctrl-better-pgs.mt'
)

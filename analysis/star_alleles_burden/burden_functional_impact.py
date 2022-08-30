"""This module examine utilisation of burden to assess the functional impact
of pharmacogenetic variants (task #0)"""
import os

import pandas as pd
import hail as hl
from hail.plot import show
from bokeh.layouts import gridplot

from analysis.star_alleles_burden.pharmvar import PharmVar
from analysis.utils.load_spark import hl_init
from analysis.utils.load_spark import wd, localfs_path, scratch_path
from analysis.utils.load_spark import tmpdir_path_iter


hl_init()
hl.plot.output_notebook()
tmpdir_path = tmpdir_path_iter()


vcf_rows = PharmVar.read_table()
cadd = hl.read_table(os.path.join(
    '/net/archive/groups/plggneuromol/imdik-zekanowski-gts',
    'data/external-data/cadd-full.ht'
))

# ex. time: 6m
vcf_rows = vcf_rows.annotate(
    cadd_score=cadd[vcf_rows.key].cadd_score
)
vcf_rows = vcf_rows.checkpoint(next(tmpdir_path))

vcf_rows.to_pandas().to_csv('results/star-cadd-burden.csv', index=False)

pg_cadd_burden = (
    vcf_rows
    .group_by(vcf_rows.gene, vcf_rows.function)
    .aggregate(cadd_sum=hl.int(hl.agg.mean(vcf_rows.cadd_score)))
)
pg_cadd_burden = pg_cadd_burden.checkpoint(next(tmpdir_path))

pg_cadd_burden = pg_cadd_burden.filter(
    ~hl.set(["function not assigned"]).contains(pg_cadd_burden.function)
)
genes = set(pg_cadd_burden.gene.collect())
for g in genes:
    pg_cadd_burden.filter(pg_cadd_burden.gene == g).show()

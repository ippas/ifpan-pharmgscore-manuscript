"""Convert pVCF into hail matrix table with non-ref only, then annotate with
scores. Filter for pgx genes only."""
# %%
import os
import re
from functools import reduce

import hail as hl

from analysis.utils.pathogenicity_scores import annotate_vcf_with_scores


def vcf_to_non_ref_mt(vcf_path, out_path, filter_pgx=False):
    mt = hl.import_vcf(
        vcf_path,
        force_bgz=True,
        array_elements_required=False,
    )

    if filter_pgx:
        bed_file_path = 'analysis/polygenic/pgx.bed'
        bed_ht = hl.import_bed(bed_file_path)
        max_ranges = hl.literal(
            {
                'chr10': 133797422,
                'chrM': 16569,
            }
        )
        intervals = bed_ht.aggregate(
            hl.agg.collect(
                hl.locus_interval(
                    bed_ht.interval.start.contig,
                    hl.max(bed_ht.interval.start.position - 1_000, 1),
                    hl.min(
                        bed_ht.interval.end.position + 1_000,
                        max_ranges.get(bed_ht.interval.start.contig)
                    ),
                )
            )
        )
        mt = mt.filter_rows(
            hl.any(lambda interval: interval.contains(mt.locus), intervals)
        )

    mt = hl.split_multi_hts(mt, permit_shuffle=True)

    mt = mt.filter_entries(
        mt.GT.is_non_ref()
    )
    mt.write(out_path)


# %%
chrom = 
chrom = {'23': 'X', '24': 'Y'}.get(chrom, chrom)
print(f'Chromosome: {chrom}')

vcf_to_non_ref_mt(
    vcf_path=vcf_path,  # * expansion for chr
    filter_pgx=True,
    out_path=f'{out_path}/vcf-mt/chr{chrom}.mt'
)

mt = hl.read_matrix_table(f'{out_path}/vcf-mt/chr{chrom}.mt')
mt_rows_annot = annotate_vcf_with_scores(mt.rows())
mt_annot = mt.annotate_rows(
    **mt_rows_annot[mt.row_key]
)
mt_annot.write(f'{out_path}/vcf-mt/chr{chrom}-annot.mt')


# %% all chrs
print('Writing...', flush=True)
vcf_mt_path = out_path + 'vcf-mt/'
mt_names = [
    f for f in os.listdir(vcf_mt_path) if re.match(r'chr\d+(_\d)?-annot\.mt', f)
]
mts = [hl.read_matrix_table(vcf_mt_path + b) for b in mt_names]
mt = hl.MatrixTable.union_rows(*mts)
mt.write(f'{out_path}/ukb200-pgx.mt')

# %%
import os
from itertools import islice
from pathlib import Path
import json
from collections import defaultdict

import hail as hl
import pandas as pd
from tqdm import tqdm
import polars as pl


def match(a, b):
    b_dict = {x: i for i, x in enumerate(b)}
    return [b_dict.get(x, None) for x in a]


def batched(iterable, n):
    """Yield successive n-sized chunks from iterable."""
    it = iter(iterable)
    while chunk := list(islice(it, n)):
        yield chunk


# %%
hl_init()
hl.plot.output_notebook()
tmpdir_path = tmpdir_path_iter()

pgx_mt_path = 'data/polygenic/pgx-vcfs.mt'
pgx_filtered_split_mt_path = 'data/polygenic/pgx-vcfs-filtered-split.mt'
n_batch = 4000


# %% filter relevant variants +/- 100kb and split multiallelic sites
pgx_mt_path  # imported raw vcf 
pgx_mt = hl.read_matrix_table(pgx_mt_path)
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

filtered_mt = pgx_mt.filter_rows(
    hl.any(lambda interval: interval.contains(pgx_mt.locus), intervals)
)
filtered_split_mt = hl.split_multi_hts(filtered_mt, permit_shuffle=True)
filtered_split_mt.write(pgx_filtered_split_mt_path)


# %% write down all samples in VCFs
pgx_mt = hl.read_matrix_table(pgx_filtered_split_mt_path)
pats = pgx_mt.s.collect()
with open('analysis/polygenic/all-eids.txt', 'w') as f:
    f.write('\n'.join(pats))  # TODO: add a new line at the end! (wc -l)


# %% split matrix table into smaller one (faster export to VCFs with one sample)
pgx_mt = hl.read_matrix_table(pgx_filtered_split_mt_path)
with open('analysis/polygenic/all-eids.txt') as f:
    pats = f.read().strip().split('\n')

for part, sub_eids in enumerate(batched(pats, 67_000)):
    part_path = f'/tmp/pgx-sub-eids-{part}.mt'
    sub_pgx_mt = pgx_mt.filter_cols(hl.literal(sub_eids).contains(pgx_mt.s))
    sub_pgx_mt = sub_pgx_mt.repartition(24)
    sub_pgx_mt.write(part_path)


# %% export each patient to a separate VCF file
vcf_metadata = hl.get_vcf_metadata(example_vcf_path)

with open('analysis/polygenic/all-eids.txt') as f:
    pats = f.read().strip().split('\n')

cached = False
for part, sub_eids in enumerate(batched(pats, n_batch)):
    print('part:', part)
    pgx_vcf_pats_dir = f'tmp/pgx-patients/part-{part}'
    os.makedirs(pgx_vcf_pats_dir, exist_ok=True)
    sub_mt = hl.read_matrix_table(f'tmp/pgx-sub-eids-{part}.mt')

    pgx_files = set(Path(pgx_vcf_pats_dir).iterdir())
    for i, sample_id in enumerate(sub_eids):
        sample_path = os.path.join(pgx_vcf_pats_dir, f'{sample_id}.vcf.bgz')
        if Path(sample_path) not in pgx_files:
            if not cached:
                sub_mt = sub_mt.cache()
                cached = True
            sample_mt = sub_mt.filter_cols(sub_mt.s == sample_id)
            hl.export_vcf(
                sample_mt,
                sample_path,
                metadata=vcf_metadata,
                tabix=True
            )
    else:
        cached = False



# %% run polygenic
polygenic.sh

# %% prepare for PharmCAT outside call format
not_supported_by_PharmCAT = {'CYP2C8', 'CYP2A6', 'CYP2A13'}
genes = {'CYP2A13', 'CYP2A6', 'CYP2B6', 'CYP2C19', 'CYP2C8', 'CYP2C9', 'CYP2D6',
         'CYP3A4', 'CYP3A5', 'CYP4F2', 'DPYD', 'NUDT15', 'SLCO1B1'} - not_supported_by_PharmCAT
out_dir = Path('data/polygenic/cromwell-output/')


eid_phenotypes_dict = defaultdict(list)
for part_path in out_dir.iterdir():
    print(part_path.name)

    phenotype_list = list(part_path.glob('**/*.phenotype.json'))
    assert len(phenotype_list) == 4000 or part_path.name == 'part-50'
    assert len(phenotype_list) == 642 or part_path.name != 'part-50'
    for ph_path in tqdm(phenotype_list):
        pat = ph_path.name.rsplit('-', maxsplit=1)[0]

        # sample_id-openpgx.tsv is not empty
        with open(ph_path.with_suffix('').with_suffix('').with_suffix('.tsv')) as f:
            if not f.read():
                print(pat, ph_path)
                break  # tsv won't be complete

        with open(ph_path) as f:
            pat_phenotype = json.load(f)

        cpic = pat_phenotype['geneReports']['CPIC']
        dpwg = pat_phenotype['geneReports']['DPWG']

        for gene in genes:
            if cpic[gene]['sourceDiplotypes'] != cpic[gene]['recommendationDiplotypes']:
                if gene == 'DPYD':
                    for a, b in zip(cpic[gene]['sourceDiplotypes'], cpic[gene]['recommendationDiplotypes']):
                        aa = a.copy()
                        aa['allele2'] = a['allele1']
                        aa['allele1'] = a['allele2']
                        if (aa != b) & (a != b):
                            print('not equal cpic', pat, gene)
                else:
                    print('not equal cpic', pat, gene)
            if dpwg[gene]['sourceDiplotypes'] != dpwg[gene]['recommendationDiplotypes']:
                if gene == 'DPYD':
                    for a, b in zip(dpwg[gene]['sourceDiplotypes'], dpwg[gene]['recommendationDiplotypes']):
                        aa = a.copy()
                        aa['allele2'] = a['allele1']
                        aa['allele1'] = a['allele2']
                        if (aa != b) & (a != b):
                            print('not equal dpwg', pat, gene)
                else:
                    print('not equal dpwg', pat, gene)
            if len(cpic[gene]['sourceDiplotypes']) != 1:
                print('lenght cpic != 1', pat, gene)
            if len(dpwg[gene]['sourceDiplotypes']) != 1:
                print('lenght dpwg != 1', pat, gene)

            eid_phenotypes_dict['eid'].append(pat)
            item = cpic[gene]['sourceDiplotypes'][0]
            for key, val in item.items():
                if isinstance(val, dict):
                    for sub_key, sub_val in val.items():
                        eid_phenotypes_dict[f'{key}_{sub_key}'].append(sub_val)
                elif isinstance(val, list):
                    if len(val) > 1:
                        print(ph_path, gene, val)
                    else:
                        eid_phenotypes_dict[key].append(val[0] if val else '')
                else:
                    eid_phenotypes_dict[key].append(val)

    cols = [
        'eid', 'gene', 'label', 'phenotypes', 'activityScore', 'lookupKey',
        'allele1_function', 'allele1_activityValue',
        'allele2_function', 'allele2_activityValue',
       'phenotypeDataSource'
    ]
    df = pd.DataFrame(eid_phenotypes_dict)
    df.to_csv(f'data/polygenic/pheno-{part_path.name}.tsv', sep='\t', index=False)

# %%
polygenic_results_list = list()
polygenic_results_paths = list(Path('data/polygenic/').glob('pheno-part*.tsv'))
for polygenic_result_path in polygenic_results_paths:
    df = pl.read_csv(
        polygenic_result_path,
        separator='\t',
        null_values='n/a',
        schema_overrides={'eid': pl.String}
    )
    polygenic_results_list.append(df)

polygenic_results = pl.concat(polygenic_results_list)

(
    polygenic_results
    .rename({'phenotypes': 'poly_phenotype', 'label': 'poly_label'})
    .write_csv('data/polygenic/polygenic-results.csv')
)
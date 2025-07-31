"""UKB phenotype annotation"""
# %%
from typing import Counter

import pandas as pd
import hail as hl
from hail.plot import show
from bokeh.io import show, output_notebook
from bokeh.layouts import gridplot
from IPython.display import display
output_notebook()

from analysis.utils.annotations import load_annotations


def codes_to_set(codes_array, first_n=4):
    codes_set = hl.set(hl.map(lambda x: x[:first_n], codes_array))
    codes_set = codes_set.remove(hl.missing(hl.tstr))
    return codes_set


# %%
CYPS = [
    'CYP1A2',
    'CYP2B6',
    'CYP2C19',
    'CYP2C8',
    'CYP2C9',
    'CYP2D6',
    'CYP2E1',
    'CYP3A4',
    'CYP3A5',
]

list_genes = (
    hl.import_table(
        'data/genes-list.csv',
        key='Gene',
        delimiter=','
    )
    .rename({'Gene': 'gene_name', 'Group': 'group'})
    .key_by('gene_name')
)
n_codes = 4


# %% controls phenotyping-adr-groups
fafb = hl.read_table('data/full-annots-for-burden.ht')
fafb2 = fafb.select(
    fafb.f_31,  # sex
    fafb.f_189,  # townsend
    fafb.f_21022,  # age at the recruitment
    fafb.group,
    fafb.intention,
)
fafb2 = fafb2.transmute(
    sex=hl.if_else(fafb2.f_31 == 'Male', 1, 0),
    townsend=hl.float(fafb2.f_189),
    age=hl.int(fafb2.f_21022),
)

# from mdd directory
dd = pd.read_csv('misc/therapies/therapies_stats_exploded.tsv', sep='\t', dtype={'eid': 'str'})
dd = dd.rename({'eid': 's'}, axis=1)
dd = dd[[
    's', 'therapies.average_dose', 'therapies.drug_name',
    'therapies.total_duration', 'therapies.start_year',
    'therapies.start_month', 'therapies.start_day',
    'changes_in_dose', 'changes_in_drug', 'no_changes',
    'number_of_therapies', 'total_duration', 'drug_names', 'unique_drugs'
]]
dd_hl = hl.Table.from_pandas(dd).key_by('s')
dd_hl = dd_hl.persist()

dd = dd_hl.transmute(
    therapies=hl.Struct(
        average_dose=hl.int(dd_hl['therapies.average_dose']),
        drug_name=dd_hl['therapies.drug_name'],
        total_duration=dd_hl['therapies.total_duration'],
        start_year=hl.int(dd_hl['therapies.start_year']),
        start_month=hl.int(dd_hl['therapies.start_month']),
        start_day=hl.int(dd_hl['therapies.start_day']),
    )
)
dd = dd.group_by('s').aggregate(
    therapies=hl.agg.collect_as_set(dd.therapies),
)
dd_hl_dist = dd_hl.distinct().select(
    'changes_in_dose', 'changes_in_drug', 'no_changes',
    'number_of_therapies', 'total_duration', 'drug_names', 'unique_drugs'
)
dd_hl_dist = dd_hl_dist.annotate(
    **dd[dd_hl_dist.key]
)
dd_hl_dist = dd_hl_dist.annotate(
    drug_names=hl.set(
        (
            dd_hl_dist.drug_names
            .replace('\{', '').replace('\}', '').replace("'", '').replace(" ", '')
            .split(',')
            .map(lambda x: hl.str(x))
        )
    )
)

therapy_ordering = lambda mt: ((
    -mt.therapies.start_year, -mt.therapies.start_month, -mt.therapies.start_day,
    -mt.therapies.average_dose
))
bc = dd_hl_dist
bc = bc.explode(bc.therapies)
bc = bc.persist()
bcg = bc.group_by('s')
bc2 = bcg.aggregate(
    start_last_therapy=hl.array([
        hl.agg.take(bc.therapies.start_year, 1, ordering=therapy_ordering(bc))[0],
        hl.agg.take(bc.therapies.start_month, 1, ordering=therapy_ordering(bc))[0],
        hl.agg.take(bc.therapies.start_day, 1, ordering=therapy_ordering(bc))[0],
    ]),
)

bc3 = bc.annotate(
    **bc2[bc.key]
)
bc3 = bc3.filter(
    (bc3.therapies.start_year == bc3.start_last_therapy[0])
    & (bc3.therapies.start_month == bc3.start_last_therapy[1])
    & (bc3.therapies.start_day == bc3.start_last_therapy[2])
)
bc3g = bc3.group_by('s')
bc3 = bc3g.aggregate(
    last_drugs=hl.agg.collect_as_set(bc3.therapies.drug_name),
    last_durations=hl.agg.collect_as_set(bc3.therapies.total_duration),
    last_average_doses=hl.agg.collect_as_set(bc3.therapies.average_dose),
)

dd2 = dd_hl_dist.annotate(
    **bc2[dd_hl_dist.key],
    **bc3[dd_hl_dist.key]
)


burden = hl.read_matrix_table('tmp/pgx/burden.mt')
burden = burden.annotate_cols(
    **fafb2[burden.col_key],
    **dd2[burden.col_key],
)
burden.write('tmp/pgx/burden-annot.mt', overwrite=True)

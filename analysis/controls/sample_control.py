"""Picking control group"""
import os
from operator import itemgetter

import hail as hl
from hail.plot import show
from IPython.display import display
import pandas as pd

from analysis.utils.load_spark import hl_init
from analysis.utils.load_spark import wd, localfs_path, scratch_path
from analysis.utils.annotations import (
    find_annot_desc, load_annotations, annotate_adr_patients,
    flatten_map_agg, dc_recoding, #icd10_match, annotate_genes
)

def convert_annot_types(annots, f_map):
    """Convert columns of a table to be easier to work with: flattening,
    type conversion, aggregation across time points.

    annots: annotations table
    f_map: dictionary: {'friendly_name': 'actual_colname'}
    """
    depressive_episodes_code = 'F32'
    recurrent_depressive_disorder_code = 'F33'
    d_trans = {
        #f_map['sex']: hl.if_else(annots[f_map['sex']] == "Male", 1, 0),
        f_map['townsend']: hl.float32(annots[f_map['townsend']]),
        f_map['bmi']: flatten_map_agg(
            annots[f_map['bmi']],
            map_fun=hl.float,
            agg_fun=hl.mean
        ),
        f_map['ever_depressed_week']: flatten_map_agg(
            annots[f_map['ever_depressed_week']],
            map_fun=lambda x: dc_recoding(x, update_dict={"_missing": 0}),
            agg_fun=hl.max
        ),
        f_map['seen_psychiatrist']: flatten_map_agg(
            annots[f_map['seen_psychiatrist']],
            map_fun=lambda x: dc_recoding(x, update_dict={"_missing": 0}),
            agg_fun=hl.max
        ),
        f_map['prolonged_feelings_of_sadness']: dc_recoding(
            annots[f_map['prolonged_feelings_of_sadness']],
            update_dict={"_missing": 0}
        ),
    }
    annots = annots.transmute(**d_trans)

    #d_ann = {
    #    f_map['depressive_episodes']: icd10_match(
    #        annots[f_map['icd_10_all']], depressive_episodes_code
    #    ),
    #    f_map['recurrent_depressive_disorder']: icd10_match(
    #        annots[f_map['icd_10_all']], recurrent_depressive_disorder_code
    #    )
    #}
    annots = annots.annotate(**d_ann)
    return annots



#----------------------------------------------------------------------

def fit_control(annots, f_map, method='exact'):
    """Create all combinations of groups based on values in f_map dictionary.
    In each group take randomly as many controls as adr patients.

    annots: annotations table
    f_map: dictrionary with column names to be fit
    """
    # make bins for continuos variables
    # convertin to str, because of bug in pandas
    qcut_col = pd.qcut(annots[f_map['townsend']], 5).astype(str)
    qcut_col = qcut_col.where(qcut_col != 'nan', None)
    annots[f_map['townsend']] = qcut_col

    qcut_col = pd.qcut(annots[f_map['bmi']], 5).astype(str)
    qcut_col = qcut_col.where(qcut_col != 'nan', None)
    annots[f_map['bmi']] = qcut_col

    # creat combinations of all groups
    groups = [
        (intervals, sub_df)
        for intervals, sub_df in annots.groupby(list(f_map.values()), dropna=False)
    ]

    control = []
    for group, sub_group in groups:
        adr = sub_group.groupby('adr').adr.count().get(True, 0)
        non_adr = sub_group.groupby('adr').adr.count().get(False, 0)
        if adr > non_adr:
            print(f"WARN: insufficient control (adr/non-adr): {adr}/{non_adr}")
            print(group)
        idx = sub_group.loc[sub_group.adr == False].sample(min(adr, non_adr))
        control.append(idx)
    df_control = pd.concat(control)
    return df_control.index


if __name__ == '__main__':
    hl_init()


    annot_raw = load_annotations(wd + 'data/ukb-annotations.ht')
    coding19 = hl.import_table(
        wd + 'raw/dataset/data-codings/coding19.tsv', key='coding'
    )

    an = annotate_adr_patients(annot_raw)

    depressive_episodes_code = 'F32'
    recurrent_depressive_disorder_code = 'F33'
    f_tofit = {
        'sex': 'f_31',
        'townsend': 'f_189',
        'bmi': 'f_21001',
        'ever_depressed_week': 'f_4598',
        'prolonged_feelings_of_sadness': 'f_20446', #!
        'seen_psychiatrist': 'f_2100',
        'icd_10_all': 'f_41270',
        'depressive_episodes': 'f_41270_' + depressive_episodes_code,
        'recurrent_depressive_disorder': 'f_41270_' + recurrent_depressive_disorder_code,
    }

    anc = convert_annot_types(an, f_map=f_tofit)

    df = anc.select(*f_tofit.values(), 'adr').to_pandas()
    df = df.set_index('f_eid')

    f_ = f_tofit.copy()
    del f_['icd_10_all']
    del f_['ever_depressed_week']
    del f_['seen_psychiatrist']
    control_index = fit_control(df, f_map=f_)

    # final set
    an_final = an.filter(
        an.adr | hl.literal(set(control_index)).contains(an.f_eid)
    )
    df_final = an_final.to_pandas()

    df_final.to_csv(
        wd + 'data/adr-control-patients.csv',
        columns=['f_eid', 'adr'],
        index=False
    )

    df_final = pd.read_csv(
        wd + 'data/adr-control-patients.csv',
        dtype={'f_eid': str}
    )

    # first set: chr1-chr19b31; second: chr19b32-chrY
    print("Creating VCF set0...")
    mt = hl.import_vcf(
        [
            wd + 'raw/pvcf-200k/ukb23156_c[1-9]_b*_v1.vcf.gz',
            wd + 'raw/pvcf-200k/ukb23156_c1[0-8]_b*_v1.vcf.gz',
            wd + 'raw/pvcf-200k/ukb23156_c19_b?_v1.vcf.gz',
            wd + 'raw/pvcf-200k/ukb23156_c19_b[1-2]?_v1.vcf.gz',
            wd + 'raw/pvcf-200k/ukb23156_c19_b3[0-1]_v1.vcf.gz',
        ],
        force_bgz=True,
        array_elements_required=False,
    )
    mt = mt.filter_cols(hl.literal(set(df_final.f_eid)).contains(mt.s))
    mt.write(wd + 'data/pvcf/c_vcf_adr_ctl_set0.mt')

    print("Creating VCF set1...")
    mt = hl.import_vcf(
        [
            wd + 'raw/pvcf-200k/ukb23156_c19_b3[2-9]_v1.vcf.gz',
            wd + 'raw/pvcf-200k/ukb23156_c19_b[4-9]?_v1.vcf.gz',
            wd + 'raw/pvcf-200k/ukb23156_c2[0-2]_b*_v1.vcf.gz',
            wd + 'raw/pvcf-200k/ukb23156_c{X,Y}_b*_v1.vcf.gz',
        ],
        force_bgz=True,
        array_elements_required=False,
    )
    mt = mt.filter_cols(hl.literal(set(df_final.f_eid)).contains(mt.s))
    mt.write(wd + 'data/pvcf/c_vcf_adr_ctl_set1.mt')


    vcf0 = hl.read_matrix_table(wd + 'data/pvcf/c_vcf_adr_ctl_set0.mt')
    vcf0_g = annotate_genes(vcf0, wd + 'data/pvcf/c_vcf_adr_ctl_set0_genes.mt')
    vcf1 = hl.read_matrix_table(wd + 'data/pvcf/c_vcf_adr_ctl_set1.mt')
    vcf1_g = annotate_genes(vcf1, wd + 'data/pvcf/c_vcf_adr_ctl_set1_genes.mt')

    mta = hl.MatrixTable.union_rows(vcf0_g, vcf1_g)
    mta.write(wd + 'data/pvcf/vcf-adr-ctrl-better.mt')

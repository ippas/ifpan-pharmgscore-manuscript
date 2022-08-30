"""this module creates a summary table of important statistics"""

import hail as hl
import pandas as pd

def counter_stat(ht, expr, group):
    """get counter of expr in group"""
    results = ht.aggregate(
        hl.agg.group_by(
            (group),
            hl.agg.counter(ht[expr])
        )
    )[True]
    return(results)

def field_stat(ht, expr, group):
    """get statistics of expr in group"""
    results = ht.aggregate(
        hl.agg.group_by(
            (group),
            hl.agg.stats(ht[expr])
        )
    )[True]
    return(results)

def make_summary_table(ht):
    """describe each of the analysed groups in terms of important statistics"""

    col_names = ['n',
                 'percentage of females',
                 'Townsend mean',
                 'Townsend SD',
                 'BMI mean',
                 'BMI SD',
                 'age when at assesment centre mean',
                 'age when at assesment centre SD',
                 'depression hospital code',
                 'anxiety hospital code',
                 'both depression and anxiety hospital code',
                 'self harm with drugs',
                 'self harm without drugs',
                 'drug abuse history',
                 'number of non-cancer illness reported mean',
                 'number of non-cancer illness reported sd',
                 'ever contemplated self harm',
                 'taken any medication for depression']
    
    table_dict = dict.fromkeys(col_names, None)
    
    f_map = {
        'sex' : 'f_31',
        'townsend' : 'f_189',
        'bmi' : 'f_21001',
        'age' : 'f_21003',
        'illness_no' : 'f_135',
        'contemplated' : 'f_20485',
        'medicated' : 'f_20546'
    }
    

    anx_or_md = ['anxiety', 'depression', 'both']
    self_harm_list = ['Yes, more than once', 'Yes, once']
    medicated_list = ['Drugs or alcohol (more than once)', 'Medication prescribed to you (for at least two weeks)', 'Unprescribed medication (more than once)']

    # add a relevant annotation field and convert others
    ht = ht.annotate(
        inpatient_anx_or_md = hl.if_else(
        (hl.literal(
                anx_or_md).contains(ht.mental_health_inpatient)
            ),
             'anx_or_md',
             'no_anx_or_md'),
         f_189 = hl.float32(ht[f_map['townsend']]),
         f_21001 = hl.float32(ht[f_map['bmi']][0][0]),
         f_21003 = hl.float32(ht[f_map['age']][0][0]),
         f_135 = hl.float32(ht[f_map['illness_no']][0][0]),
         contemplated_self_harm = hl.literal(
                     self_harm_list
                     ).contains(ht[f_map['contemplated']]),
         medicated_depression = hl.literal(
                         medicated_list
                     ).contains(
                         ht[f_map['medicated']][0][0])
              )

    groups = {
            'adr':(ht['group'] == 'adr'),
            'adr_intentional':(ht['intention'] == 'intentional'),
            'adr_unintentional':(
                (ht['intention'] == 'therapeutic') | (ht['intention'] == 'accidental')
             ),
            'all_controls':(ht['group'] == 'control'),
            'self_harm_drugs_controls':((ht['group'] == 'control') & (ht['self_harm'] == 'with_drugs')),
            'anx_md_controls':((ht['group'] == 'control') & (ht['inpatient_anx_or_md'] == 'anx_or_md'))
        }

    n_s=[]
    m_f_ratios = []
    tow_means = []
    tow_sds = []
    bmi_means = []
    bmi_sds = []
    age_means = []
    age_sds = []
    deps = []
    anxs = []
    depsanxs = []
    drugs = []
    nodrugs = []
    drug_hist = []
    ill_means = []
    ill_sds = []
    contemplated = []
    medicated = []
    
    for group in groups.values():
        
        n_s.append(
            ht.filter(group).count()
        )
    
        sex_dict = counter_stat(
            ht, f_map['sex'], group
        ) 
        
        m_f_ratios.append(
            (sex_dict['Female']*100)/(sex_dict['Female']+sex_dict['Male'])
        )
        
        tow_res = field_stat(
            ht, f_map['townsend'], group
        )
        
        tow_means.append(tow_res['mean'])
        tow_sds.append(tow_res['stdev'])
        
        bmi_res = field_stat(
            ht, f_map['bmi'], group
        )
        
        bmi_means.append(bmi_res['mean'])
        bmi_sds.append(bmi_res['stdev'])
        
        age_res = field_stat(
            ht, f_map['age'], group
        )
        
        age_means.append(age_res['mean'])
        age_sds.append(age_res['stdev'])
        
        depression_dict = counter_stat(
            ht, 'mental_health_inpatient', group
        ) 
        
        deps.append(depression_dict['depression'])
        anxs.append(depression_dict['anxiety'])
        depsanxs.append(depression_dict['both'])
        
        self_harm_dict = counter_stat(
            ht, 'self_harm', group
        )
        
        drugs.append(self_harm_dict.get('with_drugs', 0))       
        nodrugs.append(self_harm_dict.get('without_drugs', 0))
        
        drug_abuse_dict = counter_stat(
            ht, 'drug_abuse', group
        )
        
        drug_hist.append(drug_abuse_dict['yes'])
        
        ill_res = field_stat(
            ht, f_map['illness_no'], group
        )
        
        ill_means.append(ill_res['mean'])
        ill_sds.append(ill_res['stdev'])
        
        contemplated_count = counter_stat(
            ht, 'contemplated_self_harm', group
        ).get(True, 0)
        
        contemplated.append(contemplated_count)
        
        medicated_count = counter_stat(
            ht, 'medicated_depression', group
        ).get(True, 0)
        
        medicated.append(medicated_count)
  
    table_dict['n'] = n_s
    table_dict['percentage of females'] = m_f_ratios
    table_dict['Townsend mean'] = tow_means
    table_dict['Townsend SD'] = tow_sds
    table_dict['BMI mean'] = bmi_means
    table_dict['BMI SD'] = bmi_sds
    table_dict['age when at assesment centre mean'] = age_means
    table_dict['age when at assesment centre SD'] = age_sds
    table_dict['depression hospital code'] = deps
    table_dict['anxiety hospital code'] = anxs
    table_dict['both depression and anxiety hospital code'] = depsanxs
    table_dict['self harm with drugs'] = drugs
    table_dict['self harm without drugs'] = nodrugs
    table_dict['drug abuse history'] = drug_hist
    table_dict['number of non-cancer illness reported mean'] = ill_means
    table_dict['number of non-cancer illness reported sd'] = ill_sds
    table_dict['ever contemplated self harm'] = contemplated
    table_dict['taken any medication for depression'] = medicated

    df = pd.DataFrame.from_dict(table_dict)
    return df
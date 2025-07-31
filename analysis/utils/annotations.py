import os

import hail as hl


def annotate_adr_patients(annots, overwrite=False):
    annots = load_annotations('data/ukb-annotations.ht')
    """Annotate table with ADR patients based on ICD10 codes"""
    path = 'annots_200k.ht'
    if not os.path.exists(path) or overwrite:
        # filter 200k
        f_pvcf = 'f_23156'
        annots_200k = annots.filter(annots[f_pvcf] == "977")

        # filter out withdraws
        with open('raw/withdraw.csv') as f:
            withdraw = f.read().splitlines()
        annots_200k = annots_200k.filter(
            ~hl.set(withdraw).contains(annots_200k.f_eid)
        )

        # annotate ADR patients
        icd = {
                'main': 'f_41202',
                'secondary': 'f_41204',
                'all': 'f_41270',
        }
        
        _a = annots_200k.select(icd['all'], 'f_40001', 'f_40002').persist()
        _adr = (hl.set(_a[icd['all']][0]).intersection(
            hl.literal({'T430', 'T431', 'T432'})).length() > 0)
        _combined_death = flatten_map_agg(
            _a.f_40001
            ).extend(
                flatten_map_agg(
                    _a.f_40002
                )
        )
        _adr_death_cond = (hl.set(_combined_death)
                    .intersection(hl.literal({'T430', 'T431', 'T432'}))
                    .length() > 0
                    )
        sum((_adr | _adr_death_cond).collect())

        adr = (hl.set(annots_200k[icd['all']][0]).intersection(
            hl.literal({'Y490', 'Y491', 'Y492', 'T430', 'T431', 'T432'})).length() > 0)
        
        combined_death = flatten_map_agg(
            annots_200k.f_40001
            ).extend(
                flatten_map_agg(
                    annots_200k.f_40002
                )
        )

        adr_death_cond = (hl.set(combined_death)
                    .intersection(hl.literal({'Y490', 'Y491', 'Y492', 'T430', 'T431', 'T432'}))
                    .length() > 0
                    )
                

        annots_200k = annots_200k.annotate(
            group = hl.if_else(
                (adr | adr_death_cond),
                'adr',
                'control'
            )
        )

        annots_200k = annots_200k.persist()

    return hl.read_table(path)


def flatten_map_agg(expr, map_fun=None, agg_fun=None):
    """Useful for fields with many instances and indexes: flattens the field,
    converts data types and aggregates elements."""
    if map_fun is None:
        map_fun = lambda x: x
    if agg_fun is None:
        agg_fun = lambda x: x
    flatten_map = hl.flatten(expr).map(lambda x: map_fun(x))
    return agg_fun(flatten_map)


def annotate_date_codes_dict(ht):
    """Annotate table with a dictionary containg a date 
    with a hospital episode and all codes from that date"""
    ht = ht.annotate(
        f_41280 = ht['f_41280'][0],
        f_41270 = ht['f_41270'][0]
    )

    ht = ht.annotate(
        date_code_zip=(
            hl.zip(ht.f_41280, ht.f_41270)
            .filter(lambda x: hl.is_defined(x[0]))
        )
    )

    ht = ht.annotate(
        date_code_list= ht.date_code_zip.map(
            lambda dc: tuple([
                dc[0], 
                (
                    ht
                    .date_code_zip
                    .filter(lambda x: dc[0] == x[0])
                    .map(lambda y: y[1])
                )
            ])
        )
    )

    ht = ht.annotate(
        date_code_dict=hl.dict(ht.date_code_list)
    )
    
    return ht


def annotate_intent(ht):
    """define relevant phenotypes to ADRs and annotate table"""
    icd = {
        'main': 'f_41202',
        'secondary': 'f_41204',
        'all': 'f_41270',
    }
    
    self_harm_with_drugs = hl.str(ht[icd['all']]).contains('X6')
    accident_with_drugs = hl.str(ht[icd['all']]).contains('X4')

    therapeutic_dose = (
        (hl.set(ht[icd['all']])
         .intersection(hl.literal({'Y490', 'Y491', 'Y492'}))
         .length() > 0) &
        (hl.set(ht[icd['all']])
         .intersection(hl.literal({'T430', 'T431', 'T432'}))
         .length() == 0)
    )

    toxicity = (hl.set(ht[icd['all']])
                .intersection(hl.literal({'T430', 'T431', 'T432'}))
                .length() > 0
    )
    
    any_self_harm =  (
        hl.str(ht[icd['all']]).contains('Z915') |
        hl.str(ht[icd['all']]).contains('X6') |
        hl.str(ht[icd['all']]).contains('X7') |
        hl.str(ht[icd['all']]).contains('X80') |
        hl.str(ht[icd['all']]).contains('X81') |
        hl.str(ht[icd['all']]).contains('X82') |
        hl.str(ht[icd['all']]).contains('X83') |
        hl.str(ht[icd['all']]).contains('X84')
   )
    
    depression_codes = hl.literal({
    'F329', 'F339', 'F322', 'F321', 'F323', 'F331',
    'F334', 'F251', 'F315', 'F338', 'F328', 'F313',
    'F330', 'F316', 'F332'
    })
    
    anxiety_codes = hl.literal({'F419', 'F410', 'F411', 'F418', 'F606', 'F408', 'F409', 'F413'})

    depression = (
        hl.set(ht[icd['all']])
        .intersection(depression_codes)
        .length() > 0
    )
    
    anxiety = (
        hl.set(ht[icd['all']])
        .intersection(anxiety_codes)
        .length() > 0
    )
    
    combined_death = flatten_map_agg(
        ht.f_40001
    ).extend(
        flatten_map_agg(
            ht.f_40002
        )
    )

    is_dead = hl.is_defined(
        combined_death[0]
    )

    adr_death = hl.if_else(
        (hl.set(combined_death)
                .intersection(hl.literal({'Y490', 'Y491', 'Y492', 'T430', 'T431', 'T432'}))
                .length() > 0), 'adr', 'other')
    
    ht = ht.annotate(
        intention = hl.if_else(
            (therapeutic_dose & ~self_harm_with_drugs),
            'therapeutic',
            hl.if_else(
                (toxicity & ~self_harm_with_drugs),
                'accidental',
                hl.if_else(
                    (toxicity & self_harm_with_drugs),
                    'intentional',
                    hl.if_else(
                        (therapeutic_dose & self_harm_with_drugs),
                        'intentional',
                        hl.if_else(
                            (ht.group == 'adr'),
                            'unspecified','control'
                        )
                    )
                )
            )
        ),

        self_harm = hl.if_else(
            any_self_harm,
            hl.if_else(
                self_harm_with_drugs,
                'with_drugs',
                'without_drugs'
            ),
            'no_self_harm'
        ),

        drug_abuse = hl.if_else(
            ht[icd['all']].contains('Z864'),
            'yes',
            'no'
        ),

        mental_health_inpatient = hl.if_else(
            ((ht[icd['all']].contains('F412')) | (anxiety & depression)),
            'both',
            hl.if_else(
                depression,
                'depression',
                hl.if_else(
                    anxiety,
                    'anxiety',
                    'no_mental_health_inpatient')
            )
        ),
        death = hl.if_else(
            is_dead,
            adr_death,
            'alive'
        ),
    )

    return(ht)

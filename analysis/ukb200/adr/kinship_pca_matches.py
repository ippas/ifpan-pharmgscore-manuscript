"""This module removes related individuals, keeping cases preferentiably and calculates controls closest on PCA"""

import hail as hl
import pandas as pd
import numpy as np
import scipy

def annotate_relatednes_dict(ht):
    """Annotate table with a dictionary containg relatedness pair and kinship"""
    ht = ht.annotate(
        f_22011 = ht['f_22011'][0],
        f_22012 = ht['f_22012'][0]
    )

    ht = ht.annotate(
        pairs_zip=(
            hl.zip(ht.f_22011, ht.f_22012)
            .filter(lambda x: hl.is_defined(x[0]))
        )
    )

    ht = ht.annotate(
        pairs_dict=hl.dict(ht.pairs_zip)
    )
    
    return ht

def tie_breaker(l, r):
     return hl.if_else(l.is_case & ~r.is_case, -1,
                       hl.if_else(~l.is_case & r.is_case, 1, 0))

def remove_related(ht, kin_cut=0.125):
    """remove individuals with kinship above kin_cut on scale 0-5 keeping cases preferentially with tie braker"""
    ht = annotate_relatednes_dict(ht)
    kin = ht.explode(ht.pairs_zip)

    kin = kin.select(kin.pairs_zip)

    kin = kin.transmute(
        pair = kin.pairs_zip[0],
        kinship = hl.float64(kin.pairs_zip[1])
    )

    kin = kin.filter(kin.kinship < kin_cut)

    pair_numbers = kin.pair.collect()
    seen = set()
    pairs = [x for x in pair_numbers if x in seen or seen.add(x)]  

    kin = kin.filter(hl.literal(pairs).contains(kin.pair))
    
    kin = (kin.group_by(kin.pair)
           .aggregate(pairs = hl.agg.collect(kin.f_eid)))
    
    kin = kin.transmute(f_eid_1 = kin.pairs[0],
                        f_eid_2 = kin.pairs[1])
    
    kin = kin.annotate(
        group_1 = ht[kin.f_eid_1].group,
        group_2 = ht[kin.f_eid_2].group
    )

    pairs_with_case = kin.key_by(
         i=hl.struct(id=kin.f_eid_1, is_case=(kin.group_1 == 'adr')),
         j=hl.struct(id=kin.f_eid_2, is_case=(kin.group_2 == 'adr'))
    )
   
    related_samples_to_remove = hl.maximal_independent_set(
        pairs_with_case.i, pairs_with_case.j, False, tie_breaker)

    ht = ht.filter(hl.is_defined(
         related_samples_to_remove.key_by(
            s = related_samples_to_remove.node.id)[ht.key]), keep=False)
    
    return ht

import os
import re
import io
from uuid import uuid4
from collections import defaultdict

import pandas as pd
import hail as hl
from biomart import BiomartServer

from analysis.utils.load_spark import wd, localfs_path


def load_annotations(path):
    """Read annotations for all ukb patients (converted by adr.Rmd
    from the main dataset).

    path: annotations path (as ht) to read (if exists) or to convert, merge
        and write to
    """
    if not os.path.exists(path):
        print('converting annotations...')
        annot = hl.import_table(wd + 'raw/dataset/ukb48815.tsv', key='f.eid')
        annot2 = hl.import_table(wd + 'raw/dataset/ukb49742.tsv', key='f.eid')
        # don't merge removed patients at a newer annotations (how='left')
        annot = annot2.join(annot, how='left')
        # annot = annot.repartition(24)
        annot = annot.checkpoint(localfs_path + 'a.ht', overwrite=True)

        field_name_re = re.compile(
            r'^f\.(?P<field_id>\d+)' +
            r'\.(?P<instance_index>\d+)\.' +
            r'(?P<array_index>\d+)$'
        )
        fields = defaultdict(lambda: defaultdict(list))

        for col in annot.row_value.keys():
            match = field_name_re.match(col)
            if match:
                array_index = int(match.group('array_index'))
                instance_index = int(match.group('instance_index'))
                fields[match.group('field_id')][instance_index].append(array_index)

        for i, (field, sub_fields) in enumerate(fields.items()):
            if len(sub_fields) > 1 or len(sub_fields[0]) > 1:
                field_kwargs = {}
                instance_kwargs = {}
                for instance, elements in sub_fields.items():
                    cols = [f'f.{field}.{instance}.{idx}' for idx in elements]
                    instance_kwargs[f'f.{field}.{instance}'] = hl.array(
                        [annot[c] for c in cols]
                    )
                annot = annot.transmute(**instance_kwargs)
                cols = [f'f.{field}.{inst}' for inst in sub_fields]
                field_kwargs[f'f_{field}'] = hl.array([annot[c] for c in cols])
                annot = annot.transmute(**field_kwargs)
            else:
                annot = annot.transmute(**{f'f_{field}': annot[f'f.{field}.0.0']})
            if i % 5 == 0 and i > 0:  # out of memory java error
                print(i)
                annot = annot.checkpoint(localfs_path + f'annot-{str(uuid4())}.ht')

        # drop duplicatd columns
        annot = annot.drop(*[
            col for col in annot.row_value.keys() if re.match(r'^f\..*_\d$', col)
        ])
        annot = annot.rename({'f.eid': 'f_eid'})
        annot.write(path)
    return hl.read_table(path)


def load_annotations_description(path):
    if not os.path.exists(path):
        df_left = pd.read_html('raw/dataset/ukb48815.html')
        df_left = df_left[1]
        df_left = df_left.drop('Column', axis=1)
        df_left = df_left.set_index('UDI')

        df_right = pd.read_html('raw/dataset/ukb49742.html')
        df_right = df_right[1]
        df_right = df_right.drop('Column', axis=1)
        df_right = df_right.set_index('UDI')

        df = df_left.merge(
            df_right,
            how='outer',
            left_on=['UDI', 'Type', 'Description'],
            right_on=['UDI', 'Type', 'Description'],
            suffixes=('_48815', '_49742')
        )

        df = df[['Type', 'Description', 'Count_48815', 'Count_49742']]
        df.to_csv(path)
        
    return hl.Table.from_pandas(pd.read_csv(path))


data = load_annotations_description(path='data/ukb-annotations-description.csv')
def find_annot_desc(query='', data=data, collapse=False):
    """Search descriptions of fields of annotation table"""
    desc_hits = data.Description.str.lower().str.contains(query.lower())
    udi_hits = data.index.str.contains(query)
    out = data[desc_hits | udi_hits]
    if collapse:
        out = out[out.index.str.endswith('0.0')]
    return out


def annotate_adr_patients(annots, overwrite=False):
    """Annotate table with ADR patients based on ICD10 codes"""
    path = localfs_path + 'annots_200k.ht'
    if not os.path.exists(path) or overwrite:
        # filter 200k
        f_pvcf = 'f_23156'
        annots_200k = annots.filter(annots[f_pvcf] == "977")

        # filter out withdraws
        with open(wd + 'raw/w62979_20210809.csv') as f:
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

        annots_200k = annots_200k.checkpoint(path, overwrite=overwrite)

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


def dc_recoding(expr, map_dict=None, update_dict={}):
    """Recode some categorical answers. Broad version."""
    if map_dict is None:
        map_dict = {
            "Yes": 1,
            "No": 0,
            "Do not know": 0,
            "Prefer not to answer": 0,
            "_missing": 0,
        }
    map_dict.update(update_dict)
    sb = hl.switch(expr)
    for key, val in map_dict.items():
        if key[0] != '_':
            sb = sb.when(key, val)

    if map_dict['_missing'] is None:
        return sb.or_missing()
    else:
        sb = sb.when_missing(map_dict['_missing'])
        return sb.or_error('SwithBulider Error')


def icd10_match(expr, code):
    """Find ICD10 diagnosis"""
    flat = hl.flatten(expr)
    flat = flat.map(
        lambda x: hl.if_else(hl.is_missing(x), False, x.startswith(code))
    )
    return hl.int(hl.any(flat))


def translate_icd10_codes(icd10_codes):
    """Translate ICD10 codes to diagnoses"""
    coding19 = hl.import_table(
        wd + 'raw/dataset/data-codings/coding19.tsv',
        key='coding'
    )
    codes = [c for code in icd10_codes for c in (code, code[:3])]
    coding19 = coding19.filter(
       hl.literal(codes).contains(coding19.coding)
    )
    return coding19


def annotate_genes(vcf, path):
    """Annotate matrix table with UKB dataset and gene names"""
    if not os.path.exists(path):
        genes = hl.read_table(
            '/net/archive/groups/plggneuromol/imdik-zekanowski-gts/' +
            '/data/external-data/genecode_v32.ht'
        )
        genes = genes.filter(
            hl.is_valid_contig(
                genes['hg38.knownGene.chrom'],
                reference_genome='GRCh38'
            )
        )
        start = genes['hg38.knownGene.txStart']
        stop =  genes['hg38.knownGene.txEnd']
        genes = genes.transmute(
            interval=hl.locus_interval(
                genes['hg38.knownGene.chrom'],
                start,
                stop,
                reference_genome='GRCh38',
                includes_start=False
            )
        )
        genes = genes.key_by(genes.interval)

        vcf = vcf.annotate_rows(
            within_gene=hl.array(
                hl.set(
                    genes.index(vcf.locus, all_matches=True)['hg38.kgXref.geneSymbol']
                )
            )
        )
        vcf = vcf.write(path)

    return hl.read_matrix_table(path)


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


def query_biomart(
        dataset, filters, attributes, url='http://ensembl.org/biomart'):
    """Make a query to a BioMart database
    dataset: str
    filters: dict
    attributes: list
    url: str, to a ensembl server (grch37: http://grch37.ensembl.org/biomart/)
        mirrors:
        http://uswest.ensembl.org/biomart/
        http://useast.ensembl.org/biomart/
        http://asia.ensembl.org/biomart/
    return: pandas dataframe
    """

    server = BiomartServer(url)
    # server.verbose = True  # set verbose to True to get some messages
    # server.show_databases()
    # server.show_datasets()

    dataset = server.datasets[dataset]
    # dataset.show_filters()
    # dataset.show_attributes()

    response = dataset.search(
        {'filters': filters, 'attributes': attributes},
        header=1
    )

    return pd.read_table(io.StringIO(response.text), sep='\t', low_memory=False)


def reorder_cols(mt, template_mt):
    """Reorder columns of a matrix table to match the template. Useful before
    binding matrix table rows."""
    mt = mt.add_col_index('tmp_col_idx')
    new_col_order = mt.index_cols(template_mt.col_key).tmp_col_idx.collect()
    mt = mt.choose_cols(new_col_order)
    mt = mt.drop('tmp_col_idx')
    return mt

import os
import requests
import time
import pickle
from pathlib import Path
from dataclasses import dataclass

import pandas as pd
import hail as hl
from tqdm import tqdm

from analysis.utils.liftover import liftover_vcf
from analysis.utils.pharmvar import is_sub, dist_to_feature
from analysis.utils.score_colnames import (
    info,
    functional_predictions_score, functional_predictions_rankscore,
    conservation_scores, conservation_rankscores
)


SC_TRANSLATE = {
    'adme': 'APF',
    'adme2': 'APF2',
    'pgs': 'PharmGScore',
    'cadd': 'CADD',
    'fathmm_xf': 'FATHMM-XF',
    'ma': 'MutationAssessor',
    'provean': 'PROVEAN',
    'pmls_pred': 'PharmMLScore',
    'sift': 'SIFT',
    'phylop100': 'PhyloP100'
}


@dataclass
class ScoreInfo:
    col_name: str
    limits: tuple
    deleterious_on_right: bool
    threshold: float = 0



class CADD:
    """Combined Annotation-Dependent Depletion (CADD) - a method for
    objectively integrating many diverse annotations into a single
    measure (C score) for each variant. CADD is implemented as a support
    vector machine trained to differentiate 14.7 million high-frequency
    human-derived alleles from 14.7 million simulated variants. C scores
    for all 8.6 billion possible human single-nucleotide variants was
    precomputed and scoring of short insertions-deletions was enabled.
    (http://dx.doi.org/10.1038/ng.2892)

    https://cadd.gs.washington.edu/download
    """

    table_path = os.path.join(RESOURCES_PATH, 'cadd.ht')

    def __init__(self, version='GRCh38-v1.6'):
        self.url_api = 'https://cadd.gs.washington.edu/api/v1.0'
        self.version = version

    def get_scores(self, vcf, sleep=0):
        scores = {'locus': [], 'allele': [], 'CADD_raw': [], 'CADD_phred': []}
        if vcf.count() > 1000 and sleep == 0:
            raise ValueError("ERROR: VCF count > 1000; define sleep")
        loci = vcf.locus.collect()
        alleles = vcf.alleles.collect()
        for locus, allele in zip(tqdm(loci), alleles):
            if locus.reference_genome.name != 'GRCh38':
                raise ValueError
            chrom, pos = locus.contig[3:], locus.position
            ref, alt = allele
            uri = f'{self.url_api}/{self.version}/{chrom}:{pos}_{ref}_{alt}'
            r = requests.get(uri)
            if r.status_code == 200:
                scores['locus'].append(str(locus))
                scores['allele'].append(allele)
                r_json = r.json()
                if len(r_json) == 1 and isinstance(r_json, list):
                    snv = r_json[0]
                    scores['CADD_raw'].append(float(snv['RawScore']))
                    scores['CADD_phred'].append(float(snv['PHRED']))
                else:
                    scores['CADD_raw'].append(hl.missing('float'))
                    scores['CADD_phred'].append(hl.missing('float'))
            else:
                raise ValueError(r.status_code)
            time.sleep(sleep)
        scores_pd = pd.DataFrame(scores)
        scores_pd['to_key'] = (
            scores_pd['locus'] + ':' + scores_pd['allele'].agg(':'.join)
        )
        scores_pd = scores_pd.drop(['locus', 'allele'], axis=1)
        scores_hl = hl.Table.from_pandas(scores_pd)
        scores_hl = scores_hl.key_by(
            **hl.parse_variant(scores_hl['to_key'], reference_genome='GRCh38')
        )
        scores_hl = scores_hl.drop(scores_hl['to_key'])

        return scores_hl

    @classmethod
    def _write_as_hailtable(cls, bgz_file):
        cadd = hl.import_table(
            bgz_file,
            force_bgz=True,
            min_partitions=2000,
            comment=r'^##.*',
            types={
                'RawScore': hl.tfloat32,
                'PHRED': hl.tfloat32,
            }
        )
        cadd = cadd.key_by(
            **hl.parse_variant(
                'chr' + cadd['#Chrom'] + ':'
                    + cadd['Pos'] + ':'
                    + cadd['Ref'] + ':'
                    + cadd['Alt'],
                reference_genome='GRCh38'
            )
        )
        cadd = cadd.select(
            score_raw=cadd['RawScore'],
            score_phred=cadd['PHRED'],
        )

        cadd.write(cls.table_path)

    @classmethod
    def read_table(cls):
        return hl.read_table(cls.table_path)


class FathmmXF:
    """Predictions are given as p-values in the range [0, 1]: values above 0.5
    are predicted to be deleterious, while those below 0.5 are predicted to be
    neutral or benign. P-values close to the extremes (0 or 1) are
    the highest-confidence predictions that yield the highest accuracy.

    We use distinct predictors for positions either in coding regions
    (positions within coding-sequence exons) or non-coding regions (positions
    in intergenic regions, introns or non-coding genes). The coding predictor
    is based on six groups of features representing sequence conservation,
    nucleotide sequence characteristics, genomic features (codons, splice
    sites, etc.), amino acid features and expression levels in different
    tissues. The non-coding predictor uses five feature groups that encompass
    nearly the same kinds of data, the primary exception being evidence for
    open chromatin.

    (http://fathmm.biocompute.org.uk/fathmm-xf/)

    """

    table_path_phred = os.path.join(RESOURCES_PATH, 'fathmm-xf-phred.ht')

    def __init__(self):
        pass

    def get_scores(self, vcf=None):
        pass

    @classmethod
    def read_table(cls):
        return hl.read_table(cls.table_path_phred)

    @classmethod
    def _write_as_hailtable(cls, bgz_files):
        """Downloaded from http://fathmm.biocompute.org.uk/fathmm-xf/"""

        coding_bgz_file = bgz_files['coding']
        fathmm_xf_coding = hl.import_table(
            coding_bgz_file,
            no_header=True,
            force_bgz=True,
            min_partitions=1000,
            types={
                'f4': hl.tfloat64,
            }
        )
        fathmm_xf_coding = fathmm_xf_coding.select(
            chr=fathmm_xf_coding.f0,
            position=fathmm_xf_coding.f1,
            ref=fathmm_xf_coding.f2,
            alt=fathmm_xf_coding.f3,
            coding_score=fathmm_xf_coding.f4,
        )
        fathmm_xf_coding = fathmm_xf_coding.key_by(
            **hl.parse_variant(
                'chr' + fathmm_xf_coding['chr'] + ':'
                    + fathmm_xf_coding['position'] + ':'
                    + fathmm_xf_coding['ref'] + ':'
                    + fathmm_xf_coding['alt'],
                reference_genome='GRCh38'
            )
        )
        fathmm_xf_coding = fathmm_xf_coding.select(
            'coding_score'
        )
        fathmm_xf_coding = fathmm_xf_coding.distinct()

        noncoding_bgz_file = bgz_files['noncoding']
        fathmm_xf_noncoding = hl.import_table(
            noncoding_bgz_file,
            no_header=True,
            force_bgz=True,
            min_partitions=5000,
            types={
                'f4': hl.tfloat64,
            }
        )
        fathmm_xf_noncoding = fathmm_xf_noncoding.select(
            chr=fathmm_xf_noncoding.f0,
            position=fathmm_xf_noncoding.f1,
            ref=fathmm_xf_noncoding.f2,
            alt=fathmm_xf_noncoding.f3,
            noncoding_score=fathmm_xf_noncoding.f4,
        )
        fathmm_xf_noncoding = fathmm_xf_noncoding.key_by(
            **hl.parse_variant(
                'chr' + fathmm_xf_noncoding['chr'] + ':'
                    + fathmm_xf_noncoding['position'] + ':'
                    + fathmm_xf_noncoding['ref'] + ':'
                    + fathmm_xf_noncoding['alt'],
                reference_genome='GRCh38'
            )
        )
        fathmm_xf_noncoding = fathmm_xf_noncoding.select(
            'noncoding_score'
        )
        fathmm_xf_noncoding = fathmm_xf_noncoding.distinct()

        fathmm_xf = hl.Table.union(
            fathmm_xf_coding,
            fathmm_xf_noncoding,
            unify=True
        )
        fathmm_xf = fathmm_xf.annotate(
            score=hl.if_else(
                hl.is_defined(fathmm_xf.coding_score),
                fathmm_xf.coding_score,
                fathmm_xf.noncoding_score
            )
        )

        fathmm_xf.write(cls.table_path)


    def prepare_input_for_online_tool(self, vcf, ref_gen='GRCh38'):
        """Return lines to paste into an online tool to get FATHMM-XF scores
        (http://fathmm.biocompute.org.uk/fathmmMKL.htm)."""

        # check first locus
        locus = vcf.take(1)[0].locus
        if locus.reference_genome.name == ref_gen:
            vcf_input = vcf
        else:
            print("liftover...")
            vcf_input = liftover_vcf(
                vcf,
                from_rg=locus.reference_genome.name,
                to_rg=ref_gen,
                chain_file=os.path.join(
                    RESOURCES_PATH,
                    'chain-files-for-liftover/grch38_to_grch37.over.chain.gz'
                )
            )
            vcf_input = vcf_input.persist()

        loci = vcf_input.locus.collect()
        alleles = vcf_input.alleles.collect()
        lines = [
            f'{locus.contig},{locus.position},{allele[0]},{allele[1]}'.lstrip('chr')
            for locus, allele in zip(loci, alleles)
        ]
        return lines


class dbNSFP:
    """dbNSFP is a database developed for functional prediction and annotation
    of all potential non-synonymous single-nucleotide variants (nsSNVs) in
    the human genome (https://sites.google.com/site/jpopgen/dbNSFP).

    More info: dbNSFP4.3a.readme.txt"""

    table_path = os.path.join(RESOURCES_PATH, 'dbNSFP4.3a', 'dbNSFP4.3a.ht')
    table_path_single = os.path.join(
        RESOURCES_PATH, 'dbNSFP4.3a', 'dbNSFP4.3a-single.ht'
    )

    @staticmethod
    def _to_float32_with_missing(score_str):
        score_str = hl.if_else(
            score_str.matches('^\.|-$|^$'),  # MutPred - no values indicated by '-'
            hl.missing('str'),
            hl.if_else(
                score_str.matches('\|'),
                hl.str(hl.max(score_str.split('\|').map(hl.float32))),
                score_str
            )
        )
        return hl.float32(score_str)

    @classmethod
    def _write_as_hailtable(cls, dbnsf_path):
        """dbNSFP downloaded from: https://sites.google.com/site/jpopgen/dbNSFP
        (googledrive)
        """
        dbnsfp = hl.import_table(
            dbnsf_path,
            min_partitions=3600,
            missing='.',
            impute=False,
        )
        dbnsfp = dbnsfp.key_by(
            **hl.parse_variant(
                'chr' + dbnsfp['#chr'] + ':'
                + dbnsfp['pos(1-based)'] + ':'
                + dbnsfp['ref'] + ':'
                + dbnsfp['alt']
            )
        )
        dbnsfp.write(cls.table_path)

    @classmethod
    def convert_to_single_score(cls):
        dbnsfp = cls.read_table(single=False)

        info = info + ['MutationTaster_pred']
        all_score_cols = (info
        + functional_predictions_score
        + functional_predictions_rankscore
        + conservation_scores
        + conservation_rankscores)

        dbnsfp = dbnsfp.select(
            *all_score_cols
        )

        dbnsfp = dbnsfp.annotate(
            **{
                col: dbnsfp[col].split(';') for col in info
            },
            **{
                col: dbnsfp[col].split(';').map(_to_float32_with_missing)
                for col in all_score_cols if col not in info
            },
        )
        dbnsfp = dbnsfp.annotate(
            VEP_canonical_index=dbnsfp['VEP_canonical'].index('YES'),
        )
        dbnsfp = dbnsfp.annotate(
            single={
                col:
                hl.if_else(
                    hl.is_defined(dbnsfp['VEP_canonical_index']),
                    hl.if_else(
                        (dbnsfp[col].length() == 1) | (hl.len(dbnsfp['VEP_canonical']) > hl.len(dbnsfp[col])),
                        dbnsfp[col][0],
                        dbnsfp[col][dbnsfp['VEP_canonical_index']]
                    ),
                    dbnsfp[col][0] if col in info else hl.mean(dbnsfp[col])
                )
                for col in all_score_cols
            }
        )
        dbnsfp.write(cls.table_path_single)

    @classmethod
    def read_table(cls, single=True):
        if single:
            return hl.read_table(cls.table_path_single)
        else:
            return hl.read_table(cls.table_path)


class PharmGScore:
    """PharmGScore table contains normalized scores for computing the
    PharmGScore. This score is an ensemble of four other scores: CADD,
    FATHMM-XF, PROVEAN and MutationAssessor.

    preparation: pharm_g_score/make_pgs_score.py
    """

    table_path = os.path.join(RESOURCES_PATH, 'pharm-g-score-60.ht')

    def __init__(self):
        pass

    def get_scores(self, vcf=None):
        pass

    @classmethod
    def read_table(cls):
        return hl.read_table(cls.table_path)


class ADME2:
    """ADME2 optimized predcition framework
    https://www.nature.com/articles/s41397-018-0044-2"""

    table_path = None
    adme_score_info = {
        'AlphaMissense': ScoreInfo('AM', (0, 1), True, 0.152),
        'PROVEAN': ScoreInfo('PROVEAN', (-15, 5), False, -3.2925),
        'MutationAssessor': ScoreInfo('MutationAssessor', (-4, 6), True, 2.08),
        'Polyphen2_HVAR': ScoreInfo('Polyphen2_HVAR', (0, 1), True, 0.42),
        'VEST4': ScoreInfo('VEST4', (0, 1), True, 0.2905),
    }

    def __init__(self):
        pass

    @classmethod
    def read_table(cls):
        pass

    @classmethod
    def _write_as_hailtable(cls, bgz_file):
        pass

    @classmethod
    def compute_score(cls, vcf):
        """Annotate vcf with burden scores and then compute PGS score"""
        am = AlphaMissense.read_table()
        dbnsfp = dbNSFP()
        db = dbnsfp.read_table()

        if isinstance(vcf, hl.Table):
            key = vcf.key
        else:
            key = vcf.row_key

        scores = hl.struct(
            am=am[key].am_pathogenicity,
            provean=db[key].single.PROVEAN_score,
            mutation_assessor=db[key].single.MutationAssessor_score,
            polyphen=db[key].single.Polyphen2_HVAR_score,
            vest=db[key].single.VEST4_score,
        )

        adme_components = [
            scores.am >= cls.adme_score_info['AlphaMissense'].threshold,
            scores.provean < cls.adme_score_info['PROVEAN'].threshold,
            scores.mutation_assessor > cls.adme_score_info['MutationAssessor'].threshold,
            scores.polyphen > cls.adme_score_info['Polyphen2_HVAR'].threshold,
            scores.vest > cls.adme_score_info['VEST4'].threshold,
        ]

        adme = hl.mean(adme_components)

        return hl.struct(
            adme_component_scores2=scores,
            adme_components2=adme_components,
            adme2=hl.if_else(hl.is_nan(adme), hl.missing(hl.tfloat64), adme),
        )


class ADME:
    """ADME optimized predcition framework
    https://www.nature.com/articles/s41397-018-0044-2"""

    table_path = None
    adme_score_info = {
        'LRT': ScoreInfo('LRT', (0, 1), False, 0.0025),
        'MutationAssessor': ScoreInfo('MutationAssessor', (-4, 6), True, 2.0566),
        'PROVEAN': ScoreInfo('PROVEAN', (-15, 5), False, -3.286),
        'VEST': ScoreInfo('VEST4', (0, 1), True, 0.4534),
        'CADD': ScoreInfo('CADD', (0, 60), True, 19.19),
    }

    def __init__(self):
        pass

    @classmethod
    def read_table(cls):
        pass

    @classmethod
    def _write_as_hailtable(cls, bgz_file):
        pass

    @classmethod
    def compute_score(cls, vcf):
        """Annotate vcf with burden scores and then compute PGS score"""
        cadd = CADD.read_table()
        dbnsfp = dbNSFP()
        db = dbnsfp.read_table()

        if isinstance(vcf, hl.Table):
            key = vcf.key
        else:
            key = vcf.row_key

        scores = hl.struct(
            lrt=db[key].single.LRT_score,
            mutation_assessor=db[key].single.MutationAssessor_score,
            provean=db[key].single.PROVEAN_score,
            vest=db[key].single.VEST4_score,
            cadd=cadd[key].score_phred,
        )

        adme_components = [
            scores.lrt < cls.adme_score_info['LRT'].threshold,
            scores.mutation_assessor > cls.adme_score_info['MutationAssessor'].threshold,
            scores.provean < cls.adme_score_info['PROVEAN'].threshold,
            scores.vest > cls.adme_score_info['VEST'].threshold,
            scores.cadd > cls.adme_score_info['CADD'].threshold,
        ]

        adme = hl.mean(adme_components)

        return hl.struct(
            adme_component_scores=scores,
            adme_components=adme_components,
            adme=hl.if_else(hl.is_nan(adme), hl.missing(hl.tfloat64), adme),
        )


class AlphaMissense:
    """AlphaMissense is a machine learning-based tool that uses unsupervised
    protein language modeling and structural context from AlphaFold to predict
    the pathogenicity of missense genetic variants in proteins.
    https://www.science.org/doi/10.1126/science.adg7492"""

    table_path = os.path.join(
        RESOURCES_PATH,
        'alphamissense',
        'alphamissense.ht'
    )

    def __init__(self):
        pass

    @classmethod
    def read_table(cls):
        return hl.read_table(cls.table_path)

    @classmethod
    def _write_as_hailtable(cls, bgz_file):
        bgz_file = os.path.join(RESOURCES_PATH, 'alphamissense', 'AlphaMissense_hg38.tsv.gz')
        am = hl.import_table(
            bgz_file,
            force_bgz=True,
            min_partitions=48,
            comment=[r'^#\s.*', r'^#$'],
            types={
                'am_pathogenicity': hl.tfloat32,
            }
        )
        am = am.key_by(
            **hl.parse_variant(
                am['#CHROM'] + ':'
                    + am['POS'] + ':'
                    + am['REF'] + ':'
                    + am['ALT'],
                reference_genome='GRCh38'
            )
        )
        am = am.select('transcript_id', 'am_pathogenicity')

        am.write(cls.table_path)


def annotate_vcf_with_scores(vcf, scores=None, af=True):
    # annotation with basic informations
    genes_positions = pd.read_csv(
        'data/pharmacogenetic-score/genes-position.tsv',
        sep='\t',
        low_memory=False
    )
    cadd = CADD.read_table()
    fathmm_xf = FathmmXF.read_table()
    pgs = PharmGScore.read_table()

    adme = ADME()
    adme2 = ADME2()
    dbnsfp = dbNSFP.read_table(single=False)
    nsfp = dbNSFP.read_table(single=True)

    if 'star_allele' in vcf.row:
        vcf = vcf.annotate(
            is_sub=is_sub(vcf)
        )
    if 'gene_name' not in vcf.row:
        vcf = vcf.annotate(
            gene_name=pgs[vcf.key].gene_name[0]
        )
    vcf = vcf.annotate(
        dist_to_exon=dist_to_feature(vcf, genes_positions, feature='exon')
    )
    vcf = vcf.persist()

    # splitted because of performance issues
    vcf = vcf.annotate(scores=hl.struct(
        pgs=pgs[vcf.key].pgs[0],
        cadd=cadd[vcf.key].score_raw,
        cadd_phred=cadd[vcf.key].score_phred,
        fathmm_xf=fathmm_xf[vcf.key].score,
        fathmm_xf_phred=fathmm_xf[vcf.key].phred,
        ma=nsfp[vcf.key].single.MutationAssessor_score,
        af=hl.float(dbnsfp[vcf.key].gnomAD_genomes_AF),
    ))
    vcf = vcf.persist()

    vcf = vcf.annotate(scores=vcf.scores.annotate(
        sift_phred=phred_rankscore_col(nsfp[vcf.key].single.SIFT_converted_rankscore),
        phylop100_phred=phred_rankscore_col(nsfp[vcf.key].single.phyloP100way_vertebrate_rankscore),
        ma_phred=phred_rankscore_col(nsfp[vcf.key].single.MutationAssessor_rankscore),

        sift=-nsfp[vcf.key].single.SIFT_score,
    ))
    # vcf = vcf.persist()

    vcf = vcf.annotate(scores=vcf.scores.annotate(
        provean=-nsfp[vcf.key].single.PROVEAN_score,
        phylop100=nsfp[vcf.key].single.phyloP100way_vertebrate,
        provean_phred=phred_rankscore_col(nsfp[vcf.key].single.PROVEAN_converted_rankscore),

        dann_phred=phred_rankscore_col(nsfp[vcf.key].single.DANN_rankscore),
        eigen_pc_phred=phred_rankscore_col(nsfp[vcf.key].single['Eigen-PC-raw_coding_rankscore' ]),

        bstatistic_phred=phred_rankscore_col(nsfp[vcf.key].single.bStatistic_converted_rankscore),
        clinpred_phred=phred_rankscore_col(nsfp[vcf.key].single.ClinPred_rankscore),
        fathmm_phred=phred_rankscore_col(nsfp[vcf.key].single.FATHMM_converted_rankscore),
        genocanyon_phred=phred_rankscore_col(nsfp[vcf.key].single.GenoCanyon_rankscore),
        gerp_phred=phred_rankscore_col(nsfp[vcf.key].single['GERP++_RS_rankscore' ]),
        gm_phred=phred_rankscore_col(nsfp[vcf.key].single.GM12878_fitCons_rankscore),
        h1_phred=phred_rankscore_col(nsfp[vcf.key].single['H1-hESC_fitCons_rankscore' ]),
        huvec_phred=phred_rankscore_col(nsfp[vcf.key].single.HUVEC_fitCons_rankscore),
        integraged_phred=phred_rankscore_col(nsfp[vcf.key].single.integrated_fitCons_rankscore),
    ))
    vcf = vcf.persist()

    vcf = vcf.annotate(scores=vcf.scores.annotate(
        list_phred=phred_rankscore_col(nsfp[vcf.key].single['LIST-S2_rankscore' ]),
        lrt_phred=phred_rankscore_col(nsfp[vcf.key].single.LRT_converted_rankscore),
        metasvm_phred=phred_rankscore_col(nsfp[vcf.key].single.MetaSVM_rankscore),
        metalr_phred=phred_rankscore_col(nsfp[vcf.key].single.MetaLR_rankscore),
        metarnn_phred=phred_rankscore_col(nsfp[vcf.key].single.MetaRNN_rankscore),
        mpc_phred=phred_rankscore_col(nsfp[vcf.key].single.MPC_rankscore),
        mvp_phred=phred_rankscore_col(nsfp[vcf.key].single.MVP_rankscore),
        phylop30_phred=phred_rankscore_col(nsfp[vcf.key].single.phyloP30way_mammalian_rankscore),
        phylop17_phred=phred_rankscore_col(nsfp[vcf.key].single.phyloP17way_primate_rankscore),
        phastcons100_phred=phred_rankscore_col(nsfp[vcf.key].single.phastCons100way_vertebrate_rankscore),
        phastcons30_phred=phred_rankscore_col(nsfp[vcf.key].single.phastCons30way_mammalian_rankscore),
        phastcons17_phred=phred_rankscore_col(nsfp[vcf.key].single.phastCons17way_primate_rankscore),
        primateai_phred=phred_rankscore_col(nsfp[vcf.key].single.PrimateAI_rankscore),
        sciphy_phred=phred_rankscore_col(nsfp[vcf.key].single.SiPhy_29way_logOdds_rankscore),
    ))

    vcf = vcf.persist()

    adme_score = adme.compute_score(vcf)
    adme2_score = adme2.compute_score(vcf)
    vcf = vcf.annotate(scores=vcf.scores.annotate(
        **adme_score,
        **adme2_score,
    ))
    vcf = vcf.persist()

    # normalize phreds to be 0 - 1
    def get_stats(ht, name, expr):
        path = Path(f'tmp/{name}-stats.pkl')
        if path.is_file():
            with open(path, 'rb') as f:
                stats = pickle.load(f)
        else:
            if name in ['cadd', 'fathmm_xf']:
                stats = ht.aggregate(hl.agg.stats(expr))
            else:
                stats = ht.aggregate(hl.agg.stats(phred_rankscore_col(expr)))
            with open(path, 'wb') as f:
                pickle.dump(stats, f)
        return stats

    cadd_stats = get_stats(cadd, 'cadd', cadd.score_phred)
    fathmm_xf_stats = get_stats(fathmm_xf, 'fathmm_xf', fathmm_xf.phred)
    ma_stats = get_stats(nsfp, 'ma', nsfp.single.MutationAssessor_rankscore)
    provean_stats = get_stats(nsfp, 'provean', nsfp.single.PROVEAN_converted_rankscore)
    phylop_stats = get_stats(nsfp, 'phylop', nsfp.single.phyloP100way_vertebrate_rankscore)

    dann_stats = get_stats(nsfp, 'dann', nsfp.single.DANN_rankscore)
    bayes_noaf_stats = get_stats(nsfp, 'bayes_noaf', nsfp.single.BayesDel_noAF_rankscore)
    cadd_stats = get_stats(cadd, 'cadd', cadd.score_phred)
    sift4g_stats = get_stats(nsfp, 'sift4g', nsfp.single.SIFT4G_converted_rankscore)
    sift_stats = get_stats(nsfp, 'sift', nsfp.single.SIFT_converted_rankscore)
    fathmm_mkl_stats = get_stats(nsfp, 'fathmm_mkl', nsfp.single['fathmm-MKL_coding_rankscore'])
    mt_stats = get_stats(nsfp, 'mt', nsfp.single.MutationTaster_converted_rankscore)
    vest4_stats = get_stats(nsfp, 'vest4', nsfp.single.VEST4_rankscore)
    lrt_stats = get_stats(nsfp, 'lrt', nsfp.single.LRT_converted_rankscore)
    bayes_addaf_stats = get_stats(nsfp, 'bayes_addaf', nsfp.single.BayesDel_addAF_rankscore)
    revel_stats = get_stats(nsfp, 'revel', nsfp.single.REVEL_rankscore)
    eigen_stats = get_stats(nsfp, 'eigen', nsfp.single['Eigen-raw_coding_rankscore'])
    eigen_pc_stats = get_stats(nsfp, 'eigen_pc', nsfp.single['Eigen-PC-raw_coding_rankscore'])


    bstatistic_stats = get_stats(nsfp, 'bstatistic', nsfp.single.bStatistic_converted_rankscore)
    clinpred_stats = get_stats(nsfp, 'clinpred', nsfp.single.ClinPred_rankscore)
    fathmm_stats = get_stats(nsfp, 'fathmm', nsfp.single.FATHMM_converted_rankscore)
    genocanyon_stats = get_stats(nsfp, 'genocanyon', nsfp.single.GenoCanyon_rankscore)
    gerp_stats = get_stats(nsfp, 'gerp', nsfp.single['GERP++_RS_rankscore' ])
    gm_stats = get_stats(nsfp, 'gm', nsfp.single.GM12878_fitCons_rankscore)
    h1_stats = get_stats(nsfp, 'h1', nsfp.single['H1-hESC_fitCons_rankscore' ])
    huvec_stats = get_stats(nsfp, 'huvec', nsfp.single.HUVEC_fitCons_rankscore)
    integraged_stats = get_stats(nsfp, 'integraged', nsfp.single.integrated_fitCons_rankscore)
    list_stats = get_stats(nsfp, 'list', nsfp.single['LIST-S2_rankscore' ])
    lrt_stats = get_stats(nsfp, 'lrt', nsfp.single.LRT_converted_rankscore)
    metasvm_stats = get_stats(nsfp, 'metasvm', nsfp.single.MetaSVM_rankscore)
    metalr_stats = get_stats(nsfp, 'metalr', nsfp.single.MetaLR_rankscore)
    metarnn_stats = get_stats(nsfp, 'metarnn', nsfp.single.MetaRNN_rankscore)
    mpc_stats = get_stats(nsfp, 'mpc', nsfp.single.MPC_rankscore)
    mvp_stats = get_stats(nsfp, 'mvp', nsfp.single.MVP_rankscore)
    phylop30_stats = get_stats(nsfp, 'phylop30', nsfp.single.phyloP30way_mammalian_rankscore)
    phylop17_stats = get_stats(nsfp, 'phylop17', nsfp.single.phyloP17way_primate_rankscore)
    phastcons100_stats = get_stats(nsfp, 'phastcons100', nsfp.single.phastCons100way_vertebrate_rankscore)
    phastcons30_stats = get_stats(nsfp, 'phastcons30', nsfp.single.phastCons30way_mammalian_rankscore)
    phastcons17_stats = get_stats(nsfp, 'phastcons17', nsfp.single.phastCons17way_primate_rankscore)
    primateai_stats = get_stats(nsfp, 'primateai', nsfp.single.PrimateAI_rankscore)
    sciphy_stats = get_stats(nsfp, 'sciphy', nsfp.single.SiPhy_29way_logOdds_rankscore)

    def norm_score(expr, clip_max=40, max_val=None):
        return hl.min(expr, clip_max, filter_missing=False) / hl.min(max_val, clip_max)

    vcf = vcf.annotate(scores=vcf.scores.annotate(
        cadd_norm=norm_score(vcf.scores.cadd_phred, max_val=cadd_stats['max']),
        fathmm_xf_norm=norm_score(vcf.scores.fathmm_xf_phred, max_val=fathmm_xf_stats['max']),
        sift_norm=norm_score(vcf.scores.sift_phred, max_val=sift_stats['max']),
        ma_norm=norm_score(vcf.scores.ma_phred, max_val=ma_stats['max']),
        provean_norm=norm_score(vcf.scores.provean_phred, max_val=provean_stats['max']),
        phylop100_norm=norm_score(vcf.scores.phylop100_phred, max_val=phylop_stats['max']),

        dann_norm=norm_score(vcf.scores.dann_phred, max_val=dann_stats['max']),
        eigen_pc_norm=norm_score(vcf.scores.eigen_pc_phred, max_val=eigen_pc_stats['max']),

        bstatistic_norm=norm_score(vcf.scores.bstatistic_phred, max_val=bstatistic_stats['max']),
        clinpred_norm=norm_score(vcf.scores.clinpred_phred, max_val=clinpred_stats['max']),
        fathmm_norm=norm_score(vcf.scores.fathmm_phred, max_val=fathmm_stats['max']),
        genocanyon_norm=norm_score(vcf.scores.genocanyon_phred, max_val=genocanyon_stats['max']),
        gerp_norm=norm_score(vcf.scores.gerp_phred, max_val=gerp_stats['max']),
        gm_norm=norm_score(vcf.scores.gm_phred, max_val=gm_stats['max']),
        h1_norm=norm_score(vcf.scores.h1_phred, max_val=h1_stats['max']),
        huvec_norm=norm_score(vcf.scores.huvec_phred, max_val=huvec_stats['max']),
        integraged_norm=norm_score(vcf.scores.integraged_phred, max_val=integraged_stats['max']),
    ))
    vcf = vcf.persist()

    vcf = vcf.annotate(scores=vcf.scores.annotate(
        list_norm=norm_score(vcf.scores.list_phred, max_val=list_stats['max']),
        lrt_norm=norm_score(vcf.scores.lrt_phred, max_val=lrt_stats['max']),
        metasvm_norm=norm_score(vcf.scores.metasvm_phred, max_val=metasvm_stats['max']),
        metalr_norm=norm_score(vcf.scores.metalr_phred, max_val=metalr_stats['max']),
        metarnn_norm=norm_score(vcf.scores.metarnn_phred, max_val=metarnn_stats['max']),
        mpc_norm=norm_score(vcf.scores.mpc_phred, max_val=mpc_stats['max']),
        mvp_norm=norm_score(vcf.scores.mvp_phred, max_val=mvp_stats['max']),
        phylop30_norm=norm_score(vcf.scores.phylop30_phred, max_val=phylop30_stats['max']),
        phylop17_norm=norm_score(vcf.scores.phylop17_phred, max_val=phylop17_stats['max']),
        phastcons100_norm=norm_score(vcf.scores.phastcons100_phred, max_val=phastcons100_stats['max']),
        phastcons30_norm=norm_score(vcf.scores.phastcons30_phred, max_val=phastcons30_stats['max']),
        phastcons17_norm=norm_score(vcf.scores.phastcons17_phred, max_val=phastcons17_stats['max']),
        primateai_norm=norm_score(vcf.scores.primateai_phred, max_val=primateai_stats['max']),
        sciphy_norm=norm_score(vcf.scores.sciphy_phred, max_val=sciphy_stats['max']),
    ))
    vcf = vcf.persist()

    return vcf


def phred_rankscore_col(column_expr):
    return -10 * hl.log10(1 - hl.float32(column_expr))


def phred_score(ht, score_col_name, out_col_name=None, smaller_more_damaging=False):
    print(f'{datetime.now()}: start', flush=True)
    if not out_col_name:
        out_col_name = f'{score_col_name}_phred'

    undefined_missing_expr = (
        ~hl.is_finite(ht[score_col_name])
        | hl.is_missing(ht[score_col_name])
    )
    n_undefined_missing = ht.filter(undefined_missing_expr).count()
    print(n_undefined_missing)

    if smaller_more_damaging:
        ht = ht.order_by(
            ~undefined_missing_expr,
            ht[score_col_name]
        )
    else:
        ht = ht.order_by(
            ~undefined_missing_expr,
            hl.desc(ht[score_col_name])
        )

    ht = ht.persist()
    print(f'{datetime.now()}: persist 1', flush=True)

    ht = ht.add_index('score_index')
    ht = ht.persist()
    print(f'{datetime.now()}: persist 2', flush=True)

    ht_ranks = ht.filter(hl.is_defined(ht[score_col_name]))
    ht_ranks = (
        ht_ranks
        .group_by(ht_ranks[score_col_name])
        .aggregate(
            score_min_rank=hl.agg.min(ht_ranks.score_index) + 1 - n_undefined_missing
        )
    )
    n = ht_ranks.aggregate(hl.agg.max(ht_ranks.score_min_rank))
    ht = ht.annotate(
        percent_rank=ht_ranks[ht[score_col_name]].score_min_rank / n,
    )
    ht = ht.transmute(**{
        out_col_name: -10 * hl.log10(ht.percent_rank)
    })
    ht = ht.drop('score_index')
    ht = ht.persist()
    print(f'{datetime.now()}: persist 3', flush=True)

    ht = ht.key_by('locus', 'alleles')
    ht = ht.persist()
    print(f'{datetime.now()}: persist 4', flush=True)

    return ht

import os
import requests
import time

import pandas as pd
import hail as hl
from tqdm import tqdm

from analysis.utils.liftover import liftover_vcf


RESOURCES_PATH = '/net/archive/groups/plggneuromol/resources/'


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


    This base has a duplicated variants (about 0.3%). Comparing around 1000
    duplicated scores from the database with the ones retrived online it seems
    it's enough to do `hl.distinct()`. Example variants:
        * chr8:144138786 G/A/C/T
        * chr1:146989574 G/A/C/T - 7 scores for each variant
    """


    table_path = os.path.join(RESOURCES_PATH, 'fathmm-xf.ht')

    def __init__(self):
        pass

    def get_scores(self, vcf=None):
        pass

    @classmethod
    def read_table(cls):
        return hl.read_table(cls.table_path)

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
            vcf_input = vcf_input.checkpoint(next(tmpdir_path))

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

    @staticmethod
    def _to_float32_with_missing(score_str):
        score_str = hl.if_else(
            score_str.matches('^\.$'),
            hl.missing('str'),
            score_str
        )
        return hl.float32(score_str)

    def convert_subset(self, dbnsfp_subset_path):
        to_do_nothing = [
        ]
        to_split = [
        ]
        to_split_arrayfloat = [
            'MutationAssessor_score',
            'PROVEAN_score',
        ]
        to_float = [
            'MutationAssessor_rankscore',
            'PROVEAN_converted_rankscore',
        ]
        SCORE_COLS = to_do_nothing + to_split + to_split_arrayfloat + to_float
        annot_cols = ['genename', 'Ensembl_geneid', 'Ensembl_transcriptid']

        if os.path.exists(dbnsfp_subset_path):
            return hl.read_table(dbnsfp_subset_path)

        dbnsfp = self.read_table()
        dbnsfp = dbnsfp.select(
            *annot_cols, *SCORE_COLS,
        )
        dbnsfp = dbnsfp.annotate(
            **{
                col: dbnsfp[col].split(';') for col in annot_cols + to_split
            },
            **{
                col: hl.float32(dbnsfp[col]) for col in to_float
            },
            **{
                col: hl.map(self._to_float32_with_missing, dbnsfp[col].split(';'))
                for col in to_split_arrayfloat
            }
        )
        dbnsfp = dbnsfp.annotate(
            **{
                'mutation_assessor': hl.max(dbnsfp['MutationAssessor_score']),  # -5.17-6.49; >0.65 H(igh)/M/L/N(eutral) 3.5/1.935/0.8???
                'provean': hl.min(dbnsfp['PROVEAN_score']),  # -14-14; <=-2.5 D
                }
            )
        dbnsfp = dbnsfp.repartition(1000)  # can't read table later without it
        dbnsfp.write(dbnsfp_subset_path)

        return hl.read_table(dbnsfp_subset_path)

    @classmethod
    def _write_as_hailtable(cls, dbnsf_path):
        """dbNSFP downloaded from: https://sites.google.com/site/jpopgen/dbNSFP
        (googledrive)
        Then variants converted from *.gz to *.bgz with bgzip (htslib/samtools)
        """
        dbnsfp = hl.import_table(
            dbnsf_path,
            min_partitions=3600,
            missing='.',
            impute=False,  # NotImplementedError 'str' + 'int32'
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
    def read_table(cls):
        return hl.read_table(cls.table_path)


class DANN:
    """DANN is a functional prediction score retrained based on the
    training data of CADD using deep neural network. Scores range from
    0 to 1. A larger number indicate a higher probability to be damaging.

    DANN scores were downloaded from
    https://cbcl.ics.uci.edu/public_data/DANN/data/ and liftovered from
    GRCh37 to GRCh38.
    """

    table_path = os.path.join(RESOURCES_PATH, 'dann.ht')

    def __init__(self):
        pass

    def get_scores(self, vcf=None):
        pass

    @classmethod
    def read_table(cls):
        return hl.read_table(cls.table_path)

    @classmethod
    def _write_as_hailtable(cls, bgz_file):
        dann = hl.import_table(
            bgz_file,
            no_header=True,
            min_partitions=5000,
            types={
                'f4': hl.tfloat32,
            }
        )
        dann = dann.select(
            chr=dann.f0,
            position=dann.f1,
            ref=dann.f2,
            alt=dann.f3,
            score=dann.f4,
        )
        dann = dann.key_by(
            **hl.parse_variant(
                dann['chr'] + ':'
                    + dann['position'] + ':'
                    + dann['ref'] + ':'
                    + dann['alt'],
                reference_genome='GRCh37'
            )
        )

        dann_38 = liftover_vcf(
            dann,
            'GRCh37',
            'GRCh38',
            os.path.join(
                RESOURCES_PATH,
                'chain-files-for-liftover',
                'grch37_to_grch38.over.chain.gz'
            )
        )
        dann_38 = dann_38.select(
            'score'
        )

        dann_38.write(cls.table_path)


class PharmGScore:
    """PharmGScore table contains normalized scores for computing the
    PharmGScore. This score is an ensemble of four other scores: CADD,
    FATHMM XF, PROVEAN and MutationAssessor. DANN score in the table
    is NOT normalized.

    PharmGScore scores bgz file was produced by merging files:
        cat gene-our-scores/*.tsv | head -n 1 | bgzip > our-scores.tsv.bgz && \
        cat gene-our-scores/*.tsv | tail -n +2 | bgzip --threads 7 >> \
        our-scores.tsv.bgz

    with normalized scores genereated by the following script:
        analysis/star_alleles_burden/burden_function_plots.Rmd
    """

    table_path = os.path.join(RESOURCES_PATH, 'pharm-g-score-scores.ht')

    def __init__(self):
        pass

    def get_scores(self, vcf=None):
        pass

    @classmethod
    def read_table(cls):
        return hl.read_table(cls.table_path)

    @classmethod
    def _write_as_hailtable(cls, bgz_file):
        pgs = hl.import_table(
            bgz_file,
            types={
                'locus': hl.tstr,
                'allele_1': hl.tstr,
                'allele_2': hl.tstr,
                'gene_name': hl.tstr,
                'dist_to_exon': hl.tint32,
                'dist_to_gene': hl.tint32,
                'gnomad_nfe_AF': hl.tfloat32,
                'in_gnomad': hl.tbool,
                'cadd_raw': hl.tfloat32,
                'mutation_assessor': hl.tfloat32,
                'provean': hl.tfloat32,
                'fathmm_xf': hl.tfloat32,
                'dann': hl.tfloat32,  # it's oryginal score, not normalized
                'cadd_raw_oryg': hl.tfloat32,
                'fathmm_xf_oryg': hl.tfloat32,
                'provean_oryg': hl.tfloat32,
                'mutation_assessor_oryg': hl.tfloat32,
            },
        )
        pgs = pgs.annotate(
            scores=hl.array([
                pgs.cadd_raw,
                pgs.mutation_assessor,
                pgs.provean,
                pgs.fathmm_xf
            ])
        )
        pgs = pgs.transmute(
            pgs=hl.if_else(
                hl.all(hl.map(hl.is_missing, pgs.scores)),
                0,
                hl.mean(pgs.scores)
            )
        )

        pgs = pgs.annotate(
            **hl.parse_variant(
                pgs.locus + ':' + pgs.allele_1 + ':' + pgs.allele_2
            )
        )
        pgs = pgs.drop('allele_1', 'allele_2')

        key_cols = ['locus', 'alleles', 'gene_name']
        pgs = pgs.key_by(*key_cols)
        # reorder columns
        pgs = pgs.select(
            *[col for col in pgs._fields.keys() if col not in key_cols]
        )
        pgs = pgs.repartition(12000)

        pgs.write(cls.table_path)

    @staticmethod
    def _score_aggregator(expr_score, expr_gt):
        return hl.if_else(
            hl.agg.all(hl.is_missing(expr_score) | expr_gt.is_hom_ref()),
            hl.missing(hl.tfloat),
            hl.agg.sum(expr_score * expr_gt.n_alt_alleles())
        )

    @classmethod
    def compute_burden(cls, vcf):
        """Annotate vcf with burden scores and then compute PGS score"""
        vcf = hl.split_multi_hts(vcf, permit_shuffle=True)
        # TODO: remove non-variant sites?
        vcf = vcf.semi_join_rows(cls.read_table())
        vcf = vcf.annotate_rows(
            **cls.read_table()[vcf.row_key]
        )

        vcfe = vcf.entries()
        vg = vcfe.group_by('s', 'gene_name').aggregate(
            n=hl.agg.sum(vcfe.GT.n_alt_alleles()),
            cadd_raw=cls._score_aggregator(vcfe.cadd_raw, vcfe.GT),
            fathmm_xf=cls._score_aggregator(vcfe.fathmm_xf, vcfe.GT),
            provean=cls._score_aggregator(vcfe.provean, vcfe.GT),
            mutation_assessor=cls._score_aggregator(vcfe.mutation_assessor, vcfe.GT),
            score_1=hl.agg.sum(vcfe.cadd_raw_oryg * vcfe.GT.n_alt_alleles()),
            score_2=hl.agg.sum(vcfe.fathmm_xf_oryg * vcfe.GT.n_alt_alleles()),
            score_3=hl.agg.sum(vcfe.provean_oryg * vcfe.GT.n_alt_alleles()),
            score_4=hl.agg.sum(vcfe.mutation_assessor_oryg * vcfe.GT.n_alt_alleles()),
        )
        vg = vg.annotate(
            scores=hl.array(
                [vg.cadd_raw, vg.fathmm_xf, vg.provean, vg.mutation_assessor]
            ),
        )
        vg = vg.annotate(
            burden=hl.if_else(
                vg.scores.all(hl.is_missing),
                0,
                hl.mean(vg.scores),
            )
        )
        return vg


if __name__ == '__main__':
    from analysis.utils.load_spark import hl_init, wd
    from analysis.utils.load_spark import tmpdir_path_iter

    hl_init()
    hl.plot.output_notebook()
    tmpdir_path = tmpdir_path_iter()

    f = FathmmXF._write_as_hailtable({
        'coding': wd + 'raw/scores/fathmm-xf/fathmm_xf_coding_hg38.vcf.gz',
        'noncoding': wd + 'raw/scores/fathmm-xf/fathmm_xf_noncoding_hg38.vcf.gz'
    })

    # dbNSFP exports
    dbNSFP._write_as_hailtable(
        wd + 'raw/scores/dbNSFP4.3a/dbNSFP4.3a_variant.chr*.bgz',
    )

    dbnsfp = dbNSFP()
    dbnsfp.convert_subset(
        dbnsfp_subset_path=os.path.join(
            '/net/people/plgjacekh/projects/ifpan-gosborcz-ukb/',
            'data/scores/dbNSFP4.3a-subset.ht'
        )
    )

    DANN._write_as_hailtable(
        os.path.join(RESOURCES_PATH, 'dann', 'DANN_whole_genome_SNVs.tsv.bgz')
    )

    CADD._write_as_hailtable(
        os.path.join(RESOURCES_PATH, 'cadd', 'whole_genome_SNVs.tsv.gz')
    )

    PharmGScore._write_as_hailtable(
        os.path.join(wd, 'data', 'pharmacogenetic-score', 'our-scores.tsv.bgz')
    )

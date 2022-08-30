import os
import re
import requests
from typing import Dict, List
from enum import Enum, IntEnum, auto
from dataclasses import dataclass, field

import hail as hl

from analysis.utils.load_spark import wd, localfs_path, scratch_path
from analysis.utils.load_spark import tmpdir_path_iter


tmpdir_path = tmpdir_path_iter()
RESOURCES_PATH = '/net/archive/groups/plggneuromol/resources/'


def starall_natsort_key(dit, key='alleleName'):
    all_name = dit[key]
    try:
        n = re.search(r'(?:rs|\*)(\d+)\.?(\d+)?', all_name)
        return tuple(int(i) for i in n.groups() if i is not None)
    except AttributeError:
        return all_name


class EvidenceLevel(IntEnum):
    NA0 = -1
    LIMITED = 1
    MODERATE = 2
    DEFINITIVE = 3


class Function(Enum):
    FUNCTION_NOT_ASSIGNED = auto()
    NO_FUNCTION = auto()
    UNCERTAIN_FUNCTION = auto()
    UNKNOWN_FUNCTION = auto()
    NORMAL_FUNCTION = auto()
    DECREASED_FUNCTION = auto()
    INCREASED_FUNCTION = auto()


@dataclass
class Variant:
    hgvs: str


@dataclass
class Allele:
    allele_function: Function
    evidence_level: str
    variants: List[Variant]
    vcf_path: str


@dataclass
class Gene:
    gene_symbol: str
    gene_name: str
    alleles: Dict[str, Allele] = field(default_factory=dict)


class PharmVar:
    table_path = os.path.join(RESOURCES_PATH, 'pharmvar', 'pharmvar-vcf.ht')

    def __init__(self):
        self.base_url = 'https://www.pharmvar.org/'
        self.reference_collection = 'GRCh38'
        self.evidence_level_dict = {
            '0': 'NA0',
            'L': 'LIMITED',
            'M': 'MODERATE',
            'D': 'DEFINITIVE'
        }
        self. all_pharmvar_genes = {
            'database': [
                'CYP2A13', 'CYP2B6', 'CYP2C8', 'CYP2C9', 'CYP2C19', 'CYP2D6',
                'CYP2F1', 'CYP2J2', 'CYP2R1', 'CYP2S1', 'CYP2W1',
                'CYP3A4', 'CYP3A5', 'CYP3A7', 'CYP3A43',
                'CYP4F2',
                'SLCO1B1',
                'DPYD', 'NUDT15',
            ],
            'off_database': [
                'CYP1A1', 'CYP1A2', 'CYP1B1',
                'CYP2A6', 'CYP2E1',
            ],
            'legacy': [
                'CYP26A1',
            ]
        }
        self.genes = self.populate_genes()

    def request_genes(self):
        api_url = f'{self.base_url}/api-service/genes'
        query = {
            'reference-collection': self.reference_collection
        }
        response = requests.get(api_url, params=query)
        if response.status_code == 200:
            return response.json()
        else:
            raise ValueError(f"Response status code: {response.status_code}")

    def populate_genes(self):
        response = self.request_genes()
        genes = {}
        for gene_json in response:
            gene_symbol = gene_json['geneSymbol']
            alleles_dict = {}
            for allele_json in gene_json['alleles']:
                allele_name = allele_json['alleleName']
                variants_list = []
                for variants_json in allele_json['variants']:
                    vcf_path = os.path.join(
                        'raw',
                        'pharmvar-5.1.6',
                        gene_symbol,
                        self.reference_collection,
                        allele_name.replace('*', '_').replace(' ', '-') + '.vcf'
                    )
                    if not os.path.exists(vcf_path):
                        raise FileNotFoundError(wd + vcf_path)
                    variant = Variant(
                        hgvs=variants_json['hgvs'],
                    )
                    variants_list.append(variant)
                allele = Allele(
                    allele_function=allele_json['function'],
                    evidence_level=allele_json['evidenceLevel'],
                    variants=variants_list,
                    vcf_path=vcf_path
                )
                alleles_dict[allele_name] = allele
            gene = Gene(
                gene_symbol=gene_symbol,
                gene_name=gene_json['geneName'],
                alleles=alleles_dict
            )
            genes[gene_json['geneSymbol']] =  gene
        return genes

    def union_starallele_vcfs(self):
        """Import vcf file for an each allel in pharmvar and union them.
        Execution time: about 12 minutes """
        vcf_list = []
        for gene_symbol, pv_gene, in self.genes.items():
            for allele_name, pv_allele in pv_gene.alleles.items():
                vcf = hl.import_vcf(pv_allele.vcf_path)
                vcf = vcf.annotate_rows(
                    gene=gene_symbol,
                    star_allele=allele_name,
                    function=pv_allele.allele_function,
                    evidence_level=pv_allele.evidence_level,
                )
                vcf_list.append(vcf)

        # splitting list because of java error
        from math import sqrt
        m = int(sqrt(len(vcf_list))) + 1
        vcf_sub_lists = [vcf_list[i:i+m] for i in range(0, len(vcf_list), m)]
        sub_vcfs = []
        for i, vcf_sub_list in enumerate(vcf_sub_lists):
            sub_vcf = hl.MatrixTable.union_rows(*vcf_sub_list)
            sub_vcf = sub_vcf.checkpoint(next(tmpdir_path))
            sub_vcfs.append(sub_vcf)

        vcf_merged = hl.MatrixTable.union_rows(*sub_vcfs)
        vcf_merged = vcf_merged.checkpoint(next(tmpdir_path))
        return vcf_merged

    @classmethod
    def write_staralleles_table(cls, vcf):
        vcf.rows().write(cls.table_path)

    @classmethod
    def read_table(cls):
        return hl.read_table(cls.table_path)

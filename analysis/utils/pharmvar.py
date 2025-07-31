import os
import re
import json
import requests
from datetime import datetime
from typing import Dict, List
from collections import defaultdict
from enum import Enum, IntEnum, auto
from dataclasses import dataclass, field

import hail as hl


ALL_PV_GENES = {
    # 'CYP1A2',  # not present in 6.0.10
    'CYP2A6',
    'CYP2A13',
    'CYP2B6',
    'CYP2C8',
    'CYP2C9',
    'CYP2C19',
    'CYP2D6',
    'CYP3A4',
    'CYP3A5',
    'CYP4F2',
    'DPYD',
    # 'NAT2',  # not present in 6.0.10
    'NUDT15',
    'SLCO1B1',
}

PV_GENES = {
    'CYP2B6',
    'CYP2C9',
    'CYP2C19',
    'CYP2D6',
    'CYP3A5',
#    'CYP4F2',
    'DPYD',
    'NUDT15',
    'SLCO1B1',
}

CPIC_GENES = {
    'CYP2B6',
    'CYP2C9',
    'CYP2C19',
    'CYP2D6',
    'CYP3A5',
#    'CYP4F2',  # no recommendations
    'DPYD',
    'NUDT15',
    'SLCO1B1',
}

FUN_MAP = {
    'decreased function': 'decreased function',
    'function not assigned': 'undefined function',
    'increased function': 'increased function',
    'no function': 'no function',
    'normal function': 'normal function',
    'uncertain function': 'undefined function',
    'unknown function': 'undefined function',
}
FUN_MAP_ORDER = [
    'increased function', 'normal function', 'decreased function',
    'no function', 'undefined function'
]


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
    def __init__(self, version):
        self.base_url = 'https://www.pharmvar.org/'
        self.reference_collection = 'GRCh38'
        self.version = version
        self.table_path = os.path.join(
            RESOURCES_PATH,
            'pharmvar',
            f'pharmvar-vcf-{self.version}.ht'
        )
        self.table_filtered_path = os.path.join(
            'data',
            f'pharmvar-{self.version}-annotated-filtered.ht'
        )
        self.genes = None

    def request_genes(self, dump_json=False):
        api_url = f'{self.base_url}/api-service/genes'
        query = {
            'reference-collection': self.reference_collection
        }
        response = requests.get(api_url, params=query)
        if response.status_code == 200:
            genes = response.json()
            if dump_json:
                now = datetime.now().strftime("%Y-%m-%d-%H%M%S")
                json_path = os.path.join(
                    RESOURCES_PATH,
                    'pharmvar',
                    f'genes-{now}-v{self.version}.json'  # API doesn't control version, be careful
                )
                with open(json_path, 'w') as f:
                    json.dump(genes, f)
            return genes
        else:
            raise ValueError(f"Response status code: {response.status_code}")

    def populate_genes(self):
        response = self.request_genes(dump_json=True)
        genes = {}
        for gene_json in response:
            gene_symbol = gene_json['geneSymbol']
            alleles_dict = {}
            for allele_json in gene_json['alleles']:
                allele_name = allele_json['alleleName']
                variants_list = []
                for variants_json in allele_json['variants']:
                    vcf_path = os.path.join(
                        RESOURCES_PATH,
                        'pharmvar',
                        f'pharmvar-{self.version}',
                        gene_symbol,
                        self.reference_collection,
                        allele_name.replace('*', '_').replace(' ', '-') + '.vcf'
                    )
                    if not os.path.exists(vcf_path):
                        raise FileNotFoundError(vcf_path)
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
        self.genes = genes

    def union_starallele_vcfs(self):
        """Import vcf file for an each allel in pharmvar and union them.
        Execution time: about 12 minutes """
        vcf_list = []
        for gene_symbol, pv_gene, in self.genes.items():
            for allele_name, pv_allele in pv_gene.alleles.items():
                vcf = hl.import_vcf(pv_allele.vcf_path)
                vcf = vcf.annotate_rows(
                    gene_name=gene_symbol,
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
            sub_vcf = sub_vcf.persist()
            sub_vcfs.append(sub_vcf)

        vcf_merged = hl.MatrixTable.union_rows(*sub_vcfs)
        vcf_merged = vcf_merged.persist()
        return vcf_merged

    def write_staralleles_table(self, vcf):
        vcf.rows().write(self.table_path)

    def read_table(self):
        return hl.read_table(self.table_path)
    
    def read_table_filtered(self):
        return hl.read_table(self.table_filtered_path)


def dist_to_feature(ht, positions, feature):
    gene_starts = positions.groupby('Gene name')['Gene start (bp)'].nunique()
    gene_starts.name = 'n_start'
    positions_n_starts = positions.merge(
        gene_starts,
        how='outer',
        left_on='Gene name',
        right_index=True,
    )
    genes_multiple_start = set(
        positions_n_starts['Gene name'][positions_n_starts.n_start != 1]
    )

    gene_names = set(ht.gene_name.collect()) - {None}
    if feature == 'gene':
        feature_col_start = 'Gene start (bp)'
        feature_col_end = 'Gene end (bp)'
    elif feature == 'exon':
        feature_col_start = 'Exon region start (bp)'
        feature_col_end = 'Exon region end (bp)'
    feature_intervals, feature_intervals_hl = defaultdict(set), defaultdict(list)
    feature_breaks = defaultdict(set)
    positions = positions.set_index('Gene name')
    for i, gene_name in enumerate(gene_names):
        if gene_name in genes_multiple_start:
            print(f'WARNING: {gene_name} has multiple start')
        feature_positions = positions.loc[[gene_name]]
        for i, row in feature_positions.iterrows():
            contig = row['Chromosome/scaffold name']
            feature_start = row[feature_col_start]
            feature_end = row[feature_col_end]
            feature_breaks[gene_name].add(feature_start)
            feature_breaks[gene_name].add(feature_end)
            feature_interval = f'chr{contig}:{feature_start}-{feature_end + 1}'
            if feature_interval not in feature_intervals[gene_name]:
                feature_intervals_hl[gene_name].append(hl.parse_locus_interval(feature_interval))
                feature_intervals[gene_name].add(feature_interval)

    feature_breaks[hl.missing(hl.tstr)] = hl.missing(hl.tset(hl.tint))
    feature_intervals_hl[hl.missing(hl.tstr)] = hl.missing(
        hl.tarray(hl.tinterval(hl.tlocus()))
    )

    feature_breaks = hl.dict(feature_breaks)
    feature_intervals_hl = hl.dict(feature_intervals_hl)

    dist_to_feature_expr = hl.if_else(
        feature_intervals_hl[ht.gene_name].any(
            lambda inter: inter.contains(ht.locus)
        ),
        0,
        hl.min(
            feature_breaks[ht.gene_name].map(
                lambda x: hl.abs(ht.locus.position - x)
            )
        )
    )

    return dist_to_feature_expr


def is_sub(pv):
    sub_sa_re = r'^[A-Z0-9]+\*\d+\.\d+$'
    sa_re = r'^[A-Z0-9]+\*\d+$'

    pv_sub_sa = pv.filter(pv.star_allele.matches(sub_sa_re))
    pv_sa = pv.filter(pv.star_allele.matches(sa_re))
    pv2 = pv_sa.annotate(
        sub_in_sa=pv_sub_sa[pv_sa.key].star_allele
    )

    pv2 = pv.filter((pv.gene_name == 'DPYD') != (pv.star_allele.matches(r'^rs')))  # xor
    assert pv2.count() == 0, 'DPYD not only if star allele starts with "rs"'

    return pv.star_allele.matches(r'^rs|' + sub_sa_re)

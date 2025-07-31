# %%
import os
import numpy as np
import pandas as pd
import hail as hl
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score
from sklearn.metrics import roc_curve, auc
from collections import Counter

from analysis.utils.pathogenicity_scores import CADD, SC_TRANSLATE, annotate_vcf_with_scores
from analysis.pharm_ml_score.pmls import pmls_pl2

GENETIC_CODE = {  # GENETIC_CODE_NNK
    'F': ['TTT'],
    'L': ['TTG', 'CTT', 'CTG'],
    'I': ['ATT'],
    'M': ['ATG'],
    'V': ['GTT', 'GTG'],
    'S': ['TCT', 'TCG', 'AGT'],
    'P': ['CCT', 'CCG'],
    'T': ['ACT', 'ACG'],
    'A': ['GCT', 'GCG'],
    'Y': ['TAT'],
    'X': ['TAG'],  # = *
    'H': ['CAT'],
    'Q': ['CAG'],
    'N': ['AAT'],
    'K': ['AAG'],
    'D': ['GAT'],
    'E': ['GAG'],
    'C': ['TGT'],
    'W': ['TGG'],
    'R': ['CGT', 'CGG', 'AGG'],
    'G': ['GGT', 'GGG']
}


ZERO_NEUTRAL = 15

def auc_bootstrap_ci(y_true, y_score, B=1000, alpha=0.05, random_state=None):
    rng = np.random.RandomState(random_state)
    n = len(y_true)
    auc = roc_auc_score(y_true, y_score)
    
    bootstrapped_aucs = []
    for i in range(B):
        idx = rng.randint(0, n, n)
        if len(np.unique(y_true[idx])) < 2:
            continue
        score_b = roc_auc_score(y_true[idx], y_score[idx])
        bootstrapped_aucs.append(score_b)
    
    bootstrapped_aucs = np.array(bootstrapped_aucs)
    lower = np.percentile(bootstrapped_aucs, 100 * (alpha/2))
    upper = np.percentile(bootstrapped_aucs, 100 * (1 - alpha/2))
    
    return auc, lower, upper


def create_template(exon_regions_path):
    template = CADD().read_table()
    def interval_it(path):
        exons = pd.read_csv(path, sep='\t')
        for _, row in exons.iterrows():
            yield row[['Genomic coding start', 'Genomic coding end']]
    interval_dict = hl.literal({
        hl.parse_locus_interval(f'chr10:{start}-{end + 1}'): i + 1
        for i, (start, end) in enumerate(interval_it(exon_regions_path))
    })
    template = hl.filter_intervals(template, interval_dict.keys())
    assert template.count() / 3 / 3 == 491  # aminoacids

    template = template.annotate(
        exon_nr=interval_dict[
            interval_dict.keys().find(
                lambda interval: interval.contains(template.locus)
            )
        ],
        exon_start=interval_dict.keys().find(
            lambda interval: interval.contains(template.locus)
        ).start.position,
        exon_end=interval_dict.keys().find(
            lambda interval: interval.contains(template.locus)
        ).end.position
    )
    template = template.annotate(
        exon_lengths=interval_dict.keys().map(
            lambda interval: interval.end.position - interval.start.position
        )
    )
    template = template.annotate(
        allele_1=template.alleles[0],
        allele_2=template.alleles[1],
        position_seq=(
            template.locus.position
            - template.exon_start
            + hl.sum(template.exon_lengths[:template.exon_nr - 1])
        )
    )
    template = template.transmute(
        position=template.position_seq // 3 + 1,
        pos=template.position_seq % 3,
        aux=hl.struct(
            exon_nr=template.exon_nr,
            exon_start=template.exon_start,
            exon_end=template.exon_end,
            exon_lengths=template.exon_lengths
        )
    )
    return template.drop('score_raw', 'score_phred')


def compute_variant_burden(ht, score, agg_fun=hl.agg.sum):
    ht = ht.annotate(
        **{
            aa: hl.struct(
                **{
                    seq: (hl.literal(seq)[ht.pos] == ht.allele_2) * ht[score]
                    for seq in sequences
                }
            )
            for aa, sequences in GENETIC_CODE.items()
        },
    )
    ht_grouped = ht.group_by(ht.position)
    ht = ht_grouped.aggregate(
        **{
            aa: hl.struct(
                **{
                    seq: agg_fun(ht[aa][seq]) for seq in sequences
                },
                **{
                    f'{seq}_n': hl.agg.sum(hl.literal(seq)[ht.pos] == ht.allele_2)
                    for seq in sequences
                },
                **{
                    f'{seq}_af': hl.agg.sum((hl.literal(seq)[ht.pos] == ht.allele_2) * ht.af)
                    for seq in sequences
                },
            )
            for aa, sequences in GENETIC_CODE.items()
        }
    )

    df = ht.to_pandas()
    fa = df.melt(id_vars=['position'], var_name='end', value_vars=df.filter(regex=r'[A-Z]\.[A-Z]{3}$').columns, value_name=score)
    fn = df.melt(id_vars=['position'], var_name='end', value_vars=df.filter(regex=r'[A-Z]\.[A-Z]{3}_n$').columns, value_name='n')
    faf = df.melt(id_vars=['position'], var_name='end', value_vars=df.filter(regex=r'[A-Z]\.[A-Z]{3}_af$').columns, value_name='af')
    fn['end'] = fn['end'].str.replace('_n', '', regex=True)
    faf['end'] = faf['end'].str.replace('_af', '', regex=True)
    fff = pd.merge(fa, fn, on=['position', 'end'])
    fff = pd.merge(fff, faf, on=['position', 'end'])
    fff[['end', 'codon']] = fff['end'].str.split('.', expand=True)
    return fff


def plot_rocs(variants, aa_scores, score_name, score_tools, n=1, cyp=None,
    title='Score comparison'):

    print(score_name)
    cutoffs = pd.read_csv('results/manuscript-supplementary/pharmvar/cut-offs.csv')
    plt.rcParams.update({'font.size': 14})
    plt.figure(figsize=(9, 9))
    plot_dfs = {}
    for sc in score_tools:
        if (sc_min := min(variants.filter(hl.is_defined(variants[sc]))[sc].collect())) < 0:
            print(f'WARN: We assume that 0 is neutral score (is: {sc_min})')
        agg_fun = hl.agg.sum
        variants_burden = compute_variant_burden(
            variants, score=sc, agg_fun=agg_fun
        )
        df = pd.merge(
            aa_scores, variants_burden, how='left', on=['position', 'end']
        )
        df = df[df.n == n].copy()

        # if more codons for an amino acid substitution
        group_cols = [col for col in df.columns if col not in df.loc[:, sc:'codon'].columns]
        df = (
            df
            .groupby(group_cols, dropna=False)
            .agg(**{
                sc: (sc, lambda x: x.mean()),
                'codons': ('codon', lambda x: set(x.dropna())),
                'max_af': ('af', lambda x: x.max()),
            })
            .reset_index()
        )
        print(sc)

        plot_dfs[sc] = df.copy()

        q1, q3 = 0.8, 0.2
        df_auc = df[(df[score_name] > q1) | (df[score_name] < q3)]
        df_auc = df_auc[~(df_auc[score_name].isna())]
        th = cutoffs[(cutoffs.score == SC_TRANSLATE[sc]) & (cutoffs.gene_name == cyp)].J.values[0] + ZERO_NEUTRAL
        print('acc:', ((df_auc[score_name] < 0.5) == (df_auc[sc] < th)).mean())

        print('n variants ROC:', df_auc.shape)
        aa_y = df_auc[score_name] < 0.5
        print(sc, roc_auc_score(aa_y, df_auc[sc]), end='\n\n')
        _auc, auc_ci_lower, auc_ci_upper = auc_bootstrap_ci(np.array(aa_y), np.array(df_auc[sc]))

        fpr, tpr, thresholds = roc_curve(aa_y, df_auc[sc])
        roc_auc = auc(fpr, tpr)
        idx = np.argmax(tpr - fpr)
        youdens_index = thresholds[idx]
        print(sc, roc_auc_score(aa_y, df_auc[sc]), youdens_index-ZERO_NEUTRAL, end='\n\n')
        plt.plot(
            fpr, tpr,
            label=(
                f"{SC_TRANSLATE.get(sc, sc) + ': '}"
                f"AUC={roc_auc:.3f}"
            ),
            lw =2
        )
    plt.plot([0, 1], [0, 1], color='navy', lw=1.5, linestyle='--')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(title)
    plt.legend(loc="lower right")
    return plot_dfs


# %% annotate
# Downloaded from Biomart and checked with UCSC genome browser
cyp2c9_exon_regions_path = os.path.join(
    'data',
    'cyp2c9-c19-activity-abundance',
    'CYP2C9-genomic-coding-regions.txt'
)
template = create_template(cyp2c9_exon_regions_path)
cyp2c9_variants = template.annotate(
    gene_name='CYP2C9'
)
cyp2c9_variants = annotate_vcf_with_scores(cyp2c9_variants)
cyp2c9_variants = cyp2c9_variants.checkpoint(
    'data/activity-abundance/CYP2C9.ht',
    overwrite=True
)

# Downloaded from Biomart and checked with UCSC genome browser
cyp2c19_exon_regions_path = os.path.join(
    'data',
    'cyp2c9-c19-activity-abundance',
    'CYP2C19-genomic-coding-regions.txt'
)
cyp2c19_variants = create_template(cyp2c19_exon_regions_path)
cyp2c19_variants = cyp2c19_variants.annotate(
    gene_name='CYP2C19'
)
cyp2c19_variants = annotate_vcf_with_scores(cyp2c19_variants)
cyp2c19_variants = cyp2c19_variants.checkpoint(
    'data/activity-abundance/CYP2C19.ht',
    overwrite=True
)


#################### CYP2C9 ####################
# %% CYP2C9 burden scores
cyp = 'CYP2C9'
ht = hl.read_table(f'data/activity-abundance/{cyp}.ht')
dataset = ht.transmute(
    cadd_nn=ht.scores.cadd_norm,
    fathmm_xf_nn=ht.scores.fathmm_xf_norm,
    pair0_nn=ht.scores.sift_norm,
    pair1_nn=ht.scores.phylop100_norm,
)
dataset_df = dataset.to_pandas()
dataset_df['locus'] = dataset_df.locus.astype(str)
dataset_df['alleles'] = dataset_df.alleles.apply(lambda x: tuple(x))
out = pmls_pl2(
    dataset_df,
    id_columns=['gene_name', 'locus', 'alleles'],
    models_dir_prefix='data/pmls-models/sift-phylop100-',
    genes=set(dataset_df.gene_name)
)
out[['a1', 'a2']] = pd.DataFrame(out.alleles.to_list())
out = out.drop('alleles', axis=1)
out_hl = hl.Table.from_pandas(out)
out_hl = out_hl.annotate(
    locus=hl.parse_locus(out_hl.locus),
    alleles=hl.array([out_hl.a1, out_hl.a2]),
)
out_hl = out_hl.key_by('locus', 'alleles')
out_hl = out_hl.select('pmls_pred')

ht = hl.read_table(f'data/activity-abundance/{cyp}.ht')
cyp_preds = ht.transmute(
    af=ht.scores.af,
    pgs=ht.scores.pgs,
    adme=ht.scores.adme,
    adme2=ht.scores.adme2,
    cadd=ht.scores.cadd,
    fathmm_xf=ht.scores.fathmm_xf,
    ma=ht.scores.ma,
    sift=ht.scores.sift,
    phylop100=ht.scores.phylop100,
    provean=ht.scores.provean,

    pmls_pred=out_hl[ht.key].pmls_pred
)
cyp_preds = cyp_preds.persist()


# %% combine annotations with activity and abundance scores and plot

cyp_aa_scores = pd.read_csv(
    f'data/cyp2c9-c19-activity-abundance/{cyp}_activity_abundance_scores.csv'
)
assert len(cyp_aa_scores[cyp_aa_scores['class'] == 'wt']) == 1
cyp_aa_scores = cyp_aa_scores[cyp_aa_scores['class'] != 'wt']

print(cyp_aa_scores.shape)
print(Counter(cyp_aa_scores['class']))
print(cyp_aa_scores[cyp_aa_scores.activity_score.notna()].shape)
print(cyp_aa_scores[cyp_aa_scores.abundance_score.notna()].shape)

cyp_preds2 = cyp_preds.annotate(
    pgs=cyp_preds.pgs + ZERO_NEUTRAL,
    cadd=cyp_preds.cadd + ZERO_NEUTRAL,
    pmls_pred=cyp_preds.pmls_pred + ZERO_NEUTRAL,
    adme=cyp_preds.adme + ZERO_NEUTRAL,
    adme2=cyp_preds.adme2 + ZERO_NEUTRAL,
    phylop100=cyp_preds.phylop100 + ZERO_NEUTRAL,
    ma=cyp_preds.ma + ZERO_NEUTRAL,
    fathmm_xf=cyp_preds.fathmm_xf + ZERO_NEUTRAL,
    sift=cyp_preds.sift + ZERO_NEUTRAL,
    provean=cyp_preds.provean + ZERO_NEUTRAL,
)

cyp_burdens = plot_rocs(
    variants=cyp_preds2,
    aa_scores=cyp_aa_scores,
    score_name='abundance_score',
    score_tools=[
        'pgs',
        'cadd',
        'pmls_pred',
        'adme2',
        'adme',
        'phylop100',
        'ma',
        'sift',
        'fathmm_xf',
        'provean',
    ],
    n=1,
    cyp=cyp,
    title=f'ROC curves - {cyp} protein abundance',
) 

cyp_burdens_c9 = cyp_burdens
plt.savefig(f'results/manuscript-pgs/activity-abundance/abundance-{cyp}.svg')
plt.close()

cyp_burdens = plot_rocs(
    variants=cyp_preds2,
    aa_scores=cyp_aa_scores,
    score_name='activity_score',
    score_tools=[
        'pgs',
        'pmls_pred',
        'adme2',
        'adme',
        'cadd',
        'ma',
        'sift',
        'phylop100',
        'fathmm_xf',
        'provean',
    ],
    n=1,
    cyp=cyp,
    title=f'ROC curves - {cyp} enzyme activity',
)
plt.savefig(f'results/manuscript-pgs/activity-abundance/activity-{cyp}.svg')
plt.close()


#################### CYP2C19 ####################
# %% CYP2C19 ROCs & Abundance
cyp = 'CYP2C19'
ht = hl.read_table(f'data/activity-abundance/{cyp}.ht')
dataset = ht.transmute(
    cadd_nn=ht.scores.cadd_norm,
    fathmm_xf_nn=ht.scores.fathmm_xf_norm,
    pair0_nn=ht.scores.sift_norm,
    pair1_nn=ht.scores.phylop100_norm,
)
####
# adding score
dataset_df = dataset.to_pandas()
dataset_df['locus'] = dataset_df.locus.astype(str)
dataset_df['alleles'] = dataset_df.alleles.apply(lambda x: tuple(x))
out = pmls_pl2(
    dataset_df,
    id_columns=['gene_name', 'locus', 'alleles'],
    models_dir_prefix='data/pmls-models/sift-phylop100-',
    genes=set(dataset_df.gene_name)
    )
out[['a1', 'a2']] = pd.DataFrame(out.alleles.to_list())
out = out.drop('alleles', axis=1)
out_hl = hl.Table.from_pandas(out)
out_hl = out_hl.annotate(
    locus=hl.parse_locus(out_hl.locus),
    alleles=hl.array([out_hl.a1, out_hl.a2]),
)
out_hl = out_hl.key_by('locus', 'alleles')
out_hl = out_hl.select('pmls_pred')

# adding score
ht = hl.read_table(f'data/activity-abundance/{cyp}.ht')
cyp_preds = ht.transmute(
    af=ht.scores.af,
    pgs=ht.scores.pgs,
    adme=ht.scores.adme,
    adme2=ht.scores.adme2,
    cadd=ht.scores.cadd,
    fathmm_xf=ht.scores.fathmm_xf,
    ma=ht.scores.ma,
    sift=ht.scores.sift,
    phylop100=ht.scores.phylop100,
    provean=ht.scores.provean,
    pmls_pred=out_hl[ht.key].pmls_pred
)
cyp_preds = cyp_preds.persist()

cyp_aa_scores = pd.read_csv(
    f'data/cyp2c9-c19-activity-abundance/{cyp}_data.csv'
)


# %%
cyp_preds2 = cyp_preds.annotate(
    pgs=cyp_preds.pgs + ZERO_NEUTRAL,
    cadd=cyp_preds.cadd + ZERO_NEUTRAL,
    pmls_pred=cyp_preds.pmls_pred + ZERO_NEUTRAL,
    adme2=cyp_preds.adme2 + ZERO_NEUTRAL,
    adme=cyp_preds.adme + ZERO_NEUTRAL,
    phylop100=cyp_preds.phylop100 + ZERO_NEUTRAL,
    ma=cyp_preds.ma + ZERO_NEUTRAL,
    fathmm_xf=cyp_preds.fathmm_xf + ZERO_NEUTRAL,
    sift=cyp_preds.sift + ZERO_NEUTRAL,
    provean=cyp_preds.provean + ZERO_NEUTRAL,
)

j_c19={'abundance': {}, 'activity': {}}
cyp_burdens = plot_rocs(
    variants=cyp_preds2,
    aa_scores=cyp_aa_scores,
    score_name='abundance_score',
    score_tools=[
        'pgs',
        'cadd',
        'adme2',
        'pmls_pred',
        'adme',
        'phylop100',
        'fathmm_xf',
        'ma',
        'sift',
        'provean'
    ],
    n=1,
    cyp=cyp,
    title=f'ROC curves - {cyp} protein abundance',
)
cyp_burdens_c19 = cyp_burdens
plt.savefig(f'results/manuscript-pgs/activity-abundance/abundance-{cyp}.svg', bbox_inches='tight')
plt.close()

"""Odds ratio"""
# %%
import hail as hl
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
from scipy.stats.contingency import odds_ratio
from scipy.stats import chi2_contingency


def odds_ratio_stats(tb):
    se = np.sqrt((1 / tb.flatten()).sum())
    l = np.log(odds_ratio(tb).statistic)
    return (l, se), (np.exp(l), np.exp(l - 1.96 * se), np.exp(l + 1.96 * se))


# %% load burden
mt = hl.read_matrix_table('burden-annot.mt')
burden = mt.entries()

# %%
genes = ['CYP2C19']
ht = burden.filter(hl.is_defined(burden.number_of_therapies) & hl.literal(genes).contains(burden.gene_name))

group_ht = np.array(ht.number_of_therapies.collect())
pgs_ht = np.array(ht.pgs.collect())
poly_ht = np.array(ht.poly.collect())
or_scores = [
    ('PharmGScore', 10.468, np.array(ht.pgs.collect())),
    ('PharmMLScore', -0.980, np.array(ht.pmls_pred.collect())),
    ('APF', 0.200, np.array(ht.burden_adme.collect())),
    ('APF2', 0.250, np.array(ht.burden_adme2.collect())),
    ('CADD', -0.985, np.array(ht.burden_cadd.collect())),
    ('FATHMM-XF', 0.250, np.array(ht.burden_fathmm_xf.collect())),
    ('MutationAssessor', 1.500, np.array(ht.burden_ma.collect())),
    ('Phylop100', -0.793, np.array(ht.burden_phylop100.collect())),
    ('PROVEAN', 1.970, np.array(ht.burden_provean.collect())),
    ('SIFT', -0.418, np.array(ht.burden_sift.collect()))
]
g_th = 2
tbs = [
    np.array([
        [
            nonfun_case := np.sum(((poly_ht == 'No') | (poly_ht == 'Decreased')) & (group_ht > g_th)),
            nonfun_noncase := np.sum(((poly_ht == 'No') | (poly_ht == 'Decreased')) & (group_ht <= g_th))
        ],
        [
            fun_case := np.sum(((poly_ht == 'Normal')) & (group_ht > g_th)),
            fun_noncase := np.sum(((poly_ht == 'Normal')) & (group_ht <= g_th))
        ]
    ])
]
tbs = tbs + [
    np.array([
        [
            nonfun_case := np.sum((score_ht >= t_score) & (group_ht > g_th)),
            nonfun_noncase := np.sum((score_ht >= t_score) & (group_ht <= g_th))
        ],
        [
            fun_case := np.sum((score_ht < t_score) & (group_ht > g_th)),
            fun_noncase := np.sum((score_ht < t_score) & (group_ht <= g_th))
        ]
    ])
    for _, t_score, score_ht in or_scores
]

# %% draw dot plot with errors
lcc = [odds_ratio_stats(tb)[1] for tb in tbs]
means =  [l[0] for l in lcc]
stds = [[-(l[1] - l[0]), l[2] - l[0]] for l in lcc]
group_labels = ['Adverse drug reaction', 'Number of therapies']

# order
fixed_idx = [0, 1, 2]
new_order = fixed_idx + list(range(3, len(means)))
categories = list(range(len(new_order)))
y_vals = [means[i] for i in new_order]
y_errs = [stds[i] for i in new_order]
lower_err = [err[0] for err in y_errs]
upper_err = [err[1] for err in y_errs]
errors = [lower_err, upper_err]
#
positions = [0, 1, 2] + list(range(4, 4 + (len(new_order) - 3)))
positions = list(range(len(new_order)))

fig, ax = plt.subplots(figsize=(10, 7))
plt.title(f'Odds ratios of genetic predictors for two therapy switches', fontsize=18, pad=20)
ax.errorbar(
    positions, y_vals, yerr=errors, fmt='o', capsize=5, capthick=1, markersize=8,
    label=f'Number of therapies > {g_th} (with 95% CI)'
)
p_values = []
for tb, (name, *_) in zip(tbs, [('sa', None, None)] + or_scores):
    chi2, p, dof, expected = chi2_contingency(tb, correction=False)
    p_values.append(p)

for x, y, err_up, p in zip(positions, y_vals, errors[1], p_values):
    if p < 0.05:
        ax.text(x, y + err_up + 0.08, '*', ha='center', va='bottom', fontsize=16)
        ax.text(x, y + err_up + 0.03, f'({p:.3f})', ha='center', va='bottom', fontsize=10)

ax.set_ylabel('Odds ratio', fontsize=16, labelpad=20)
ax.set_yticks([i / 10 for i in range(8, 16)])
plt.yticks(fontsize=14)
ax.set_xticks(positions)
xlabels = ['Star alleles', 'PharmGScore', 'PharmMLScore'] + [or_scores[i-1][0] for i in new_order[3:]]
print(or_df := pd.DataFrame({
    'Score': xlabels,
    'or': y_vals,
    'lower': (y - l for y, l in zip(y_vals, lower_err)),
    'upper': (y + u for y, u in zip(y_vals, upper_err)),
    'pval': p_values,
}))

or_df['Odds ratio'] = or_df['or'].round(2)
or_df['p-value'] = or_df['pval'].round(2)
or_df['95% confidence interval'] = (
    or_df['lower'].round(2).astype(str)
    + ' - '
    + or_df['upper'].round(2).astype(str)
)
or_df[['Score', 'Odds ratio', '95% confidence interval', 'p-value']]

ax.set_xticklabels(xlabels, rotation=60, fontsize=16, ha='right')
ax.yaxis.grid(True)
ax.axhline(y=1, color='red', linewidth=2, linestyle='--', label='Statistical significance')
ax.legend(loc='upper right', fontsize=14)
plt.tight_layout()
plt.savefig(f'results/manuscript-pgs/ukb200/odds-ratios-n-therapies.svg')
plt.close()


# %%
genes = ['CYP2C19']
mode = 'all'  # 'accidental'

ht = burden.filter(hl.is_defined(burden.group) & hl.literal(genes).contains(burden.gene_name))
if mode == 'all':
    ht = ht.annotate(
        group=hl.if_else(ht.group == 'adr', 1, 0),
    )
elif mode == 'accidental':
    ht = ht.annotate(
        group=hl.if_else(
            (hl.literal(['accidental', 'therapeutic']).contains(ht.intention)),
            1,
            0
        ),
    )
else:
    raise ValueError

ht.aggregate(hl.agg.counter(ht.group))

group_ht = np.array(ht.group.collect())
pgs_ht = np.array(ht.pgs.collect())
poly_ht = np.array(ht.poly.collect())
or_scores = [
    ('PharmGScore', 10.468, np.array(ht.pgs.collect())),
    ('PharmMLScore', -0.980, np.array(ht.pmls_pred.collect())),
    ('APF', 0.200, np.array(ht.burden_adme.collect())),
    ('APF2', 0.250, np.array(ht.burden_adme2.collect())),
    ('CADD', -0.985, np.array(ht.burden_cadd.collect())),
    ('FATHMM-XF', 0.250, np.array(ht.burden_fathmm_xf.collect())),
    ('MutationAssessor', 1.500, np.array(ht.burden_ma.collect())),
    ('Phylop100', -0.793, np.array(ht.burden_phylop100.collect())),
    ('PROVEAN', 1.970, np.array(ht.burden_provean.collect())),
    ('SIFT', -0.418, np.array(ht.burden_sift.collect()))
]
g_th = 1
tbs = [
    np.array([
        [
            nonfun_case := np.sum(((poly_ht == 'No') | (poly_ht == 'Decreased')) & (group_ht >= g_th)),
            nonfun_noncase := np.sum(((poly_ht == 'No') | (poly_ht == 'Decreased')) & (group_ht < g_th))
        ],
        [
            fun_case := np.sum(((poly_ht == 'Normal')) & (group_ht >= g_th)),
            fun_noncase := np.sum(((poly_ht == 'Normal')) & (group_ht < g_th))
        ]
    ])
]
tbs = tbs + [
    np.array([
        [
            nonfun_case := np.sum((score_ht >= t_score) & (group_ht >= g_th)),
            nonfun_noncase := np.sum((score_ht >= t_score) & (group_ht < g_th))
        ],
        [
            fun_case := np.sum((score_ht < t_score) & (group_ht >= g_th)),
            fun_noncase := np.sum((score_ht < t_score) & (group_ht < g_th))
        ]
    ])
    for _, t_score, score_ht in or_scores
]

# %%

lcc = [odds_ratio_stats(tb)[1] for tb in tbs]
means =  [l[0] for l in lcc]
stds = [[-(l[1] - l[0]), l[2] - l[0]] for l in lcc]
group_labels = ['Adverse drug reaction', 'Number of therapies']

categories = list(range(len(means)))

y_vals = [means[cat_idx] for cat_idx in range(len(categories))]
y_errs = [stds[cat_idx] for cat_idx in range(len(categories))]
lower_err = [err[0] for err in y_errs]
upper_err = [err[1] for err in y_errs]
errors = [lower_err, upper_err]

# order
fixed_idx = [0, 1, 2]
new_order = fixed_idx + list(range(3, len(means)))
positions = [0, 1, 2] + list(range(4, 4 + (len(new_order) - 3)))
positions = list(range(len(new_order)))
#


fig, ax = plt.subplots(figsize=(10, 7))
plt.title(f'Odds ratios of genetic predictors for adverse drug reactions', fontsize=18, pad=20)
ax.errorbar(
    positions, y_vals, yerr=errors, fmt='o', capsize=5, capthick=1, markersize=8,c='#ff7f0e',
    label='Adverse drug reaction (with 95% CI)'
)
p_values = []
for tb, (name, *_) in zip(tbs, [('sa', None, None)] + or_scores):
    chi2, p, dof, expected = chi2_contingency(tb, correction=False)
    p_values.append(p)

for x, y, err_up, p in zip(positions, y_vals, errors[1], p_values):
    if p < 0.05:
        ax.text(x, y + err_up + 0.1, '*', ha='center', va='bottom', fontsize=16)
        ax.text(x, y + err_up + 0.03, f'({p:.3f})', ha='center', va='bottom', fontsize=10)

ax.set_ylabel('Odds ratio', fontsize=16, labelpad=20)
if mode == 'all':
    ax.set_yticks([i / 10 for i in range(4, 24, 2)])
elif mode == 'accidental':
    ax.set_yticks([i / 10 for i in range(0, 56, 5)])
plt.yticks(fontsize=14)
ax.set_xticks(positions)

xlabels = ['Star alleles', 'PharmGScore', 'PharmMLScore'] + [or_scores[i-1][0] for i in new_order[3:]]
print(or_df := pd.DataFrame({
    'Score': xlabels,
    'or': y_vals,
    'lower': (y - l for y, l in zip(y_vals, lower_err)),
    'upper': (y + u for y, u in zip(y_vals, upper_err)),
    'pval': p_values,
}))
or_df['Odds ratio'] = or_df['or'].round(2)
or_df['p-value'] = or_df['pval'].round(2)
or_df['95% confidence interval'] = (
    or_df['lower'].round(2).astype(str)
    + ' - '
    + or_df['upper'].round(2).astype(str)
)
or_df[['Score', 'Odds ratio', '95% confidence interval', 'p-value']]
ax.set_xticklabels(xlabels, rotation=60, fontsize=16, ha='right')
ax.yaxis.grid(True)
ax.axhline(y=1, color='red', linewidth=2, linestyle='--', label='Statistical significance')
ax.legend(loc='upper right', fontsize=14)
plt.tight_layout()
if mode == 'all':
    plt.savefig(f'results/manuscript-pgs/ukb200/odds-ratios-adr.svg')
elif mode == 'accidental':
    plt.savefig(f'results/manuscript-supplementary/ukb200/odds-ratios-adr-accidental.svg')
plt.close()

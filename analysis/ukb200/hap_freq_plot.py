"""Haplotype, phenotype frequencies"""
# %%
import hail as hl
import pandas as pd

from analysis.utils.ukb_const import PHENOTYPE_MAP


# %% read polygenic resutls
poly = pd.read_csv('data/polygenic/polygenic-results.csv', low_memory=False)


# %% read pgxpop results
pop = pd.read_csv('data/pgxpop/pgxpop-results.csv', low_memory=False)
mapping = {
    # SLCO1B1
    '*1A': '*1',
    '*1B': '*37'
}
def hap_to_set(df_col):
    return (
        df_col
        .str.split(';', n=1).str[0]
        .str.findall(r'\*\w+')
        .map(set)
        .map(lambda s: {mapping.get(allele, allele) for allele in s})
    )

pop['hap_1_set'] = hap_to_set(pop['hap_1'])
pop['hap_2_set'] = hap_to_set(pop['hap_2'])

# %% compare haplotypes
poly_h = poly[poly.gene.isin(['CYP2B6', 'CYP2C19', 'CYP2C9', 'CYP2D6', 'CYP3A5', 'CYP4F2', 'NUDT15', 'SLCO1B1'])]
pop_h = pop[pop.gene.isin(['CYP2B6', 'CYP2C19', 'CYP2C9', 'CYP2D6', 'CYP3A5', 'CYP4F2', 'NUDT15', 'SLCO1B1'])]

pp_hap = poly_h.merge(pop_h, on=['eid', 'gene'])
pp_hap[['eid', 'gene', 'allele1_name', 'allele2_name', 'hap_1_set', 'hap_2_set']]  # poly, pop

same_haplotypes = (
    (
        pp_hap.apply(lambda row: row['allele1_name'] in row['hap_1_set'], axis=1) 
        & pp_hap.apply(lambda row: row['allele2_name'] in row['hap_2_set'], axis=1) 
    )
    | (
        pp_hap.apply(lambda row: row['allele2_name'] in row['hap_1_set'], axis=1) 
        & pp_hap.apply(lambda row: row['allele1_name'] in row['hap_2_set'], axis=1) 
    )
)
pp_hap[same_haplotypes][['eid', 'gene', 'allele1_name', 'allele2_name', 'hap_1_set', 'hap_2_set']]  # poly, pop
pp_hap[~same_haplotypes][['eid', 'gene', 'allele1_name', 'allele2_name', 'hap_1_set', 'hap_2_set']]  # poly, pop

pp_hap[same_haplotypes].shape[0] / pp_hap.shape[0]


# %% compare phenotypes

# CYP2B6, CYP2C19, CYP2C9, CYP2D6, CYP3A5, CYP4F2 (no recommendations), DPYD, NUDT15, SLCO1B1
# brak DPYD z powodu innego kodowania:
poly_p = poly[poly.gene.isin(['CYP2B6', 'CYP2C19', 'CYP2C9', 'CYP2D6', 'CYP3A5', 'NUDT15', 'SLCO1B1'])]
pop_p = pop[pop.gene.isin(['CYP2B6', 'CYP2C19', 'CYP2C9', 'CYP2D6', 'CYP3A5', 'NUDT15', 'SLCO1B1'])]

pp_p = poly_p.merge(pop_p, on=['eid', 'gene'])
pp_p['poly_pheno'] = pp_p['poly_phenotype'].map(PHENOTYPE_MAP)
pp_p['pop_pheno'] = pp_p['pop_phenotype'].map(PHENOTYPE_MAP)
pp_p[['eid', 'gene', 'poly_pheno', 'pop_pheno']]  # poly, pop

same_phenotypes = pp_p['poly_pheno'] == pp_p['pop_pheno']
pp_p[same_phenotypes][['eid', 'gene', 'poly_pheno', 'pop_pheno']]
pp_p[same_phenotypes].shape[0] / pp_p.shape[0]

# %% merge - histogram frequencies

df = poly.merge(pop, on=['eid', 'gene'])
df['pop_label2'] = df.apply(lambda row: '/'.join(row['hap_1_set'] | row['hap_2_set']), axis=1)

th_freq = 0.001  # 0.5%
colors = ['#d496a7', '#2b061e', '#725752']
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import PercentFormatter

for gene in genes:
    print(f"Processing {gene}")
    df_gene = df[df.gene == gene]

    freq_dict = {}
    for name, col in [('Polygenic', 'poly_label'),
                      ('PGxPOP',   'pop_label2')]:
        alleles = (
            df_gene[col]
              .str.split('/', expand=True)
              .stack()
        )
        freq_dict[name] = alleles.value_counts(normalize=True)

    freq_dict['PharmGKB'] = gkb[gkb.gene_name == gene].set_index('allele')['european']

    alleles_to_plot = set()
    for freq in freq_dict.values():
        alleles_to_plot |= set(freq[freq >= th_freq].index)

    freq_df = pd.DataFrame({
        name: (freq.reindex(alleles_to_plot, fill_value=0) * 100)
        for name, freq in freq_dict.items()
    }).sort_index(key=lambda x: x.str.strip('*').astype(int))
    print(freq_df)

    plt.rcParams.update({'font.size': 14})
    ax = freq_df.plot(
        kind='bar',
        figsize=(10,6),
        width=0.8,
        color=colors
    )
    ax.yaxis.set_major_locator(MultipleLocator(5))
    ax.grid(axis='y', color='lightgray', linestyle='--', linewidth=0.5)
    ax.set_axisbelow(True)

    ax.yaxis.set_major_formatter(PercentFormatter(xmax=100, decimals=0))
    ax.set_xlabel('Star Allele')
    ax.set_ylabel('Frequency')
    ax.set_title(f'Star Allele Frequency in {gene} (≥ {th_freq*100:.1f}% in any source)')
    ax.legend(title='Source')
    plt.tight_layout()

    # save & show
    outpath = f'results/{gene}-star-allele-frequencies.svg'
    plt.savefig(outpath)
    plt.show()


# %
import plotly.graph_objects as go
import plotly.io as pio
pio.renderers.default = "notebook"

df = poly.merge(pop, on=['eid', 'gene'])
df['poly_pheno'] = df['poly_phenotype'].map(PHENOTYPE_MAP)
df['pop_pheno'] = df['pop_phenotype'].map(PHENOTYPE_MAP)


palette = {
    'No (polygenic)': '#b27069',
    'No (PGxPOP)': '#b27069',
    'Decreased (polygenic)': '#e1bb80',
    'Decreased (PGxPOP)': '#e1bb80',
    'Normal (polygenic)': '#e1edbd',
    'Normal (PGxPOP)': '#e1edbd',
    'Increased (polygenic)': '#5b7d8d',
    'Increased (PGxPOP)': '#5b7d8d',
}

# fallback colour if a label isn't in your palette:
default_color = 'lightgray'

for gene in genes:
    print(gene)
    sankey = df[df.gene == gene][['poly_pheno', 'pop_pheno']]
    sankey = sankey.groupby(sankey.columns.tolist(), dropna=False).size().reset_index(name='value').fillna('Undefined')
    source_list = sankey['poly_pheno'].to_list()
    target_list = sankey['pop_pheno'].to_list()

    source_labels = [f"{label} (polygenic)" for label in source_list]
    target_labels = [f"{label} (PGxPOP)" for label in target_list]
    labels = list(set(source_labels + target_labels))
    labels = [
    'Increased (polygenic)',
    'Normal (polygenic)',
    'Decreased (polygenic)',
    'No (polygenic)',
    'Undefined (polygenic)',
    'Increased (PGxPOP)',
    'Normal (PGxPOP)',
    'Decreased (PGxPOP)',
    'No (PGxPOP)',
    'Undefined (PGxPOP)',
    ]
    source_indices = [labels.index(src) for src in source_labels]
    target_indices = [labels.index(tgt) for tgt in target_labels]
    values = sankey['value'].to_list()

    # Calculate incoming and outgoing sums for each node
    node_values = pd.Series(0, index=range(len(labels)))
    for i, value in enumerate(values):
        node_values[source_indices[i]] += value
        node_values[target_indices[i]] += value

    node_colors = [
        palette.get(lbl, default_color)
        for lbl in labels
    ]
    # Append counts to the labels
    updated_labels = [f"{label} ({node_values[i]})" for i, label in enumerate(labels)]

    import plotly.graph_objects as go
    fig = go.Figure(data=[go.Sankey(
        node=dict(
            pad=20,
            thickness=50,
            line=dict(color="black", width=0.0),
            label=updated_labels,
            color=node_colors,
        ),
        link=dict(
            source=source_indices,
            target=target_indices,
            value=sankey['value'].to_list(),
        )
    )])
    fig.update_layout(
        title={
            'text': f"Phenotype concordance for {gene}",
            'font': dict(size=16),
        },
        font=dict(size=14)
    )
    fig.write_image(f'results/{gene}-sankey.svg')
    fig.show()

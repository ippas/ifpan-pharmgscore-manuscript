# %%
import os
import json
from random import randint
from itertools import combinations
import numpy as np
import pandas as pd
import hail as hl
import copy
import torch
import torch.nn as nn
import torch.optim as optim
from sklearn.metrics import accuracy_score, roc_auc_score

from analysis.utils.pharmvar import PharmVar

from analysis.pharm_ml_score.pmls import SAFunction
from analysis.utils.pathogenicity_scores import annotate_vcf_with_scores

pv_version = '6.0.10'

df = pd.read_csv('data/dbnsfp_select_scores.tsv', sep='\t')
df = df[(df.ratio > 0.9) & (df.increased_normal_decreased_no > 0.75)]
print(len(scs := set(df.score_name2)), scs)

# %% training
pv_data = PharmVar(pv_version).read_table_filtered()

pv = PharmVar(pv_version).read_table()
pv = pv.transmute(gene_name=pv.gene)
pv = annotate_vcf_with_scores(pv).persist()
pv_data_all_sub = pv.filter(pv.is_sub).persist()

test_scores = [
    'dann',
    'bayes_noaf',
    'sift4g',
    'sift',
    'fathmm_mkl',
    'mt',
    'vest4',
    'lrt',
    'bayes_addaf',
    'revel',
    'eigen',
    'eigen_pc',
    'provean',
    'phylop100',

]
nsfp = dbNSFP.read_table(single=True)
combins = list(combinations(test_scores, 2))

# %%
for pair in combins:
    k_pair = '-'.join(pair)
    print(f'Pair: {k_pair}')

    res_path = f'data/pmls-results/{k_pair}.json'
    res = dict()

    pv_data = PharmVar(pv_version).read_table_filtered()
    scores = [
        'cadd_norm',
        'fathmm_xf_norm',
        'pair0_norm',
        'pair1_norm',
    ]
    pair0, pair1 = [p + '_norm' for p in pair]
    dataset = pv_data.transmute(
        cadd_norm=pv_data.scores.cadd_norm,
        fathmm_xf_norm=pv_data.scores.fathmm_xf_norm,
        pair0_norm=pv_data.scores[pair0],
        pair1_norm=pv_data.scores[pair1],
        dist_to_exon=pv_data.dist_to_exon,
    )
    dataset = dataset.select(
        'gene_name', 'star_allele', 'function', 'is_sub', 'dist_to_exon',
        *scores
    )
    dataset_df = dataset.to_pandas()

    code_dict = {
        'increased function': 0,
        'normal function': 0,
        'decreased function': 1,
        'no function': 1,
    }
    dataset_df['y'] = dataset_df['function'].apply(lambda x: code_dict.get(x))
    code_dict = {
        'increased function': 0.5,
        'normal function': 1,
        'decreased function': 0.5,
        'no function': 1,
    }
    dataset_df['loss_weight'] = dataset_df['function'].apply(lambda x: code_dict.get(x))
    dataset_df = dataset_df[~dataset_df.y.isna()]
    dataset_df[scores] = dataset_df[scores].apply(lambda row: row.fillna(row.mean()), axis=1)
    
    dataset_df['n'] = None
    dataset_grouped = dataset_df.groupby(['gene_name', 'star_allele', 'y', 'loss_weight'])
    dataset_df = dataset_grouped.agg({'n': lambda x: x.size - 1} | {k: 'sum' for k in scores})
    dataset_df = dataset_df.reset_index()

    summary_accs = []
    nn_state_dicts = {}

    n_epoch = 25_000
    OO = pd.DataFrame()
    genes = dataset_df.groupby('gene_name').agg(n=('gene_name', 'size')).index.to_list()
    for gene in genes:
        print('-----------------', gene, '-----------------')
        gene_filter = dataset_df['gene_name'] == gene
        train_df = dataset_df[~gene_filter]
        test_df = dataset_df[gene_filter]
        model_path =  f'data/pmls-models/{k_pair}-{gene}.pt'

        X = torch.tensor(train_df[['n', *scores]].values.astype(float), dtype=torch.float32)
        Y = torch.tensor(train_df['y'].values.astype(float)[:, None], dtype=torch.float32)
        loss_weight = torch.tensor(train_df['loss_weight'].values.astype(float)[:, None], dtype=torch.float32)
        train_y = Y.numpy().flatten()

        nn_state_dicts[gene] = {}
        net = SAFunction()
        criterion = nn.BCEWithLogitsLoss(weight=loss_weight)
        optimizer = optim.Adam(net.parameters(), lr=1e-4, weight_decay=0.0)

        train_losses = []
        train_accs = []
        min_loss = np.inf
        for epoch in range(1, n_epoch + 1):
            output = net(X)
            loss = criterion(output, Y)
            loss.backward()
            optimizer.step()
            train_losses.append(loss.item())

            train_acc = accuracy_score(train_y, output > 0.5)
            train_accs.append(train_acc)
            if train_losses[-1] < min_loss:
                min_loss = train_losses[-1]
                best_epoch = epoch
                nn_state_dicts[gene]['best'] = copy.deepcopy(net.state_dict())

        torch.save(nn_state_dicts[gene]['best'], model_path)

        test_input = torch.tensor(test_df[['n', *scores]].values.astype(float), dtype=torch.float32)
        test_output = torch.tensor(test_df['y'].values.astype(float)[:, None], dtype=torch.float32)
        test_y = test_output.numpy().flatten()

        net = SAFunction()
        net.load_state_dict(torch.load(model_path))
        net.eval()

        with torch.no_grad():
            test_pred = net(test_input).numpy().flatten()
            train_pred = net(X).numpy().flatten()
        oo2 = pd.DataFrame({
            'test_pred': test_pred,
            'test_y': test_y,
            'star_allele': test_df.star_allele
        })

        OO = pd.concat([OO, oo2])
        train_acc = accuracy_score(train_y, train_pred > 0.5)
        if test_y.shape[0]:
            test_acc = accuracy_score(test_y, test_pred > 0.5)
            test_auc =  roc_auc_score(test_y, test_pred)
        else:
            test_acc, test_auc = -1, -1
        summary_acc = f'{gene} (train/test/test auc):\t{train_acc:.5f}\t{test_acc:.5f}\t{test_auc:.5f}\n'
        summary_accs.append(summary_acc)

    print(*summary_accs, sep='', flush=True)
    c_dict = {
        'increased function': 0,
        'normal function': 0,
        'decreased function': 1,
        'no function': 1,
    }
    all_sub_df = pv_data_all_sub.key_by('star_allele', 'function').select().distinct().to_pandas()
    all_sub_df['test_y'] = all_sub_df['function'].apply(lambda x: c_dict.get(x))
    all_sub_df = all_sub_df[~all_sub_df.test_y.isnull()]
    all_sub_df['test_pred'] = -99
    all_sub_df = all_sub_df[['star_allele', 'test_y', 'test_pred']]
    OO = pd.concat([OO, all_sub_df]).drop_duplicates(subset='star_allele', keep='first')

    res[k_pair] = roc_auc_score(OO.test_y, OO.test_pred), accuracy_score(OO.test_y, OO.test_pred > 0)
    print(res)

    with open(res_path, 'w') as f:
        json.dump(res, f)

# %%

res_out = dict()
out_dir = 'data/pmls-results/'
for file in os.listdir(out_dir):
    with open(out_dir + '/' + file, 'r') as f:
        res = json.load(f)
        res_out.update(res)

with open('ml-results.json', 'w') as f:
    json.dump(res_out, f)

print(f"{'Pair':25}{'AUC':10}{'Accuracy'}")
for n, s in sorted(res_out.items(), key=lambda x: x[1][0], reverse=True):
    print(f'{n:25}{s[0]:<10.5f}{s[1]:<10.5f}')

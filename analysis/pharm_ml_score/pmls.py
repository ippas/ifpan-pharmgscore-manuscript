import pandas as pd
import polars as pl
import numpy as np
import torch
import torch.nn as nn


class SAFunction(nn.Module):
    def __init__(self):
        super(SAFunction, self).__init__()
        self.sigmoid = nn.Sigmoid()
        self.fc1 = nn.Linear(5, 4)
        self.fc2 = nn.Linear(4, 4)
        self.fco = nn.Linear(4, 1)

    def forward(self, x):
        x = self.fc1(x)
        x = self.sigmoid(x)
        x = self.fc2(x)
        x = self.sigmoid(x)
        return self.fco(x)


def fill_na_with_row_mean(df, scores):
    arr = df[scores].to_numpy(na_value=None).astype(float)
    row_means = np.nanmean(arr, axis=1)
    inds = np.where(np.isnan(arr))
    arr[inds] = np.take(row_means, inds[0])
    df[scores] = arr
    return df


def pmls_pl2(df, id_columns, models_dir_prefix, genes):
    scores = ['cadd_nn', 'fathmm_xf_nn', 'pair0_nn', 'pair1_nn']

    df = fill_na_with_row_mean(df, scores)
    df = pl.from_pandas(df)
    df = df.with_columns(pl.lit(None).alias('n'))
    agg_exprs = [pl.count().alias('n')] + [pl.col(k).sum().alias(k) for k in scores]
    df = df.group_by(id_columns).agg(agg_exprs)
    df = df.with_columns((pl.col('n') - 1).alias('n'))
    df = df.to_pandas()

    preds_df = []
    print('predicting...')
    for gene in genes:
        print(gene)
        input_df = df[df.gene_name == gene].reset_index(drop=True)
        x = torch.tensor(
            input_df[['n', *scores]].values.astype(float),
            dtype=torch.float32
        )

        net = SAFunction()
        if gene in ['CYP2B6', 'CYP2C19', 'CYP2C9', 'CYP2D6', 'CYP3A5', 'DPYD', 'NUDT15', 'SLCO1B1']:
            net.load_state_dict(torch.load(f'{models_dir_prefix}{gene}.pt'))
        else:
            net.load_state_dict(torch.load(f'{models_dir_prefix}OTHER.pt'))
        net.eval()

        with torch.no_grad():
            pred = pd.Series(
                net(x).numpy().flatten(), name='pmls_pred'
            ).reset_index(drop=True)

        preds_df.append(pd.concat([input_df, pred], axis=1))

    return pd.concat(preds_df)

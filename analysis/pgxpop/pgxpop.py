import hail as hl

import polars as pl

# %%
# batch vcf into 67000 participants
# run pgxpop-run.sh in PGxPOP directory

# %%

pgxpop_results_list = list()
pgxpop_results_paths = list(Path('data/pgxpop/').glob('pgxpop-part*.txt'))
for pgxpop_result_path in pgxpop_results_paths:
    df = pl.read_csv(
        pgxpop_result_path,
        null_values='None',
        schema_overrides={'sample_id': pl.String, 'activity_score': pl.Float32}
    )
    pgxpop_results_list.append(df)

pgxpop_results = pl.concat(pgxpop_results_list)

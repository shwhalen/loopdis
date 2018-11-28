#!/usr/bin/env python

import pandas as pd

chroms = set(map(str, range(1, 23)))

gwas_df = pd.read_csv(
    'gwas_catalog_v1.0.1-associations_e91_r2018-01-31.tsv',
    sep = '\t',
    header = 0,
    usecols = ['DISEASE/TRAIT', 'CHR_ID', 'CHR_POS'],
)

gwas_df.columns = ['gwas_phenotype', 'gwas_chrom', 'gwas_start']
gwas_df = gwas_df.query('gwas_chrom in @chroms')
gwas_df['gwas_start'] = (
    gwas_df['gwas_start']
    .str.split(';', expand = True)
    .iloc[:, 0]
    .astype(int)
)
gwas_df['gwas_end'] = gwas_df['gwas_start'] + 1
gwas_df.info()

(
    gwas_df
    [['gwas_chrom', 'gwas_start', 'gwas_end', 'gwas_phenotype']]
    .sort_values(['gwas_chrom', 'gwas_start', 'gwas_phenotype'])
    .reset_index(drop = True)
    .to_feather('gwas_catalog.feather')
)

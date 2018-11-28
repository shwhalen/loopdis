#!/usr/bin/env python

import pandas as pd

from chromatics import *
from common import *

output_fn = 'output/gencode19_lpa_cell_idr0.1.feather'

tss_df = pd.read_csv(
    f'{db_dir}/gencode/tss/gencode.v19.TSS.notlow.gff.gz',
    sep = '\t',
    header = None
)
start_sites = (
    tss_df
    [3]
    .astype(int)
)
gene_ids = (
    tss_df
    [8]
    .str.extract(r'(ENSG[\d\.]+)', expand = False)
)
transcript_ids = (
    tss_df
    [8]
    .str.extract(r'(ENST[\d\.]+)', expand = False)
)
tss_df = pd.concat([start_sites, gene_ids, transcript_ids], axis = 1)
tss_df.columns = ['gene_tss', 'gene_id', 'transcript_id']

expression_df = pd.read_csv(
    f'{db_dir}/gencode/expression/gencodev19_genes_with_RPKM_and_npIDR_oct2014.txt.gz',
    sep = ' '
)
expression_df = pd.melt(
    expression_df,
    id_vars = 'gene_id'
)

variables_df = (
    expression_df
    ['variable']
    .str.split('[:.]', expand = True)
)
variables_df.columns = ['lab_ids', 'rna_extract', 'cell_type', 'localization']

values_df = (
    expression_df
    ['value']
    .str.split('[:]', expand = True)
    .apply(
        lambda x: pd.to_numeric(x, 'coerce')
    )
)
values_df.columns = ['rpkm1', 'rpkm2', 'idr']

expression_df = (
    pd.concat([expression_df['gene_id'], variables_df, values_df], axis = 1)
    .query('rna_extract == "longPolyA" and localization == "cell" and idr <= 0.1')
    .eval('rpkm = (rpkm1 + rpkm2) / 2')
    .drop(['lab_ids', 'rpkm1', 'rpkm2', 'idr'], axis = 1)
    .merge(
        tss_df,
        on = 'gene_id',
        how = 'inner'
    )
)

expression_df = add_transcript_coords_to_transcript_ids(
    add_gene_coords_to_gene_names(
        add_gene_names_to_gene_ids(
            expression_df
        )
    )
    .reset_index(drop = True)
)

print(expression_df.head())
expression_df.info()

expression_df.to_feather(output_fn)

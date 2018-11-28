#!/usr/bin/env python

import feather
import json
import pandas as pd

from chromatics import *
from common import *

def write_config(example):
    print(example)

    locus_start = example[['pir_start', 'bait_start']].min()
    locus_end = example[['pir_end', 'bait_end']].max()
    margin_bp = (locus_end - locus_start) * 0.15
    locus_start = ((locus_start - margin_bp) // contact_bin_size) * contact_bin_size
    locus_end = ((locus_end + margin_bp) // contact_bin_size) * contact_bin_size

    config = {
        'eqtl_gene_cell_type': example['eqtl_gene_cell_type'],
        'ld_bin_size': 1000,
        'locus_chrom': int(example['pir_chr'][3:]),
        'locus_start': int(locus_start),
        'locus_end': int(locus_end),
        'pir_start': int(example['pir_start']),
        'pir_end': int(example['pir_end']),
        'bait_start': int(example['bait_start']),
        'bait_end': int(example['bait_end']),
        'annotations': [
            (int(example['eqtl_gene_tss_start']), example['eqtl_gene_name'] + ' (eQTL Target)'),
            (int(example['eqtl_closest_gene_tss_start']), example['eqtl_closest_gene_name'] + ' (Closest Gene to SNP)'),
            (int(example['lead_eqtl_snp_start']), example['lead_eqtl_snp_id'] + ' (eQTL SNP)')
        ]
    }

    json.dump(
        config,
        open(f'config/{example["eqtl_gene_name"]}.json', 'w'),
        indent = 2
    )


expression_intersect_columns = [
    'gene_chrom',
    'gene_tss_start',
    'gene_tss_end',
    'gene_start',
    'gene_end',
    'gene_name',
    'gene_id',
    'transcript_id',
    'gene_cell_type',
    'rpkm'
]
lead_eqtl_intersect_columns = [
    'lead_eqtl_snp_chr',
    'lead_eqtl_snp_start',
    'lead_eqtl_snp_end'
]

print('loading expression')
expression_df = (
    feather.read_dataframe('output/gencode19_lpa_cell_idr0.1.feather')
    .query('cell_type in @encode_cell_types')
    .drop(['rna_extract', 'localization'], axis = 1)
    .rename(
        columns = {
            'cell_type': 'gene_cell_type',
            'gene_tss': 'gene_tss_start'
        }
    )
    .assign(gene_tss_end = lambda x: x['gene_tss_start'] + 1)
)
high_expression_threshold = (
    expression_df
    ['rpkm']
    .quantile(0.25)
)
expression_df = expression_df.query('rpkm > @high_expression_threshold')

eqtls_df = get_eqtls()
eqtls_df['pir_chr'] = 'chr' + eqtls_df['pir_chr'].astype(str)
eqtls_df['bait_chr'] = 'chr' + eqtls_df['bait_chr'].astype(str)
eqtls_df['lead_eqtl_snp_chr'] = 'chr' + eqtls_df['lead_eqtl_snp_chr'].astype(str)
eqtls_df['lead_eqtl_snp_end'] = eqtls_df['lead_eqtl_snp_start'] + 1

print('merging expression')
eqtls_df = pd.merge(
    eqtls_df,
    expression_df.rename(columns = lambda x: 'eqtl_' + x),
    on = 'eqtl_gene_name',
    how = 'inner'
)

print('limiting interaction range')
eqtls_df = eqtls_df.query(
    'abs(lead_eqtl_snp_start - eqtl_gene_start) > 1e5 and \
    abs(lead_eqtl_snp_start - eqtl_gene_start) < 1e7'
)

print('loading rao loops')
loops_df = pd.concat(
    [
        get_rao_loops(cell_type)
        .assign(loop_cell_type = cell_type)
        for cell_type in encode_cell_types
    ],
    ignore_index = True
)
loops_df['f1_chr'] = 'chr' + loops_df['f1_chr'].astype(str)
loops_df['f2_chr'] = 'chr' + loops_df['f2_chr'].astype(str)

print('requiring eqtl overlap')
eqtl_rao_loops_df = bedtools(
    'pairtopair',
    loops_df[juicer_columns],
    (
        eqtls_df
        [['eqtl_gene_chrom', 'eqtl_gene_tss_start', 'eqtl_gene_tss_end'] + lead_eqtl_intersect_columns]
        .drop_duplicates()
    )
)

print('restoring eqtl data')
eqtls_df = pd.merge(
    eqtls_df,
    eqtl_rao_loops_df,
    on = lead_eqtl_intersect_columns + ['eqtl_gene_chrom', 'eqtl_gene_tss_start', 'eqtl_gene_tss_end'],
    how = 'inner'
)

print('restoring loop data')
eqtls_df = pd.merge(
    eqtls_df,
    loops_df,
    on = juicer_columns,
    how = 'left'
)

print('requiring loop and eqtl target cell type match')
eqtls_df = eqtls_df.query('loop_cell_type == eqtl_gene_cell_type')

print('computing closest gene to eqtl snp')
closest_genes_df = bedtools(
    'closest -t first',
    eqtls_df[lead_eqtl_intersect_columns],
    (
        expression_df
        [expression_intersect_columns]
        .rename(columns = lambda x: 'eqtl_closest_' + x)
        .query('eqtl_closest_gene_name not in @eqtls_df.eqtl_gene_name')
    )
)

print('merging closest genes')
eqtls_df = pd.merge(
    eqtls_df,
    closest_genes_df,
    on = lead_eqtl_intersect_columns,
    how = 'left'
)

# don't count trascripts separately
eqtls_df.drop_duplicates(['lead_eqtl_snp_start', 'eqtl_gene_name', 'eqtl_closest_gene_name'], inplace = True)

for i in range(len(eqtls_df)):
    print('writing config')
    write_config(eqtls_df.iloc[i])

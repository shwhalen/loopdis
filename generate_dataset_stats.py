#!/usr/bin/env python

import pandas as pd

from common import *

def get_snp_pair_count(chrom, super_pop):
    ld_fn = f'{output_dir}/plink/{chrom}-{super_pop}.ld.h5'
    return (
        pd.HDFStore(ld_fn, mode = 'r')
        .get_storer('ld')
        .nrows
    )


cell_type_stats = [
    ['Hi-C', 'Lymphoblastoid cells', 'GM12878', 9, 4907147001, 9448],
    ['Hi-C', 'Fetal lung fibroblasts', 'IMR90', 2, 1136673290, 8040],
    ['Hi-C', 'Epidermal skin keratinocytes', 'NHEK', 1, 664899299, 4929],
    ['Hi-C', 'Erythroleukemia cells', 'K562', 2, 932208867, 6057],
    ['Hi-C', 'Umbilical vein endothelial cells', 'HUVEC', 1, 460393495, 3865],
    ['PCHi-C', 'Megakaryocytes', 'MK', 4, 653848788, 150779],
    ['PCHi-C', 'Erythroblasts', 'Ery', 3, 588786672, 151215],
    ['PCHi-C', 'Neutrophils', 'Neu', 3, 736055569, 142435],
    ['PCHi-C', 'Monocytes', 'Mon', 3, 572357387, 165947],
    ['PCHi-C', 'Macrophages M0', 'Mac0', 3, 668675248, 180190],
    ['PCHi-C', 'Macrophages M1', 'Mac1', 3, 497683496, 171031],
    ['PCHi-C', 'Macrophages M2', 'Mac2', 3, 523561551, 186172],
    ['PCHi-C', 'Endothelial precursors', 'EP', 3, 420536621, 145888],
    ['PCHi-C', 'Naive B cells', 'nB', 3, 629928642, 189720],
    ['PCHi-C', 'Total B cells', 'tB', 3, 702533922, 213539],
    ['PCHi-C', 'Fetal thymus', 'FoeT', 3, 776491344, 166743],
    ['PCHi-C', 'Naive CD4$^+$ T cells', 'nCD4', 4, 844697853, 210074],
    ['PCHi-C', 'Total CD4$^+$ T cells', 'tCD4', 3, 836974777, 199525],
    ['PCHi-C', 'Non-activated total CD4$^+$ T cells', 'naCD4', 3, 721030702, 211720],
    ['PCHi-C', 'Activated total CD4$^+$ T cells', 'aCD4', 3, 749720649, 213235],
    ['PCHi-C', 'Naive CD8$^+$ T cells', 'nCD8', 3, 747834572, 216232],
    ['PCHi-C', 'Total CD8$^+$ T cells', 'tCD8', 3, 628771947, 204382]
]
cell_type_stats_df = pd.DataFrame(
    cell_type_stats,
    columns = ['Assay', 'Cell Type', 'Acronym', 'Reps.', 'Unique Read Pairs', 'Sig. Interactions']
)
cell_type_stats_df.to_latex(
    'output/cell_type_stats-table.tex',
    index = False,
    escape = False,
    formatters = {
        'Unique Read Pairs': comma_formatter,
        'Sig. Interactions': comma_formatter
    }
)

ld_blocks_df = pd.concat(
    [
        get_plink_ld_blocks(None, super_pop)
        .eval('distance = ld_block_end - ld_block_start')
        .drop(['ld_block_start', 'ld_block_end'], axis = 1)
        .assign(super_pop = super_pop)
        for super_pop in super_pops
    ],
    ignore_index = True
)

ld_block_quantiles_df = (
    ld_blocks_df
    .groupby('super_pop')
    .quantile([0.05, 0.5, 0.95])
    .reset_index(level = 1)
    .rename(columns = {'level_1': 'quantile'})
    .pivot_table(index = 'super_pop', columns = 'quantile', values = ['distance', 'ld_block_nsnps'])
)

ld_block_counts_df = (
    ld_blocks_df
    .groupby('super_pop')
    .size()
    .rename('total_blocks')
)

ld_counts = [
    (get_snp_pair_count(chrom, super_pop), chrom, super_pop)
    for chrom in chroms
    for super_pop in super_pops
]
ld_counts_df = pd.DataFrame(
    ld_counts,
    columns = ['total_snps', 'chrom', 'super_pop']
)
ld_counts_df = (
    ld_counts_df
    .groupby('super_pop')
    ['total_snps']
    .sum()
)

ld_stats_df = (
    pd.concat([ld_counts_df, ld_block_counts_df, ld_block_quantiles_df], axis = 1)
    .astype(int)
    .reset_index()
)
ld_stats_df.columns = [
    'Superpop.',
    '# SNP Pairs',
    '# LD Blocks',
    'BL .05',
    'BL .5',
    'BL .95',
    'BP .05',
    'BP .5',
    'BP .95'
]
formatters = dict(
    zip(
        ld_stats_df.columns[1:],
        [comma_formatter] * (ld_stats_df.shape[1] - 1)
    )
)
ld_stats_df.to_latex(
    'output/ld_stats-table.tex',
    index = False,
    formatters = formatters
)

#!/usr/bin/env python

import pandas as pd

from chromatics import *
from common import *
from sklearn.externals.joblib import Parallel, delayed

def permute_distances(ld_blocks_df, domain_boundaries_df, permutation):
    if permutation > 0:
        ld_blocks_df = bedtools(
            f'shuffle -chrom -seed {permutation}',
            ld_blocks_df,
            genome = 'GRCh37',
            sort = False
        )
    ld_blocks_df['ld_block_name'] = concat_coords(ld_blocks_df)

    # convert block enpoints into unique boundary coordinates for intersecting
    a = (
        ld_blocks_df
        [['ld_block_chr', 'ld_block_start']]
        .rename(columns = {'ld_block_start': 'ld_block_boundary_start'})
    )
    b = (
        ld_blocks_df
        [['ld_block_chr', 'ld_block_end']]
        .rename(columns = {'ld_block_end': 'ld_block_boundary_start'})
    )
    ld_block_boundaries_df = pd.concat(
        [a, b],
        ignore_index = True,
        sort = False
    )
    ld_block_boundaries_df['ld_block_boundary_end'] = ld_block_boundaries_df['ld_block_boundary_start'] + 1
    ld_block_boundaries_df.drop_duplicates(['ld_block_chr', 'ld_block_boundary_start'], inplace = True)

    domain_boundaries_closest_ld_boundaries_df = (
        bedtools(
            'closest -d -t first',
            domain_boundaries_df,
            ld_block_boundaries_df,
            genome = 'GRCh37'
        )
        .query('distance != -1') # filter failed intersections
    )

    return pd.DataFrame({
        'cell_line': cell_line,
        'super_pop': super_pop,
        'permutation': permutation,
        'distance': domain_boundaries_closest_ld_boundaries_df['distance']
    })


cell_line = 'GM12878'
super_pop = 'AFR'
n_permutations = 1000

domains_df = get_rao_domains(cell_line)

# convert contact domain boundaries into unique boundary coordinates for intersecting
domain_boundaries_df = (
    domains_df
    .melt(
        id_vars = 'domain_chr',
        value_vars = ['domain_start', 'domain_end'],
        value_name = 'domain_boundary_start'
    )
    .drop('variable', axis = 1)
    .rename(columns = {'domain_chr': 'domain_boundary_chr'})
    .assign(domain_boundary_end = lambda x: x['domain_boundary_start'] + 1)
)
domain_boundaries_df.drop_duplicates(['domain_boundary_chr', 'domain_boundary_start'], inplace = True)

ld_blocks_df = (
    get_plink_ld_blocks(None, super_pop)
    [['ld_block_chr', 'ld_block_start', 'ld_block_end']]
)

ld_blocks_closest_domain_boundaries_df = bedtools(
    'closest -d -t first',
    ld_blocks_df,
    domain_boundaries_df,
    genome = 'GRCh37'
)
ld_blocks_closest_domain_boundaries_df['cell_line'] = cell_line
ld_blocks_closest_domain_boundaries_df['super_pop'] = super_pop
ld_blocks_closest_domain_boundaries_df.to_feather(f'output/boundary_distances-ld_blocks_to_domains-{cell_line}-{super_pop}.feather')

print(ld_blocks_df.eval('ld_block_end - ld_block_start').describe())
print(
    bedtools('intersect -wa -wb', domain_boundaries_df, ld_blocks_df, genome = 'GRCh37')
    .eval('abs(ld_block_end - ld_block_start)')
    .describe()
)
print(
    bedtools('closest -d -t first -io', domain_boundaries_df, ld_blocks_df, genome = 'GRCh37')
    .eval('abs(ld_block_end - ld_block_start)')
    .describe()
)

# distances = Parallel(8)(
#     delayed(permute_distances)(ld_blocks_df, domain_boundaries_df, _)
#     for _ in tqdm(list(range(n_permutations)))
# )

# (
#     pd.concat(distances, ignore_index = True, sort = False)
#     .to_feather(f'output/boundary_distances-domains_to_ld_blocks-{cell_line}-{super_pop}.feather')
# )

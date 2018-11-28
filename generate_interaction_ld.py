#!/usr/bin/env python

import itertools
import os
import pandas as pd
import sys

from chromatics import *
from common import *
from tqdm import tqdm

def get_pairwise_ld(group_df):
    f1_mask = group_df.eval('f1_start <= snp_start <= f1_end')
    if group_df.iloc[0]['f2_start'] < group_df.iloc[0]['f1_start']:
        f1_mask = ~f1_mask
    f1_snps = group_df.loc[f1_mask, 'snp_start'].tolist()
    f2_snps = group_df.loc[~f1_mask, 'snp_start'].tolist()
    pairs = itertools.product(f1_snps, f2_snps)
    return ld.reindex(pairs).tolist()


def split_javierre_contacts():
    contacts_df = (
        get_javierre_contacts()
        .query('bait_chr == @chrom and abs(dist) < @max_interaction_distance')
    )
    significant_pchic_mask = contacts_df[cell_type] >= significant_pchic_threshold

    loops_df = contacts_df.loc[significant_pchic_mask, bait_columns[:3] + oe_columns[:3]]
    loops_df.columns = juicer_columns

    nonloops_df = contacts_df.loc[~significant_pchic_mask, bait_columns[:3] + oe_columns[:3]]
    nonloops_df.columns = juicer_columns

    return distance_match(loops_df, nonloops_df)


def split_rao_contacts():
    loops_df = (
        get_rao_loops(cell_type)
        .query('f1_chr == @chrom and (f2_start - f1_end) < @max_interaction_distance')
    )

    chrom_sizes = pd.read_csv(
        f'{genomes_dir}/GRCh37/GRCh37.chrom.sizes',
        sep = '\t',
        header = None,
        names = ['chrom', 'size'],
        index_col = 'chrom',
        squeeze = True
    )
    chrom_size = chrom_sizes.loc[chrom]

    nonloops = []
    np.random.seed(0)
    for index, row in loops_df.iterrows():
        while row['f1_start'] in frozenset(loops_df['f1_start']):
            interaction_distance = row['f2_start'] - row['f1_start']
            fragment_length = row['f1_end'] - row['f1_start']
            new_start = np.random.randint(chrom_size) % (chrom_size - interaction_distance - fragment_length)

            row['f1_start'] = new_start
            row['f1_end'] = row['f1_start'] + fragment_length

            row['f2_start'] = new_start + interaction_distance
            row['f2_end'] = row['f2_start'] + fragment_length
        nonloops.append(row)
    nonloops_df = pd.DataFrame(nonloops)
    assert not loops_df.equals(nonloops_df)
    assert len(loops_df) == len(nonloops_df)

    return loops_df, nonloops_df


chrom = sys.argv[1]
super_pop = sys.argv[2]
cell_type = sys.argv[3]
dataset = 'javierre' if cell_type in blood_cell_types else 'rao'

output_fn = f'{output_dir}/{dataset}/fragment_snps_ld-{chrom}-{super_pop}-{cell_type}.h5'
if os.path.exists(output_fn):
    print(f'{output_fn} exists, exiting...')
    sys.exit()

tqdm.pandas()

# use full list of snps for intersection so that snps below plink threshold will return nans, then nans converted to 0 (or some other approach) so they can still influence average.  excluding snps below threshold would inflate average.
snps_df = pd.read_csv(
    f'{output_dir}/plink/{chrom}-{super_pop}.bim',
    sep = '\t',
    header = None,
    usecols = [0, 3],
    names = ['snp_chr', 'snp_start']
)
snps_df['snp_end'] = snps_df['snp_start']

print('reading, splitting, and distance matching fragments')
if dataset == 'javierre':
    loop_fragments_df, nonloop_fragments_df = split_javierre_contacts()
else:
    loop_fragments_df, nonloop_fragments_df = split_rao_contacts()

print('intersecting looping fragments with snps')
loop_fragment_snps_df = bedtools(
    'pairtobed -type both',
    loop_fragments_df,
    snps_df,
    genome = 'GRCh37',
    sort = False
)
loop_fragment_snps_df['is_significant'] = 1

print('intersecting nonlooping fragments with snps')
nonloop_fragment_snps_df = bedtools(
    'pairtobed -type both',
    nonloop_fragments_df,
    snps_df,
    genome = 'GRCh37',
    sort = False
)
nonloop_fragment_snps_df['is_significant'] = 0

print('combining fragment snps')
fragment_snps_df = (
    pd.concat([loop_fragment_snps_df, nonloop_fragment_snps_df])
    .drop(['snp_chr', 'snp_end'], axis = 1)
)
assert fragment_snps_df.eval('f1_end > f1_start and f2_end > f2_start').all()

print('reading and filtering ld')
# separating these into fragment-specific snps is tricky due to dataset differences.
# if modified to do so, don't use global snp sets inside get_pairwise_ld as it expects disjoint sets which is not true outside the groupby.
all_snps = frozenset(fragment_snps_df['snp_start'])
ld = (
    get_plink_ld(chrom, super_pop)
    .query('(snp1_start in @all_snps) and (snp2_start in @all_snps)')
    .squeeze()
)

print('computing ld between fragment snps')
fragment_snps_ld_df = (
    fragment_snps_df
    .groupby(['f1_start', 'f1_end', 'f2_start', 'f2_end', 'is_significant'])
    .progress_apply(get_pairwise_ld)
    .rename('ld')
    .reset_index()
)

import warnings
warnings.filterwarnings('ignore', category = pd.io.pytables.PerformanceWarning)

fragment_snps_ld_df.to_hdf(
    output_fn,
    'ld',
    mode = 'w',
    complevel = 1,
    complib = 'lzo'
)

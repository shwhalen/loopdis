#!/usr/bin/env python

import os
import pandas as pd
import sklearn.externals.joblib as joblib
import sys

from chromatics import *
from common import *
from tqdm import tqdm

def get_loop_superpop_crossings(loops_df, super_pop):
    ld_blocks_df = get_plink_ld_blocks(None, super_pop)
    crossings_df = bedtools(
        'intersect -c',
        loops_df,
        ld_blocks_df,
        genome = 'GRCh37'
    )
    crossings_df['super_pop'] = super_pop
    return crossings_df[loops_df.columns.tolist() + ['super_pop', 'count']]


def get_loop_crossings(loops_df):
    crossings = joblib.Parallel(-1)(
        joblib.delayed(get_loop_superpop_crossings)(loops_df, _)
        for _ in super_pops
    )
    return pd.concat(crossings, ignore_index = True)


def get_domain_superpop_crossings(domains_df, super_pop):
    ld_blocks_df = get_plink_ld_blocks(None, super_pop)
    crossings_df = bedtools(
        'intersect -c',
        ld_blocks_df,
        domains_df,
        genome = 'GRCh37'
    )
    crossings_df['super_pop'] = super_pop
    return crossings_df[ld_blocks_df.columns.tolist() + ['super_pop', 'count']]


def get_domain_crossings(domains_df):
    crossings = joblib.Parallel(-1)(
        joblib.delayed(get_domain_superpop_crossings)(domains_df, _)
        for _ in super_pops
    )
    return pd.concat(crossings, ignore_index = True)


# rao loops
output_fn = f'{output_dir}/interaction_crossings-ld_blocks-rao.feather'
crossings = []
for cell_type in tqdm(encode_cell_types, 'cell line'):
    loops_df = get_rao_loops(cell_type)
    loops_df['window_chr'] = loops_df['f1_chr']
    loops_df['window_start'] = loops_df[['f1_start', 'f2_start']].min(axis = 1)
    loops_df['window_end'] = loops_df[['f1_end', 'f2_end']].max(axis = 1)

    crossings_df = get_loop_crossings(loops_df[['window_chr', 'window_start', 'window_end']])
    crossings_df['cell_type'] = cell_type
    crossings.append(crossings_df)
crossings_df = pd.concat(crossings, ignore_index = True)
crossings_df['window_chr'] = crossings_df['window_chr'].astype('str')
crossings_df['super_pop'] = crossings_df['super_pop'].astype('category')
crossings_df['cell_type'] = crossings_df['cell_type'].astype('category')
crossings_df.to_feather(output_fn)

# javierre loops
output_fn = f'{output_dir}/interaction_crossings-ld_blocks-javierre.feather'
crossings = []
for cell_type in tqdm(blood_cell_types, 'cell line'):
    loops_df = get_javierre_loops().query(f'{cell_type} >= @significant_pchic_threshold')
    loops_df['window_chr'] = loops_df['bait_chr']
    loops_df['window_start'] = loops_df[['bait_start', 'oe_start']].min(axis = 1)
    loops_df['window_end'] = loops_df[['bait_end', 'oe_end']].max(axis = 1)

    crossings_df = get_loop_crossings(loops_df[['window_chr', 'window_start', 'window_end']])
    crossings_df['cell_type'] = cell_type
    crossings.append(crossings_df)
crossings_df = pd.concat(crossings, ignore_index = True)
crossings_df['window_chr'] = crossings_df['window_chr'].astype('str')
crossings_df['super_pop'] = crossings_df['super_pop'].astype('category')
crossings_df['cell_type'] = crossings_df['cell_type'].astype('category')
crossings_df.to_feather(output_fn)

# rao domains only, there are no pre-called javierre domains
output_fn = f'{output_dir}/interaction_crossings-domains-rao.feather'
crossings = []
for cell_type in tqdm(encode_cell_types, 'cell_line'):
    domains_df = get_rao_domains(cell_type)
    crossings_df = get_domain_crossings(domains_df)
    crossings_df['cell_type'] = cell_type
    crossings.append(crossings_df)
crossings_df = pd.concat(crossings, ignore_index = True)
crossings_df.to_feather(output_fn)

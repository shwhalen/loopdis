#!/usr/bin/env python

import itertools
import numpy as np
import os
import pandas as pd
import sys

from chromatics import *
from common import *
from tqdm import tqdm

def bin(df, bin_size):
    return (
        df
        .rename_axis(['start_bin', 'end_bin'])
        .reset_index()
        .assign(
            start_bin = lambda x: x['start_bin'] // bin_size,
            end_bin = lambda x: x['end_bin'] // bin_size
        )
    )


def get_bin_extremes(df1, df2):
    min_bin = min(
        df1.index.get_level_values('start_bin').min(),
        df2.index.get_level_values('start_bin').min()
    )
    max_bin = max(
        df1.index.get_level_values('end_bin').max(),
        df2.index.get_level_values('end_bin').max()
    )
    return min_bin, max_bin


def get_windowed_indices(window_size, min_bin, max_bin):
    row_indices, col_indices = np.triu_indices(window_size)
    paired_indices = [
        zip(row_indices + i, col_indices + i)
        for i in range(min_bin, max_bin, window_size)
    ]
    return list(itertools.chain.from_iterable(paired_indices))


def get_windowed_concordance(aggregate_ld_df, contacts_df, window_size, strong_contact_threshold):
    windowed_indices = get_windowed_indices(
        window_size,
        *get_bin_extremes(aggregate_ld_df, contacts_df)
    )
    windowed_aggregate_ld_df = aggregate_ld_df.loc[windowed_indices]
    windowed_contacts_df = contacts_df.loc[windowed_indices]
    windowed_aggregate_ld_contacts_df = pd.concat([windowed_aggregate_ld_df, windowed_contacts_df], axis = 1)

    table = pd.crosstab(
        windowed_aggregate_ld_contacts_df['r2'] > strong_ld_threshold,
        windowed_aggregate_ld_contacts_df['contact_value'] > strong_contact_threshold,
        margins = True,
        normalize = True
    )
    return table.loc[True, True], table.loc[True, 'All'] * table.loc['All', True]


def get_hic_ld_concordance(super_pop, chrom):
    output_fn = f'{output_dir}/hic_ld_concordance/hic_ld_concordance-{chrom}-{super_pop}-{cell_type}-{oe_suffix}-strong_ld_q{strong_ld_quantile}-strong_ld_threshold_{strong_ld_threshold}-strong_contact_q{strong_contact_quantile}_from_window_size_{window_sizes[-1] * contact_bin_size}.feather'
    if os.path.exists(output_fn):
        print(f'{output_fn} exists, skipping...')
        return

    # faster than a groupby, bottleneck is np.percentile
    aggregate_ld_df = pd.pivot_table(
        bin(get_plink_ld(chrom, super_pop), ld_bin_size),
        index = ['start_bin', 'end_bin'],
        values = 'r2',
        aggfunc = lambda x: np.percentile(x, 100 * strong_ld_quantile)
    )

    contacts_df = bin(get_rao_contacts(chrom, cell_type, observed_over_expected), contact_bin_size)
    contacts_df['contact_value'] = contacts_df['contact_value'].astype(np.float32)
    contacts_df.set_index(['start_bin', 'end_bin'], inplace = True)

    paired_indices = get_windowed_indices(
        window_sizes[-1],
        *get_bin_extremes(aggregate_ld_df, contacts_df)
    )
    strong_contact_threshold = (
        contacts_df
        .loc[paired_indices, 'contact_value']
        .quantile(strong_contact_quantile)
    )

    results = [
        get_windowed_concordance(aggregate_ld_df, contacts_df, window_size, strong_contact_threshold)
        for window_size in window_sizes
    ]
    percent_both_high, percent_both_high_baseline = zip(*results)

    stats_df = pd.DataFrame({
        'window_size': window_sizes,
        'percent_both_high': percent_both_high,
        'baseline': False,
        'cell_type': cell_type,
        'chrom': chrom,
        'super_pop': super_pop
    })
    stats_baseline_df = stats_df.copy()
    stats_baseline_df['percent_both_high'] = percent_both_high_baseline
    stats_baseline_df['baseline'] = True
    (
        pd.concat([stats_df, stats_baseline_df], ignore_index = True)
        .to_feather(output_fn)
    )


cell_type = sys.argv[1]

window_sizes = [1, 2, 4, 8, 16, 32, 64, 128, 256]
ld_bin_size = contact_bin_size
observed_over_expected = False # expectation correlates with distance, don't apply expectation to contacts but not ld
oe_suffix = 'oe' if observed_over_expected else 'o'
strong_ld_quantile = 0.75
strong_contact_quantile = 0.75

params = list(itertools.product(super_pops, chroms))
for super_pop, chrom in tqdm(params):
    get_hic_ld_concordance(super_pop, chrom)

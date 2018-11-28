#!/usr/bin/env python

import itertools
import numpy as np
import pandas as pd
import sklearn.externals.joblib as joblib

from common import *
from tqdm import tqdm

def get_max_ld(chrom, super_pop, cell_type):
    dataset = 'javierre' if cell_type in blood_cell_types else 'rao'
    return (
        pd.read_hdf(f'{output_dir}/{dataset}/fragment_snps_ld-{chrom}-{super_pop}-{cell_type}.h5', 'ld')
        .assign(f1_chr = chrom, f2_chr = chrom, cell_type = cell_type, super_pop = super_pop)
        .set_index(['cell_type' ,'super_pop'] + juicer_columns + ['is_significant'])
        .squeeze()
        .apply(lambda x: [0] if len(x) == 0 else x)
        .apply(np.nan_to_num)
        .apply(np.max)
        .rename('max_ld')
    )


params = list(itertools.product(chroms, super_pops, encode_cell_types + blood_cell_types))
max_lds = joblib.Parallel(-1)(
    joblib.delayed(get_max_ld)(chrom, super_pop, cell_type)
    for chrom, super_pop, cell_type in tqdm(params)
)

max_ld_df = (
    pd.concat(max_lds)
    .reset_index()
)
max_ld_df['is_significant'] = max_ld_df['is_significant'].astype('category')
max_ld_df['super_pop'] = max_ld_df['super_pop'].astype('category')
max_ld_df['cell_type'] = max_ld_df['cell_type'].astype('category')
max_ld_df.to_feather(f'{output_dir}/max_interaction_ld.feather')

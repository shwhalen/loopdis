#!/usr/bin/env python

import pandas as pd

from common import *
from tqdm import tqdm

chrom = 20
csv_output_fn = f'{output_dir}/plink/{chrom}.ld.csv'
feather_output_fn = f'{output_dir}/plink/{chrom}.ld.feather'

for super_pop in tqdm(super_pops):
    mode = 'a'
    header = False
    if super_pop == super_pops[0]:
        mode = 'w'
        header = True

    x = pd.read_csv(
        f'{output_dir}/plink/{chrom}-{super_pop}.ld.gz',
        sep = r'\s+',
        header = 0,
        names = ['chr_a', 'bp_a', 'snp_a', 'chr_b', 'bp_b', 'snp_b', 'r2', 'dprime'],
        usecols = ['bp_a', 'bp_b', 'r2', 'dprime']
    )

    x.eval('distance = abs(bp_b - bp_a)', inplace = True)
    x['distance'] = x['distance'].astype(int)
    del x['bp_a']
    del x['bp_b']

    x['super_pop'] = super_pop

    (
        x
        [['super_pop', 'distance', 'r2', 'dprime']]
        .to_csv(csv_output_fn, mode = mode, header = header, index = False)
    )

# writing feather will fail if super_pop is a string instead of a categorical
x = pd.read_csv(
    csv_output_fn,
    dtype = {'super_pop': 'category', 'distance': int, 'r2': float, 'dprime': float}
)
x.to_feather(feather_output_fn)

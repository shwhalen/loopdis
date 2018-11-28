#!/usr/bin/env python

import numpy as np
import pandas as pd
import sklearn.externals.joblib as joblib
import sys

from chromatics import *
from common import *
from tqdm import tqdm

def get_domain_snps(chrom, shuffle):
    ld_fn = f'{output_dir}/plink/{chrom}-{super_pop}.ld.h5'
    store = pd.HDFStore(ld_fn, mode = 'r')
    nrows = store.get_storer('ld').nrows

    np.random.seed(0)
    subsampled_rows = np.random.randint(0, nrows, size = int(5e6))
    ld_df = (
        pd.read_hdf(ld_fn, 'ld', mode = 'r', where = pd.Index(subsampled_rows))
        .reset_index()
    )
    ld_df['r2'] = ld_df['r2'].round(3) # reduces memory footprint considerably when dumping to stdin of bedtools
    ld_df.insert(0, 'chrom', chrom)

    domains_df = (
        get_rao_domains(cell_type)
        .query(f'domain_chr == "{chrom}"')
    )
    merged_domains_df = bedtools(
        'merge',
        domains_df,
        genome = 'GRCh37'
    )
    if len(merged_domains_df) == 0:
        print(f'no contact domains for chrom {chrom} in cell type {cell_type}')
        return

    if shuffle:
        merged_domains_df = bedtools(
            'shuffle -chrom -noOverlapping',
            merged_domains_df,
            genome = 'GRCh37'
        )

    domain_snps_df = bedtools(
        'intersect -sorted -wo',
        ld_df,
        merged_domains_df,
        genome = 'GRCh37'
    )
    del ld_df
    assert domain_snps_df.eval('snp2_start > snp1_start').all()

    return (
        domain_snps_df
        .drop_duplicates(subset = ['snp1_start', 'snp2_start'])
        .eval('distance = snp2_start - snp1_start')
        .eval('intra = distance == overlap') # snp pair is intra-domain if overlap spans full distance between snps
        [['distance', 'r2', 'intra']]
    )


super_pop = sys.argv[1]
cell_type = sys.argv[2]
shuffle = int(sys.argv[3])

shuffle_suffix = 'shuffled' if shuffle else 'unshuffled'
output_fn = f'{output_dir}/domain_snps-{cell_type}-{super_pop}-{shuffle_suffix}.feather'

if os.path.exists(output_fn):
    print(f'{output_fn} exists, exiting...')
    sys.exit()

domain_snps = joblib.Parallel(-1)(
    joblib.delayed(get_domain_snps)(_, shuffle)
    for _ in tqdm(chroms)
)
domain_snps_df = pd.concat(domain_snps, ignore_index = True)
domain_snps_df.to_feather(output_fn)

#!/usr/bin/env python

from common import *
from chromatics import *

pchic_columns = ['oe_chr', 'oe_start', 'oe_end', 'bait_chr', 'bait_start', 'bait_end']

ld_blocks_df = (
    get_plink_ld_blocks(None, 'EAS')
    [['ld_block_chr', 'ld_block_start', 'ld_block_end']]
)
ld_blocks_df['ld_block_chr'] = ensembl_to_ucsc_chroms(ld_blocks_df['ld_block_chr'])

eqtls_df = (
    pd.read_excel('fairfax/ng.2205-S2.xls', sheet_name = 2)
    [['CHR', 'BP', 'Gene']]
)
eqtls_df.columns = ['eqtl_chrom', 'eqtl_start', 'gene_name']

eqtls_df['eqtl_chrom'] = ensembl_to_ucsc_chroms(eqtls_df['eqtl_chrom'])
eqtls_df['eqtl_end'] = eqtls_df['eqtl_start'] + 1
eqtls_df = add_gene_data(eqtls_df)
eqtls_df = eqtls_df[['eqtl_chrom', 'eqtl_start', 'eqtl_end', 'gene_chrom', 'gene_start', 'gene_end', 'gene_name']]
eqtls_df['eqtl_name'] = concat_coords(
    eqtls_df,
    ['eqtl_chrom', 'eqtl_start', 'gene_name']
)

eqtls_df = eqtls_df.query('abs(eqtl_start - gene_start) < 2e6 or abs(eqtl_start - gene_end) < 2e6')
eqtl_to_gene_start = eqtls_df.eval('abs(eqtl_start - gene_start)')
eqtl_to_gene_end = eqtls_df.eval('abs(eqtl_start - gene_end)')
eqtl_distances = pd.DataFrame([eqtl_to_gene_start, eqtl_to_gene_end]).min(axis = 0)
eqtls_df['distance_bin'] = pd.qcut(eqtl_distances, q = 4)

contacts_df = (
    get_javierre_loops()
    [pchic_columns]
)
contacts_df['oe_chr'] = ensembl_to_ucsc_chroms(contacts_df['oe_chr'])
contacts_df['bait_chr'] = ensembl_to_ucsc_chroms(contacts_df['bait_chr'])

results = []
for distance_bin, binned_eqtls_df in eqtls_df.groupby('distance_bin'):
    binned_eqtls_df = binned_eqtls_df.drop('distance_bin', axis = 1)
    contacts_eqtls_df = bedtools(
        'pairtopair -type both',
        contacts_df,
        binned_eqtls_df
    )
    contacts_eqtls_df['eqtl_name'] = concat_coords(
        contacts_eqtls_df,
        ['eqtl_chrom', 'eqtl_start', 'gene_name']
    )
    eqtls_with_contacts_percent = binned_eqtls_df.eval('eqtl_name in @contacts_eqtls_df.eqtl_name').mean()

    ld_blocks_eqtls_df = bedtools(
        'pairtobed -type both',
        binned_eqtls_df,
        ld_blocks_df
    )
    ld_blocks_eqtls_df['ld_block_name'] = concat_coords(
        ld_blocks_eqtls_df,
        ['ld_block_chr', 'ld_block_start', 'ld_block_end']
    )
    ld_blocks_eqtls_df = (
        ld_blocks_eqtls_df
        .groupby(['eqtl_chrom', 'eqtl_start', 'gene_name'])
        ['ld_block_name']
        .nunique()
        .rename('n_unique_ld_blocks')
        .reset_index()
        .query('n_unique_ld_blocks == 1')
    )
    ld_blocks_eqtls_df['eqtl_name'] = concat_coords(
        ld_blocks_eqtls_df,
        ['eqtl_chrom', 'eqtl_start', 'gene_name']
    )
    eqtls_with_same_ld_block_percent = binned_eqtls_df.eval('eqtl_name in @ld_blocks_eqtls_df.eqtl_name').mean()

    results.append((
        str(distance_bin),
        eqtls_with_contacts_percent,
        eqtls_with_same_ld_block_percent
    ))
pd.DataFrame(results, columns = ['distance_bin', 'contacts_percent', 'same_ld_block_percent']).to_feather('output/eqtl_support.feather')

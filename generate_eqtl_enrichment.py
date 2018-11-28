#!/usr/bin/env python

import feather
import pandas as pd

from common import *

# fairfax eQTLs are for b-cells
eqtls_df = get_eqtls()
eqtl_cell_type = 'nB'

print('\nreading ld data')
max_interaction_ld_df = (
    feather.read_dataframe(f'{output_dir}/max_interaction_ld.feather')
    .query(f'cell_type == "{eqtl_cell_type}"')
    .rename(columns = {
        'f1_chr': 'bait_chr',
        'f1_start': 'bait_start',
        'f1_end': 'bait_end',
        'f2_chr': 'oe_chr',
        'f2_start': 'oe_start',
        'f2_end': 'oe_end'
    })
)

stacked_index_cols = max_interaction_ld_df.columns.tolist()[:-1]
stacked_index_cols.remove('super_pop')
super_pop_lds = [
    max_interaction_ld_df
    .query('super_pop == @super_pop')
    .rename(columns = {'max_ld': 'max_ld_' + super_pop})
    .drop('super_pop', axis = 1)
    .set_index(stacked_index_cols)
    for super_pop in super_pops
]
max_interaction_ld_df = (
    pd.concat(super_pop_lds, axis = 1)
    .fillna(0)
    .reset_index()
    .assign(
        bait_chr = lambda x: x['bait_chr'].astype(str),
        oe_chr = lambda x: x['oe_chr'].astype(str)
    )
)

print('\nreading contacts')
contacts_df = get_javierre_contacts()

contacts_df['bait_gene_names'] = (
    contacts_df['bait_gene_names']
    .fillna('')
    .str.split(';', expand = False)
    .apply(
        lambda x: set(x)
    )
)

contacts_df['oe_gene_names'] = (
    contacts_df['oe_gene_names']
    .str.replace('^.', '')
    .str.split(';', expand = False)
    .apply(
        lambda x: set(x)
    )
)

contacts_df = get_closest_genes(contacts_df, 'GRCh37')

contacts_df = pd.merge(
    contacts_df,
    max_interaction_ld_df,
    on = bait_columns[:3] + oe_columns[:3],
    how = 'inner'
)

print('\ncomputing interaction ids')
eqtl_ids = set(
    eqtls_df['pir_id'].astype(str) +
    ':' +
    eqtls_df['bait_id'].astype(str)
)
contacts_df['interaction_id'] = (
    contacts_df['oe_id'].astype(str) +
    ':' +
    contacts_df['bait_id'].astype(str)
)

print('\ncomputing if oe contains an eqtl for baited promoter')
contacts_df['oe_contains_eqtl'] = (
    contacts_df
    ['interaction_id']
    .isin(eqtl_ids)
)

print('\ncomputing if bait is closest gene to oe')
contacts_df['oe_closest_gene_name'] = (
    contacts_df
    ['oe_closest_gene_name']
    .apply(
        lambda x: set([x])
    )
)

# element-wise set subtraction, if length of resulting set is 0 then the closest gene was found in the bait genes. much faster than set subtraction with apply.
contacts_df['is_bait_closest_gene_to_oe'] = (
    (
        contacts_df['oe_closest_gene_name'] -
        contacts_df['bait_gene_names']
    )
    .str.len() == 0
)

print('\ncomputing if interaction is significant')
contacts_df['is_significant'] = contacts_df[eqtl_cell_type] >= significant_pchic_threshold

eqtl_stats_df = contacts_df[
    bait_columns[:3] +
    oe_columns[:3] +
    ['oe_contains_eqtl', 'is_bait_closest_gene_to_oe', 'is_significant'] +
    ['max_ld_' + _ for _ in super_pops]
]
eqtl_stats_df.to_feather(f'{output_dir}/eqtl_enrichment.feather')

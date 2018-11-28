#!/usr/bin/env python

import pandas as pd
import sklearn.externals.joblib as joblib

from chromatics import *
from common import *
from tqdm import tqdm

def get_accuracies(cell_type):
    cell_type_contacts_df = get_closest_genes(
        contacts_df.query(f'{cell_type} >= @significant_pchic_threshold'),
        genome = 'GRCh37'
    )

    cell_type_contacts_df['oe_closest_gene_name'] = (
        cell_type_contacts_df
        ['oe_closest_gene_name']
        .apply(
            lambda x: set([x])
        )
    )

    cell_type_contacts_df['is_bait_closest_gene_to_oe'] = (
        (
            cell_type_contacts_df['oe_closest_gene_name'] -
            cell_type_contacts_df['bait_gene_names']
        )
        .str.len() == 0
    )

    closest_accuracy = (
        cell_type_contacts_df
        ['is_bait_closest_gene_to_oe']
        .mean()
    )

    super_pop_accuracies = [
        get_contact_ld_blocks(cell_type_contacts_df, super_pop)
        ['oe_shares_ld_block_with_bait']
        .mean()
        for super_pop in super_pops
    ]

    return [cell_type, closest_accuracy] + super_pop_accuracies


def get_contact_ld_blocks(contacts_df, super_pop):
    def get_ld_block_set(fragment_columns, prefix):
        fragment_contacts_df = (
            contacts_df
            [fragment_columns[:3]]
            .drop_duplicates()
        )
        fragment_ld_blocks_df = ld_blocks_df.rename(columns = lambda x: prefix + x)
        fragment_ld_block_columns = fragment_ld_blocks_df.columns.tolist()
        fragment_ld_blocks_df = bedtools(
            'intersect -sorted -wa -wb',
            fragment_contacts_df,
            fragment_ld_blocks_df,
            genome = 'GRCh37',
        )

        fragment_ld_blocks_df[prefix + 'ld_block_name'] = concat_coords(
            fragment_ld_blocks_df,
            fragment_ld_block_columns
        )

        return (
            fragment_ld_blocks_df
            .groupby(fragment_columns[:3])
            .apply(
                lambda x: set(
                    x
                    [prefix + 'ld_block_name']
                    .unique()
                )
            )
            .rename(prefix + 'ld_block_names')
            .reset_index()
        )

    ld_blocks_df = get_plink_ld_blocks(None, super_pop)

    oe_ld_blocks_df = get_ld_block_set(oe_columns[:3], 'oe_')
    bait_ld_blocks_df = get_ld_block_set(bait_columns[:3], 'bait_')

    contacts_df = pd.merge(contacts_df, oe_ld_blocks_df, on = oe_columns[:3])
    contacts_df = pd.merge(contacts_df, bait_ld_blocks_df, on = bait_columns[:3])

    # note the negation operator
    contacts_df['oe_shares_ld_block_with_bait'] = ~(
        contacts_df
        .apply(
            lambda x: x['oe_ld_block_names'].isdisjoint(x['bait_ld_block_names']),
            axis = 1
        )
    )

    return contacts_df


contacts_df = get_javierre_contacts()
contacts_df['bait_gene_names'] = (
    contacts_df
    ['bait_gene_names']
    .fillna('')
    .str.split(';', expand = False)
    .apply(set)
)
del contacts_df['oe_gene_names']

results = joblib.Parallel(-1)(
    joblib.delayed(get_accuracies)(_)
    for _ in tqdm(blood_cell_types, 'cell line')
)

stats_df = pd.DataFrame(
    results,
    columns = ['Cell Type', 'Closest Gene'] + [f'LD ({_})' for _ in super_pops]
)

stats_df.to_latex(
    'output/target_stats-table.tex',
    index = False,
    float_format = percent_formatter
)

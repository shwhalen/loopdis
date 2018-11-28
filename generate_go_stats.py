#!/usr/bin/env python

import feather
import itertools
import numpy as np
import os
import pandas as pd
import scipy.stats as stats
import sklearn.externals.joblib as joblib

from chromatics import *
from common import *
from tqdm import tqdm

"""
universe is all genes
2 tables per go term, one for closest gene, one for pir of baited gene

loop over all phenotypes
    loop over all go terms
        find genes closest to gwas snps
        find bait genes with gwas snps in pir

        is current gene closest to any gwas snp? y/n
        current gene annotated with current go term? y/n

        is current gene the bait of a pir (no really means NA, but treat as no)?  does pir contain a gwas snp? y/n
        current (bait) gene annotated with current go term? y/n

        sum of 2x2 table should be total # genes
"""
def get_pvalues(method, phenotype, go_term_chunk, genes, gene_hits):
    assert not (isinstance(genes, set) or isinstance(gene_hits, set))
    gene_misses = np.setdiff1d(genes, gene_hits)
    pvalues = []
    for go_term in go_term_chunk:
        go_term_genes = list(go_term_to_genes[go_term])
        gene_hit_go_hit_mask = np.in1d(gene_hits, go_term_genes)
        gene_miss_go_hit_mask = np.in1d(gene_misses, go_term_genes)

        table = [
            [sum(~gene_miss_go_hit_mask), sum(gene_miss_go_hit_mask)],
            [sum(~gene_hit_go_hit_mask), sum(gene_hit_go_hit_mask)]
        ]
        odds_ratio, pvalue = stats.fisher_exact(table)
        pvalues.append((method, phenotype, go_term, pvalue))
    return pvalues


def get_phenotype_stats(phenotype, distance):
    output_fn = f'output/go_stats/go_stats-{phenotype}-{distance}.feather'
    if os.path.exists(output_fn) and not overwrite:
        return

    current_gwas_df = (
        gwas_df
        .query('gwas_phenotype == @phenotype')
        .drop('gwas_phenotype', axis = 1)
    )

    # closest
    gwas_closest_genes_df = bedtools(
        'closest -d -t all',
        current_gwas_df,
        all_genes_df,
        genome = 'GRCh37'
    )
    if distance == 'proximal':
        gwas_closest_genes_df = gwas_closest_genes_df.query('distance < @proximal_cutoff_bp')
    else:
        gwas_closest_genes_df = gwas_closest_genes_df.query('distance >= @proximal_cutoff_bp')
    gwas_closest_genes = gwas_closest_genes_df['gene_name'].unique()

    # percent of snps where closest gene agrees with pchic bait
    gwas_closest_pchic_df = bedtools(
        'intersect -wa -wb',
        gwas_closest_genes_df,
        pchic_loops_df[oe_columns[:3] + ['bait_gene_names']],
        genome = 'GRCh37'
    )
    gwas_closest_pchic_df['bait_gene_names'] = (
        gwas_closest_pchic_df
        ['bait_gene_names']
        .str.split(';')
    )
    closest_is_bait = (
        gwas_closest_pchic_df
        .apply(lambda x: x['gene_name'] in x['bait_gene_names'], axis = 1)
        .sum()
    )
    if debug:
        print(f'{closest_is_bait / len(gwas_closest_pchic_df):.2f}')

    closest_pvalues = joblib.Parallel(n_jobs)(
        joblib.delayed(get_pvalues)('closest', phenotype, go_term_chunk, all_genes, gwas_closest_genes)
        for go_term_chunk in tqdm(go_term_chunks, 'go term (closest)')
    )

    # pchic
    gwas_oe_hits_df = bedtools(
        'intersect -wa -u',
        pchic_loops_df[oe_columns[:3] + ['bait_gene_names']],
        current_gwas_df,
        genome = 'GRCh37'
    )
    gwas_bait_genes = list(
        set(
            itertools.chain.from_iterable(
                gwas_oe_hits_df
                ['bait_gene_names']
                .str.split(';')
            )
        )
    )
    gwas_bait_genes = np.intersect1d(all_genes, gwas_bait_genes)

    pchic_pvalues = joblib.Parallel(n_jobs)(
        joblib.delayed(get_pvalues)('pchic', phenotype, go_term_chunk, all_genes, gwas_bait_genes)
        for go_term_chunk in tqdm(go_term_chunks, 'go term (pchic)')
    )

    # same ld block
    gwas_ld_df = bedtools(
        'intersect -wa -u',
        ld_blocks_df,
        current_gwas_df,
        genome = 'GRCh37'
    )
    gwas_ld_genes_df = bedtools(
        'intersect -wa -wb',
        gwas_ld_df,
        all_genes_df,
        genome = 'GRCh37'
    )
    gwas_ld_genes = gwas_ld_genes_df['gene_name'].unique()

    ld_block_pvalues = joblib.Parallel(n_jobs)(
        joblib.delayed(get_pvalues)('ld_block', phenotype, go_term_chunk, all_genes, gwas_ld_genes)
        for go_term_chunk in tqdm(go_term_chunks, 'go term (ld block)')
    )

    stats_df = pd.DataFrame(
        list(
            itertools.chain.from_iterable(
                closest_pvalues + ld_block_pvalues + pchic_pvalues
            )
        ),
        columns = ['method', 'phenotype', 'go_term', 'pvalue']
    )
    stats_df.to_feather(output_fn)


n_jobs = 12
overwrite = True
go_terms_fn = 'output/go_stats/go_terms.h5'
debug = False
proximal_cutoff_bp = 5000

if not os.path.exists(go_terms_fn):
    generate_go_terms(go_terms_fn)

gene_to_go_terms = pd.read_hdf(go_terms_fn, 'gene_to_go_terms')
go_term_to_genes = pd.read_hdf(go_terms_fn, 'go_term_to_genes')
go_term_chunks = np.array_split(go_term_to_genes.index.tolist(), n_jobs * 10)
assert go_term_to_genes.index.str.startswith('GO:').all()

all_genes_df = pd.read_hdf(go_terms_fn, 'genes')
all_genes = all_genes_df['gene_name'].unique()

ld_blocks_df = (
    get_plink_ld_blocks(None, 'EUR')
    .drop('ld_block_nsnps', axis = 1)
)

pchic_loops_df = (
    get_javierre_loops()
    .dropna(subset = ['bait_gene_names'])
    [oe_columns[:3] + bait_columns[:3] + ['bait_gene_names']]
)

gwas_df = feather.read_dataframe('gwas/gwas_catalog.feather')
gwas_df['gwas_phenotype'] = (
    gwas_df
    ['gwas_phenotype']
    .str.lower()
    .str.replace(' ', '_')
    .str.replace('/', '_over_')
)

phenotypes = (
    gwas_df
    ['gwas_phenotype']
    .value_counts()
    .head(30)
    .index
    .tolist()
)
for distance in tqdm(['proximal', 'distal'], 'distance'):
    for phenotype in tqdm(phenotypes, 'phenotype'):
        get_phenotype_stats(phenotype, distance)

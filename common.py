import chromatics
import feather
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import scipy.stats as stats
import seaborn as sns

from glob import glob
from sklearn.externals.joblib import Parallel, delayed
from tqdm import tqdm

def distance_match(loops_df, nonloops_df, distance_bin_count = 5, negatives_per_positive = 3):
    # bin positives
    loops_df['distance_bin'], distance_bins = pd.qcut(
        get_interaction_distance(loops_df),
        q = distance_bin_count,
        retbins = True
    )

    # re-use positive bins for negatives
    nonloops_df['distance_bin'] = pd.cut(
        get_interaction_distance(nonloops_df),
        bins = distance_bins,
        include_lowest = True
    )

    # subsample negatives within each distance bin
    loop_bin_counts = loops_df['distance_bin'].value_counts()
    distance_matched_negatives = []
    for distance_bin, group in nonloops_df.groupby('distance_bin', as_index = False):
        n_samples = negatives_per_positive * loop_bin_counts.loc[distance_bin]
        distance_matched_negatives.append(
            group.sample(n_samples, random_state = 0)
        )
    nonloops_df = pd.concat(distance_matched_negatives, ignore_index = True)

    positive_counts = loops_df['distance_bin'].value_counts()
    negative_counts = nonloops_df['distance_bin'].value_counts()
    assert (positive_counts * negatives_per_positive).equals(negative_counts)
    del loops_df['distance_bin']
    del nonloops_df['distance_bin']

    return loops_df, nonloops_df


def fast_sort(df):
    order = np.lexsort(
        [df[col].values for col in reversed(list(df.columns))]
    )
    for col in list(df.columns):
        df[col] = df[col].values[order]


def generate_go_terms(output_fn):
    from pybiomart import Server

    dataset = (
        Server('http://grch37.ensembl.org')
        .marts['ENSEMBL_MART_ENSEMBL']
        .datasets['hsapiens_gene_ensembl']
    )

    results = (
        dataset
        .query(attributes = [
            'chromosome_name',
            'start_position',
            'end_position',
            'external_gene_name',
            'go_id'
        ])
        .rename(columns = {
            'Chromosome/scaffold name': 'gene_chrom',
            'Gene start (bp)': 'gene_start',
            'Gene end (bp)': 'gene_end',
            'Gene name': 'gene_name',
            'GO term accession': 'go_terms'
        })
        .dropna()
    )

    gene_to_go_terms = (
        results
        .groupby('gene_name')
        ['go_terms']
        .apply(set)
    )

    gene_to_go_terms.to_hdf(
        output_fn,
        'gene_to_go_terms',
        mode = 'w',
        complevel = 1
    )

    go_term_to_genes = (
        results
        .groupby('go_terms')
        ['gene_name']
        .apply(set)
        .rename_axis('go_term')
        .rename('gene_names')
    )

    go_term_to_genes.to_hdf(
        output_fn,
        'go_term_to_genes',
        mode = 'a',
        complevel = 1
    )

    str_chroms = list(map(str, chroms))
    (
        results
        [['gene_chrom', 'gene_start', 'gene_end', 'gene_name']]
        .drop_duplicates('gene_name')
        .query('gene_chrom in @str_chroms')
        .reset_index(drop = True)
        .to_hdf(
            output_fn,
            'genes',
            mode = 'a',
            complevel = 1
        )
    )


def get_closest_genes(contacts_df, genome):
    oe_closest_genes_df = chromatics.bedtools(
        'closest -t first',
        (
            contacts_df
            [oe_columns[:3]]
            .drop_duplicates()
        ),
        (
            get_genes()
            .rename(columns = lambda x: 'oe_closest_' + x)
        ),
        genome = genome
    )

    return pd.merge(
        contacts_df,
        oe_closest_genes_df,
        on = oe_columns[:3],
        how = 'left'
    )


# lead eqtl may not overlap pir, but one of the pir_overlapping_snp_ids will be.
# lead eqtl may overlap bait, but the hit is due to ld with pir_overlapping_snp_ids.
def get_eqtls():
    print('reading eqtls')

    bc_sheetname = 'B cells (Fairfax et al. 2012)'
    wb_sheetname = 'Whole Blood (Westra et al. 2013' # missing closing paren in excel file
    eqtl_columns = ['eqtl_gene_name', 'lead_eqtl_snp_id', 'lead_eqtl_snp_chr', 'lead_eqtl_snp_start', 'pir_overlapping_snp_ids', 'r2_with_lead_eqtl_snp'] + pir_columns + ['pir_contains_promoter', 'array_probe_id'] + bait_columns[:4]
    eqtls_df = pd.read_excel(
        'javierre/Supplemental/mmc2.xlsx',
        sheet_name = bc_sheetname,
        skiprows = 3,
        names = eqtl_columns
    )

    if eqtls_df.dtypes['bait_chr'] == int:
        print('\nfixing mislabeled westra eqtl bait columns')
        eqtls_df.rename(
            columns = {
                'bait_chr': 'bait_id',
                'bait_start': 'bait_chr',
                'bait_end': 'bait_start',
                'bait_id': 'bait_end'
            },
            inplace = True
        )

    eqtls_df['lead_eqtl_snp_chr'] = eqtls_df['lead_eqtl_snp_chr'].str.extract(r'(\d+|X|Y)', expand = False)
    eqtls_df['pir_chr'] = eqtls_df['pir_chr'].str.extract(r'(\d+|X|Y)', expand = False)
    eqtls_df['bait_chr'] = eqtls_df['bait_chr'].str.extract(r'(\d+|X|Y)', expand = False)
    eqtls_df['lead_eqtl_snp_end'] = eqtls_df['lead_eqtl_snp_start']
    return eqtls_df


def get_genes():
    genes_df = pd.read_csv(
        f'{db_dir}/biomart/hg19_pc_genes.tsv',
        sep = '\t',
        header = 0,
        names = ['gene_name', 'gene_start', 'gene_end', 'gene_chr', 'gene_id'],
        dtype = {'gene_chr': str}
    )
    genes_df = genes_df[['gene_chr', 'gene_start', 'gene_end', 'gene_name']]
    return genes_df.sort_values(genes_df.columns.tolist())


def get_interaction_distance(contacts_df):
    # needs to include fragments not just window, as rare javierre interactions have overlapping fragments
    # In [1]: x.eval('bait_start < oe_start and bait_end > oe_start').sum()
    # Out[1]: 3
    return contacts_df[['f1_end', 'f2_end']].max(axis = 1) - contacts_df[['f1_start', 'f2_start']].min(axis = 1)


def get_javierre_contacts():
    output_fn = f'javierre/Detected interactions/PCHi-C/PCHiC_peak_matrix_cutoff0.feather'
    if os.path.exists(output_fn):
        return feather.read_dataframe(output_fn)

    contacts_df = pd.read_csv(
        'javierre/Detected interactions/PCHi-C/PCHiC_peak_matrix_cutoff0.txt.gz',
        sep = '\t',
        header = 0,
        names = bait_columns + oe_columns + ['dist'] + blood_cell_types,
        dtype = {
            'bait_chr': str,
            'bait_start': int,
            'bait_end': int,
            'bait_id': int,
            'bait_gene_names': str,
            'oe_chr': str,
            'oe_start': int,
            'oe_end': int,
            'oe_id': int,
            'oe_gene_names': str,
            'dist': np.float32,
            blood_cell_types[0]: np.float32,
            blood_cell_types[1]: np.float32,
            blood_cell_types[2]: np.float32,
            blood_cell_types[3]: np.float32,
            blood_cell_types[4]: np.float32,
            blood_cell_types[5]: np.float32,
            blood_cell_types[6]: np.float32,
            blood_cell_types[7]: np.float32,
            blood_cell_types[8]: np.float32,
            blood_cell_types[9]: np.float32,
            blood_cell_types[10]: np.float32,
            blood_cell_types[11]: np.float32,
            blood_cell_types[12]: np.float32,
            blood_cell_types[13]: np.float32,
            blood_cell_types[14]: np.float32,
            blood_cell_types[15]: np.float32,
            blood_cell_types[16]: np.float32
        }
    )
    contacts_df.to_feather(output_fn)
    return contacts_df


def get_javierre_loops():
    return pd.read_csv(
        'javierre/Detected interactions/PCHi-C/PCHiC_peak_matrix_cutoff5.txt.gz',
        sep = '\t',
        header = 0,
        names = bait_columns + oe_columns + ['dist'] + blood_cell_types + ['cluster_id', 'cluster_post_prob'],
        dtype = {
            'bait_chr': str,
            'oe_chr': str
        }
    )


def get_plink_ld(chrom, super_pop):
    if chrom is None:
        return pd.concat([get_plink_ld(_, super_pop) for _ in chroms])

    output_fn = f'{output_dir}/plink/{chrom}-{super_pop}.ld.h5'
    if os.path.exists(output_fn):
        return pd.read_hdf(output_fn, 'ld')

    ld_df = pd.read_csv(
        f'{output_dir}/plink/{chrom}-{super_pop}.ld.gz',
        sep = r'\s+',
        header = 0,
        usecols = [1, 4, 6],
        names = [
            'snp1_start',
            'snp2_start',
            'r2'
        ],
        dtype = {
            'snp1_start': int,
            'snp2_start': int,
            'r2': np.float32
        }
    )

    fast_sort(ld_df)
    ld_df.set_index(['snp1_start', 'snp2_start'], inplace = True)
    ld_df = ld_df[~ld_df.index.duplicated()]
    ld_df.to_hdf(
        output_fn,
        'ld',
        format = 'table',
        mode = 'w',
        complevel = 1,
        complib = 'lzo'
    )
    return ld_df


def get_plink_ld_blocks(chrom, super_pop):
    if chrom is None:
        return pd.concat([get_plink_ld_blocks(_, super_pop) for _ in chroms], ignore_index = True)

    return pd.read_csv(
        f'{output_dir}/plink/{chrom}-{super_pop}.blocks.det',
        sep = r'\s+',
        header = 0,
        usecols = [0, 1, 2, 4, 5],
        names = [
            'ld_block_chr',
            'ld_block_start',
            'ld_block_end',
            'ld_block_nsnps',
            'ld_block_snps'
        ],
        dtype = {
            'ld_block_chr': str,
            'ld_block_start': int,
            'ld_block_end': int,
            'ld_block_nsnps': int,
            'ld_block_snps': str
        }
    )


def get_rao_contacts(chrom, cell_type, observed_over_expected = True):
    if chrom is None:
        return pd.concat([get_rao_contacts(_, cell_type, observed_over_expected) for _ in chroms])

    oe_flag = 'oe' if observed_over_expected else 'o'
    output_fn = f'{output_dir}/rao/contacts-{chrom}-{cell_type}-vcnorm-{oe_flag}.feather'
    if os.path.exists(output_fn):
        return feather.read_dataframe(output_fn)

    rao_prefix = f'{db_dir}/chromatin_organization/rao2014/contacts/{cell_type}/{contact_bin_size // 1000}kb_resolution_intrachromosomal/chr{chrom}/MAPQGE30/chr{chrom}_{contact_bin_size // 1000}kb'

    contacts_df = pd.read_csv(
        f'{rao_prefix}.RAWobserved.gz',
        sep = '\t',
        header = None,
        names = ['f1_start', 'f2_start', 'contact_value'],
        dtype = {'contact_value': np.float32}
    )

    normalization_fn = f'{rao_prefix}.VCnorm'
    normalization = list(map(float, open(normalization_fn).readlines()))
    contacts_df['contact_value'] = [
        count / (normalization[f1_start // contact_bin_size] * normalization[f2_start // contact_bin_size])
        for f1_start, f2_start, count in contacts_df.itertuples(index = False)
    ]

    if observed_over_expected:
        expected_fn = f'{rao_prefix}.VCexpected'
        expected = list(map(float, open(expected_fn).readlines()))
        contacts_df['contact_value'] = [
            count / expected[(f2_start - f1_start) // contact_bin_size]
            for f1_start, f2_start, count in contacts_df.itertuples(index = False)
        ]

    fast_sort(contacts_df)
    contacts_df.to_feather(output_fn)
    return contacts_df


def get_rao_domains(cell_type):
    return pd.read_csv(
        glob(f'{db_dir}/chromatin_organization/rao2014/domains/*{cell_type}*domainlist.txt.gz')[0],
        sep = '\t',
        header = 0,
        usecols = range(3),
        names = domain_columns,
        dtype = {'domain_chr': str}
    )


def get_rao_loops(cell_type):
    return pd.read_csv(
        glob(f'{db_dir}/chromatin_organization/rao2014/loops/*{cell_type}*looplist.txt.gz')[0],
        sep = '\t',
        header = 0,
        usecols = range(6),
        names = juicer_columns,
        dtype = {'f1_chr': str, 'f2_chr': str}
    )


def spearman(x, y):
    return stats.spearmanr(x, y)[0]
spearman.__name__ = r'$\rho$'


db_dir = os.environ['DB_DIR']
genomes_dir = os.environ['GENOMES_DIR']
output_dir = os.environ['LD_OUTPUT_DIR']

lead_eqtl_columns = [
    'lead_eqtl_snp_chr',
    'lead_eqtl_snp_start',
    'lead_eqtl_snp_end',
    'lead_eqtl_snp_id',
    'eqtl_gene_name'
]
pir_columns = [
    'pir_chr',
    'pir_start',
    'pir_end',
    'pir_id'
]
bait_columns = [
    'bait_chr',
    'bait_start',
    'bait_end',
    'bait_id',
    'bait_gene_names'
]
oe_columns = [
    'oe_chr',
    'oe_start',
    'oe_end',
    'oe_id',
    'oe_gene_names'
]
juicer_columns = [
    'f1_chr',
    'f1_start',
    'f1_end',
    'f2_chr',
    'f2_start',
    'f2_end'
]
domain_columns = [
    'domain_chr',
    'domain_start',
    'domain_end'
]

chroms = list(range(1, 23))
super_pops = [
    'AFR',
    'AMR',
    'EAS',
    'EUR',
    'SAS'
]
encode_cell_types = [
    'K562',
    'GM12878',
    'IMR-90',
    'NHEK',
    'HUVEC'
]
blood_cell_types = [
    'Mon',
    'Mac0',
    'Mac1',
    'Mac2',
    'Neu',
    'MK',
    'EP',
    'Ery',
    'FoeT',
    'nCD4',
    'tCD4',
    'aCD4',
    'naCD4',
    'nCD8',
    'tCD8',
    'nB',
    'tB'
]

max_interaction_distance = 2e6
maf = 0.05
contact_bin_size = 5000
strong_ld_threshold = 0.8
significant_pchic_threshold = 5

comma_formatter = lambda x: '{:,}'.format(x)
percent_formatter = lambda x: '{:.1%}'.format(x)

plt.rc('font', **{'family': 'sans-serif', 'sans-serif': ['Arial']})
sns.set(style = 'white')

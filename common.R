library(feather)
library(readr)

library(dplyr)

get_javierre_contacts <- function() {
  output_fn <- 'javierre/Detected interactions/PCHi-C/PCHiC_peak_matrix_cutoff0.feather'

  if (!file.exists(output_fn)) {
    contacts_df <- read_tsv(
      'javierre/Detected interactions/PCHi-C/PCHiC_peak_matrix_cutoff0.txt.gz',
      col_names = c(bait_columns, oe_columns, 'dist', blood_cell_types),
      skip = 1
    )
    write_feather(
      contacts_df,
      output_fn
    )
  }

  read_feather(output_fn)
}

get_javierre_loops <- function() {
  output_fn <- 'javierre/Detected interactions/PCHi-C/PCHiC_peak_matrix_cutoff5.feather'

  if (!file.exists(output_fn)) {
    loops_df <- read_tsv(
      'javierre/Detected interactions/PCHi-C/PCHiC_peak_matrix_cutoff5.txt.gz',
      col_names = c(bait_columns, oe_columns, 'dist', blood_cell_types, 'cluster_id', 'cluster_post_prob'),
      skip = 1
    )
    write_feather(
      loops_df,
      output_fn
    )
  }

  read_feather(output_fn)
}

get_plink_ld_blocks <- function(chrom, super_pop) {
  if (is.na(chrom)) {
    return (
      do.call(rbind, lapply(1:22, get_plink_ld_blocks, super_pop = super_pop))
    )
  }

  read_table(
    sprintf('%s/plink/%s-%s.blocks.det', output_dir, chrom, super_pop),
    col_names = c('ld_block_chr', 'ld_block_start', 'ld_block_end', 'ld_block_kb', 'ld_block_nsnps', 'ld_block_snps'),
    col_types = c('ciidic'),
    skip = 1
  ) %>%
  mutate(super_pop = super_pop)
}

get_rao_loops <- function(cell_line) {
  fns <- list.files(
    sprintf('%s/chromatin_organization/rao2014/loops', db_dir),
    glob2rx(sprintf('GSE63525_%s*.bed.gz', cell_line)),
    full.names = T
  )
  read_tsv(
    fns[1],
    col_names = c('chr1', 'start1', 'end1', 'chr2', 'start2', 'end2')
  )
}

theme_loopdis <- function(base_size = 14) {
  theme_minimal(base_size) +
  theme(
    axis.text = element_text(size = base_size),
    axis.title = element_text(size= base_size),
    axis.title.x = element_text(margin = margin(15, 0, 0, 0)),
    axis.title.y = element_text(margin = margin(0, 15, 0, 0)),
    legend.text = element_text(size = base_size),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(1.5, 'lines')
  )
}

db_dir <- Sys.getenv('DB_DIR')
output_dir <- Sys.getenv('LD_OUTPUT_DIR')

bait_columns <- c(
  'bait_chr',
  'bait_start',
  'bait_end',
  'bait_id',
  'bait_gene_names'
)
oe_columns <- c(
  'oe_chr',
  'oe_start',
  'oe_end',
  'oe_id',
  'oe_gene_names'
)
encode_cell_types = c(
  'K562',
  'GM12878',
  'IMR-90',
  'NHEK',
  'HUVEC'
)
blood_cell_types <- c(
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
)
rao_columns <- c(
  'chr1',
  'x1',
  'x2',
  'chr2',
  'y1',
  'y2',
  'color',
  'o',
  'e_bl',
  'e_donut',
  'e_h',
  'e_v',
  'fdr_bl',
  'fdr_donut',
  'fdr_h',
  'fdr_v',
  'num_collapsed',
  'centroid1',
  'centroid2',
  'radius'
)
gtex_columns <- c(
  'bin',
  'chrom',
  'start',
  'end',
  'name',
  'score',
  'target_id',
  'target',
  'distance',
  'max_effect',
  'effect_type',
  'max_pvalue',
  'exp_count',
  'exp_names',
  'exp_scores',
  'exp_pvalues',
  'exp_probs'
)
super_pops <- c(
  'AFR',
  'AMR',
  'EAS',
  'EUR',
  'SAS'
)

contact_bin_bp <- 5000
contact_bin_kb <- contact_bin_bp / 1000
max_interaction_bp <- 2000000
max_interaction_kb <- max_interaction_bp / 1000
strong_ld_threshold <- 0.8

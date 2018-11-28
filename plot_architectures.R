#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)
library(ggstance)
library(readr)
library(readxl)
library(scales)

source('common.R')

get_fairfax_eqtls <- function(fn, cell_label) {
  read_csv(fn) %>%
  mutate(
    probe_midpoint = (probe_start + probe_end) / 2,
    distance = BP - probe_midpoint,
    source = sprintf('%s eQTLs\n(Fairfax et al. 2012)', cell_label),
    architecture = 'Genetic'
  ) %>%
  select(
    distance,
    source,
    architecture
  )
}

get_gtex_eqtls <- function() {
  read_tsv(
    sprintf('%s/gtex/gtexEqtlCluster.txt.gz', db_dir),
    col_names = gtex_columns
  ) %>%
  mutate(
    distance = distance,
    source = 'Combined eQTLs\n(GTEx Consortium 2017)',
    architecture = 'Genetic'
  ) %>%
  select(
    distance,
    source,
    architecture
  )
}

get_pairwise_ld <- function(super_pop) {
  read_table(sprintf('%s/plink/%s-%s.ld.gz', output_dir, ld_chrom, super_pop)) %>%
    filter(R2 >= strong_ld_threshold) %>%
    mutate(
      distance = BP_A - BP_B,
      source = sprintf('Strong LD Pairs\n(1KG Phase 3, chr %s)', ld_chrom),
      architecture = 'Genetic',
    ) %>%
    select(
      distance,
      source,
      architecture
    )
}

process_javierre_loops <- function() {
  get_javierre_loops() %>%
  mutate(
    distance = pmax(bait_end, oe_end) - pmin(bait_start, oe_start),
    source = 'Chromatin Loops\n(Javierre et al. 2016)',
    architecture = 'Physical'
  ) %>%
  select(
    distance,
    source,
    architecture
  )
}

process_rao_loops <- function(cell_line) {
  get_rao_loops(cell_line) %>%
    mutate(
      distance = pmax(end1, end2) - pmin(start1, start2),
      source = sprintf('Chromatin Loops\n(Rao et al. 2014)', cell_line),
      architecture = 'Physical'
    ) %>%
    select(
      distance,
      source,
      architecture
    )
}

get_ld_blocks <- function() {
  fns <- list.files(
    sprintf('%s/plink', output_dir),
    glob2rx('*.blocks.det'),
    full.names = T
  )
  do.call(rbind, lapply(fns, read_table)) %>%
    mutate(
      distance = KB * 1000,
      source = 'LD Blocks\n(1KG Phase 3)',
      architecture = 'Genetic'
    ) %>%
    select(
      distance,
      source,
      architecture
    )
}

get_contact_domains <- function() {
  domain_col_types <- list(
    chr1 = 'c',
    x1 = 'i',
    x2 = 'i',
    chr2 = 'c',
    y1 = 'i',
    y2 = 'i',
    color = 'c',
    f1 = 'd',
    f2 = 'd',
    f3 = 'd',
    f4 = 'd',
    f5 = 'd'
  )
  fns <- list.files(
    sprintf('%s/chromatin_organization/rao2014/domains', db_dir),
    glob2rx('*domainlist.txt.gz'),
    full.names = T
  )
  do.call(rbind, lapply(fns, read_tsv, col_types = domain_col_types)) %>%
    mutate(
      distance = x2 - x1,
      source = 'Contact Domains\n(Rao et al. 2014)',
      architecture = 'Physical'
    ) %>%
    select(
      distance,
      source,
      architecture
    )
}

get_tads <- function() {
  fns <- list.files(
    sprintf('%s/chromatin_organization/dixon2012/domains', db_dir),
    glob2rx('*Combined.bed'),
    full.names = T
  )
  do.call(rbind, lapply(fns, read_tsv, col_names = c('chr', 'start', 'end'))) %>%
    mutate(
      distance = end - start,
      source = 'Topological Domains\n(Dixon et al. 2012)',
      architecture = 'Physical'
    ) %>%
    select(
      distance,
      source,
      architecture
    )
}

pdf(NULL)

ld_chrom <- 22
strong_ld_fn <- sprintf('%s/strong_ld-%s.feather', output_dir, ld_chrom)
if (!file.exists(strong_ld_fn)) {
  ld_df <- do.call(rbind, lapply(super_pops, get_pairwise_ld))
  write_feather(ld_df, strong_ld_fn)
}
ld_df <- read_feather(strong_ld_fn)

ld_blocks_df <- get_ld_blocks()
gtex_df <- get_gtex_eqtls()
fairfax_bcells_df <- get_fairfax_eqtls('fairfax/ng.2205-S2-bcell_cis.csv', 'B cell & Monocyte')
fairfax_monocytes_df <- get_fairfax_eqtls('fairfax/ng.2205-S2-mono_cis.csv', 'B cell & Monocyte')
fairfax_shared_df <- get_fairfax_eqtls('fairfax/ng.2205-S2-both_cis.csv', 'B cell & Monocyte')
contact_domains_df <- get_contact_domains()
rao_df <- do.call(rbind, lapply(encode_cell_types, process_rao_loops))
javierre_df <- process_javierre_loops()
tads_df <- get_tads()

all_df <- (
  rbind(
    ld_blocks_df,
    ld_df,
    gtex_df,
    fairfax_bcells_df,
    fairfax_monocytes_df,
    fairfax_shared_df,
    contact_domains_df,
    rao_df,
    javierre_df,
    tads_df
  ) %>%
  mutate(
    distance_kb = abs(distance / 1e3),
    source = factor(
      source,
      levels = c(
        unique(tads_df$source),
        unique(javierre_df$source),
        unique(rao_df$source),
        unique(contact_domains_df$source),
        unique(fairfax_bcells_df$source),
        unique(gtex_df$source),
        unique(ld_df$source),
        unique(ld_blocks_df$source)
      )
    )
  )
)

all_df %>%
  group_by(source) %>%
  summarize_at(
    vars(distance_kb),
    funs(min, median, max, sd)
  )

all_df %>%
  ggplot(
    aes(
      x = distance_kb,
      y = source,
      color = architecture
    )
  ) +
  geom_violinh(
    trim = T,
    scale = 'area',
    position = 'identity'
  ) +
  geom_boxploth(
    width = 0.08,
    outlier.shape = NA,
    show.legend = F
  ) +
  labs(
    x = expression(log[10] * '(Distance or Size (kilobases))'),
    y = 'Genomic Feature',
    color = 'Architecture'
  ) +
  scale_color_brewer(palette = 'Paired') +
  scale_fill_brewer(palette = 'Paired') +
  scale_x_log10(labels = comma) +
  scale_y_discrete(expand = expand_scale(0, 1)) +
  annotation_logticks(
    sides = 'b',
    short = unit(0.05, 'cm'),
    mid = unit(0.1, 'cm'),
    long = unit(0.2, 'cm')
  ) +
  theme_loopdis(12) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
    ),
    legend.position = 'right'
  )

ggsave('output/architectures.pdf', width = 6, height = 5)

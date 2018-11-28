#!/usr/bin/env Rscript

library(dplyr)
library(feather)
library(foreach)
library(ggplot2)
library(Matrix)
library(reshape2)

source('common.R')

plot_hic_scaling <- function(x_max, bins, output_fn) {
  pdf(NULL)
  contacts_df %>%
    ggplot(aes(distance_kb, contact_value)) +
    geom_hex(bins = bins) +
    guides(fill = F) +
    labs(
      x = 'Genomic Distance (kilobases)',
      y = expression(log[10] * '(Hi-C Contact Count)')
    ) +
    scale_fill_distiller(
      palette = 'Purples',
      na.value = '#f2f0f7',
      direction = 1
    ) +
    xlim(0, x_max) +
    theme_loopdis() +
    theme(
      axis.text.x = element_text(
        angle = 90,
        vjust = 0.5
      )
    )
  ggsave(output_fn, width = 4, height = 4)
}

cell_type <- commandArgs(T)[1]
observed_over_expected <- as.integer(commandArgs(T)[2])
oe_flag <- if (observed_over_expected) 'oe' else 'o'
normalization <- 'vcnorm'
sparse_or_dense <- 'sparse'
stopifnot(cell_type %in% encode_cell_types)

contacts_df <-
  foreach(chrom = 1:22, .combine = rbind) %do% {
    contacts_fn <- sprintf('%s/rao/contacts-%s-%s-%s-%s.feather', output_dir, chrom, cell_type, normalization, oe_flag)
    contacts_df <- read_feather(contacts_fn)
    if (sparse_or_dense == 'dense') {
      contacts_df <-
        sparseMatrix(
          contacts_df$f1_start / contact_bin_bp,
          contacts_df$f2_start / contact_bin_bp,
          x = contacts_df$contact_value
        ) %>%
        as.matrix %>%
        melt %>%
        setNames(c('f1_start', 'f2_start', 'contact_value')) %>%
        filter(f2_start > f1_start) %>%
        mutate(
          f1_start = f1_start * contact_bin_bp,
          f2_start = f2_start * contact_bin_bp,
          contact_value = contact_value + 1
        )
    }
    contacts_df %>%
      mutate(
        distance_kb = (f2_start - f1_start) / 1e3,
        contact_value = log10(contact_value)
      ) %>%
      filter(distance_kb <= max_interaction_kb) %>%
      select(-f1_start, -f2_start)
  }

plot_hic_scaling(max_interaction_kb, 20, sprintf('output/hic_scaling-%s-%s-%s-%s.pdf', cell_type, normalization, oe_flag, sparse_or_dense))
plot_hic_scaling(max_interaction_kb / 20, 20, sprintf('output/hic_scaling-%s-%s-%s-%s-zoomed.pdf', cell_type, normalization, oe_flag, sparse_or_dense))

#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)

source('common.R')

plot_hic_scaling <- function(x_max, bins, output_fn) {
  p <-
    contacts_df %>%
    ggplot(aes_string('interaction_distance', cell_type)) +
    geom_hex(bins = bins) +
    geom_hline(
      yintercept = log10(5),
      alpha = 0.5,
      linetype = 'dashed'
    ) +
    guides(fill = F) +
    labs(
      x = 'Genomic Distance (kilobases)',
      y = expression(log[10] * '(CHiCAGO Score)')
    ) +
    scale_fill_distiller(
      palette = 'Purples',
      na.value = '#f2f0f7',
      direction = 1
    ) +
    scale_y_log10(limits = c(0.05, 10)) +
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
stopifnot(cell_type %in% blood_cell_types)

contacts_df <-
  get_javierre_contacts() %>%
  mutate(
    interaction_end = pmax(bait_end, oe_end),
    interaction_start = pmin(bait_start, oe_start),
    interaction_distance = (interaction_end - interaction_start) / 1e3
  ) %>%
  filter(interaction_distance < max_interaction_kb)

plot_hic_scaling(max_interaction_kb, 20, sprintf('output/hic_scaling-%s.pdf', cell_type))

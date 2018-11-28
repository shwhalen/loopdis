#!/usr/bin/env Rscript

library(dplyr)
library(feather)
library(ggplot2)

source('common.R')

plot_ld_scaling <- function(x_max, bins, output_fn) {
  p <-
    max_ld_df %>%
    ggplot(
      aes(
        window_length,
        max_ld
      )
    ) +
    geom_hex(bins = bins) +
    guides(fill = F) +
    labs(
      x = 'Genomic Distance (kilobases)',
      y = expression('Max Pairwise LD (' ~ R^2 ~ ')')
    ) +
    scale_fill_distiller(
      palette = 'Greens',
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

plot_ld_scaling_facet <- function(x_max, bins, output_fn) {
  p <-
    max_ld_df %>%
    ggplot(
      aes(
        window_length,
        max_ld,
        fill = super_pop
      )
    ) +
    facet_wrap(
      ~ super_pop,
      nrow = 2,
      labeller = as_labeller(correlation_labels, label_parsed)
    ) +
    geom_hex(
      aes(alpha = ..count..),
      bins = 20
    ) +
    guides(alpha = F, fill = F) +
    labs(
      x = 'Genomic Distance (kilobases)',
      y = expression('Max Pairwise LD (' ~ R^2 ~ ')'),
      fill = 'Super Population'
    ) +
    scale_fill_brewer(palette = 'Set1') +
    xlim(0, 100) +
    theme_loopdis() +
    theme(
      axis.text.x = element_text(
        angle = 90,
        vjust = 0.5
      )
    )

  ggsave(output_fn, width = 8, height = 6)
}

max_ld_df <-
  read_feather(sprintf('%s/max_interaction_ld.feather', output_dir)) %>%
  mutate(window_length = (pmax(f1_start, f2_start) - pmin(f1_end, f2_end)) / 1e3)

correlation_df <-
  max_ld_df %>%
  group_by(super_pop) %>%
  summarize(correlation = round(cor(window_length, max_ld, method = 'spearman'), 2))

correlation_labels <- setNames(
  paste(correlation_df$super_pop, '~(rho==', correlation_df$correlation, ')'),
  correlation_df$super_pop
)

pdf(NULL)

plot_ld_scaling(
  2000,
  20,
  'output/interaction_ld_scaling.pdf'
)

plot_ld_scaling(
  100,
  20,
  'output/interaction_ld_scaling-zoomed.pdf'
)

plot_ld_scaling_facet(
  100,
  20,
  'output/interaction_ld_scaling-facets-zoomed.pdf'
)

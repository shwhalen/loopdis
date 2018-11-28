#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)
library(hexbin)
library(readr)
library(viridis)

source('common.R')

plot_ld_scaling <- function(x_max, bins, y_var, y_label, output_fn) {
  p <-
    ld_df %>%
    ggplot(aes_string('distance', y_var)) +
    geom_hex(bins = bins) +
    guides(fill = F) +
    labs(
      x = 'Genomic Distance (kilobases)',
      y = y_label
    ) +
    # scale_fill_viridis(option = 'C') +
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

plot_ld_scaling_facet <- function(x_max, bins, y_var, y_label, output_fn) {
  p <-
    ld_df %>%
    ggplot(aes_string('distance', y_var, fill = 'super_pop')) +
    facet_wrap(
      ~ super_pop,
      nrow = 2,
      labeller = as_labeller(correlation_labels, label_parsed)
    ) +
    geom_hex(
      aes(alpha = ..count..),
      bins = bins
    ) +
    guides(alpha = F, fill = F) +
    labs(
      x = 'Genomic Distance (kilobases)',
      y = y_label,
      fill = 'Super Population'
    ) +
    scale_fill_brewer(palette = 'Set1') +
    xlim(0, x_max) +
    theme_loopdis() +
    theme(
      axis.text.x = element_text(
        angle = 90,
        vjust = 0.5
      )
    )

  ggsave(output_fn, width = 7.5, height = 6)
}

chrom <- 14
r2_y_label <- expression('LD (' ~ R^2 ~ ')')
dprime_y_label <- expression('LD ( D\' )')

ld_df <-
  read_feather(sprintf('%s/plink/%s.ld.feather', output_dir, chrom)) %>%
  mutate(distance = distance / 1e3)
cat('loaded ld\n')

correlation_df <-
  ld_df %>%
  group_by(super_pop) %>%
  summarize(
    correlation = cor(distance, r2, method = 'spearman') %>% round(2)
  )
cat('correlations computed\n')

correlation_labels <- setNames(
  paste(correlation_df$super_pop, '~(rho==', correlation_df$correlation, ')'),
  correlation_df$super_pop
)

pdf(NULL)

plot_ld_scaling(
  2000,
  20,
  'r2',
  r2_y_label,
  sprintf('output/ld_scaling-chr%s-r2.pdf', chrom)
)

plot_ld_scaling_facet(
  2000,
  20,
  'r2',
  r2_y_label,
  sprintf('output/ld_scaling-chr%s-r2-facets.pdf', chrom)
)

plot_ld_scaling(
  100,
  20,
  'r2',
  r2_y_label,
  sprintf('output/ld_scaling-chr%s-r2-zoomed.pdf', chrom)
)

plot_ld_scaling_facet(
  100,
  20,
  'r2',
  r2_y_label,
  sprintf('output/ld_scaling-chr%s-r2-facets-zoomed.pdf', chrom)
)

# dprime
correlation_df <-
  ld_df %>%
  group_by(super_pop) %>%
  summarize(
    correlation = cor(distance, dprime, method = 'spearman') %>% round(2)
  )

correlation_labels <- setNames(
  paste(correlation_df$super_pop, '~(rho==', correlation_df$correlation, ')'),
  correlation_df$super_pop
)

plot_ld_scaling(
  2000,
  20,
  'dprime',
  dprime_y_label,
  sprintf('output/ld_scaling-chr%s-dprime.pdf', chrom)
)

#!/usr/bin/env Rscript

library(dplyr)
library(feather)
library(ggplot2)
library(scales)

source('common.R')

plot_ld_block_crossings <- function(dataset, y_breaks) {
  breaks <- c(0, 10, 50, 100, 500, 100000)
  labels <- c('0-10', '10-50', '50-100', '100-500', '> 500')

  counts_df <-
    read_feather(sprintf('%s/interaction_crossings-ld_blocks-%s.feather', output_dir, dataset)) %>%
    group_by(super_pop) %>%
    do(
      data.frame(
        count = hist(.$count, breaks = breaks, plot = F)$counts,
        bin = factor(seq(labels), labels = labels)
      )
    )

  pdf(NULL)
  counts_df %>%
    ggplot() +
    geom_col(
      aes(
        bin,
        count,
        fill = super_pop
      )
    ) +
    guides(fill = F) +
    facet_wrap(~ super_pop, scales = 'free_x') +
    labs(
      x = '# Crossed LD Blocks',
      y = expression(log[10] * '(Count)')
    ) +
    scale_fill_brewer(palette = 'Set1') +
    scale_y_continuous(
      trans = 'log10',
      breaks = y_breaks,
      labels = scientific
    ) +
    theme_loopdis() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5),
    )
  ggsave(
    sprintf('output/interaction_crossings-ld_blocks-%s.pdf', dataset),
    width = 5,
    height = 5.5
  )
}

plot_rao_domain_crossings <- function(output_fn) {
  breaks <- c(0, 1, 2, 3, 4, 5, 6, 7)
  labels <- c(1, 2, 3, 4, 5, 6, 7)

  counts_df <-
    read_feather(sprintf('%s/interaction_crossings-domains-rao.feather', output_dir)) %>%
    group_by(super_pop) %>%
    do(
      data.frame(
        count = hist(.$count, breaks = breaks, plot = F)$counts,
        bin = factor(seq(labels))
      )
    )

  pdf(NULL)
  counts_df %>%
    ggplot() +
    geom_col(
      aes(
        x = bin,
        y = count,
        fill = super_pop
      )
    ) +
    guides(fill = F) +
    facet_wrap(~ super_pop, scales = 'free_x') +
    labs(
      x = '# Crossed Contact Domains',
      y = expression(log[10] * '(Count)')
    ) +
    scale_fill_brewer(palette = 'Set1') +
    scale_y_continuous(trans = 'log10') +
    theme_loopdis()
  ggsave(
    output_fn,
    width = 5,
    height = 4
  )
}

plot_ld_block_crossings('rao', c(1e0, 1e2, 1e4))
plot_ld_block_crossings('javierre', c(1e1, 1e3, 1e5))
plot_rao_domain_crossings('output/interaction_crossings-domains-rao.pdf')

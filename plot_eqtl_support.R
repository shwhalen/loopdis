#!/usr/bin/env Rscript

library(feather)
library(ggplot2)
library(reshape2)
library(tidyr)

library(dplyr)

source('common.R')

pdf(NULL)

read_feather('output/eqtl_support.feather') %>%
  melt(id.vars = 'distance_bin') %>%
  mutate(
    variable = factor(
      variable,
      levels = c('contacts_percent', 'same_ld_block_percent'),
      labels = c('PCHi-C Loops (Blood)', 'Same LD Block (EAS)')
    ),
    value = value * 100,
  ) %>%
  extract(
    distance_bin,
    into = c('distance_bin_start', 'distance_bin_end'),
    '(\\d+\\.\\d+), (\\d+\\.\\d+)',
    convert = T
  ) %>%
  mutate(
    distance_bin_start_kb = ceiling(distance_bin_start / 1000),
    distance_bin_end_kb = ceiling(distance_bin_end / 1000),
    distance_bin = factor(
      distance_bin_start,
      labels = paste0(distance_bin_start_kb, '-', distance_bin_end_kb, sep = '') %>% unique
    )
  ) %>%
  ggplot(aes(distance_bin, y = value, color = variable, group = variable)) +
  geom_line() +
  geom_point() +
  labs(x = 'Genomic Distance\n(kilobases)', y = 'B-cell eQTLs (%)', color = 'eQTL Support Type') +
  scale_color_brewer(palette = 'Set1') +
  theme_loopdis() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5)
  )
ggsave('output/eqtl_support.pdf', width = 5.5, height = 4)

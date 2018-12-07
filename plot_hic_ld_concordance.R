#!/usr/bin/env Rscript

library(dplyr)
library(feather)
library(ggplot2)
library(stringr)

source('common.R')

fns <- list.files(
  path = sprintf('%s/hic_ld_concordance', output_dir),
  glob2rx('hic_ld_concordance*.feather'),
  full.names = T
)
stats_df <-
  do.call(rbind, lapply(fns, read_feather)) %>%
  mutate(
    window_size = factor(window_size * 5), # convert bin count to kb
    cell_type = str_replace(cell_type, 'IMR90', 'IMR-90'),
    cell_type = factor(cell_type, levels = c('K562', 'GM12878', 'IMR-90', 'NHEK', 'HUVEC')),
    super_pop = factor(super_pop, levels = c('AFR', 'AMR', 'EAS', 'EUR', 'SAS')),
    percent_both_high = percent_both_high * 100,
    baseline = factor(baseline, levels = c(F, T), labels = c('Observed', 'Expected'))
  )

p <-
  stats_df %>%
  ggplot(aes(x = window_size, y = percent_both_high, color = baseline)) +
  facet_grid(cell_type ~ super_pop) +
  geom_boxplot(size = 0.1, outlier.size = 0.1) +
  labs(
    x = 'Window Size (kilobases)',
    y = 'Strong LD & Contact Frequency Co-Occurrence (Percent)',
    color = ''
  ) +
  scale_color_brewer(palette = 'Set1') +
  theme_loopdis() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      margin = margin(15, 0, 0, 0)
    ),
    legend.position = 'top',
    panel.spacing = unit(1, 'lines')
  )

ggsave(
  'output/hic_ld_concordance.pdf',
  width = 10.5,
  height = 10
)

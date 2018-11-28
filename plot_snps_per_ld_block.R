#!/usr/bin/env Rscript

library(ggplot2)
library(Hmisc)
library(viridis)
library(xtable)

library(dplyr)

source('common.R')

pdf(NULL)

blocks_df <-
  do.call(
    rbind,
    lapply(super_pops, get_plink_ld_blocks, chrom = NA)
  ) %>%
  select(ld_block_kb, ld_block_nsnps, super_pop) %>%
  group_by(super_pop) %>%
  mutate(
    ld_block_kb_bin = cut2(ld_block_kb, g = 4)
  ) %>%
  ungroup

blocks_df %>%
  group_by(super_pop, ld_block_kb_bin) %>%
  summarize(median(ld_block_nsnps)) %>%
  xtable(
    display = c('d', 's','s','d')
  ) %>%
  setNames(c('Super Pop.', 'LD Block Size Bin (kb)', 'Median # SNPs')) %>%
  print(
    file = 'output/snps_per_ld_block-table.tex',
    booktabs = T,
    floating = F,
    include.rownames = F
  )

blocks_df %>%
  ggplot(aes(ld_block_kb, ld_block_nsnps)) +
  facet_wrap(~ super_pop) +
  geom_hex(bins = 20) +
  labs(
    x = 'LD Block Size (kilobases)',
    y = '# SNPs per LD Block',
    fill = 'Count'
  ) +
  scale_fill_viridis(
    trans = 'log10',
    option = 'A'
  ) +
  theme_loopdis() +
  theme(
    panel.background = element_rect(fill = 'black')
  )

ggsave(
  'output/snps_per_ld_block.pdf',
  width = 7.5,
  height = 5
)

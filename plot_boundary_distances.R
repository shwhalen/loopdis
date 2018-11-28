#!/usr/bin/env Rscript

library(feather)
library(ggplot2)
library(viridis)

library(dplyr)

source('common.R')

plot_ld_blocks_to_domains <- function(distances_df, plot_fn) {
  distances_df %>%
    ggplot(aes(distance_kb, ld_block_size_kb)) +
    facet_wrap(~ super_pop) +
    geom_hex(bins = 20) +
    labs(
      x = 'Genomic Distance from LD Block\nto Closest Domain Boundary (kilobases)',
      y = 'LD Block Size (kilobases)',
      fill = 'Count'
    ) +
    scale_fill_viridis(
      trans = 'log10',
      option = 'A'
    ) +
    theme_loopdis() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5),
      panel.background = element_rect(fill = 'black')
    )

  ggsave(
    plot_fn,
    width = 7,
    height = 5
  )
}

read_domain_to_ld_block_distances <- function(super_pop) {
  read_feather(sprintf('output/boundary_distances-domains_to_ld_blocks-%s-%s.feather', cell_line, super_pop))
}

read_ld_block_to_domain_distances <- function(super_pop) {
  read_feather(sprintf('output/boundary_distances-ld_blocks_to_domains-%s-%s.feather', cell_line, super_pop))
}

cell_line <- 'GM12878'

pdf(NULL)

distances_df <-
  do.call(
    rbind,
    lapply(super_pops, read_ld_block_to_domain_distances)
  ) %>%
  mutate(
    distance_kb = distance / 1000,
    ld_block_size_kb = (ld_block_end - ld_block_start) / 1000
  )

plot_ld_blocks_to_domains(
  distances_df,
  sprintf('output/boundary_distances-ld_blocks_to_domains-%s.pdf', cell_line)
)

plot_ld_blocks_to_domains(
  distances_df %>% filter(distance_kb < 1, ld_block_size_kb < 20),
  sprintf('output/boundary_distances-ld_blocks_to_domains-%s-zoomed.pdf', cell_line)
)

distances_df <-
  do.call(
    rbind,
    lapply(super_pops, read_domain_to_ld_block_distances)
  ) %>%
  group_by(cell_line, super_pop, permutation) %>%
  summarize(distance = median(distance)) %>%
  ungroup %>%
  mutate(
    linetype = 'Permuted',
    distance = distance / 1000
  )

unpermuted_distances_df <-
  distances_df %>%
  filter(permutation == 0) %>%
  select(cell_line, super_pop, distance) %>%
  mutate(linetype = 'Unpermuted')

distances_df %>%
  filter(permutation > 0) %>%
  ggplot(
    aes(
      x = distance,
      color = super_pop,
      linetype = linetype
    )
  ) +
  facet_wrap(~ super_pop) +
  geom_density() +
  geom_vline(
    data = unpermuted_distances_df,
    aes(
      xintercept = distance,
      color = super_pop,
      linetype = linetype
    )
  ) +
  guides(color = F) +
  labs(
    x = 'Median Genomic Distance from Domain Boundary\nto Closest Haplotype Breakpoint (kilobases)',
    y = 'Density',
    fill = 'Super Population',
    linetype = 'Distance Type'
  ) +
  scale_color_brewer(palette = 'Set1') +
  scale_fill_brewer(palette = 'Set1') +
  theme_loopdis()

ggsave(
  sprintf('output/boundary_distances-domains_to_ld_blocks-%s.pdf', cell_line),
  width = 7.5,
  height = 5
)

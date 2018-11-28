#!/usr/bin/env Rscript

library(feather)
library(ggplot2)
library(Hmisc)

library(dplyr)

source('common.R')

plot_domain_crossing_ld <- function(shuffle_suffix) {
  output_fn <- sprintf('output/domain_crossing_ld-%s-%s-%s.pdf', cell_type, super_pop, shuffle_suffix)

  domain_snps_df <-
    read_feather(sprintf('%s/domain_snps-%s-%s-%s.feather', output_dir, cell_type, super_pop, shuffle_suffix)) %>%
    filter(distance <= 50000) %>%
    mutate(
      bin = cut2(ceiling(distance / 1e3), g = 10)
    )

  pvals_df <-
    domain_snps_df %>%
    group_by(bin, intra) %>%
    sample_n(3300) %>%
    group_by(bin) %>%
    summarize(
      pval = wilcox.test(r2[intra], r2[!intra])$p.value
    )
  cat(shuffle_suffix, '\n')
  print(pvals_df)
  cat('\n')

  p <-
    domain_snps_df %>%
    mutate(
      intra = factor(
        intra,
        levels = c(F, T),
        labels = c('Inter-\nDomain', 'Intra-\nDomain')
      )
    ) %>%
    group_by(bin, intra) %>%
    summarize(
      mean = mean(r2),
      median = median(r2),
      se = sd(r2) / sqrt(n())
    ) %>%
    ggplot(
      aes(
        x = bin,
        color = intra,
        group = intra
      )
    ) +
    geom_pointrange(
      aes(
        y = mean,
        ymin = mean - se,
        ymax = mean + se
      ),
      size = 0.1
    ) +
    labs(
      x = 'Genomic Distance (kilobases)',
      y = expression('Mean LD (' ~ R^2 ~ ')'),
      color = 'SNP Pair Type'
    ) +
    scale_color_brewer(palette = 'Paired') +
    ylim(0, 0.45) +
    theme_loopdis() +
    theme(
      axis.text.x = element_text(angle = 90),
      legend.position = 'top'
    )

  ggsave(output_fn, width = 5, height = 5)
}

super_pop <- commandArgs(T)[1]
cell_type <- commandArgs(T)[2]

pdf(NULL)
set.seed(0)

plot_domain_crossing_ld('shuffled')
plot_domain_crossing_ld('unshuffled')

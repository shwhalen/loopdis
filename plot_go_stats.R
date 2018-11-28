#!/usr/bin/env Rscript

library(feather)
library(ggplot2)
library(stringr)

library(dplyr)

pdf(NULL)

distance <- 'distal'

fns <- list.files(
  'output/go_stats',
  glob2rx('go_stats*.feather'),
  full.names = T
)
stats_df <- do.call(
  rbind,
  lapply(fns, read_feather)
)

fdrs <- c(0.01, 0.05, 0.1)
labels_df <- data.frame(
  x = 95,
  y = -log10(fdrs) + 0.05,
  label = c('1% FDR', '5% FDR', '10% FDR')
)

for (current_phenotype in unique(stats_df$phenotype)) {
  title <- gsub('body_mass_index', 'bmi', current_phenotype)
  title <- gsub('joint_analysis_main_effects_and_physical_activity_interaction', 'joint_analysis', title)
  print(title)

  phenotype_stats_df <-
    stats_df %>%
    filter(phenotype == current_phenotype) %>%
    group_by(method) %>%
    mutate(qvalue = p.adjust(pvalue, 'BH')) %>%
    arrange(qvalue) %>%
    mutate(rank = row_number()) %>%
    ungroup() %>%
    mutate(
      method = factor(
        method,
        levels = c('closest', 'ld_block', 'pchic'),
        labels = c('Closest Gene', 'LD Block', 'PCHi-C')
      ),
      pvalue = -log10(pvalue),
      qvalue = -log10(qvalue)
    ) %>%
    filter(rank < 50)

  phenotype_stats_df %>%
    ggplot() +
    geom_line(
      aes(
        x = rank,
        y = qvalue,
        color = method,
        group = method
      )
    ) +
    geom_hline(
      yintercept = -log10(fdrs),
      color = 'black',
      alpha = 0.5,
      size = 0.5,
      linetype = 'dotted'
    ) +
    labs(
      x = 'GO Term (sorted by q-value)',
      y = '-log10(q-value)',
      color = 'Method'
    ) +
    scale_color_brewer(palette = 'Set1') +
    theme_minimal() +
    theme(
      legend.margin = margin(0, 0, 0, 0),
      legend.position = 'top',
      legend.title = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_blank()
    )

  output_fn <- sprintf('output/go_stats/go_stats-%s-qvalues.pdf', current_phenotype)
  ggsave(output_fn, width = 4.5, height = 4)
}

#!/usr/bin/env Rscript

library(feather)
library(ggplot2)
library(stringr)

library(dplyr)

pdf(NULL)

# note that proximal and distal is only stored separately for the closest gene method.
# to combine everything, keep method = closest rows from proximal files, and keep all rows from distal files.
# or vice versa.
proximal_fns <- list.files(
  'output/go_stats',
  glob2rx('go_stats*-proximal.feather'),
  full.names = T
)
proximal_stats_df <-
  do.call(
    rbind,
    lapply(proximal_fns, read_feather)
  ) %>%
  filter(method == 'closest') %>%
  mutate(method = str_replace(method, 'closest', 'closest_proximal'))

distal_fns <- list.files(
  'output/go_stats',
  glob2rx('go_stats*-distal.feather'),
  full.names = T
)
distal_stats_df <-
  do.call(
    rbind,
    lapply(distal_fns, read_feather)
  ) %>%
  mutate(method = str_replace(method, 'closest', 'closest_distal'))

stats_df <- rbind(proximal_stats_df, distal_stats_df)
print(unique(stats_df$method))

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
        levels = c('closest_proximal', 'closest_distal', 'ld_block', 'pchic'),
        labels = c('Closest Gene (< 5kb)', 'Closest Gene (> 5kb)', 'LD Block', 'PCHi-C')
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
    theme_minimal(base_size = 14) +
    theme(
      legend.margin = margin(0, 0, 0, 0),
      legend.position = 'right',
      legend.title = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 10),
      axis.text.x = element_blank()
    )

  output_fn <- sprintf('output/go_stats/go_stats-%s-proximal_vs_distal-qvalues.pdf', current_phenotype)
  ggsave(output_fn, width = 6, height = 4)
}

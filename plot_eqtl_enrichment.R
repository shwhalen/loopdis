#!/usr/bin/env Rscript

library(feather)
library(ggplot2)
library(Hmisc)
library(reshape2)
library(scales)

library(dplyr)

source('common.R')

get_log_odds_ratio <- function(x, y) {
  fit <- glm(
    y ~ x,
    data = data.frame(x = as.numeric(x), y = y),
    family = 'binomial'
  )
  fit_summary <- summary(fit)
  ci <- confint(fit)
  list(
    log_odds_ratio = fit_summary$coefficients['x', 'Estimate'],
    pvalue = fit_summary$coefficients['x', 'Pr(>|z|)'],
    confint_lower = ci[2, 1],
    confint_upper = ci[2, 2]
  )
}

pdf(NULL)

alpha <- 1e-4

eqtl_stats_df <-
  read_feather(sprintf('%s/eqtl_enrichment.feather', output_dir)) %>%
  mutate(
    interaction_length_mb = (pmax(oe_end, bait_end) - pmin(oe_start, bait_start)) / 1e6
  ) %>%
  filter(interaction_length_mb < 2) %>%
  mutate(
    interaction_bin = cut2(
      interaction_length_mb,
      cuts = c(
        min(interaction_length_mb),
        0.2,
        max(interaction_length_mb)
      ),
      digits = 1,
      minmax = F
    ),
    interaction_bin = factor(interaction_bin, labels = c('0-0.2', '0.2-2'))
  )
table(eqtl_stats_df$interaction_bin)

eqtl_stats_df %>%
  group_by(interaction_bin) %>%
  count(is_significant, oe_contains_eqtl)

eqtl_stats_df %>%
  group_by(interaction_bin) %>%
  count(is_bait_closest_gene_to_oe, oe_contains_eqtl)

eqtl_stats_df %>%
  group_by(interaction_bin) %>%
  count(
    max_ld_AFR >= strong_ld_threshold |
    max_ld_AMR >= strong_ld_threshold |
    max_ld_EAS >= strong_ld_threshold |
    max_ld_EUR >= strong_ld_threshold |
    max_ld_SAS >= strong_ld_threshold,
    oe_contains_eqtl
  )

chromatin_interactions_df <-
  eqtl_stats_df %>%
  group_by(interaction_bin) %>%
  do(
    data.frame(
      type = '(a)\nStatistically Significant\nChromatin Interactions',
      get_log_odds_ratio(.$is_significant, .$oe_contains_eqtl)
    )
  ) %>%
  ungroup

closest_gene_df <-
  eqtl_stats_df %>%
  group_by(interaction_bin) %>%
  do(
    data.frame(
      type = '(b)\nClosest Gene\nInside Bait Fragments',
      get_log_odds_ratio(.$is_bait_closest_gene_to_oe, .$oe_contains_eqtl)
    )
  ) %>%
  ungroup

strong_ld_df <-
  eqtl_stats_df %>%
  group_by(interaction_bin) %>%
  do(
    data.frame(
      type = '(c)\nStrong LD\nBetween Fragments',
      get_log_odds_ratio(
        .$max_ld_AFR >= strong_ld_threshold |
          .$max_ld_AMR >= strong_ld_threshold |
          .$max_ld_EAS >= strong_ld_threshold |
          .$max_ld_EUR >= strong_ld_threshold |
          .$max_ld_SAS >= strong_ld_threshold,
        .$oe_contains_eqtl
      )
    )
  ) %>%
  ungroup

rbind(chromatin_interactions_df, closest_gene_df, strong_ld_df) %>%
  ggplot(
    aes(
      interaction_bin,
      log_odds_ratio,
      group = type,
      shape = pvalue < alpha
    )
  ) +
  facet_wrap(~ type, scales = 'fixed') +
  geom_errorbar(
    aes(
      x = interaction_bin,
      ymin = confint_lower,
      ymax = confint_upper
    ),
    width = 0.2,
    color = 'grey'
  ) +
  geom_point(size = 3) +
  geom_hline(
    yintercept = 0,
    linetype = 'dotted'
  ) +
  labs(
    x = 'Genomic Distance (megabases)',
    y = 'Log Odds Ratio for eQTL Enrichment',
    shape = sprintf('p-value < %.0e', alpha)
  ) +
  theme_loopdis() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    legend.position = 'top'
  )

ggsave('output/eqtl_enrichment.pdf', width = 6.5, height = 5)

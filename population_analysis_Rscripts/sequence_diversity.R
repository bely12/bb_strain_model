# for measuring sequence diversity
library(tidyverse)
library(broom)

setwd('/Users/brandonely/Desktop/bb_strain_model_dev/evo_sim/results_tables/sim_output/immune_results/')
df1 <- read.table('immune_sim_run_stats.tsv', header = T, sep = '\t')

setwd('/Users/brandonely/Desktop/bb_strain_model_dev/evo_sim/results_tables/sim_output/no_selection_results/')
df2 <- read.table('neutral_sim_run_stats.tsv', header = T, sep = '\t')

setwd('/Users/brandonely/Desktop/bb_strain_model_dev/evo_sim/results_tables/sim_output/hybrid_results/')
df3 <- read.table('hybrid_ExpMod_sim_run_stats.tsv', header = T, sep = '\t')

df <- bind_rows(df1, df2, df3)
rm(df1, df2, df3)

df2 <- df %>%
  filter(year == 500 & rec_rate == 0.01)

df2 <- df2 %>% select(selection, infection_rate, avg_antigen_distance)

df2_summary <- df2 %>%
  group_by(selection) %>%
  summarise(
    mean_pw_dist = mean(avg_antigen_distance),
    dist_sd = sd(avg_antigen_distance)
  )

df2 %>%
  ggplot(aes(x = selection, y = avg_antigen_distance, fill = selection)) +
  geom_boxplot(outliers = F) +
  geom_jitter() +
  theme_bw()

setwd('/Users/brandonely/Desktop/bb_strain_model_dev/evo_sim/results_tables/new_plots/')
#ggsave('mean_pw_dists.png', height = 8, width = 6)

# ANOVA
anova_df <- df2 %>%
  aov(avg_antigen_distance ~ selection, data = .) %>%
  tidy()
#write.table(anova_df, 'anova_res_mean_pw_dist.tsv', quote = F, row.names = F, sep = '\t')


# Tukey HSD
tukey_df <- df2 %>%
  aov(avg_antigen_distance ~ selection, data = .) %>%
  TukeyHSD() %>%
  tidy()
#write.table(tukey_df, 'tukey_res_mean_pw_dist.tsv', quote = F, row.names = F, sep = '\t')






###### per site entropy

setwd('/Users/brandonely/Desktop/bb_strain_model_dev/evo_sim/results_tables/sim_output/immune_results/')
df1 <- read.table('immune_variant_freqs_sampled.tsv', header = T, sep = '\t', colClasses = c('character','numeric','numeric','character', 'character', 'character','numeric','numeric','character'))
df1 <- df1 %>% filter(rec_rate == 0.01 & freq >= 1)
df1 <- df1 %>% select(strain, run_id, selection)

setwd('/Users/brandonely/Desktop/bb_strain_model_dev/evo_sim/results_tables/sim_output/no_selection_results/')
df2 <- read.table('neutral_variant_freqs_sampled.tsv', header = T, sep = '\t', colClasses = c('character','numeric','numeric','character', 'character', 'character','numeric','numeric','character'))
df2 <- df2 %>% filter(rec_rate == 0.01 & freq >= 1)
df2 <- df2 %>% select(strain, run_id, selection)

setwd('/Users/brandonely/Desktop/bb_strain_model_dev/evo_sim/results_tables/sim_output/hybrid_results/')
df3 <- read.table('hybrid_ExpMod_variant_freqs_sampled.tsv', header = T, sep = '\t', colClasses = c('character','numeric','numeric','character', 'character', 'character','numeric','numeric','character'))
df3 <- df3 %>% filter(freq >= 1)
df3 <- df3 %>% select(strain, run_id, selection)


# define entropy function for binary bits
entropy_bin <- function(x) {
  p <- mean(x)
  if (p == 0 || p == 1) return(0)
  -p * log2(p) - (1 - p) * log2(1 - p)
}

# compute mean site-wise entropy per run_id
entropy_results <- df3 %>%
  group_by(run_id) %>%
  summarise(mean_entropy = {
    # split binary strings into bit columns
    bit_matrix <- strain %>%
      strsplit(split = "") %>%
      do.call(rbind, .) %>%
      apply(2, as.numeric) %>%
      as.data.frame()
    
    # compute entropy per site
    site_entropy <- apply(bit_matrix, 2, entropy_bin)
    
    # average across all sites
    mean(site_entropy)
  }) %>%
  ungroup()

imm_entropy <- entropy_results
mean(imm_entropy$mean_entropy)

no_entropy <- entropy_results
mean(no_entropy$mean_entropy)

hyb_entropy <- entropy_results
mean(hyb_entropy$mean_entropy)




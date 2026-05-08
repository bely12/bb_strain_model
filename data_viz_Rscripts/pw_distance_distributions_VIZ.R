library(tidyverse)
library(diptest)
library(e1071)
library(moments)
library(purrr)
library(broom)

##### gepurrr##### get pw distances #####
header <- c('run_id', 'distance', 'selection', 'gene', 'rec_rate', 'mut_rate', 'gen_fit')

# sample hybrid selection
setwd('/Users/brandonely/Desktop/bb_strain_model_dev/evo_sim/results_tables/sim_output/hybrid_results/')

dists <- read.delim('hybrid_ExpMod_pw_dists_sampled.tsv', header = T, col.names = header) 
dists$host_pop_stability <- 'stable'
#sampled_dists <- sample_n(dists, 1000)
sampled_dists <- dists %>%
  group_by(run_id) %>%
  slice_sample(n = 1000) %>%
  ungroup()
rm(dists)

# sample immune selection
setwd('/Users/brandonely/Desktop/bb_strain_model_dev/evo_sim/results_tables/sim_output/immune_results/')

dists <- read.delim('immune_pw_dists_sampled.tsv', header = T, col.names = header) 
dists <- dists %>% filter(rec_rate == 0.01)
dists$host_pop_stability <- 'stable'
#sampled_dists4 <- sample_n(dists, 1000)
sampled_dists2 <- dists %>%
  group_by(run_id) %>%
  slice_sample(n = 1000) %>%
  ungroup()
rm(dists)

# sample no selection 
setwd('/Users/brandonely/Desktop/bb_strain_model_dev/evo_sim/results_tables/sim_output/no_selection_results/')
dists <- read.delim('neutral_pw_dists_sampled.tsv', header = T, col.names = header) 
dists <- dists %>% filter(rec_rate == 0.01)
dists$host_pop_stability <- 'stable'
#sampled_dists5 <- sample_n(dists, 1000)
sampled_dists3 <- dists %>%
  group_by(run_id) %>%
  slice_sample(n = 1000) %>%
  ungroup()
rm(dists)

# concatenate
all_sampled_dists <- bind_rows(sampled_dists, sampled_dists2, sampled_dists3)
rm(list = setdiff(ls(), "all_sampled_dists"))

# add a category for distinct parameters
all_sampled_dists$parameters <- paste(all_sampled_dists$selection, all_sampled_dists$gen_fit, sep = '_')

# mean pw distance
df <- all_sampled_dists %>%
  group_by(selection, run_id) %>%
  summarise(
    mean_pw_dist = mean(distance))

# plot
df %>%
  ggplot(aes(x = selection, y = mean_pw_dist, fill = selection)) +
  geom_boxplot() +
  geom_jitter() +
  labs(title = 'Mean p-w distances', x = 'selection', y = 'mean distance per run') +
  theme_bw() +
  theme(plot.title = element_text(size = 25),
        axis.title.y = element_text(size = 25, face = 'bold'),
        axis.text = element_text(size = 25, face = 'bold'),
        axis.title.x = element_text(size = 25, face = 'bold'),
        legend.title = element_text(size = 25, face = 'bold'),
        legend.text = element_text(size = 25, face = 'bold'))
setwd('/Users/brandonely/Desktop/bb_strain_model_dev/evo_sim/results_tables/new_plots/')
#ggsave('mean_pw_antigen_distance_boxplot_DIS.png', height = 8, width = 8)  



# plot
all_sampled_dists %>%
  ggplot(aes(x = distance, fill = selection)) +
  geom_histogram(bins = 19, color = 'black') +
  #geom_density() +
  #geom_vline(xintercept = 3, linetype = 2) +
  #geom_vline(xintercept = 8, linetype = 2) +
  theme_bw() +
  labs(x = 'p-w distance') +
  facet_wrap(~selection, ncol = 1,
             labeller = as_labeller(
               c("hybrid" = "Hybrid", 
                 "immune" = "Immune", 
                 "neutral" = "No selection"))) +
  
  theme(axis.title.y = element_text(size = 25, face = 'bold'),
        axis.text = element_text(size = 25, face = 'bold'),
        axis.title.x = element_text(size = 25, face = 'bold'),
        legend.text = element_text(size = 25),
        legend.title = element_text(size = 25),
        strip.text = element_text(size = 25))

setwd('/Users/brandonely/Desktop/bb_strain_model_dev/evo_sim/results_tables/new_plots/')
#ggsave('all_pw_antigen_distance_histogram_DIS.png', height = 8, width = 6)  

df <- all_sampled_dists %>% select(distance, selection, run_id)

# define lower, intermediate, high ranges
low_range <- c(0, 3)    # distances before gap
mid_range <- c(4, 8) # intermediate / gap
high_range <- c(9,20) # high distances

mid_frac_norm <- df %>%
  group_by(selection, run_id) %>%
  summarise(
    n_low = sum(distance >= low_range[1] & distance <= low_range[2]),
    n_mid = sum(distance >= mid_range[1] & distance <= mid_range[2]),
    n_high = sum(distance >= high_range[1] & distance <= high_range[2]),
    mid_frac_rel = n_mid / (n_low + n_mid + n_high),
    .groups = "drop"
  )

#print(mid_frac_norm)

# ANOVA
anova_df <- mid_frac_norm %>%
  aov(mid_frac_rel ~ selection, data = .) %>%
  tidy()

#write.table(anova_df, 'anova_res_intermediate_gap.tsv', quote = F, row.names = F, sep = '\t')

# Tukey HSD
tukey_df <- mid_frac_norm %>%
  aov(mid_frac_rel ~ selection, data = .) %>%
  TukeyHSD() %>%
  tidy()
#write.table(tukey_df, 'tukey_res_intermediate_gap.tsv', quote = F, row.names = F, sep = '\t')

mid_frac_norm %>% 
  ggplot(aes(x = selection, y = mid_frac_rel, fill = selection)) +
  #geom_col(width = 0.6) +
  geom_boxplot(outliers = F) +
  geom_jitter() +
  #geom_text(aes(label = round(mid_frac_rel, 2)), vjust = -0.5, size = 5) +
  #scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = "Selection",
    y = "relative frequency",
    title = "Frequncy of Intermediate Genotypes"
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        plot.title = element_text(size = 25),
        axis.title.y = element_text(size = 25, face = 'bold'),
        axis.text = element_text(size = 25, face = 'bold'),
        axis.title.x = element_text(size = 25, face = 'bold'))
        #legend.text = element_text(size = 25),
        #legend.title = element_text(size = 25),
        #strip.text = element_text(size = 25)))

setwd('/Users/brandonely/Desktop/bb_strain_model_dev/evo_sim/results_tables/new_plots/')
#ggsave('intermediate_genotype_rel_freq_DIS.png', height = 8, width = 8)



######

# function to compute bimodality coefficient
compute_bc <- function(x) {
  n <- length(x)
  g <- skewness(x)
  k <- kurtosis(x)
  BC <- (g^2 + 1) / (k + 3 * ((n - 1)^2) / ((n - 2) * (n - 3)))
  return(BC)
}

# compute
bimodality_metrics <- df %>%
  group_by(selection, run_id) %>%
  summarise(
    dip_stat = dip.test(distance)$statistic,
    dip_p = dip.test(distance)$p.value,
    skewness = skewness(distance),
    kurtosis = kurtosis(distance),
    BC = compute_bc(distance),
    .groups = "drop"
  )

print(bimodality_metrics)

# ANOVA
anova_df <- bimodality_metrics %>%
  aov(BC ~ selection, data = .) %>%
  tidy()
#write.table(anova_df, 'anova_res_BimodCoeff.tsv', quote = F, row.names = F, sep = '\t')

# Tukey HSD
tukey_df <- bimodality_metrics %>%
  aov(BC ~ selection, data = .) %>%
  TukeyHSD() %>%
  tidy()
#write.table(tukey_df, 'tukey_res_BimodCoeff.tsv', quote = F, row.names = F, sep = '\t')


bimodality_metrics %>%
  ggplot(aes(x = selection, y = BC, fill = selection)) +
  geom_boxplot(outliers = F) +
  geom_jitter() +
  labs(title = 'Bimodality of p-w distance distributions',y = 'bimodiality coefficient') +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(size = 25),
        axis.title.y = element_text(size = 25, face = 'bold'),
        axis.text = element_text(size = 25, face = 'bold'),
        axis.title.x = element_text(size = 25, face = 'bold'))

setwd('/Users/brandonely/Desktop/bb_strain_model_dev/evo_sim/results_tables/new_plots/')
#ggsave('bimodality_coefficient_DIS.png', height = 8, width = 8)



  
  
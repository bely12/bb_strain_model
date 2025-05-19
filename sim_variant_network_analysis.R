# network analysis
setwd('/Users/brandonely/Desktop/bb_strain_model_dev/R_analysis_scripts/')
library(tidyverse)
library(ggpubr)
library(igraph)
source("sim_analysis_functions.R")

# upload data
setwd('/Users/brandonely/Desktop/bb_strain_model_dev/sim_results/v5000_yrs500_results/')
header <- c('variant', 'counts', 'frequency', 'run_id')
immune <- read.delim('immune_5000_20_500_25SET_sampled_variant_frequencies.tsv', header = T, col.names = header, colClasses = c('character','numeric','numeric','character'))
adaptive <- read.delim('hs_5000_20_500_25SET_sampled_variant_frequencies.tsv', header = T, col.names = header, colClasses = c('character','numeric','numeric','character'))
neutral <- read.delim('neutral_5000_20_500_25SET_sampled_variant_frequencies.tsv', header = T, col.names = header, colClasses = c('character','numeric','numeric','character'))


# look at network stats for a single run and single threshold
# df <- data %>% filter(run_id == 'run_20')
# network <- antigen_network(df, threshold = 3)

# visualize a network
plot_single <- neutral %>% filter(run_id == 'run_15')
#png("neutral_network_run20_thresh4.png", width = 800, height = 600)
antigen_network_viz(plot_single, threshold = 4)
#dev.off()

# look at network stats for a single run and range of thresholds
catcher = list()
index = 1
for (i in 1:5) {
  catcher[[index]] <- antigen_network(df, threshold = i)
  index = index + 1
network <- bind_rows(catcher)
}

# look at network stats for all runs and a range of thresholds 
catcher = list()
index = 1
for (run in unique(neutral$run_id)) {
  temp_df <- neutral %>% filter(run_id == run)
  for (i in 1:5) {
    catcher[[index]] <- antigen_network(temp_df, threshold = i)
    index = index + 1
  }
  neutral_network <- bind_rows(catcher)
}
neutral_network$selection <- 'neutral'
networks_df <- bind_rows(immune_network, adaptive_network, neutral_network)
#write.table(networks_df, file='sim_networks_analysis_not_sampled.tsv', quote=FALSE, sep='\t', row.names = F)


# network analysis results
setwd('/Users/brandonely/Desktop/bb_strain_model_dev/sim_results/v5000_yrs500_results/')
df <- read.delim('sim_networks_analysis_sampled.tsv', header = T)
df <- df %>% select(-6,-7)

# summary
df_summary <- df %>%
  group_by(selection, threshold) %>%
  summarise(edge_density = mean(edge_dens),
            edge_density_sd = sd(edge_dens),
            mean_degrees = mean(mean_degree),
            mean_degrees_sd = sd(mean_degree),
            cluster_coeff = mean(clust_coef),
            cluster_coefff_sd = sd(clust_coef),
            #diam = mean(diam),
            #mean_dist = mean(mean_dist),
            mod = mean(modul),
            mod_sd = sd(modul))


typical_runs <- df %>% filter(
  edge_dens > (df_summary$edge_density - df_summary$edge_density_sd[1]*5) & 
  edge_dens < (df_summary$edge_density[1] + df_summary$edge_density_sd[1]*5) &
  mean_degree > (df_summary$mean_degrees[1] - df_summary$mean_degrees_sd[1]*5) & 
  mean_degree < (df_summary$mean_degrees[1] + df_summary$mean_degrees_sd[1]*5) &
  clust_coef > (df_summary$cluster_coeff[1] - df_summary$cluster_coefff_sd[1]*5) &
  clust_coef < (df_summary$cluster_coeff[1] + df_summary$cluster_coefff_sd[1]*5) &
  modul > (df_summary$mod[1] - df_summary$mod_sd[1]*5) &
  modul < (df_summary$mod[1] + df_summary$mod_sd[1]*5))




# wide to long
df <- df %>% select(-6,-7)
df_long <- pivot_longer(df, cols = c(3:6), names_to = 'variable', values_to = 'values')

# plot line of means
ggplot(df_long, aes(x = threshold, y = values, color = selection, group = selection)) +
  geom_line(stat = "summary", fun = mean, size = 1.2) +     
  geom_point(stat = "summary", fun = mean, size = 2) +       
  facet_wrap(~ variable, scales = "free_y") +               
  theme_bw()

#ggsave('network_analysis_line_means.png', width = 8, height = 8)

# box plots
ggplot(df_long, aes(x = factor(threshold), y = values, fill = selection)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, position = position_dodge(width = 1.0)) + 
  geom_jitter(aes(color = selection), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), size = 1.5, alpha = 0.6) +
  stat_compare_means(aes(group = selection), method = "anova", label = "p.format") +
  facet_wrap(~ variable, scales = "free_y") +
  theme_bw()

ggsave('network_analysis_boxplots.png', width = 12, height = 8)


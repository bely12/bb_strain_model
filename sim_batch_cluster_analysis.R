setwd('/Users/brandonely/Desktop/bb_strain_model_dev/')
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(stringdist)
library(cluster)
source("sim_analysis_functions.R")

# upload data 
setwd('/Users/brandonely/Desktop/bb_strain_model_dev/new_sim_output/test_set_500yrs/')
cr_data <- read.delim('cr_len20_vec1000_yrs500_variant_frequencies.tsv', header = T, colClasses = c('character','numeric','numeric','character'))
hs_data <- read.delim('hs_len20_vec1000_yrs500_variant_frequencies.tsv', header = T, colClasses = c('character','numeric','numeric','character'))
ctrl_data <- read.delim('ctrl_len20_vec1000_yrs500_variant_frequencies.tsv', header = F, col.names = colnames(cr_data), colClasses = c('character','numeric','numeric','character'))


# compute cross reactivity values
df_catcher <- list() # initialize empty list to store df for each run
index <- 1 # set index for list collection
for (run in unique(cr_data$run_id)) {
  filtered_df <- cr_data %>% filter(run_id == run) # filter df for run
  if (nrow(filtered_df) < 5) { # min number of variants, otherwise will get error 
    next
  }
  df_catcher[[index]] <- batch_clustering_data(filtered_df) # collect data for run
  index <- index+1 # set index for next round
}

# create df
cr_df <- bind_rows(df_catcher) # bind rows to create master df
cr_df$selection <- 'cross reactivity' # label rows with sim condition 

# compute host specialization values
df_catcher <- list() # initialize empty list to store df for each run
index <- 1 # set index for list collection
for (run in unique(hs_data$run_id)) {
  filtered_df <- hs_data %>% filter(run_id == run) # filter df for run
  if (nrow(filtered_df) < 5) { # min number of variants, otherwise will get error 
    next
  }
  df_catcher[[index]] <- batch_clustering_data(filtered_df) # collect data for run
  index <- index+1 # set index for next round
}

# create df
hs_df <- bind_rows(df_catcher)
hs_df$selection <- 'host specialization'


# compute control values
df_catcher <- list() # initialize empty list to store df for each run
index <- 1 # set index for list collection
for (run in unique(ctrl_data$run_id)) {
  filtered_df <- ctrl_data %>% filter(run_id == run) # filter df for run
  if (nrow(filtered_df) < 5) { # min number of variants, otherwise will get error 
    next
  }
  df_catcher[[index]] <- batch_clustering_data(filtered_df) # collect data for run
  index <- index+1 # set index for next round
}

# create df
ctrl_df <- bind_rows(df_catcher)
ctrl_df$selection <- 'control'


# combine the 3 df's and convert to long format for plotting
combined_df <- bind_rows(cr_df, hs_df, ctrl_df)
df_long <- combined_df %>%
  select(-1) %>%
  pivot_longer(cols = -selection,  # Keep the "selection" column as it is
               names_to = "variable",  # Name of the new column with numeric column names
               values_to = "value")    # Values for the numeric columns

title.labs <- c('Mean within cluster p-w antigen distance', 
                'Mean across clusters p-w antigen distance',
                'Optimal k clusters',
                'Variants in final population',
                'Cluster sillhouette score')

names(title.labs) <- c('avg_in_dist',
                       'avg_out_dist',
                       'k_clusters',
                       'n_variants',
                       'sil_score')

# plot
ggplot(df_long, aes(x = selection, y = value)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 1.0)) +
  geom_jitter(aes(colour = selection)) +
  scale_color_manual(values = c("cross reactivity" = "skyblue3", "host specialization" = "red3", "control" = "yellow2")) +
  facet_wrap(~ variable, scales = "free_y", labeller = labeller(variable = title.labs), ncol = 3) +  
  labs(y = "", x = '') +
  theme_bw() +
  theme(legend.position = "right", panel.spacing = unit(1, "cm"), axis.text = element_text(size = 15), axis.text.x = element_blank()) +# axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  stat_compare_means(method = "t.test",
                     comparisons = list(c("cross reactivity", "host specialization"),
                                        c("cross reactivity", "control"),
                                        c("host specialization", "control")),
                     hide.ns = TRUE,
                     aes(group = selection),
                     label = "p.signif")

setwd('/Users/brandonely/Desktop/bb_strain_model_dev/new_sim_output/test_set_500yrs/plots/')
#ggsave('test_set_500_batch_cluster_boxplots.jpeg', height = 10, width = 12)

# dot plot for in/out distances
df_clust_dists <- combined_df %>% select(5:7)

ggplot(df_clust_dists, aes(x = avg_in_dist, y = avg_out_dist, color = selection)) +
  geom_point(alpha = 0.7) +
  #geom_jitter(alpha = 0.7) +
  scale_color_manual(values = c("cross reactivity" = "skyblue3", "host specialization" = "red3", "control" = "yellow2")) +
  xlim(0, 7) +
  ylim(0, 20) +  
  labs(x = "Mean in cluster distance", y = "Mean across cluster distance") +
  theme_bw() +
  theme(legend.text = element_text(size = 10),
        legend.position = 'top',
        legend.title = element_text(size = 0),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15))

setwd('/Users/brandonely/Desktop/bb_strain_model_dev/new_sim_output/test_set_500yrs/plots/')
ggsave('test_set_500_batch_cluster_Adist_DotPlot.jpeg', height = 8, width = 4.5)



# individual plots
ggplot(df_long %>% filter(variable == 'sil_score'), aes(x = selection, y = value)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 1.0)) +
  geom_jitter(aes(colour = selection)) +
  scale_color_manual(values = c("cross reactivity" = "skyblue3", "host specialization" = "red3", "control" = "yellow2")) +
  labs(y = "Silhouette score", x = '', title = 'Average silhouette score') +
  stat_compare_means(method = 't.test', aes(group = selection), label = 'p.signif', 
                     comparisons = list(c("cross reactivity", "host specialization"),
                                        c("cross reactivity", "control"),
                                        c("host specialization", "control"))) +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text = element_text(size = 15), 
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 20, hjust = 0.5))


setwd('/Users/brandonely/Desktop/bb_strain_model_dev/new_sim_output/test_set_500yrs/plots/')
#ggsave('test_set_500_batch_sil_scores.jpeg', height = 4, width = 6)

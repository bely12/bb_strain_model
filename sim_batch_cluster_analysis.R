setwd('/Users/brandonely/Desktop/bb_strain_model_dev/')
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(stringdist)
library(cluster)
source("sim_analysis_functions.R")

# upload data 
setwd('/Users/brandonely/Desktop/bb_strain_model_dev/new_sim_output/')
cr_data <- read.delim('cr_len20_vec1000_yrs250_variant_frequencies.tsv', header = T, colClasses = c('character','numeric','numeric','character'))
hs_data <- read.delim('hs_len20_vec1000_yrs250_variant_frequencies.tsv', header = T, colClasses = c('character','numeric','numeric','character'))

# first condition 
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

# second condition 
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

# combine the 2 df's and convert to long format for plotting
combined_df <- bind_rows(cr_df, hs_df)
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
  scale_color_manual(values = c("cross reactivity" = "skyblue3", "host specialization" = "red3")) +
  facet_wrap(~ variable, scales = "free_y", labeller = labeller(variable = title.labs), ncol = 3) +  
  labs(x = "Category", y = "Value") +
  theme_bw() +
  theme(legend.position = "top", panel.spacing = unit(1, "cm")) +
  stat_compare_means(method = "t.test",
                     hide.ns = TRUE,
                     aes(group = selection), 
                     label = "p.signif")


# dot plot for in/out distances
df_clust_dists <- combined_df %>% select(5:7)

ggplot(df_clust_dists, aes(x = avg_in_dist, y = avg_out_dist, color = selection)) +
  geom_point() +
  scale_color_manual(values = c("cross reactivity" = "skyblue3", "host specialization" = "red3")) +
  xlim(0, 20) +
  ylim(0, 20) +  
  labs(x = "Mean in cluster distance", y = "Mean across cluster distance") +
  theme_bw()


setwd('/Users/brandonely/Desktop/bb_strain_model_dev/')
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggExtra)
library(stringdist)
library(cluster)
library(vegan)
source("sim_analysis_functions.R")

### upload data
setwd('/Users/brandonely/Desktop/bb_strain_model_dev/sim_results/results/')

header <- c('variant', 'counts', 'frequency', 'run_id') # col names 

cr_data <- read.delim('immune_5000_20_500_variant_frequencies.tsv', header = T, col.names = header, colClasses = c('character','numeric','numeric','character'))
#cr_data <- cr_data %>% filter(frequency > 0.9)

hs_data <- read.delim('hs_5000_20_500_variant_frequencies.tsv', header = T, col.names = header, colClasses = c('character','numeric','numeric','character'))
#hs_data <- hs_data %>% filter(frequency > 0.9)

ctrl_data <- read.delim('neutral_5000_20_500_variant_frequencies.tsv', header = T, col.names = header, colClasses = c('character','numeric','numeric','character'))
#ctrl_data <- ctrl_data %>% filter(frequency > 0.9)

### clustering for each condition
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
cr_df <- bind_rows(df_catcher) # bind rows to create master df
cr_df$selection <- 'immune selection' # label rows with sim condition 

df_catcher <- list() 
index <- 1 
for (run in unique(hs_data$run_id)) {
  filtered_df <- hs_data %>% filter(run_id == run) # filter df for run
  if (nrow(filtered_df) < 5) { # min number of variants, otherwise will get error 
    next
  }
  df_catcher[[index]] <- batch_clustering_data(filtered_df) # collect data for run
  index <- index+1 # set index for next round
}
hs_df <- bind_rows(df_catcher)
hs_df$selection <- 'adaptive selection'

df_catcher <- list() 
index <- 1 
for (run in unique(ctrl_data$run_id)) {
  filtered_df <- ctrl_data %>% filter(run_id == run) # filter df for run
  if (nrow(filtered_df) < 5) { # min number of variants, otherwise will get error 
    next
  }
  df_catcher[[index]] <- batch_clustering_data(filtered_df) # collect data for run
  index <- index+1 # set index for next round
}
ctrl_df <- bind_rows(df_catcher)
ctrl_df$selection <- 'neutral selection'

### combine and format clustering data
combined_df <- bind_rows(cr_df, hs_df, ctrl_df)
combined_df_cleaned <- combined_df %>% select(1,3,5:8) # only keep what I need
df_long <- combined_df_cleaned %>%
  select(-1) %>%
  pivot_longer(cols = -selection,
               names_to = "variable", 
               values_to = "value")    

### make titles for plots in panel
title.labs <- c('Mean within cluster p-w antigen distance', 
                'Mean across clusters p-w antigen distance',
                'Cluster sillhouette score', 'shannon evenness')
names(title.labs) <- c('avg_in_dist',
                       'avg_out_dist',
                       'sil_score',
                       'ShanE')

### boxplot panel for all clustering results ###
ggplot(df_long, aes(x = selection, y = value)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 1.0)) +
  geom_jitter(aes(colour = selection)) +
  scale_color_manual(values = c("immune selection" = "skyblue3", "adaptive selection" = "red3", "neutral selection" = "yellow2")) +
  facet_wrap(~ variable, scales = "free_y", labeller = labeller(variable = title.labs), ncol = 3) +
  labs(y = "", x = '') +
  theme_bw() +
  theme(legend.position = "right", panel.spacing = unit(1, "cm"), axis.text = element_text(size = 15), axis.text.x = element_blank()) +# axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  stat_compare_means(method = "t.test",
                     comparisons = list(c("immune selection", "adaptive selection"),
                                        c("immune selection", "neutral selection"),
                                        c("adaptive selection", "neutral selection")),
                     hide.ns = F,
                     aes(group = selection),
                     label = "p.signif")

setwd('/Users/brandonely/Desktop/bb_strain_model_dev/sim_results/results/plots/')
#ggsave('v5000_20_500_batch_cluster_boxplots_UNFILTERED.jpeg', height = 10, width = 15)


### individual boxplots ###

### sil and shannon scores together ###
# make titlew 
title.labs2 <- c('Cluster sillhouette score', 'shannon evenness')
names(title.labs2) <- c('sil_score', 'ShanE')

# filter df and assign factors to put in desired order
df_long_subset1 <- df_long %>% filter(variable == 'sil_score' | variable == 'ShanE')
df_long_subset1$variable <- factor(df_long_subset1$variable, levels = c('sil_score', 'ShanE'))

# plot
ggplot(df_long_subset1, aes(x = selection, y = value)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 1.0)) +
  geom_jitter(aes(colour = selection)) +
  scale_color_manual(values = c("immune selection" = "skyblue3", "adaptive selection" = "red3", "neutral selection" = "yellow2")) +
  facet_wrap(~ variable, scales = 'free_y', labeller = labeller(variable = title.labs2)) +
  labs(y = "value") +
  stat_compare_means(method = "t.test",
                     comparisons = list(c("immune selection", "adaptive selection"),
                                        c("immune selection", "neutral selection"),
                                        c("adaptive selection", "neutral selection")),
                     hide.ns = F,
                     aes(group = selection),
                     label = "p.signif") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 20, hjust = 0.5))

#ggsave('sil_and_shannon_boxplot.jpeg', height = 6, width = 6)

### in and out distance together in one boxplot ###
# make titles
title.labs3 <- c('Mean in cluster distance', 'Mean across cluster distance')
names(title.labs3) <- c('avg_in_dist', 'avg_out_dist')

# filter df
df_long_subset2 <- df_long %>% filter(variable == 'avg_in_dist' | variable == 'avg_out_dist')

# plot
ggplot(df_long_subset2, aes(x = selection, y = value)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 1.0)) +
  geom_jitter(aes(colour = selection)) +
  scale_color_manual(values = c("immune selection" = "skyblue3", "adaptive selection" = "red3", "neutral selection" = "yellow2")) +
  facet_wrap(~ variable, labeller = labeller(variable = title.labs3)) +
  labs(y = "value") +
  stat_compare_means(method = "t.test",
                     comparisons = list(c("immune selection", "adaptive selection"),
                                        c("immune selection", "neutral selection"),
                                        c("adaptive selection", "neutral selection")),
                     hide.ns = F,
                     aes(group = selection),
                     label = "p.signif") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 20, hjust = 0.5))

#ggsave('mean_in_across_dist_boxplot.jpeg', height = 6, width = 6)

### dot plot for in/out distances with margin density histograms ###
df_clust_dists <- combined_df_cleaned %>% select(3,4,6)

p <- ggplot(df_clust_dists, aes(x = avg_in_dist, y = avg_out_dist, color = selection)) +
  geom_point(alpha = 0.5) +
  #stat_ellipse(level = 0.95, type = "norm", size = 1) +
  scale_color_manual(values = c("immune selection" = "skyblue3", "adaptive selection" = "red3", "neutral selection" = "yellow2")) +
  xlim(0, 7) +
  ylim(0, 15) +  
  labs(x = "Mean in cluster distance", y = "Mean across cluster distance") +
  theme_bw() +
  theme(legend.text = element_text(size = 10),
        legend.position = 'bottom',
        legend.title = element_text(size = 0),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15))

ggMarginal(p, type = "density", groupColour = F, groupFill = TRUE, alpha = 0.5)
p2 <- ggMarginal(p, type = "density", groupColour = F, groupFill = TRUE, alpha = 0.5)

#ggsave('v5000_20_500_batch_cluster_AntDist_DotPlot_UNFILTERED.jpeg', plot = p2, height = 8, width = 8)

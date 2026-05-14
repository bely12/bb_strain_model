setwd('/Users/brandonely/Desktop/bb_strain_model_dev/R_analysis_scripts/')
library(tidyverse)
library(purrr)
library(ggpubr)
library(ggExtra)
library(stringdist)
library(cluster)
library(ggdendro)
library(RColorBrewer)
library(pheatmap)
library(viridis)
library(vegan)
source("sim_analysis_functions.R")

##### cluster analysis on everything, then bind all results to one table for analysis #####
### upload data one res table at a time, combine results after clustering 
setwd('/Users/brandonely/Desktop/bb_strain_model_dev/evo_sim/results_tables/sim_output/adaptive_results/')

pop <- read.delim('adaptive_HighGenFit_variant_freqs_sampled.tsv', header = T, colClasses = c('character','numeric','numeric', 'character', 'character', 'character','numeric','numeric','character'))

#pop <- read.delim('hybrid_HighGenFit_variant_freqs_sampled.tsv', header = T, colClasses = c('character','numeric','numeric', 'character', 'character','character', 'character', 'character','numeric','numeric','character'))
#pop <- pop %>% filter(rec_rate == 0.01)
length(unique(pop$run_id))

#header_names <- c("strain", "count", "freq", "run_id", "selection", "gene_type", "rec_rate", "mut_rate", "gen_fitness")

pop$host_pops <- 'even'
pop$host_pop_stability <- 'stable'

# df <- pop %>%
#   group_by(run_id) %>%
#   filter(n() >= 30) %>% # only keep groups with at least 30 rows
#   slice_max(order_by = count, n = 50, with_ties = FALSE) %>% # take top 30 counts
#   ungroup()


# filter to get rid of "ghost variants"
pop <- pop %>% filter(freq >= 1)

# check how many variants are in these pops, readjust filter above if needed
x <- pop %>% group_by(run_id) %>% summarise(count = n())
mean(x$count)
x_filtered <- x %>% filter(count >= 10)

df <- pop %>% filter(run_id %in% x_filtered$run_id)
length(unique(df$run_id))

##### if using full multi strain data with batch_clustering_data_V2 #####
# df <- df %>%
#   mutate(
#     immune = substr(strain, 1, 20),
#     adaptive = substr(strain, 21, 40))
###### 

# cluster 
df_catcher <- list() # initialize empty list to store df for each run
index <- 1 # set index for list collection
for (run in unique(df$run_id)) { 
  filtered_df <- df %>% filter(run_id == run) # filter df for run
  if (nrow(filtered_df) < 5) { # min number of variants, otherwise will get error 
    next
  }
  df_catcher[[index]] <- batch_clustering_data(filtered_df) # collect data for run
  index <- index+1 # set index for next round
  print(index)
}


# bind all cluster_res together
cluster_res1 <- map_dfr(df_catcher, "cluster_res")
# bind all variant_cluster_labels together
variant_labels1 <- map_dfr(df_catcher, "variant_cluster_labels")
###
# bind all cluster_res together
cluster_res2 <- map_dfr(df_catcher, "cluster_res")
# bind all variant_cluster_labels together
variant_labels2 <- map_dfr(df_catcher, "variant_cluster_labels")
###
# bind all cluster_res together
cluster_res3 <- map_dfr(df_catcher, "cluster_res")
# bind all variant_cluster_labels together
variant_labels3 <- map_dfr(df_catcher, "variant_cluster_labels")


master_clust_res <- bind_rows(cluster_res1, cluster_res2, cluster_res3)
#write.table(cluster_res1, 'adpt_ExpMod_NoRec_cluster_results.tsv', sep = '\t', quote = F, row.names = F)
master_variant_cluster_labels <- bind_rows(variant_labels1, variant_labels2, variant_labels3)
#write.table(variant_labels1, 'adpt_ExpMod_NoRec_cluster_labels.tsv', sep = '\t', quote = F, row.names = F)

# check new dataframes before leaving
setwd('/Users/brandonely/Desktop/bb_strain_model_dev/evo_sim/results_tables/cluster_analysis_results/')
#master_clust <- read.delim('adaptive_cluster_results.tsv', sep = '\t', header = T)
#new <- bind_rows(master_clust, cluster_res1)
#write.table(new, 'adaptive_cluster_results.tsv', sep = '\t', quote = F, row.names = F)

#master_var_labels <- read.delim('hybrid_cluster_labels.tsv', sep = '\t', header = T, colClasses = c('character', 'integer', 'numeric', 'numeric', 'character', 'character', 'character', 'character', 'numeric', 'numeric', 'character', 'character'))
#new <- bind_rows(master_var_labels, variant_labels1)
#write.table(new, 'adaptive_cluster_labels.tsv', sep = '\t', quote = F, row.names = F)

####### LOOKING AT VARIANCE IN HS VALS WHEN CLUSTERING COMBINED GENES VS JUST IMMUNE GENE #####
df <- variant_labels1 %>%
  mutate(
    immune = substr(strain, 1, 20),
    adaptive = substr(strain, 21, 40)
  )

df <- variant_labels1 %>% mutate(hs_val = str_count(pattern = '1', adaptive)/20)
df <- df %>% select(run_id,strain, immune, adaptive, cluster_label, hs_val)

df_var <- df %>%
  group_by(run_id, cluster_label) %>%
  summarise(
    cluster_hs_var = var(hs_val)
  )

mean(df_var$cluster_hs_var, na.rm = TRUE)

hist(df_var$cluster_hs_var, breaks = 50)


######


##### analyze clustering results #####
setwd('/Users/brandonely/Desktop/bb_strain_model_dev/evo_sim/results_tables/cluster_analysis_results/')

# upload data
clust_res <- read.delim('hybrid_cluster_results.tsv', sep = '\t', header = T)
clust_res$category <- paste(clust_res$selection, clust_res$gen_fit)

clust_res2 <- read.delim('immune_cluster_results.tsv', sep = '\t', header = T)
clust_res2$category <- paste(clust_res2$selection, clust_res2$gen_fit)

clust_res3 <- read.delim('adaptive_cluster_results.tsv', sep = '\t', header = T)
clust_res3$category <- paste(clust_res3$selection, clust_res3$gen_fit)

clust_res4 <- read.delim('no_selection_cluster_results.tsv', sep = '\t', header = T)
clust_res4$category <- paste(clust_res4$selection, clust_res4$gen_fit)


df <- bind_rows(clust_res, clust_res2, clust_res3, clust_res4)
df2 <- bind_rows(clust_res, clust_res2, clust_res3, clust_res4)

# convert to long format
df <- df %>% 
  #filter(selection != 'hybrid' & rec_rate == 0.01) %>%
  pivot_longer(
  cols = n_variants:ShanE,
  names_to = 'variable',
  values_to = 'value'
)

# counts for each category
x <- df %>% 
  mutate(category = paste(selection,host_pops,host_pop_stability,gen_fit, sep = '_')) %>%
  filter(rec_rate == 0.01)

x <- df %>%
  group_by(category) %>%
  summarise(unique_runs = n_distinct(run_id))
unique(df$variable)
# immune and neutral selection at different recombination rates
unique(df$variable)
set_compare <- c("immune high", "neutral high", "hybrid exp_mod")

df %>% 
  filter((selection == 'immune' | selection == 'neutral'), (variable == 'avg_in_dist' | variable == 'avg_out_dist')) %>%
  #mutate(category = paste(selection,rec_rate,sep = '_')) %>%
  
  ggplot(aes(x = category, y = value, fill = category)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~variable, scales = 'free_y')+
  #facet_grid(variable ~ selection, scales = 'free_y') +
  theme_bw()

setwd('/Users/brandonely/Desktop/bb_strain_model_dev/evo_sim/results_tables/new_plots/')
ggsave('cluster_results_InOut_box_imm_no_selection.png', width = 8, height = 10)

df %>% 
  filter(category %in% set_compare & variable != 'avg_in_dist' & variable != 'avg_out_dist') %>%
  #mutate(category = paste(selection,rec_rate,sep = '_')) %>%
  
  ggplot(aes(x = category, y = value, fill = category)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~variable, scales = 'free_y')+
  #facet_grid(variable ~ selection, scales = 'free_y') +
  theme_bw()

setwd('/Users/brandonely/Desktop/bb_strain_model_dev/evo_sim/results_tables/new_plots/')
ggsave('cluster_results_SetCompare.png', width = 8, height = 10)


# everything together
df %>%
  #mutate(category = paste(selection,host_pops,host_pop_stability,gen_fit, sep = '_')) %>%
  #filter(rec_rate == 0.01 & variable != "ShanE" & variable != 'k_clusters') %>%
  
  ggplot(aes(x = category, y = value, fill = category)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7)) +
  facet_wrap(~variable, scales = 'free_y', ncol = 3)

setwd('/Users/brandonely/Desktop/bb_strain_model_dev/evo_sim/results_tables/new_plots/')
#ggsave('cluster_results_ALL.png', width = 8, height = 10)

# dot plot in/out dists with margin density histogram
p <- df2 %>%
  filter(category %in% set_compare) %>%
  #mutate(category = paste(selection,host_pops,host_pop_stability,gen_fit, sep = '_')) %>%
  select(category, avg_in_dist, avg_out_dist) %>%
  
  ggplot(aes(x = avg_in_dist, y = avg_out_dist, color = category)) +
  geom_point(alpha = 0.5) +
  #stat_ellipse(level = 0.99, type = "norm", size = 1) +
  labs(x = "Mean in cluster distance", y = "Mean across cluster distance") +
  theme_bw() +
  geom_hline(yintercept = 10, linetype = 2) +
  geom_vline(xintercept = 3, linetype = 2) +
  theme(legend.text = element_text(size = 10),
        legend.position = 'bottom',
        legend.title = element_text(size = 0),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15))
  #facet_wrap(~ category, ncol = 3)

#ggsave('InOut_ClusterDists_SetCompare.png', width = 8, height = 8)

ggMarginal(p, type = "density", groupColour = F, groupFill = TRUE, alpha = 0.5)
p2 <- ggMarginal(p, type = "density", groupColour = F, groupFill = TRUE, alpha = 0.5)
#ggsave('cluster_AntDist_DotPlot_SetCompare.jpeg', plot = p2, height = 8, width = 8)

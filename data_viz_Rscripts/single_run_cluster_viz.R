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

##### upload data to get pools of typical runs #####
setwd('/Users/brandonely/Desktop/bb_strain_model_dev/evo_sim/results_tables/cluster_analysis_results/')

clust_res <- read.delim('hybrid_cluster_results.tsv', sep = '\t', header = T)
clust_res$category <- paste(clust_res$selection, clust_res$gen_fit)
clust_res <- clust_res %>% filter(gen_fit == 'exp_mod')

clust_res2 <- read.delim('immune_cluster_results.tsv', sep = '\t', header = T)
clust_res2$category <- paste(clust_res2$selection, clust_res2$gen_fit)

clust_res3 <- read.delim('no_selection_cluster_results.tsv', sep = '\t', header = T)
clust_res3$category <- paste(clust_res3$selection, clust_res3$gen_fit)


df <- bind_rows(clust_res, clust_res2, clust_res3)
summary_df <- df %>% filter(selection == 'hybrid') # switch out with desired selection mode for each

summary_df_stats <- data.frame(
  mean_sil = mean(summary_df$sil_score),
  sil_sd = sd(summary_df$sil_score),
  mean_k = mean(summary_df$k_clusters),
  k_sd = sd(summary_df$k_clusters),
  mean_var = mean(summary_df$n_variants),
  var_sd = sd(summary_df$n_variants),
  mean_InDist = mean(summary_df$avg_in_dist),
  In_sd = sd(summary_df$avg_in_dist),
  mean_OutDist = mean(summary_df$avg_out_dist),
  Out_sd = sd(summary_df$avg_out_dist),
  mean_shan = mean(summary_df$ShanE),
  shan_sd = sd(summary_df$ShanE))

# CREATE A POOL OF RUNS THAT ARE WITHIN 1 STANDARD DEVIATIO OF THE MEAN FOR ALL RUNS
typical_runs <- summary_df %>% filter(sil_score > (summary_df_stats$mean_sil[1] - summary_df_stats$sil_sd[1]) & 
                                        sil_score < (summary_df_stats$mean_sil[1] + summary_df_stats$sil_sd[1]) &
                                        n_variants > (summary_df_stats$mean_var[1] - summary_df_stats$var_sd[1]) & 
                                        n_variants < (summary_df_stats$mean_var[1] + summary_df_stats$var_sd[1]) &
                                        k_clusters > (summary_df_stats$mean_k[1] - summary_df_stats$k_sd[1]) &
                                        k_clusters < (summary_df_stats$mean_k[1] + summary_df_stats$k_sd[1]) &
                                        avg_in_dist > (summary_df_stats$mean_InDist[1] - summary_df_stats$In_sd[1]) &
                                        avg_in_dist < (summary_df_stats$mean_InDist[1] + summary_df_stats$In_sd[1]) &
                                        avg_out_dist > (summary_df_stats$mean_OutDist[1] - summary_df_stats$Out_sd[1]) &
                                        avg_out_dist < (summary_df_stats$mean_OutDist[1] + summary_df_stats$Out_sd[1]) &
                                        ShanE > (summary_df_stats$mean_shan[1] - summary_df_stats$shan_sd[1]) &
                                        ShanE < (summary_df_stats$mean_shan[1] + summary_df_stats$shan_sd[1]))


##### RUN CLUSTER ANALYSIS ON SELECTED RUN #####
setwd('/Users/brandonely/Desktop/bb_strain_model_dev/evo_sim/results_tables/sim_output/no_selection_results/')
data <- read.delim('neutral_variant_freqs_sampled.tsv', header = T, colClasses = c('character','numeric','numeric', 'character', 'character', 'character','numeric','numeric','character'))

# if choosing run randomly
selected_run <- typical_runs[sample(nrow(typical_runs), 1), ]
cat('chosen run id:',selected_run$run_id)
run_pop <- data %>% filter(run_id == selected_run$run_id, freq >= 1.0) # choose a run to work

# if choosing run manually
immune_pop <- data %>% filter(run_id == 'run_94', freq >= 1.0)
neutral_pop <- data %>% filter(run_id == 'run_4', freq >= 1.0)
hybrid_pop <- data %>% filter(run_id == 'run_19', freq >= 1.0)

# run clustering on the individual run 
immune_res <- single_run_clustering(immune_pop) 
neutral_res <- single_run_clustering(neutral_pop)
hybrid_res <- single_run_clustering(hybrid_pop)


### for immune selection ### ### run_94
### PLOT DENDROGRAM
dendro_data <- dendro_data(immune_res$cluster_res)
ggplot(segment(dendro_data)) +
  geom_segment(aes(x = -y, y = x, xend = -yend, yend = xend)) +
  theme_bw() +
  labs(title = paste("Selection: immune\nk = ",immune_res$stats_df$k_clusters,
                        '\nsilhouette score = ',round(x = immune_res$stats_df$sil_score,digits = 3),
                        '\nMean in cluster distance = ',round(x = immune_res$stats_df$avg_in_dist, digits = 3),
                        '\nMean cross cluster distance = ',round(x=immune_res$stats_df$avg_out_dist, digits = 3))) +
  theme(axis.text = element_blank(),axis.title = element_blank(),panel.grid = element_blank(),
        plot.title = element_text(size = 15),
        plot.subtitle = element_text(size = 15))

setwd('/Users/brandonely/Desktop/bb_strain_model_dev/evo_sim/results_tables/new_plots/')
#ggsave('immune_cluster_dendro_run94.png', width = 8, height = 8)


### MAKE HEATMAP
# set up matrix with variants as row/col names
variant_names <- as.character(trimws(immune_res$variant_cluster_labels$strain))
mat <- as.matrix(immune_res$dist_matrix)
rownames(mat) <- variant_names
colnames(mat) <- variant_names

# set up annotation df for pheatmap
annot_rows <- immune_res$variant_cluster_labels %>% select(strain, cluster_label)
rownames(annot_rows) <- annot_rows$strain
annot_rows <- annot_rows %>% select(-1)

# set colors for annotation
clusters <- sort(unique(annot_rows$cluster_label))
#cluster_colors <- setNames(brewer.pal(length(clusters), "Set3"), clusters)
cluster_colors <- unique(c(
  brewer.pal(9, "Set1"),
  brewer.pal(8, "Set2"),
  brewer.pal(12, "Set3"),
  brewer.pal(8, "Dark2")
))[1:31]
annotation_colors <- list(cluster_label = cluster_colors)

breaks_vals <- seq(0, 20, length.out = 21)
pheatmap(
  mat,               
  cluster_rows = immune_res$cluster_res,      
  cluster_cols = immune_res$cluster_res, 
  annotation_row = annot_rows,
  annotation_col = annot_rows,
  annotation_colors = annotation_colors,
  annotation_legend = F,
  annotation_names_row = F,
  annotation_names_col = F,     
  color = viridis(20),
  breaks = breaks_vals,
  border_color = NA,
  treeheight_col = 0,
  show_rownames = F,          
  show_colnames = F,         
  filename = 'immune_cluster_heatmap_run94.jpeg',
  width = 5,
  height = 5)


### no selection ### ### run_4
dendro_data <- dendro_data(neutral_res$cluster_res)
ggplot(segment(dendro_data)) +
  geom_segment(aes(x = -y, y = x, xend = -yend, yend = xend)) +
  theme_bw() +
  labs(title = paste("Selection: no selection\nk = ",neutral_res$stats_df$k_clusters,
                     '\nsilhouette score = ',round(x = neutral_res$stats_df$sil_score,digits = 3),
                     '\nMean in cluster distance = ',round(x = neutral_res$stats_df$avg_in_dist, digits = 3),
                     '\nMean cross cluster distance = ',round(x=neutral_res$stats_df$avg_out_dist, digits = 3))) +
  theme(axis.text = element_blank(),axis.title = element_blank(),panel.grid = element_blank(),
        plot.title = element_text(size = 15),
        plot.subtitle = element_text(size = 15))

#ggsave('no_selection_cluster_dendro_run4.png', width = 8, height = 8)

### MAKE HEATMAP
# set up matrix with variants as row/col names
variant_names <- as.character(trimws(neutral_res$variant_cluster_labels$strain))
mat <- as.matrix(neutral_res$dist_matrix)
rownames(mat) <- variant_names
colnames(mat) <- variant_names

# set up annotation df for pheatmap
annot_rows <- neutral_res$variant_cluster_labels %>% select(strain, cluster_label)
rownames(annot_rows) <- annot_rows$strain
annot_rows <- annot_rows %>% select(-1)

# set colors for annotation
clusters <- sort(unique(annot_rows$cluster_label))
#cluster_colors <- setNames(brewer.pal(length(clusters), "Set3"), clusters)
cluster_colors <- unique(c(
  brewer.pal(9, "Set1"),
  brewer.pal(8, "Set2"),
  brewer.pal(12, "Set3"),
  brewer.pal(8, "Dark2")
))[1:31]
annotation_colors <- list(cluster_label = cluster_colors)

breaks_vals <- seq(0, 20, length.out = 21)
pheatmap(
  mat,               
  cluster_rows = neutral_res$cluster_res,      
  cluster_cols = neutral_res$cluster_res, 
  annotation_row = annot_rows,
  annotation_col = annot_rows,
  annotation_colors = annotation_colors,
  annotation_legend = F,
  annotation_names_row = F,
  annotation_names_col = F,     
  color = viridis(20),
  breaks = breaks_vals,
  border_color = NA,
  treeheight_col = 0,
  show_rownames = F,          
  show_colnames = F,         
  filename = 'no_selection_cluster_heatmap_run4.jpeg',
  width = 5,
  height = 5)

### hybrid selection ### ### run_19 
dendro_data <- dendro_data(hybrid_res$cluster_res)
ggplot(segment(dendro_data)) +
  geom_segment(aes(x = -y, y = x, xend = -yend, yend = xend)) +
  theme_bw() +
  labs(title = paste("Selection: hybrid\nk = ",hybrid_res$stats_df$k_clusters,
                     '\nsilhouette score = ',round(x = hybrid_res$stats_df$sil_score,digits = 3),
                     '\nMean in cluster distance = ',round(x = hybrid_res$stats_df$avg_in_dist, digits = 3),
                     '\nMean cross cluster distance = ',round(x=hybrid_res$stats_df$avg_out_dist, digits = 3))) +
  theme(axis.text = element_blank(),axis.title = element_blank(),panel.grid = element_blank(),
        plot.title = element_text(size = 15),
        plot.subtitle = element_text(size = 15))

#ggsave('hybrid_cluster_dendro_run19.png', width = 8, height = 8)

### MAKE HEATMAP
# set up matrix with variants as row/col names
variant_names <- as.character(trimws(hybrid_res$variant_cluster_labels$strain))
mat <- as.matrix(hybrid_res$dist_matrix)
rownames(mat) <- variant_names
colnames(mat) <- variant_names

# set up annotation df for pheatmap
annot_rows <- hybrid_res$variant_cluster_labels %>% select(strain, cluster_label)
rownames(annot_rows) <- annot_rows$strain
annot_rows <- annot_rows %>% select(-1)

# set colors for annotation
clusters <- sort(unique(annot_rows$cluster_label))
#cluster_colors <- setNames(brewer.pal(length(clusters), "Set3"), clusters)
cluster_colors <- unique(c(
  brewer.pal(9, "Set1"),
  brewer.pal(8, "Set2"),
  brewer.pal(12, "Set3"),
  brewer.pal(8, "Dark2")
))[1:31]
annotation_colors <- list(cluster_label = cluster_colors)

breaks_vals <- seq(0, 20, length.out = 21)
pheatmap(
  mat,               
  cluster_rows = hybrid_res$cluster_res,      
  cluster_cols = hybrid_res$cluster_res, 
  annotation_row = annot_rows,
  annotation_col = annot_rows,
  annotation_colors = annotation_colors,
  annotation_legend = F,
  annotation_names_row = F,
  annotation_names_col = F,     
  color = viridis(20),
  breaks = breaks_vals,
  border_color = NA,
  treeheight_col = 0,
  show_rownames = F,          
  show_colnames = F,         
  filename = 'hybrid_cluster_heatmap_run19.jpeg',
  width = 5,
  height = 5)

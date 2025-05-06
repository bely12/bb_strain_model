# sim cluster analysis for looking at single run
setwd('/Users/brandonely/Desktop/bb_strain_model_dev/')
library(tidyverse)
library(cluster)  
library(pheatmap)
library(viridis) 
library(stringdist)
library(ggdendro)
library(vegan)
library(RColorBrewer)
source("sim_analysis_functions.R")

### UPLOAD DATA
setwd('/Users/brandonely/Desktop/bb_strain_model_dev/sim_results/results/')
header <- c('variant', 'counts', 'frequency', 'run_id')
data <- read.delim('hs_5000_20_500_25SET_sampled_variant_frequencies.tsv', header = T, col.names = header, colClasses = c('character','numeric','numeric','character'))


### CHOOSING A SINGLE RUN TO SHOW CLUSTERING RESULTS FOR

# STEP 1: PERFOFRM CLUSTERING ON ALL RUNS
df_catcher <- list()
index <- 1 
for (run in unique(data$run_id)) {
  filtered_df <- data %>% filter(run_id == run)
  if (nrow(filtered_df) < 5) { 
    next
  }
  df_catcher[[index]] <- batch_clustering_data(filtered_df)
  index <- index+1
}
summary_df <- bind_rows(df_catcher)

# STEP 2: GET MEANS AND STANDARD DEVIATIONS ON CLUSTERING RESULTS
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

# STEP 3: CREATE A POOL OF RUNS THAT ARE WITHIN 1 STANDARD DEVIATIO OF THE MEAN FOR ALL RUNS
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


### RUN CLUSTER ANALYSIS ON SELECTED RUN
run_pop <- data %>% filter(run_id == 'run_20') # choose a run to work
cluster_results <- single_run_clustering(run_pop) # execute with homemade function 

### PLOT DENDROGRAM
dendro_data <- dendro_data(cluster_results$cluster_res)
ggplot(segment(dendro_data)) +
  geom_segment(aes(x = -y, y = x, xend = -yend, yend = xend)) +
  theme_minimal() +
  labs(title = "Hierarchical clustering of variants in final population", 
       subtitle = paste("Selection: ___\nk = ",cluster_results$stats_df$k_clusters,
         '\nsilhouette score = ',round(x = cluster_results$stats_df$sil_score,digits = 3),
         '\nMean in cluster distance = ',cluster_results$stats_df$avg_in_dist,
         '\nMean cross cluster distance = ',cluster_results$stats_df$avg_out_dist)) +
  theme(axis.text = element_blank(),axis.title = element_blank(),panel.grid = element_blank(),
        plot.title = element_text(size = 20),
        plot.subtitle = element_text(size = 20))

#setwd('/Users/brandonely/Desktop/bb_strain_model_dev/new_sim_output/test_set_500yrs/plots/')
#ggsave('test_set_500_ctrl_singleRun12_hcDendro.jpeg', height = 5, width = 10)


### MAKE HEATMAP
# set up matrix with variants as row/col names
variant_names <- as.character(trimws(cluster_results$variant_cluster_labels$variant))
mat <- as.matrix(cluster_results$dist_matrix)
rownames(mat) <- variant_names
colnames(mat) <- variant_names

# set up annotation df for pheatmap
annot_rows <- cluster_results$variant_cluster_labels %>% select(variant, cluster_label)
rownames(annot_rows) <- annot_rows$variant
annot_rows <- annot_rows %>% select(-1)

# set colors for annotation
clusters <- sort(unique(annot_rows$cluster_label))
cluster_colors <- setNames(brewer.pal(length(clusters), "Set1"), clusters)
annotation_colors <- list(cluster_label = cluster_colors)

breaks_vals <- seq(0, 20, length.out = 21)
pheatmap(
  mat,               
  cluster_rows = cluster_results$cluster_res,      
  cluster_cols = cluster_results$cluster_res, 
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
  show_colnames = F)          
  #filename = 'ctrl_singleRun12_heatmap.jpeg'


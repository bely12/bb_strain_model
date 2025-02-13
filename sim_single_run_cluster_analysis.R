# sim cluster analysis for looking at single run
setwd('/Users/brandonely/Desktop/bb_strain_model_dev/')
library(tidyverse)
library(cluster)  
library(pheatmap)
library(viridis) 
library(stringdist)
library(RColorBrewer)
library(wesanderson)
library(ggdendro)
source("sim_analysis_functions.R")

# upload data 
setwd('/Users/brandonely/Desktop/bb_strain_model_dev/new_sim_output/test_set_500yrs/')

cr_data <- read.delim('cr_len20_vec1000_yrs500_variant_frequencies.tsv', header = T, colClasses = c('character','numeric','numeric','character'))
hs_data <- read.delim('hs_len20_vec1000_yrs500_variant_frequencies.tsv', header = T, colClasses = c('character','numeric','numeric','character'))
ctrl_data <- read.delim('ctrl_len20_vec1000_yrs500_variant_frequencies.tsv', header = F, col.names = colnames(cr_data), colClasses = c('character','numeric','numeric','character'))

selected_data <- ctrl_data # set to desired data to avoid changing variable names

# compute stats from sim results to inform decision on which sim to choose
# get stats on number of variants in final pop for sim runs
result <- selected_data %>%
  group_by(run_id) %>%
  tally() %>%
  arrange(desc(n))
mean(result$n)
median(result$n)
max(result$n)
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
getmode(result$n)

# compute k and sil values for sim runs
df_catcher <- list() # initialize empty list to store df for each run
index <- 1 # set index for list collection
for (run in unique(selected_data$run_id)) {
  filtered_df <- selected_data %>% filter(run_id == run) # filter df for run
  if (nrow(filtered_df) < 5) { # min number of variants, otherwise will get error 
    next
  }
  df_catcher[[index]] <- batch_clustering_data(filtered_df) # collect data for run
  index <- index+1 # set index for next round
}
# create df
selected_summary <- bind_rows(df_catcher) # bind rows to create master df

# start analysis for single sim run
run_pop <- selected_data %>% filter(run_id == 'run_12') # choosing a run to work with

# create binary matrix and distance matrix
binary_matrix <- do.call(rbind, strsplit(run_pop$variant, split = ""))
binary_matrix_str <- apply(binary_matrix, 1, paste, collapse = "")
binary_matrix_numeric <- matrix(as.numeric(binary_matrix), nrow = nrow(binary_matrix), ncol = ncol(binary_matrix))
dist_matrix <- stringdistmatrix(binary_matrix_str, binary_matrix_str, method = "hamming")
dist_obj <- as.dist(dist_matrix) # convert to a distance object, required for hclust function 

hc_result <- hclust(dist_obj, method = "complete") # perform clustering
#plot(hc_result, main = "Hierarchical Clustering") # simple plot of clustering results

max_k <- length(binary_matrix_str) -1 # set total number of clusters allowed, set to n-1

sil_scores <- sapply(2:max_k, function(k) { # get sil scores for each k value
  clusters <- cutree(hc_result, k)
  sil <- silhouette(clusters, dist_matrix)
  mean(sil[, 3])
})

max_sil_score <- max(sil_scores) # highest sil score
optimal_k <- which.max(sil_scores) + 1 # final k value based on highest sil score
clusters <- cutree(hc_result, k = optimal_k) # cut dendrogram using optimal k
#sil <- silhouette(clusters, dist = dist_matrix) # re-compute sil score for plot
#plot(sil) # simple plot

# get average within and across cluster distances 
run_pop$cluster_label <- clusters # add cluster labels to df
in_cluster <- round(x = InGroupDistance(run_pop), digits = 3) # in cluster distance
cross_cluster <- round(x = OutGroupDistance(run_pop), digits = 3) # cross cluster distance

# plot dendrogram
dendro_data <- dendro_data(hc_result)
ggplot(segment(dendro_data)) +
  geom_segment(aes(x = -y, y = x, xend = -yend, yend = xend)) +
  theme_minimal() +
  labs(title = "Hierarchical clustering of variants in final population", 
       subtitle = paste("Selection: Neutral selection\nk = ",optimal_k,
                        '\nsilhouette score = ',round(x = max_sil_score,digits = 3),
                        '\nMean in cluster distance = ',in_cluster,
                        '\nMean cross cluster distance = ',cross_cluster)) +
  theme(axis.text = element_blank(),axis.title = element_blank(),panel.grid = element_blank(),
        plot.title = element_text(size = 20),
        plot.subtitle = element_text(size = 20))

setwd('/Users/brandonely/Desktop/bb_strain_model_dev/new_sim_output/test_set_500yrs/plots/')
ggsave('test_set_500_ctrl_singleRun12_hcDendro.jpeg', height = 5, width = 10)

# heatmap of pairwise distances with clustering
#zzz <- wes_palette("Zissou1", 20, type = "continuous")
breaks_vals <- seq(0, 20, length.out = 21)
#breaks_vals <- seq(0,16, length.out = 17)
pheatmap(
  dist_matrix,               
  cluster_rows = hc_result,      
  cluster_cols = hc_result,      
  display_numbers = FALSE,       
  color = viridis(20),
  #color = zzz,
  breaks = breaks_vals,
  border_color = NA,
  #main = "Final population p-w antigen distances\nSelection: Host specialization",
  scale = "none",               
  fontsize = 8,                 
  show_rownames = FALSE,          
  show_colnames = FALSE,          
  fontsize_row = 8,              
  fontsize_col = 8,
  filename = 'ctrl_singleRun12_heatmap.jpeg'
  #cellwidth = 3, 
  #cellheight = 3,
  #treeheight_row = 0,           
  #treeheight_col = 0
)




# sim cluster analysis for looking at single run
setwd('/Users/brandonely/Desktop/bb_strain_model_dev/')
library(tidyverse)
library(ggplot2)
library(cluster)  
library(Rtsne) 
library(pheatmap)
library(viridis) 
library(stringdist)
library(RColorBrewer)
library(wesanderson)
library(ggdendro)
source("sim_analysis_functions.R")

# upload data 
setwd('/Users/brandonely/Desktop/bb_strain_model_dev/new_sim_output/')
pop_data <- read.delim('hs_len20_vec1000_yrs250_variant_frequencies.tsv', header = T, colClasses = c('character','numeric','numeric','character')) # upload
run_pop <- pop_data %>% filter(run_id == 'run_29') # choosing a run to work with

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
# sil <- silhouette(clusters, dist = dist_matrix) # re-compute sil score for plot
# plot(sil) # simple plot

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
       subtitle = paste("Selection: Host specialization\nk = ",optimal_k,
                        '\nsilhouette score = ',round(x = max_sil_score,digits = 3),
                        '\nMean in cluster distance = ',in_cluster,
                        '\nMean cross cluster distance = ',cross_cluster)) +
  theme(axis.text = element_blank(),axis.title = element_blank(),panel.grid = element_blank())

### pca analysis ###
pca_result <- prcomp(binary_matrix_numeric, center = TRUE, scale. = F) # set to TRUE if cross reac, FALSE if host specialization --- not sure why
pca_data <- as.data.frame(pca_result$x[, 1:2])  # First 2 components
colnames(pca_data) <- c('PC1', 'PC2')
pca_data$cluster <- as.factor(clusters)

pal <- wes_palette("Zissou1", optimal_k, type = "continuous")
ggplot(pca_data, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(show.legend = FALSE, size = 3, alpha = 0.9) +
  labs(title = "PCA for variant clusters",
       subtitle = 'Selection: Host specialization',
       x = "PC1", y = "PC2") +
  theme_bw() +
  geom_hline(yintercept = 0, alpha = 0.5, linewidth = 0.2) +
  geom_vline(xintercept = 0, alpha = 0.5, linewidth = 0.2) +
  scale_fill_gradientn(colours = pal)


### t-SNE analysis ###
tsne_result <- Rtsne(dist_matrix, dims = 2, pca = TRUE, perplexity = 5, check_duplicates = FALSE)

# combine t-SNE and cluster labels
tsne_data <- data.frame(tsne_result$Y)
tsne_data$cluster <- as.factor(clusters)
pal <- wes_palette("Zissou1", optimal_k, type = "continuous")
# Plot the t-SNE results with hierarchical clusters
ggplot(tsne_data, aes(x = X1, y = X2, color = cluster)) +
  geom_point(show.legend = FALSE, size = 3, alpha = 0.9) +
  labs(title = "t-SNE for variant clusters",
       subtitle = 'Selection: Host specialization',
       x = "t-SNE Dimension 1", y = "t-SNE Dimension 2") +
  theme_bw() +
  #theme(panel.grid = element_blank(), panel.border = element_rect()) +
  geom_hline(yintercept = 0, alpha = 0.5, linewidth = 0.2) +
  geom_vline(xintercept = 0, alpha = 0.5, linewidth = 0.2) +
  scale_fill_gradientn(colours = pal)

# heatmap of pairwise distances with clustering
zzz <- wes_palette("Zissou1", 20, type = "continuous")
breaks_vals <- seq(0, 20, length.out = 21)
pheatmap(
  dist_matrix,               
  cluster_rows = hc_result,      
  cluster_cols = hc_result,      
  display_numbers = FALSE,       
  color = viridis(20),
  #color = zzz,
  breaks = breaks_vals,
  border_color = NA,
  main = "Final population p-w antigen distances\nSelection: Host specialization",
  scale = "none",               
  fontsize = 8,                 
  show_rownames = FALSE,          
  show_colnames = FALSE,          
  fontsize_row = 8,              
  fontsize_col = 8,
  #cellwidth = 3, 
  #cellheight = 3,
  #treeheight_row = 0,           
  #treeheight_col = 0
)


library(tidyverse)
library(ggplot2)
library(cluster)  
library(Rtsne) 
library(pheatmap)
library(viridis) 
library(stringdist)

setwd('/Users/brandonely/Desktop/bb_strain_model_dev/test_outputs/')

##### distribution of pairwise antigen distances #####
cr <- read.delim('cr_20_1000_300_sim_pairwise_antigen_dists.tsv', header = F)
hs <- read.delim('hs_20_1000_300_sim_pairwise_antigen_dists.tsv', header = F)
neg <- read.delim('neg_20_1000_300_sim_pairwise_antigen_dists.tsv', header = F)
cr$source <- "Cross Reactivity"
hs$source <- "Host Specialization"
neg$source <- "Random Fitness"
df_combined <- bind_rows(cr, hs, neg)

ggplot(df_combined, aes(x = V1, fill = source)) +
  geom_histogram(binwidth = 1,alpha = 0.5) +
  facet_wrap(~ source, scales = "free_y") + 
  theme_bw() +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title = "Distribution of pairwise antigen distances", x = "Distance", y = "Frequency") +
  theme(legend.position = "none", panel.grid = element_blank())

# random sampling 
# cr_sampled <- cr[sample(nrow(cr), 1000), ]
# hs_sampled <- hs[sample(nrow(hs), 1000), ]
# neg_sampled <- neg[sample(nrow(neg), 1000), ]
# sampled_combined <- bind_rows(cr_sampled, hs_sampled, neg_sampled)
# 
# ggplot(sampled_combined, aes(x = V1, fill = source)) +
#   geom_histogram(binwidth = 1,alpha = 0.5) +
#   facet_wrap(~ source, scales = "fixed") +
#   theme_bw() +
#   scale_y_continuous(expand = c(0, 0)) +
#   labs(title = "Sampling of pairwise antigen distances", x = "Distance", y = "Frequency") +
#   theme(legend.position = "none", panel.grid = element_blank())

##### cluster analysis #####
data <- read.delim('cr_20_1000_300_sim_final_antigens.tsv', header = FALSE, colClasses = 'character')
#head(data)

# create binary matrix
binary_matrix <- do.call(rbind, strsplit(data$V1, split = ""))
binary_matrix_str <- apply(binary_matrix, 1, paste, collapse = "")
# create distance matrix
dist_matrix <- stringdistmatrix(binary_matrix_str, binary_matrix_str, method = "hamming")
# convert to a distance object, required for hclust function 
dist_obj <- as.dist(dist_matrix)

# perform clustering and print dendrogram
hc_result <- hclust(dist_obj, method = "complete")
plot(hc_result, main = "Hierarchical Clustering")

# defining the range of clusters to calculate sil scores for - not sure about this???
if (length(binary_matrix_str) < 20) {
  z <- length(binary_matrix_str) - 1
} else {
  z <- 0.1 * length(binary_matrix_str)
}

# get sil scores for clusters; does for clusters k=2 to k= z, as defined in lines above
sil_scores <- sapply(2:z, function(k) {
  clusters <- cutree(hc_result, k)
  sil <- silhouette(clusters, dist_matrix)
  mean(sil[, 3])  # average silhouette score
})

# define silhouette score for sim based on chosen number of clusters
#max_sil_score <- max(sil_scores)
# define the optimal number of clusters 
optimal_k <- which.max(sil_scores) + 1

# cut dendrogram to define clusters
clusters <- cutree(hc_result, k = optimal_k) # change k for cluster number
# compute the silhouette score
sil <- silhouette(clusters, dist = dist_matrix)
# plot the silhouette scores
plot(sil)

##### t-SNE analysis #####
tsne_result <- Rtsne(dist_matrix, dims = 2, pca = TRUE, perplexity = 5, check_duplicates = FALSE)

# combine t-SNE and cluster labels
tsne_data <- data.frame(tsne_result$Y)
tsne_data$cluster <- as.factor(clusters)

#colors <- brewer.pal(min(optimal_k, 12), "Set1")
# Plot the t-SNE results with hierarchical clusters
ggplot(tsne_data, aes(x = X1, y = X2, color = cluster)) +
  #geom_point(size = 4) +
  geom_point(show.legend = FALSE) +
  labs(title = "t-SNE visualization with Hierarchical Clustering",
       subtitle = 'Host speciliazation, antigen len=20, vector pop=1000, years=50, mut=0.01',
       x = "t-SNE Dimension 1", y = "t-SNE Dimension 2") +
  theme_minimal() +
  scale_color_viridis(discrete = TRUE)
  #scale_color_manual(values = colors)


### pairwise distance heatmap ###
pheatmap(
  dist_matrix,               
  cluster_rows = hc_result,      
  cluster_cols = hc_result,      
  display_numbers = FALSE,       
  color = viridis(100), 
  border_color = NA,
  main = "Pairwise Antigen Distances",  
  scale = "none",               
  fontsize = 10,                 
  show_rownames = FALSE,          
  show_colnames = FALSE,          
  fontsize_row = 8,              
  fontsize_col = 8,
  #cellwidth = 3, 
  #cellheight = 3,
  #treeheight_row = 0,           
  treeheight_col = 0
)





library(tidyverse)
library(ggplot2)
library(cluster)  
library(Rtsne) 

setwd('/Users/brandonely/Desktop/ospC_strain_model/test_outputs/')
data <- read.delim('test2_sim_final_antigens.tsv', header = FALSE, colClasses = 'character')
#head(data)

# tranform binary strings to numeric vectors
binary_matrix <- do.call(rbind, strsplit(data$V1, split = ""))
binary_matrix <- apply(binary_matrix, 2, as.numeric)

# k-means clustering 
kmeans_result <- kmeans(binary_matrix, centers = 2)  # change centers for num of clusters
# t-SNE analysis
tsne_result <- Rtsne(binary_matrix, dims = 2, pca = FALSE, perplexity = 12, check_duplicates = FALSE)

# combine t-SNE and cluster labels
tsne_data <- data.frame(tsne_result$Y)
tsne_data$cluster <- as.factor(kmeans_result$cluster)

# plot
ggplot(tsne_data, aes(x = X1, y = X2, color = cluster)) +
  geom_point(size = 4) +
  labs(title = "t-SNE with K-means clusters",
       x = "t-SNE Dimension 1", y = "t-SNE Dimension 2") +
  theme_minimal() +
  scale_color_manual(values = c("purple", "cyan3", "orange", 'red', 'yellow', 'green'))


# heirarchical clustering

# perform hierarchical clustering; euclidean distance, complete linkage
dist_matrix <- dist(binary_matrix, method = "euclidian")  
hclust_result <- hclust(dist_matrix, method = "ward")

# cut dendrogram to create clusters
clusters <- cutree(hclust_result, k = 4) # change k for cluster number

# t-SNE analysis
tsne_result <- Rtsne(binary_matrix, dims = 2, pca = TRUE, perplexity = 8, check_duplicates = FALSE)

# combine t-SNE and cluster labels
tsne_data <- data.frame(tsne_result$Y)
tsne_data$cluster <- as.factor(clusters)

# Plot the t-SNE results with hierarchical clusters
ggplot(tsne_data, aes(x = X1, y = X2, color = cluster)) +
  geom_point(size = 4) +
  labs(title = "t-SNE visualization with Hierarchical Clustering",
       subtitle = 'Host speciliazation, antigen len=20, vector pop=1000, years=50, mut=0.01',
       x = "t-SNE Dimension 1", y = "t-SNE Dimension 2") +
  theme_minimal() +
  scale_color_manual(values = c("purple", "cyan3", "orange", 'red', 'yellow', 'green','blue'))


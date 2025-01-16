# functions for sim analysis

# hamming distance, self explanatory
hamming_distance <- function(x, y) {
  sum(strsplit(x, "")[[1]] != strsplit(y, "")[[1]])
}

# avg within group dist function; data needs 2 cols: variant, cluster_label; nk is optimal k
InGroupDistance <- function(data) { 
  nk <- length(unique(data$cluster_label))
  in_dist <- list()
  index <- 1
  for (k in 1:nk) {
    group <- data %>%
      filter(cluster_label == k) %>%
      pull(variant)
    for (item in group) {
      for (n in 1:length(group)) {
        if (item != group[n]) {
          in_dist[index] <- hamming_distance(item, group[n])
          index <- index + 1
        }
      }
    }
  }
  return(mean(unlist(in_dist)))
}

# avg across group dist function
OutGroupDistance <- function(data) {
  nk <- length(unique(data$cluster_label))
  out_dist <- list() 
  index <- 1  
  for (k1 in 1:(nk - 1)) {
    for (k2 in (k1 + 1):nk) {  
      group1 <- data %>%
        filter(cluster_label == k1) %>%
        pull(variant) 
      group2 <- data %>%
        filter(cluster_label == k2) %>%
        pull(variant) 
      
      for (item1 in group1) {
        for (item2 in group2) {
          out_dist[[index]] <- hamming_distance(item1, item2)
          index <- index + 1
        }
      }
    }
  }
  return(mean(unlist(out_dist)))
}

# get pairwise hamming distances for final population 
# input - final pop counts/frequency output from single run of sim, returns df for making histogram
final_antigen_distances <- function(data) {
  # part 1 - get all non-unique variants based on counts into a list
  variants <- list()
  i <- 1
  for (item in data$variant) {
    multiplied <- replicate(data$counts[i], item)
    variants[[i]] <- multiplied
    i <- i + 1
  }
  
  # part 2 - calculate pairwise hamming distances for all
  final_pop <- unlist(variants, recursive = F)
  pairwise_distances <- list()
  n = length(final_pop)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      # Calculate the Hamming distance between each pair
      dist <- hamming_distance(final_pop[[i]], final_pop[[j]])
      # Store the pairwise distance
      pairwise_distances[[length(pairwise_distances) + 1]] <- dist
    }
  }
  pairwise_distances_vector <- unlist(pairwise_distances)
  pairwise_distances_df <- data.frame(distance = pairwise_distances_vector)
  return(pairwise_distances_df)
}

# collects sil score, optimal k for clustering, avg wtihin and cross cluster distances
batch_clustering_data <- function(data) {
  
  # create binary matrix and distance matrix
  binary_matrix <- do.call(rbind, strsplit(data$variant, split = ""))
  binary_matrix_str <- apply(binary_matrix, 1, paste, collapse = "")
  binary_matrix_numeric <- matrix(as.numeric(binary_matrix), nrow = nrow(binary_matrix), ncol = ncol(binary_matrix))
  dist_matrix <- stringdistmatrix(binary_matrix_str, binary_matrix_str, method = "hamming")
  dist_obj <- as.dist(dist_matrix) # convert to a distance object, required for hclust function 
  
  hc_result <- hclust(dist_obj, method = "complete") # perform clustering
  max_k <- length(binary_matrix_str) -1 # set total number of clusters allowed, set to n-1
  
  sil_scores <- sapply(2:max_k, function(k) { # get sil scores for each k value
    clusters <- cutree(hc_result, k)
    sil <- silhouette(clusters, dist_matrix)
    mean(sil[, 3])
  })
  
  # define all data to be collected and return as a df
  max_sil_score <- max(sil_scores) # highest sil score
  optimal_k <- which.max(sil_scores) + 1 # final k value based on highest sil score
  clusters <- cutree(hc_result, k = optimal_k) # cut dendrogram using optimal k
  data$cluster_label <- clusters # add cluster labels to df
  in_dist <- InGroupDistance(data) # in cluster distance
  out_dist <- OutGroupDistance(data) # cross cluster distance
  
  return(data.frame(
    rund_id=data$run_id[1],
    n_variants= nrow(data),
    sil_score=max_sil_score, 
    k_clusters=optimal_k, 
    avg_in_dist=in_dist, 
    avg_out_dist=out_dist))
}
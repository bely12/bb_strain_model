# functions for sim analysis

# hamming distance, self explanatory
hamming_distance <- function(x, y) {
  sum(strsplit(x, "")[[1]] != strsplit(y, "")[[1]])
}

# avg within group dist function; data needs 2 cols: strain, cluster_label; nk is optimal k
InGroupDistance <- function(data) { 
  nk <- length(unique(data$cluster_label))
  in_dist <- list()
  index <- 1
  for (k in 1:nk) {
    group <- data %>%
      filter(cluster_label == k) %>%
      pull(strain)
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
        pull(strain) 
      group2 <- data %>%
        filter(cluster_label == k2) %>%
        pull(strain) 
      
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
  for (item in data$strain) {
    multiplied <- replicate(data$count[i], item)
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
  binary_matrix <- do.call(rbind, strsplit(data$strain, split = ""))
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
  
  #data$cluster_label <- clusters
  x <- data %>%
    group_by(cluster_label) %>%
    summarise(
      totals = sum(count),
      run_id = run_id
    )
  x$freq <- x$totals/sum(x$totals) * 100

  H <- vegan::diversity(x$freq, index = "shannon", base = exp(1))  
  N <- sum(x$freq > 0)
  E <- H / log(N)
  
  w <- data %>% select(strain, cluster_label, count, freq, run_id, selection, gene_type, gen_fitness, rec_rate, mut_rate, host_pops, host_pop_stability) #added, might need to take out
  
  df = data.frame(
    run_id=data$run_id[1],
    n_variants= nrow(data),
    sil_score=max_sil_score, 
    k_clusters=optimal_k, 
    avg_in_dist=in_dist, 
    avg_out_dist=out_dist,
    ShanE = E,
    selection = unique(data$selection),
    rec_rate = unique(data$rec_rate),
    mut_rate = unique(data$mut_rate),
    gen_fit = unique(data$gen_fitness),
    host_pops = unique(data$host_pops),
    host_pop_stability = unique(data$host_pop_stability)
  )
  
  # return(data.frame(
  #   run_id=data$run_id[1],
  #   n_variants= nrow(data),
  #   sil_score=max_sil_score, 
  #   k_clusters=optimal_k, 
  #   avg_in_dist=in_dist, 
  #   avg_out_dist=out_dist,
  #   ShanE = E,
  #   selection = unique(data$selection),
  #   rec_rate = unique(data$rec_rate),
  #   mut_rate = unique(data$mut_rate),
  #   gen_fit = unique(data$gen_fitness),
  #   host_pops = unique(data$host_pops),
  #   host_pop_stability = unique(data$host_pop_stability)
  #   ))
  
  return(list('cluster_res' = df, 'freq_table' = x, 'variant_cluster_labels' = w))
}


batch_clustering_data_V2 <- function(data) {
  
  # create binary matrix and distance matrix
  binary_matrix <- do.call(rbind, strsplit(data$immune, split = ""))
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
  
  #data$cluster_label <- clusters
  x <- data %>%
    group_by(cluster_label) %>%
    summarise(
      totals = sum(count),
      run_id = run_id
    )
  x$freq <- x$totals/sum(x$totals) * 100
  
  H <- vegan::diversity(x$freq, index = "shannon", base = exp(1))  
  N <- sum(x$freq > 0)
  E <- H / log(N)
  
  w <- data %>% select(immune, adaptive, strain, cluster_label, count, freq, run_id, selection, gene_type, gen_fitness, rec_rate, mut_rate, host_pops, host_pop_stability) #added, might need to take out
  
  df = data.frame(
    run_id=data$run_id[1],
    n_variants= nrow(data),
    sil_score=max_sil_score, 
    k_clusters=optimal_k, 
    avg_in_dist=in_dist, 
    avg_out_dist=out_dist,
    ShanE = E,
    selection = unique(data$selection),
    rec_rate = unique(data$rec_rate),
    mut_rate = unique(data$mut_rate),
    gen_fit = unique(data$gen_fitness),
    host_pops = unique(data$host_pops),
    host_pop_stability = unique(data$host_pop_stability)
  )
  
  # return(data.frame(
  #   run_id=data$run_id[1],
  #   n_variants= nrow(data),
  #   sil_score=max_sil_score, 
  #   k_clusters=optimal_k, 
  #   avg_in_dist=in_dist, 
  #   avg_out_dist=out_dist,
  #   ShanE = E,
  #   selection = unique(data$selection),
  #   rec_rate = unique(data$rec_rate),
  #   mut_rate = unique(data$mut_rate),
  #   gen_fit = unique(data$gen_fitness),
  #   host_pops = unique(data$host_pops),
  #   host_pop_stability = unique(data$host_pop_stability)
  #   ))
  
  return(list('cluster_res' = df, 'freq_table' = x, 'variant_cluster_labels' = w))
}


# single run clustering; returns a list of 3 items: 1-df of cluster summary stats, 2-dist object 3-cluster res, 4-freq table for clusters 
single_run_clustering <- function(data, strict = 'no', threshold = 2) {
  
  # create binary matrix and distance matrix
  binary_matrix <- do.call(rbind, strsplit(data$strain, split = ""))
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
  
  ### just added - not sure I would use this it doesn't seem to work the intended way
  if (strict == 'yes') {
    clusters <- cutree(hc_result, h = threshold)
  } else {
    clusters <- cutree(hc_result, k = optimal_k)
  }
  ###
  
  #clusters <- cutree(hc_result, k = optimal_k) # cut dendrogram using optimal k
  data$cluster_label <- clusters # add cluster labels to df
  in_dist <- InGroupDistance(data) # in cluster distance
  out_dist <- OutGroupDistance(data) # cross cluster distance
  
  #data$cluster_label <- clusters
  w <- data %>% select(strain, cluster_label,count, freq)
  
  x <- data %>%
    group_by(cluster_label) %>%
    summarise(
      totals = sum(count)
    )
  x$freq <- x$totals/sum(x$totals) * 100
  H <- diversity(x$freq, index = "shannon", base = exp(1))  
  N <- sum(x$freq > 0)
  E <- H / log(N)
  
  df <- data.frame(
    run_id=data$run_id[1],
    n_variants= nrow(data),
    sil_score=max_sil_score, 
    k_clusters=optimal_k, 
    avg_in_dist=in_dist, 
    avg_out_dist=out_dist,
    ShanE = E)
  
  return(list('stats_df' = df, 'dist_matrix' = dist_obj, 'cluster_res' = hc_result, 'freq_table' = x, 'variant_cluster_labels' = w))
}



# function for getting largest connected component for each igraph object
get_lcc_fraction <- function(g) {
  comps <- components(g)
  largest_size <- max(comps$csize)
  return(largest_size / vcount(g))  # normalized for differing n nodes
}

# network analysis
antigen_network <- function(data, threshold) {
  # assign edges 
  edges <- c()
  for (i in 1:(nrow(data)- 1)) {
    for (j in (i +1):nrow(data)) {
      dist <- hamming_distance(data$strain[i], data$strain[j])
      if (dist <= threshold*nchar(data$strain[i])) {
        edges <- c(edges, data$strain[i], data$strain[j])
      }
    }
  }
  
  # contstruct the network
  g <- make_graph(edges, directed = FALSE)
  
  # compute metrics and store in df
  comm <- cluster_louvain(g)
  large_comp <- get_lcc_fraction(g)
  return(data.frame(run_id = data$run_id[1],
                    threshold = threshold,
                    lrg_comp = large_comp, 
                    edge_dens = edge_density(g),
                    mean_degree = mean(degree(g)),
                    #clust_coef = transitivity(g, type = "global"),
                    modul = modularity(comm),
                    selection = unique(data$selection),
                    rec_rate = unique(data$rec_rate),
                    mut_rate = unique(data$mut_rate),
                    #host_pops = unique(data$host_pops),
                    #host_pop_stability = unique(data$host_pop_stability),
                    gen_fitness = unique(data$gen_fitness)))
}

# visualize antigen network
antigen_network_viz <- function(data, threshold) {
  # assign edges 
  edges <- c()
  for (i in 1:(nrow(data)- 1)) {
    for (j in (i +1):nrow(data)) {
      dist <- hamming_distance(data$strain[i], data$strain[j])
      if (dist <= threshold) {
        edges <- c(edges, data$strain[i], data$strain[j])
      }
    }
  }
  
  # contstruct the network
  g <- make_graph(edges, directed = FALSE)
  
  # set size and color of nodes according to frequency of variants
  V(g)$frequency <- data$freq[match(V(g)$name, data$strain)]
  V(g)$size <- log1p(V(g)$frequency) * 5 
  freqs <- V(g)$frequency
  palette_func <- colorRampPalette(c("lightblue", "darkred"))
  colors <- palette_func(100)
  freq_scaled <- as.integer( (freqs - min(freqs)) / (max(freqs) - min(freqs)) * 99 ) + 1
  V(g)$color <- colors[freq_scaled]
  
  # plot
  return(plot(g,
              vertex.label = NA,
              vertex.frame.color = "gray30",
              layout = layout_with_fr))
}



cluster_host_assoc <- function(hi, labels, low_thresh = 0.33, high_thresh = 0.67) {
  
  # ensure no duplicates in labels, but keep run_id to avoid cross-run contamination
  labels_clean <- labels %>%
    distinct(run_id, strain, cluster_label)
  
  # collapse duplicates in hi per run_id + immune
  hi_collapsed <- hi %>%
    group_by(gen_fitness, run_id, host_pop_stability, selection, immune) %>%
    summarise(
      rodent_count   = sum(rodent_count, na.rm = TRUE),
      bird_count     = sum(bird_count, na.rm = TRUE),
      combined_count = sum(combined_count, na.rm = TRUE),
      .groups = "drop"
    )
  
  # join with labels by run_id and immune strain
  joined <- hi_collapsed %>%
    inner_join(labels_clean, by = c("run_id", "immune" = "strain"))
  
  # aggregate to cluster level per run
  df_clustered <- joined %>%
    group_by(gen_fitness, run_id, host_pop_stability, cluster_label, selection) %>%
    summarise(
      rodent_count   = sum(rodent_count),
      bird_count     = sum(bird_count),
      combined_count = sum(combined_count),
      .groups = "drop_last"
    ) %>%
    mutate(
      host_assoc = bird_count / combined_count,
      classification = case_when(
        host_assoc <= low_thresh ~ "rodent_specialist",
        host_assoc >= high_thresh ~ "bird_specialist",
        TRUE ~ "generalist"
      )
    )
  
  return(df_clustered)
}
library(Biostrings)
library(igraph)
library(tidyverse)


##### functions #####

# mutate sequences
mutate_seq <- function(seq, mut_rate = 0.02) {
  s <- strsplit(seq, "")[[1]]
  n_mut <- round(length(s) * mut_rate)
  if (n_mut > 0) {
    pos <- sample(1:length(s), n_mut)
    bases <- c("A","T","C","G")
    s[pos] <- sapply(s[pos], function(b) sample(setdiff(bases, b), 1))
  }
  paste0(s, collapse="")
}

# hamming distance for ospC nucleotide seqs
hamming_nt <- function(s1, s2) {
  sum(strsplit(s1,"")[[1]] != strsplit(s2,"")[[1]])
}

# construct network graph
antigen_network_viz_ospC <- function(data, allowed_muts = 4) {
  
  edges_list <- data.frame(from=character(),
                           to=character(),
                           stringsAsFactors=FALSE)
  
  n <- nrow(data)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      h <- hamming_nt(data$seq[i], data$seq[j])
      if (h <= allowed_muts) {
        edges_list <- rbind(edges_list,
                            data.frame(from = data$node_name[i],
                                       to   = data$node_name[j],
                                       stringsAsFactors=FALSE))
      }
    }
  }
  
  g <- graph_from_data_frame(edges_list,
                             directed=FALSE,
                             vertices=data.frame(name=data$node_name))
  
  V(g)$frequency <- data$freq[match(V(g)$name, data$node_name)]
  V(g)$size <- log1p(V(g)$frequency) * 6
  
  freqs <- V(g)$frequency
  pal <- colorRampPalette(c("lightblue","darkred"))
  cols <- pal(100)
  freq_scaled <- as.integer((freqs - min(freqs)) / (max(freqs) - min(freqs)) * 99) + 1
  V(g)$color <- cols[freq_scaled]
  
  plot(g,
       vertex.label=NA,
       vertex.frame.color="gray30",
       layout=layout_with_fr)
  
  return(g)
}

##### import ospC seqs into df #####
setwd('/Users/brandonely/Desktop')

# import ospC multifasta
ospC <- readDNAStringSet("ospC-ref-nuc.fas")

# format ospC pool into dataframe
ospC_df <- data.frame(
  allele = names(ospC),
  seq    = as.character(ospC),
  stringsAsFactors = FALSE)

##### sample and mutate ospC seqs to create natural pop #####
# parameters for sample + mutate to create realistic natural population 
set.seed(123)
n_sample <- 50
mut_rate <- 0.02

# sample from ospC pool
sampled <- data.frame(
  allele = sample(ospC_df$allele, n_sample, replace=TRUE),
  seq = NA,
  stringsAsFactors=FALSE)

# mutate sampled seqs
sampled$seq <- mapply(function(a) {
  orig <- ospC_df$seq[ospC_df$allele==a]
  mutate_seq(orig, mut_rate)
}, sampled$allele)

# compute frequencies of sampled alleles (igraph expects this)
sampled_summary <- sampled %>%
  count(allele, seq, name="freq")

# unique node name for each variant (avoids problems with igraph)
sampled_summary$node_name <- paste0(sampled_summary$allele, "_", seq_len(nrow(sampled_summary)))

##### build network #####
# scaled threshold parameters
model_threshold <- 0.5  
ospC_length <- 635      
allowed_muts <- round(model_threshold * ospC_length)

# initiate file to save to
png("ospC_network_50.png", width=800, height=600)

# build network
antigen_network_viz_ospC(sampled_summary, allowed_muts = allowed_muts)

# save/reset
dev.off()

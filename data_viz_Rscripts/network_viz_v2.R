library(tidyverse)
library(ggpubr)
library(igraph)
library(broom)
setwd('/Users/brandonely/Desktop/bb_strain_model_dev/R_analysis_scripts/')
source("sim_analysis_functions.R")


# upload raw data to build new network
setwd('/Users/brandonely/Desktop/bb_strain_model_dev/evo_sim/results_tables/sim_output/adaptive_results/')
df <- read.delim('adaptive_ExpMod_variant_freqs_sampled.tsv', header = T, colClasses = c('character','numeric','numeric','character', 'character', 'character','numeric','numeric','character'))

# get rid of ghost variants 
df <- df %>% filter(freq >= 1)

# choose specific sim run to build network for
df <- df %>% filter(run_id == 'run_25')

# visualize the network
setwd('/Users/brandonely/Desktop/bb_strain_model_dev/evo_sim/results_tables/new_plots/')
png("adaptive_network_run25_thresh50_xxx.png", width = 800, height = 600)
antigen_network_viz(df, threshold = 10)
dev.off()

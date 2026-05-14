# network analysis
setwd('/Users/brandonely/Desktop/bb_strain_model_dev/R_analysis_scripts/')
library(tidyverse)
library(ggpubr)
library(igraph)
library(broom)
source("sim_analysis_functions.R")

##### network analysis on all runs #####

#header_names <- c("strain", "count", "freq", "run_id", "selection", "gene_type", "rec_rate", "mut_rate", "gen_fitness")

setwd('/Users/brandonely/Desktop/bb_strain_model_dev/evo_sim/results_tables/final_results/adaptive_results/')
df <- read.delim('adaptive_LowGenFit_variant_freqs_sampled.tsv', header = T, colClasses = c('character','numeric','numeric','character', 'character', 'character','numeric','numeric','character'))
df <- df %>% filter(rec_rate == 0.01)
df$host_pop_stability <- 'stable'

# filter to remove "ghost variants"
df <- df %>% filter(freq >= 1)

# check to see how many variants I have per run 
x <- df %>% group_by(run_id) %>% summarise(counts = n())
x <- x %>% filter(counts >=10) # keep only runs with at least 15 variants

# to cut down on the number of runs I do network analysis on, only if neccessary
#run_filter <- sample(x$run_id, 25, replace = FALSE)
df <- df %>% filter(run_id %in% x$run_id)

# look at network stats for all runs and a range of thresholds 
catcher = list()
index = 1
proporotional_thresholds <- seq(0.05, 0.5, by = 0.05) # on a 20 bit strings, gives 1:10 
for (run in unique(df$run_id)) {
  temp_df <- df %>% filter(run_id == run)
  for (i in proporotional_thresholds) {
    catcher[[index]] <- antigen_network(temp_df, threshold = i)
    index = index + 1
    print(index)
  }
}

# collect results for all conditions
no_net <- bind_rows(catcher)
no_net <- no_net %>% filter(!if_any(everything(), ~ is.na(.) | is.nan(.) | is.infinite(.)))
no_net$host_pop_stability <- 'stable'

imm_net <- bind_rows(catcher)
imm_net <- imm_net %>% filter(!if_any(everything(), ~ is.na(.) | is.nan(.) | is.infinite(.)))
imm_net$host_pop_stability <- 'stable'

hyb_ExpMod_net <- bind_rows(catcher)
hyb_ExpMod_net <- hyb_ExpMod_net %>% filter(!if_any(everything(), ~ is.na(.) | is.nan(.) | is.infinite(.)))

hyb_Linear_net <- bind_rows(catcher)
hyb_Linear_net <- hyb_Linear_net %>% filter(!if_any(everything(), ~ is.na(.) | is.nan(.) | is.infinite(.)))

hyb_Low_net <- bind_rows(catcher)
hyb_Low_net <- hyb_Low_net %>% filter(!if_any(everything(), ~ is.na(.) | is.nan(.) | is.infinite(.)))

adpt_ExpMod_net <- bind_rows(catcher)
adpt_ExpMod_net <- adpt_ExpMod_net %>% filter(!if_any(everything(), ~ is.na(.) | is.nan(.) | is.infinite(.)))

adpt_Linear_net <- bind_rows(catcher)
adpt_Linear_net <- adpt_Linear_net %>% filter(!if_any(everything(), ~ is.na(.) | is.nan(.) | is.infinite(.)))

adpt_Low_net <- bind_rows(catcher)
adpt_Low_net <- adpt_Low_net %>% filter(!if_any(everything(), ~ is.na(.) | is.nan(.) | is.infinite(.)))

# concat results to single table and save
networks_df <- bind_rows(adpt_ExpMod_net, adpt_Linear_net, adpt_Low_net)
networks_df$host_pop_stability <- 'stable'
setwd('/Users/brandonely/Desktop/bb_strain_model_dev/evo_sim/results_tables/network_analysis_results/')
#write.table(networks_df, file='adaptive_selection_network_results_freq1_10min_cutoffs.tsv', quote=FALSE, sep='\t', row.names = F)





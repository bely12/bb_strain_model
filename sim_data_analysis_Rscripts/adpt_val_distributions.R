setwd('/Users/brandonely/Desktop/bb_strain_model_dev/R_analysis_scripts/')
library(tidyverse)
library(ggplot2)
library(ggpubr)
source("sim_analysis_functions.R")

header <- c('run_id', 'hs_val', 'selection', 'gene', 'rec_rate', 'mut_rate', 'gen_fit')

# hybrid selection vals
setwd('/Users/brandonely/Desktop/bb_strain_model_dev/evo_sim/results_tables/sim_output/hybrid_results/')

vals <- read.delim('hybrid_ExpMod_adpt_vals_sampled.tsv', header = F, col.names = header)
vals$host_pop_stability <- 'stable'
sampled_vals <- sample_n(tbl = vals, size = 1000, replace = F)
rm(vals)

vals <- read.delim('hybrid_LinearGenFit_adpt_vals_sampled.tsv', header = F, col.names = header)
vals$host_pop_stability <- 'stable'
sampled_vals2 <- sample_n(tbl = vals, size = 1000, replace = F)
rm(vals)

vals <- read.delim('hybrid_LowGenFit_adpt_vals_sampled.tsv', header = F, col.names = header)
vals$host_pop_stability <- 'stable'
sampled_vals3 <- sample_n(tbl = vals, size = 1000, replace = F)
rm(vals)

vals <- read.delim('hybrid_HighGenFit_adpt_vals_sampled.tsv', header = F, col.names = header)
vals$host_pop_stability <- 'stable'
sampled_vals4 <- sample_n(tbl = vals, size = 1000, replace = F)
rm(vals)

# adaptive selection vals
setwd('/Users/brandonely/Desktop/bb_strain_model_dev/evo_sim/results_tables/sim_output/adaptive_results/')

vals <- read.delim('adaptive_ExpMod_adpt_vals_sampled.tsv', header = F, col.names = header)
vals$host_pop_stability <- 'stable'
sampled_vals5 <- sample_n(tbl = vals, size = 1000, replace = F)
rm(vals)

vals <- read.delim('adaptive_LinearGenFit_adpt_vals_sampled.tsv', header = F, col.names = header)
vals$host_pop_stability <- 'stable'
sampled_vals6 <- sample_n(tbl = vals, size = 1000, replace = F)
rm(vals)

vals <- read.delim('adaptive_LowGenFit_adpt_vals_sampled.tsv', header = F, col.names = header)
vals$host_pop_stability <- 'stable'
sampled_vals7 <- sample_n(tbl = vals, size = 1000, replace = F)
rm(vals)

vals <- read.delim('adaptive_HighGenFit_adpt_vals_sampled.tsv', header = F, col.names = header)
vals$host_pop_stability <- 'stable'
sampled_vals8 <- sample_n(tbl = vals, size = 1000, replace = F)
rm(vals)



# concatenate 
vals <- bind_rows(sampled_vals, sampled_vals2, sampled_vals3, sampled_vals4, sampled_vals5, sampled_vals6, sampled_vals7, sampled_vals8)

setwd('/Users/brandonely/Desktop/bb_strain_model_dev/evo_sim/results_tables/host_association_results/')
write.table(vals, 'sampled_adaptive_trait_vals_ALL.tsv', quote = F, row.names = F, sep = '\t')
# add category for distinct parameters
vals$category <- paste(vals$selection, vals$gen_fit)

# plot
vals %>%
  ggplot(aes(x = hs_val, fill = category)) +
  geom_histogram(bins = 6, color = 'black') +
  theme_bw() +
  facet_wrap(~category, ncol = 4)

setwd('/Users/brandonely/Desktop/bb_strain_model_dev/evo_sim/results_tables/new_plots/')
#ggsave('hs_trait_val_distribution.png', width = 6, height = 6)
  
  

# sim antigen distance distribution analysis 
setwd('/Users/brandonely/Desktop/bb_strain_model_dev/')
library(tidyverse)
library(ggplot2)
source("sim_analysis_functions.R")

# upload data from sim output, format, and sample 
setwd('/Users/brandonely/Desktop/bb_strain_model_dev/new_sim_output/')
cr_ad <- read.delim('cr_len20_vec1000_yrs250_sampled_pairwise_antigen_dists.tsv', header = F) # upload
cr_sampled <- cr_ad[sample(nrow(cr_ad), 5000), ] # sample
cr_sampled$source <- 'cross reactivity' # label 
hs_ad <- read.delim('hs_len20_vec1000_yrs250_sampled_pairwise_antigen_dists.tsv', header = F) # upload
hs_sampled <- hs_ad[sample(nrow(hs_ad), 5000), ] # sample
hs_sampled$source <- 'host specialization' # label 
df_combined <- bind_rows(cr_sampled, hs_sampled) # combine

# plot single histogram 
# ggplot(cr_sampled, aes(x = V2)) +
#   geom_histogram(binwidth = 1, alpha = 0.5, fill = 'blue')

# plot multi-panel histogram
ggplot(df_combined, aes(x = V2, fill = source)) +
  geom_histogram(binwidth = 1,alpha = 0.5) +
  facet_wrap(~ source, scales = 'free_y', nrow = 2) +
  theme_bw() +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title = "Distribution of pairwise antigen distances", x = "Distance", y = "Frequency") +
  theme(legend.position = "none", panel.grid = element_blank())

# get antigen distance distributions for final population 
pop_data <- read.delim('cr_len20_vec1000_yrs250_variant_frequencies.tsv', header = T, colClasses = c('character','numeric','numeric','character')) #upload data 

# look at final number of variants in pop, i think i have this data elsewhere though...
# result <- pop_data %>%
#   group_by(run_id) %>%
#   tally() %>%
#   arrange(desc(n)) 
# mean(result$n)
# max(result$n)

run_pop <- pop_data %>% filter(run_id == 'run_79') # choosing a run to work with
ant_dists <- final_antigen_distances(run_pop) # custom function for calculating distances
ant_dists$selection <- 'cross reactivity' # label for multi-panel hist or just to keep as df when sampling
sampled_dists <- ant_dists[sample(nrow(ant_dists), 1000), ] # sample the distances 

# plot
ggplot(sampled_dists, aes(x = distance)) +
  geom_histogram(binwidth = 1, fill = "blue", alpha = 0.7) +
  labs(title = "Pairwise antigen distances",
       x = "Distance",
       y = "Frequency") +
  theme_minimal()


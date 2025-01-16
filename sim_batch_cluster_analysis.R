# like sim_cluster_analysis.R but for doing all runs and computing/plotting averages

setwd('/Users/brandonely/Desktop/bb_strain_model_dev/')
library(tidyverse)
library(ggplot2)
library(ggpubr)
source("sim_analysis_functions.R")

# upload data 
setwd('/Users/brandonely/Desktop/bb_strain_model_dev/new_sim_output/')
cr_data <- read.delim('cr_len20_vec1000_yrs250_variant_frequencies.tsv', header = T, colClasses = c('character','numeric','numeric','character'))
hs_data <- read.delim('hs_len20_vec1000_yrs250_variant_frequencies.tsv', header = T, colClasses = c('character','numeric','numeric','character'))

# first condition 
df_catcher <- list() # initialize empty list to store df for each run
index <- 1 # set index for list collection
for (run in unique(cr_data$run_id)) {
  filtered_df <- cr_data %>% filter(run_id == run) # filter df for run
  if (nrow(filtered_df) < 5) { # min number of variants, otherwise will get error 
    next
  }
  df_catcher[[index]] <- batch_clustering_data(filtered_df) # collect data for run
  index <- index+1 # set index for next round
}

# create df
cr_df <- bind_rows(df_catcher) # bind rows to create master df
cr_df$selection <- 'cross_reactivity' # label rows with sim condition 

# second condition 
df_catcher <- list() # initialize empty list to store df for each run
index <- 1 # set index for list collection
for (run in unique(hs_data$run_id)) {
  filtered_df <- hs_data %>% filter(run_id == run) # filter df for run
  if (nrow(filtered_df) < 5) { # min number of variants, otherwise will get error 
    next
  }
  df_catcher[[index]] <- batch_clustering_data(filtered_df) # collect data for run
  index <- index+1 # set index for next round
}

# create df
hs_df <- bind_rows(df_catcher)
hs_df$selection <- 'host_specialization'

# combine the 2 df's and convert to long format for plotting
combined_df <- bind_rows(cr_df, hs_df)
df_long <- combined_df %>%
  select(-1) %>%
  pivot_longer(cols = -selection,  # Keep the "selection" column as it is
               names_to = "variable",  # Name of the new column with numeric column names
               values_to = "value")    # Values for the numeric columns

# plot
ggplot(df_long, aes(x = selection, y = value, fill = selection)) +
  geom_boxplot() +
  facet_wrap(~ variable, scales = "free_y", ncol = 3) +  
  labs(x = "Category", y = "Value") +
  theme_bw() +
  theme(legend.position = "right") +
  stat_compare_means(method = "t.test",
                     hide.ns = TRUE,
                     aes(group = selection), 
                     label = "p.signif")

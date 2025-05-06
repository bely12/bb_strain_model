setwd('/Users/brandonely/Desktop/bb_strain_model_dev/')
library(tidyverse)
library(ggplot2)
library(ggpubr)
source("sim_analysis_functions.R")

### upload data 
setwd('/Users/brandonely/Desktop/bb_strain_model_dev/new_sim_output/set_larger_pops/')

header <- c('run_tag', 'year', 'infection_rate', 'diversity', 'active_strain_count', 'avg_strains_carried', 'unique_lineages')

cr_data <- read.delim('immune_5000_20_500_sim_data.tsv', header = T, col.names = header)
hs_data <- read.delim('hs_5000_20_500_sim_data.tsv', header = T, col.names = header)
ctrl_data <- read.delim('neutral_5000_20_500_sim_data.tsv', header = T, col.names = header)

### filter to remove failed sims - likely not neccessary with large pop sizes
# x_runs_to_remove <- x_data %>% filter(active_strain_count == 0)
# x_data <- x_data %>% filter(!run_tag %in% c(unique(x_runs_to_remove$run_tag)))

### get means by year for each condition
x <- cr_data %>%
  select(-1) %>%
  group_by(year) %>%
  summarise(
    avg_infection_rate = mean(infection_rate),
    #avg_diversity = mean(diversity),
    #avg_antigen_distance = mean(avg_antigen_distance),
    avg_active_strain_count = mean(active_strain_count),
    avg_avg_strains_carried = mean(avg_strains_carried),
    #avg_unique_lineages = mean(unique_lineages)
  )
x$selection <- 'cross reactivity'

y <- hs_data %>%
  select(-1) %>%
  group_by(year) %>%
  summarise(
    avg_infection_rate = mean(infection_rate),
    #avg_diversity = mean(diversity),
    #avg_antigen_distance = mean(avg_antigen_distance),
    avg_active_strain_count = mean(active_strain_count),
    avg_avg_strains_carried = mean(avg_strains_carried),
    #avg_unique_lineages = mean(unique_lineages)
  )
y$selection <- 'host specialization'

z <- ctrl_data %>%
  select(-1) %>%
  group_by(year) %>%
  summarise(
    avg_infection_rate = mean(infection_rate),
    #avg_antigen_distance = mean(avg_antigen_distance),
    avg_active_strain_count = mean(active_strain_count),
    avg_strains_carried = mean(avg_strains_carried),
    #avg_unique_lineages = mean(unique_lineages)
  )
z$selection <- 'control'

df <- bind_rows(x, y, z) # combine df's

df_long <- df %>%
  pivot_longer(cols = -c(year, selection), 
               names_to = "variable", 
               values_to = "mean_value")

### create names for plots in panel
title.labs <- c('Active strains in population',
                'Number of strains carried per tick',
                'Infection rate')
names(title.labs) <- c('avg_active_strain_count',
                       'avg_strains_carried',
                       'avg_infection_rate')

### plot 
ggplot(df_long, aes(x = year, y = mean_value, group = interaction(variable, selection), color = selection)) +
  geom_line() +
  scale_color_manual(values = c("cross reactivity" = "skyblue3", "host specialization" = "red3", "control" = "yellow2")) +
  #geom_point() +
  facet_wrap(~variable, scales = "free_y", ncol = 1, labeller = labeller(variable = title.labs)) +
  theme_bw() +
  labs(x = "Year", y = "Mean Value", title = "Means for simulation data") +
  theme(legend.position = "bottom")

#setwd('/Users/brandonely/Desktop/bb_strain_model_dev/new_sim_output/set_larger_pops/plots/')
#ggsave('v5000_20_500_gen_data.jpeg', height = 8, width = 6)


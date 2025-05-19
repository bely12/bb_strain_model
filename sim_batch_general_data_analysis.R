setwd('/Users/brandonely/Desktop/bb_strain_model_dev/')
library(tidyverse)
library(ggplot2)
library(ggpubr)
source("sim_analysis_functions.R")

### upload data 
setwd('/Users/brandonely/Desktop/bb_strain_model_dev/sim_results/results/')

#header <- c('run_tag', 'year', 'infection_rate', 'diversity', 'active_strain_count', 'avg_strains_carried', 'unique_lineages')

cr_data <- read.delim('immune_5000_20_500_25SET_sim_data.tsv', header = T)
hs_data <- read.delim('hs_5000_20_500_25SET_sim_data.tsv', header = T)
ctrl_data <- read.delim('neutral_5000_20_500_25SET_sim_data.tsv', header = T)

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
    avg_antigen_distance = mean(avg_antigen_distance),
    avg_active_strain_count = mean(active_strain_count),
    avg_strains_carried = mean(avg_strains_carried),
    #avg_unique_lineages = mean(unique_lineages)
  )
x$selection <- 'immune selection'

y <- hs_data %>%
  select(-1) %>%
  group_by(year) %>%
  summarise(
    avg_infection_rate = mean(infection_rate),
    #avg_diversity = mean(diversity),
    avg_antigen_distance = mean(avg_antigen_distance),
    avg_active_strain_count = mean(active_strain_count),
    avg_strains_carried = mean(avg_strains_carried),
    #avg_unique_lineages = mean(unique_lineages)
  )
y$selection <- 'adaptive selection'

z <- ctrl_data %>%
  select(-1) %>%
  group_by(year) %>%
  summarise(
    avg_infection_rate = mean(infection_rate),
    #avg_diversity = mean(diversity),
    avg_antigen_distance = mean(avg_antigen_distance),
    avg_active_strain_count = mean(active_strain_count),
    avg_strains_carried = mean(avg_strains_carried),
    #avg_unique_lineages = mean(unique_lineages)
  )
z$selection <- 'neutral selection'

df <- bind_rows(x, y, z) # combine df's

df_long <- df %>%
  pivot_longer(cols = -c(year, selection), 
               names_to = "variable", 
               values_to = "mean_value")

# assign factor levels to put things in specific
df_long$variable <- factor(df_long$variable, levels = c(
  "avg_infection_rate",
  "avg_antigen_distance",
  "avg_active_strain_count",
  "avg_strains_carried"))

### create names for plots in panel
title.labs <- c('Mean infection rate',
                'Mean antigen distance',
                'Mean variants in population',
                'Mean strains carried by vector')
names(title.labs) <- c(
  "avg_infection_rate",
  "avg_antigen_distance",
  "avg_active_strain_count",
  "avg_strains_carried")

### plot 
ggplot(df_long, aes(x = year, y = mean_value, group = interaction(variable, selection), color = selection)) +
  geom_line() +
  scale_color_manual(values = c("immune selection" = "skyblue3", "adaptive selection" = "red3", "neutral selection" = "green4")) +
  #geom_point() +
  facet_wrap(~variable, scales = "free_y", ncol = 1, labeller = labeller(variable = title.labs)) +
  theme_bw() +
  labs(x = "Year", y = "Mean Value", title = "Means for simulation data") +
  theme(legend.position = "bottom")

#setwd('/Users/brandonely/Desktop/bb_strain_model_dev/sim_results/results/plots/')
#ggsave('v5000_20_500_gen_data.jpeg', height = 8, width = 8)


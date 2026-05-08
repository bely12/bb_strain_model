setwd('/Users/brandonely/Desktop/bb_strain_model_dev/R_analysis_scripts/')
library(tidyverse)
library(ggplot2)
library(ggpubr)
source("sim_analysis_functions.R")

setwd('/Users/brandonely/Desktop/bb_strain_model_dev/evo_sim/results_tables/')

# upload data
df <- read.table('general_stats_all_sims.tsv', header = T, sep = '\t')

# filter and label groups
df2 <- df %>% 
  filter(rec_rate == 0.01 & !(host_pops == "sixty_forty" & host_pop_stability == "stable")) %>%
  mutate(category = paste(selection, host_pops, host_pop_stability, gen_fitness, sep = '_'))

# get yearly means for each group 
df_summary <- df2 %>%
  group_by(category, year) %>%
  summarise(
    avg_infection_rate = mean(infection_rate),
    avg_antigen_distance = mean(avg_antigen_distance),
    avg_strains_carried = mean(avg_strains_carried),
    avg_spec_weight = mean(spec_weight)
  )

# convert to long format for plotting
df_long <- df_summary %>%
  pivot_longer(cols = -c(year, category), 
               names_to = "variable", 
               values_to = "mean_value")

# assign factor levels to put things in specific order
df_long$variable <- factor(df_long$variable, levels = c(
  "avg_infection_rate",
  "avg_strains_carried",
  "avg_antigen_distance",
  "avg_spec_weight"))

# create names for plots in panel
title.labs <- c('Mean infection rate',
                'Mean strains carried by vector',
                'Mean antigen distance',
                'Mean specialization weight')
names(title.labs) <- c(
  "avg_infection_rate",
  "avg_strains_carried",
  "avg_antigen_distance",
  "avg_spec_weight")

# plot 
df_long %>%
  ggplot(aes(x = year, y = mean_value, group = interaction(variable, category), color = category)) +
  geom_line() +
  facet_wrap(~variable, scales = "free_y", ncol = 1, labeller = labeller(variable = title.labs)) +
  #facet_wrap(~interaction(variable, rec_rate), scales = "free_y", ncol = 1, labeller = labeller(variable = title.labs)) +
  #facet_grid(variable ~ rec_rate, scales = "free_y", labeller = label_both) +
  theme_bw() +
  labs(x = "Year", y = "Mean Value", title = "Means for simulation data") +
  theme(legend.position = "bottom")

#ggsave('run_data.jpeg', height = 8, width = 8)

df2 %>%
  filter(year == 500) %>%
  pivot_longer(cols = 2:6, names_to = 'variable', values_to = 'value') %>%
  ggplot(aes(x = category, y = value, fill = selection)) +
  geom_boxplot() +
  facet_wrap(~variable, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8), legend.position = "bottom")

#ggsave('RunData_Year500_boxplots.jpeg', height = 8, width = 12)

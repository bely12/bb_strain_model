setwd('/Users/brandonely/Desktop/bb_strain_model_dev/')
library(tidyverse)
library(ggplot2)
library(ggpubr)
source("sim_analysis_functions.R")

# upload data 
setwd('/Users/brandonely/Desktop/bb_strain_model_dev/new_sim_output/')
cr_data <- read.delim('cr_len20_vec1000_yrs250_sim_data.tsv', header = T)
hs_data <- read.delim('hs_len20_vec1000_yrs250_sim_data.tsv', header = T)

# filter to remove failed sims
runs_to_remove <- hs_data %>% filter(active_strain_count == 0)
hs_data <- hs_data %>% filter(!run_tag %in% c(unique(runs_to_remove$run_tag)))

# get means by year for condition 1
x <- cr_data %>%
  select(-1) %>%
  group_by(year) %>%
  summarise(
    avg_infection_rate = mean(infection_rate),
    avg_diversity = mean(diversity),
    avg_antigen_distance = mean(avg_antigen_distance),
    avg_active_strain_count = mean(active_strain_count),
    avg_avg_strains_carried = mean(avg_strains_carried),
    avg_unique_lineages = mean(unique_lineages)
  )
x$selection <- 'cross reactivity'

# get means by year for condition 2
y <- hs_data %>%
  select(-1) %>%
  group_by(year) %>%
  summarise(
    avg_infection_rate = mean(infection_rate),
    avg_diversity = mean(diversity),
    avg_antigen_distance = mean(avg_antigen_distance),
    avg_active_strain_count = mean(active_strain_count),
    avg_avg_strains_carried = mean(avg_strains_carried),
    avg_unique_lineages = mean(unique_lineages)
  )
y$selection <- 'host specialization'

# combine df's
df <- bind_rows(x, y)
df_long <- df %>%
  pivot_longer(cols = -c(year, selection), 
               names_to = "variable", 
               values_to = "mean_value")

title.labs <- c('Active strains in population', 
                'Mean pairwise antigen distance',
                'Number of strains carried per tick',
                'Simpson diversity index',
                'Infection rate',
                'Surviving lineages in population')

names(title.labs) <- c('avg_active_strain_count',
                       'avg_antigen_distance',
                       'avg_avg_strains_carried',
                       'avg_diversity',
                       'avg_infection_rate',
                       'avg_unique_lineages')
# plot 
ggplot(df_long, aes(x = year, y = mean_value, group = interaction(variable, selection), color = selection, linetype = selection)) +
  geom_line() +
  scale_color_manual(values = c("cross reactivity" = "skyblue3", "host specialization" = "red3")) +
  #geom_point() +
  facet_wrap(~variable, scales = "free_y", ncol = 1, labeller = labeller(variable = title.labs)) +
  theme_bw() +
  labs(x = "Year", y = "Mean Value", title = "Means for simulation data", subtitle = 'N = 100 simulations per condition') +
  theme(legend.position = "bottom")



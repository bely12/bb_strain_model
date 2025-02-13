setwd('/Users/brandonely/Desktop/bb_strain_model_dev/')
library(tidyverse)
library(ggplot2)
library(ggpubr)
source("sim_analysis_functions.R")

# upload data 
setwd('/Users/brandonely/Desktop/bb_strain_model_dev/new_sim_output/test_set_500yrs/')
cr_data <- read.delim('cr_len20_vec1000_yrs500_sim_data.tsv', header = T)
hs_data <- read.delim('hs_len20_vec1000_yrs500_sim_data.tsv', header = T)
ctrl_data <- read.delim('ctrl_len20_vec1000_yrs500_sim_data.tsv', header = F, col.names = colnames(cr_data))

# filter to remove failed sims
cr_runs_to_remove <- cr_data %>% filter(active_strain_count == 0)
cr_data <- cr_data %>% filter(!run_tag %in% c(unique(cr_runs_to_remove$run_tag)))
hs_runs_to_remove <- hs_data %>% filter(active_strain_count == 0)
hs_data <- hs_data %>% filter(!run_tag %in% c(unique(hs_runs_to_remove$run_tag)))
ctrl_runs_to_remove <- ctrl_data %>% filter(active_strain_count == 0)
ctrl_data <- ctrl_data %>% filter(!run_tag %in% c(unique(ctrl_runs_to_remove$run_tag)))


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

# get means by year for condition 3
z <- ctrl_data %>%
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
z$selection <- 'control'


# combine df's
df <- bind_rows(x, y, z)
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
ggplot(df_long, aes(x = year, y = mean_value, group = interaction(variable, selection), color = selection)) +
  geom_line() +
  scale_color_manual(values = c("cross reactivity" = "skyblue3", "host specialization" = "red3", "control" = "yellow2")) +
  #geom_point() +
  facet_wrap(~variable, scales = "free_y", ncol = 1, labeller = labeller(variable = title.labs)) +
  theme_bw() +
  labs(x = "Year", y = "Mean Value", title = "Means for simulation data", subtitle = 'N = 100 simulations per condition') +
  theme(legend.position = "bottom")

setwd('/Users/brandonely/Desktop/bb_strain_model_dev/new_sim_output/test_set_500yrs/plots/')
#ggsave('test_set_500_gen_sim_data.jpeg', height = 8, width = 6)

# select just a few variables of interest instead

title.labs2 <- c('Active strains in population', 
                'Mean pairwise antigen distance',
                'Number of strains carried per tick',
                'Infection rate')

names(title.labs2) <- c('avg_active_strain_count',
                       'avg_antigen_distance',
                       'avg_avg_strains_carried',
                       'avg_infection_rate')

ggplot(df_long %>% filter(variable == 'avg_active_strain_count' | variable == 'avg_antigen_distance' | variable == 'avg_avg_strains_carried' | variable == 'avg_infection_rate'), 
       aes(x = year, y = mean_value, group = interaction(variable, selection), color = selection)) +
  geom_line() +
  scale_color_manual(values = c("cross reactivity" = "skyblue3", "host specialization" = "red3", "control" = "yellow2")) +
  #geom_point() +
  facet_wrap(~variable, scales = "free_y", ncol = 1, labeller = labeller(variable = title.labs2)) +
  theme_bw() +
  labs(x = "Year", y = "Mean Value", title = "Mean values for n = 100 simulations per condition") +
  theme(legend.position = "top",
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        strip.text = element_text(size = 15),
        axis.text = element_text(size = 15), 
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 20, hjust = 0.5))

setwd('/Users/brandonely/Desktop/bb_strain_model_dev/new_sim_output/test_set_500yrs/plots/')
ggsave('test_set_500_gen_sim_data_v2.jpeg', height = 8, width = 8)



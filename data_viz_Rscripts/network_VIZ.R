# network analysis

library(tidyverse)
library(ggpubr)
library(igraph)
library(broom)
setwd('/Users/brandonely/Desktop/bb_strain_model_dev/R_analysis_scripts/')
source("sim_analysis_functions.R")

##### main viz #####
# load data
setwd('/Users/brandonely/Desktop/bb_strain_model_dev/evo_sim/results_tables/network_analysis_results/')
df <- read.delim('network_analysis_res_ALL.tsv', header = T, sep = '\t')
set_compare <- c("neutral_high", "immune_high", "hybrid_exp_mod")
df <- df %>% filter(category %in% set_compare)

sim_counts <- df %>%
  group_by(category) %>%
  summarise(
    n_runs = n()/10
  )

# summary
df_summary <- df %>%
  group_by(category, threshold) %>%
  summarise(edge_density = mean(edge_dens),
            edge_density_sd = sd(edge_dens),
            mean_degrees = mean(mean_degree),
            mean_degrees_sd = sd(mean_degree),
            mod = mean(modul),
            mod_sd = sd(modul),
            lcc = mean(lrg_comp),
            lcc_sd = sd(lrg_comp))

# wide to long
df_long <- pivot_longer(df, cols = c(3:6), names_to = 'variable', values_to = 'values')

# plot line of means
df_long %>% 
  filter(variable != 'mean_degree' & category %in% set_compare) %>%
  mutate(selection = recode(selection, "neutral" = "No selection"),
         variable = factor(variable, levels = c("edge_dens", "modul", "lrg_comp"))) %>%
  ggplot(aes(x = threshold, y = values, color = selection, group = selection)) +
  geom_line(stat = "summary", fun = mean, size = 1.2) +     
  geom_point(stat = "summary", fun = mean, size = 2) + 
  labs(title = 'Mean values for network metrics across thresholds',
       subtitle = 'Connection threshold range 0.05-0.5, step = 0.05') +
  facet_wrap(~ variable, 
             scales = "free_y", ncol = 1, 
             labeller = as_labeller(c("edge_dens" = "Edge Density", 
                                      "modul" = "Modularity Score",
                                      "lrg_comp" = "LCC Fraction"))) +
  theme_bw() +
  theme(plot.title = element_text(size = 25, face = 'bold'),
        plot.subtitle = element_text(size = 20),
        axis.title.y = element_text(size = 25, face = 'bold'),
        axis.text = element_text(size = 25, face = 'bold'),
        axis.title.x = element_text(size = 25, face = 'bold'),
        legend.title = element_text(size = 25, face = 'bold'),
        legend.text = element_text(size = 25, face = 'bold'), 
        strip.text = element_text(size = 25))

setwd('/Users/brandonely/Desktop/bb_strain_model_dev/evo_sim/results_tables/new_plots/')
#ggsave('network_analysis_line_means_V2_DIS.png', width = 10, height = 8)

##### growth of network metrics via mean slope across thresholds #####

### modularity ###
#re-format for workflow
mod_df <- df_long %>% filter(variable == "modul" & threshold < 0.4 & threshold > 0.15) %>% 
  mutate(mod_score = values) %>%
  select(selection, run_id, threshold, mod_score)

# slope of lcc line vs threshold
mod_slopes <- mod_df %>%
  group_by(selection, run_id) %>%
  summarise(
    slope = coef(lm(mod_score ~ threshold))[2],
    .groups = "drop")

# anova + tukey for slope differences
mod_anova_fit <- aov(slope ~ selection, data = mod_slopes)
mod_anova_res <- tidy(mod_anova_fit)
#write.table(anova_res, 'anova_res_mod_slope.tsv', quote = F, row.names = F, sep = '\t')

mod_tukey_results <- tidy(TukeyHSD(mod_anova_fit))
#write.table(tukey_results, 'tukey_res_mod_slope.tsv', quote = F, row.names = F, sep = '\t')

# average slope with se 
mod_slopes_summary <- mod_slopes %>%
  group_by(selection) %>%
  summarise(
    mean_slope = -1*mean(slope, na.rm = TRUE),
    se_slope = sd(slope, na.rm = TRUE) / sqrt(n()),
    n_runs = n(),
    .groups = "drop")


### edge density ###
edge_df <- df_long %>% filter(variable == "edge_dens" & threshold < 0.4 & threshold > 0.15) %>% 
  mutate(edge_score = values) %>%
  select(selection, run_id, threshold, edge_score)

# slope of lcc line vs threshold
edge_slopes <- edge_df %>%
  group_by(selection, run_id) %>%
  summarise(
    slope = coef(lm(edge_score ~ threshold))[2],
    .groups = "drop")

# anova + tukey for slope differences
edge_anova_fit <- aov(slope ~ selection, data = edge_slopes)
edge_anova_res <- tidy(edge_anova_fit)
#write.table(anova_res, 'anova_res_edge_dens_slope.tsv', quote = F, row.names = F, sep = '\t')

edge_tukey_results <- tidy(TukeyHSD(edge_anova_fit))
#write.table(tukey_results, 'tukey_res_edge_dens_slope.tsv', quote = F, row.names = F, sep = '\t')

# average slope with se 
edge_slopes_summary <- edge_slopes %>%
  group_by(selection) %>%
  summarise(
    mean_slope = mean(slope, na.rm = TRUE),
    se_slope = sd(slope, na.rm = TRUE) / sqrt(n()),
    n_runs = n(),
    .groups = "drop")

### lcc ###
lcc_df <- df_long %>% filter(variable == "lrg_comp" & threshold < 0.4 & threshold > 0.15) %>% 
  mutate(lcc_fraction = values) %>%
  select(selection, run_id, threshold, lcc_fraction)

# slope of lcc line vs threshold
lcc_slopes <- lcc_df %>%
  group_by(selection, run_id) %>%
  summarise(
    slope = coef(lm(lcc_fraction ~ threshold))[2],
    .groups = "drop")

# anova + tukey for slope differences
lcc_anova_fit <- aov(slope ~ selection, data = lcc_slopes)
lcc_anova_res <- tidy(lcc_anova_fit)
#write.table(anova_res, 'anova_res_lcc_slope.tsv', quote = F, row.names = F, sep = '\t')

lcc_tukey_results <- tidy(TukeyHSD(lcc_anova_fit))
#write.table(tukey_results, 'tukey_res_lcc_slope.tsv', quote = F, row.names = F, sep = '\t')

# average slope with se 
lcc_slopes_summary <- lcc_slopes %>%
  group_by(selection) %>%
  summarise(
    mean_slope = mean(slope, na.rm = TRUE),
    se_slope = sd(slope, na.rm = TRUE) / sqrt(n()),
    n_runs = n(),
    .groups = "drop")

### plot 
mod_slopes_summary$metric <- 'Modularity Score'
edge_slopes_summary$metric <- 'Edge Density'
lcc_slopes_summary$metric <- 'LCC Fraction'
slopes <- bind_rows(mod_slopes_summary, edge_slopes_summary, lcc_slopes_summary)


slopes %>%
  mutate(selection = recode(selection, 'neutral' = 'No selection'),
         metric = factor(metric, levels = c("Edge Density", "Modularity Score", "LCC Fraction"))) %>%
  ggplot(aes(x = selection, y = mean_slope, fill = selection)) +
  geom_col(width = 0.5) +
  geom_errorbar(aes(ymin = mean_slope - se_slope, ymax = mean_slope + se_slope), width = 0.2, color = 'black') +
  labs(title = 'Slope growth for mean values\nacross intermediate thresholds', 
       subtitle = '0.15 < Connection Threshold < 0.4') +
  ylab(label = 'Mean slope') +
  facet_wrap(~ metric, ncol = 1, scales = 'free_y') +
  theme_bw() +
  theme(plot.title = element_text(size = 25, face = 'bold'),
        plot.subtitle = element_text(size = 20),
        axis.title.y = element_text(size = 25, face = 'bold'),
        axis.text = element_text(size = 25, face = 'bold'),
        axis.title.x = element_text(size = 25, face = 'bold'),
        legend.position = 'none',
        #legend.title = element_text(size = 25, face = 'bold'),
        #legend.text = element_text(size = 25, face = 'bold'), 
        strip.text = element_text(size = 25))

#ggsave('all_slopes_DIS.png', width = 7, height = 10)

##### individual network graphs #####
setwd('/Users/brandonely/Desktop/bb_strain_model_dev/evo_sim/results_tables/network_analysis_results/')
# load network data analysis for help choosing a run
df <- read.delim('network_analysis_res_ALL.tsv', header = T, sep = '\t')
df <- df %>% filter(selection == 'immune', gen_fitness == 'high', threshold == 0.25)

# upload raw data to build new network
setwd('/Users/brandonely/Desktop/bb_strain_model_dev/evo_sim/results_tables/sim_output/no_selection_results/')
df2 <- read.delim('neutral_variant_freqs_sampled.tsv', header = T, colClasses = c('character','numeric','numeric','character', 'character', 'character','numeric','numeric','character'))

# get rid of ghost variants 
df2 <- df2 %>% filter(freq >= 1)

# choose specific sim run to build network for
df2 <- df2 %>% filter(run_id == 'run_4')

# build the network 
network <- antigen_network(df2, threshold = 0.25)

# visualize the network
setwd('/Users/brandonely/Desktop/bb_strain_model_dev/evo_sim/results_tables/new_plots/')
plot_single <- df2 %>% filter(run_id == 'run_4')
png("no_selection_network_run4_thresh25.png", width = 800, height = 600)
antigen_network_viz(plot_single, threshold = 5)
dev.off()

# look at network stats for a single run and range of thresholds
# catcher = list()
# index = 1
# for (i in 1:10) {
#   catcher[[index]] <- antigen_network(df, threshold = i)
#   index = index + 1
#   network <- bind_rows(catcher)
# }

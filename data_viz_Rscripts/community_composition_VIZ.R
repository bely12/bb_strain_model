### data viz for coexistence and community compositions 
library(tidyverse)

##### community composition  #####
setwd('/Users/brandonely/Desktop/bb_strain_model_dev/evo_sim/results_tables/host_association_results/')

hyb <- read.delim('hybrid_community_composition_probabilities_p01.tsv', header = T, sep = '\t')
adp <- read.delim('adaptive_community_composition_probabilities_p01.tsv', header = T, sep = '\t')
df <- bind_rows(hyb, adp)
df$category <- sub("_stable$", "", df$category)
df[c('selection', 'model')] <- str_split_fixed(df$category, '_', 2)


plot_labels <- c(
  "g" = 'Generalists only',
  "gs" = 'Generalists and 1 specialist',
  "gss" = 'Generalists and 2 specialists',
  's' = '1 specialist',
  'ss' = '2 specialists'
)

plot_colors <- c(
  "g"   = "cyan2",
  "gs"  = "orange",     
  "gss" = "forestgreen", 
  "s"   = "dodgerblue",  
  "ss"  = "purple"       
)

# grouped bar plot
df %>%
  filter(model != 'linear') %>%
  ggplot(aes(x = category, y = probability, fill = community)) +
  geom_col(position = position_dodge()) +
  scale_fill_manual(values = plot_colors, labels = plot_labels) +
  theme_bw() +
  labs(x = "parameters", y = "probability", fill = "Community composition") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

setwd('/Users/brandonely/Desktop/bb_strain_model_dev/evo_sim/results_tables/new_plots/')
#ggsave('adaptive_coexistence_grouped_bar.png', width = 8, height = 6)

# stacked bar plot 

df %>%
  filter(model != 'linear') %>%
  mutate(model = factor(model, levels = c("low", "exp_mod", "high"))) %>%
  ggplot(aes(x = category, y = probability, fill = community)) +
  geom_col() +  
  scale_fill_manual(values = plot_colors, labels = plot_labels) +
  scale_x_discrete(labels = c(
    "adaptive_low" = "Adaptive",
    "hybrid_low" = "Hybrid",
    "adaptive_exp_mod" = "Adaptive",
    "hybrid_exp_mod" = "Hybrid",
    #'adaptive_linear' = 'Adaptvie',
    #'hybrid_linear' = 'Hybrid',
    'adaptive_high' = 'Adaptvie',
    'hybrid_high' = 'Hybrid')) +
  theme_bw() +
  labs(x = "Parameters", y = "Proportion", fill = "Community composition") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~model, 
             scales = "free_x", 
             drop = TRUE, 
             labeller = as_labeller(
               c("exp_mod" = "Moderate penalty", 
                 #"linear" = "Linear",
                 'high' = 'Low penalty',
                 "low" = "High penalty"))) +
  theme(panel.spacing = unit(1.5, "lines"),
        axis.title.y = element_text(size = 25),
        axis.text = element_text(size = 25, face = 'bold'),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 25),
        legend.title = element_text(size = 25),
        strip.text = element_text(size = 25))


#ggsave('community_composition_probabilities_StackedBar_DIS.png', width = 15, height = 6)

##### adaptive trait vals #####

setwd('/Users/brandonely/Desktop/bb_strain_model_dev/evo_sim/results_tables/host_association_results/')

vals <- read.delim('sampled_adaptive_trait_vals_ALL.tsv', header = T, sep = '\t')
vals$category <- paste(vals$selection, vals$gen_fit)
unique(vals$category)

# plot
vals %>%
  filter(gen_fit != 'linear') %>%
  mutate(category = factor(category, levels = c("adaptive low", "adaptive exp_mod", "adaptive high", 'hybrid low', 'hybrid exp_mod', 'hybrid high'))) %>%
  ggplot(aes(x = hs_val, fill = selection)) +
  geom_histogram(bins = 6, color = 'black') +
  labs(x = 'Adaptive trait value') +
  theme_bw() +
  facet_wrap(~category, ncol = 3, 
             labeller = as_labeller(c("adaptive exp_mod" = "Adaptive Moderate Penalty", 
                                      #"adaptive linear" = "Adaptive Low Penalty",
                                      'adaptive high' = 'Adaptive Low Penalty',
                                      "adaptive low" = "Adaptive High Penalty",
                                      "hybrid exp_mod" = "Hybrid Moderate Penalty", 
                                      #"hybrid linear" = "Hybrid Low Penalty",
                                      'hybrid high' = 'Hybrid Low Penalty',
                                      "hybrid low" = "Hybrid High Penalty"))) +
  theme(legend.position = 'none',
        axis.title.y = element_text(size = 25, face = 'bold'),
        axis.text = element_text(size = 25, face = 'bold'),
        axis.title.x = element_text(size = 25, face = 'bold'),
        strip.text = element_text(size = 25))
  

setwd('/Users/brandonely/Desktop/bb_strain_model_dev/evo_sim/results_tables/new_plots/')
#ggsave('hs_trait_val_distribution_DIS.png', width = 14, height = 8)

##### community specialization weights #####

setwd('/Users/brandonely/Desktop/bb_strain_model_dev/evo_sim/results_tables/sim_output/hybrid_results/')
df1 <- read.table('hybrid_ExpMod_sim_run_stats.tsv', header = T, sep = '\t')
df2 <- read.table('hybrid_LowGenFit_sim_run_stats.tsv', header = T, sep = '\t')
df3 <- read.table('hybrid_HighGenFit_sim_run_stats.tsv', header = T, sep = '\t')

setwd('/Users/brandonely/Desktop/bb_strain_model_dev/evo_sim/results_tables/sim_output/adaptive_results/')
df4 <- read.table('adaptive_ExpMod_sim_run_stats.tsv', header = T, sep = '\t')
df5 <- read.table('adaptive_LowGenFit_sim_run_stats.tsv', header = T, sep = '\t')
df6 <- read.table('adaptive_HighGenFit_sim_run_stats.tsv', header = T, sep = '\t')

df <- bind_rows(df1, df2, df3, df4, df5, df6)
rm(df1, df2, df3, df4, df5, df6)
df <- df %>% filter(year == 500)

df_summary <- df %>% 
  group_by(selection, gen_fitness) %>%
  summarise(
    mean_spec_weight = mean(spec_weight)
  )


df %>% mutate(category = paste(selection, gen_fitness, sep = '_')) %>%
  mutate(model = factor(gen_fitness, levels = c("low", "exp_mod", "high"))) %>%
  ggplot(aes(x = category, y = spec_weight, fill = selection)) +
  geom_boxplot(outliers = F) +
  geom_jitter(alpha = 0.2) +
  geom_hline(yintercept = 0.25, colour = 'red', linetype = 2, alpha = 0.3) +
  # annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0.08, alpha = 0.3, fill = "#3399ff") +
  # annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.08, ymax = 0.16, alpha = 0.25, fill = "#99ccff") +
  # annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.16, ymax = 0.25, alpha = 0.2, fill = "#cce5ff") +
  # annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.25, ymax = 0.34, alpha = 0.2, fill = "#ffcccc") +
  # annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.34, ymax = 0.42, alpha = 0.25, fill = "#ff9999") +
  # annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.42, ymax = 0.5, alpha = 0.3, fill = "#ff6666") +
  scale_x_discrete(labels = c( "adaptive_low" = "Adaptive", "hybrid_low" = "Hybrid", "adaptive_exp_mod" = "Adaptive", "hybrid_exp_mod" = "Hybrid", 'adaptive_high' = 'Adaptvie', 'hybrid_high' = 'Hybrid')) +
  theme_bw() +
  labs(x = "", y = "Community specialization weight", fill = "Selection") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~model, 
             scales = "free_x", 
             drop = TRUE, 
             labeller = as_labeller(
               c("exp_mod" = "Moderate penalty", 
                 "high" = "Low penalty", 
                 "low" = "High penalty"))) +
  theme(axis.title.y = element_text(size = 20, face = 'bold'),
        axis.text = element_text(size = 25, face = 'bold'),
        #axis.title.x = element_text(size = 25, face = 'bold'),
        legend.text = element_text(size = 25),
        legend.title = element_text(size = 25),
        strip.text = element_text(size = 25))

setwd('/Users/brandonely/Desktop/bb_strain_model_dev/evo_sim/results_tables/new_plots/')
ggsave('community_spec_weights_DIS.png', width = 12, height = 6)



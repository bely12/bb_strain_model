library(tidyverse)
library(ggExtra)
library(broom)

setwd('/Users/brandonely/Desktop/bb_strain_model_dev/evo_sim/results_tables/cluster_analysis_results/')

# upload data
clust_res <- read.delim('hybrid_cluster_results.tsv', sep = '\t', header = T)
clust_res$category <- paste(clust_res$selection, clust_res$gen_fit)
clust_res <- clust_res %>% filter(gen_fit == 'exp_mod')

clust_res2 <- read.delim('immune_cluster_results.tsv', sep = '\t', header = T)
clust_res2$category <- paste(clust_res2$selection, clust_res2$gen_fit)

clust_res3 <- read.delim('adaptive_cluster_results.tsv', sep = '\t', header = T)
clust_res3$category <- paste(clust_res3$selection, clust_res3$gen_fit)

clust_res4 <- read.delim('no_selection_cluster_results.tsv', sep = '\t', header = T)
clust_res4$category <- paste(clust_res4$selection, clust_res4$gen_fit)


df <- bind_rows(clust_res, clust_res2, clust_res4)
df2 <- bind_rows(clust_res, clust_res2, clust_res4)

# convert to long format
df <- df %>% 
  #filter(selection != 'hybrid' & rec_rate == 0.01) %>%
  pivot_longer(
    cols = n_variants:ShanE,
    names_to = 'variable',
    values_to = 'value'
  )

# counts for each category
x <- df %>% 
  mutate(category = paste(selection,host_pops,host_pop_stability,gen_fit, sep = '_')) %>%
  filter(rec_rate == 0.01)

x <- df %>%
  group_by(category) %>%
  summarise(unique_runs = n_distinct(run_id))
unique(df$variable)

# ANOVA
anova_df <- df2 %>%
  aov(avg_out_dist ~ selection, data = .) %>%
  tidy()
#write.table(anova_df, 'anova_res_avgOutDist.tsv', quote = F, row.names = F, sep = '\t')

# Tukey HSD
tukey_df <- df2 %>%
  aov(ShanE ~ selection, data = .) %>%
  TukeyHSD() %>%
  tidy()
#write.table(tukey_df, 'tukey_res_ShanE.tsv', quote = F, row.names = F, sep = '\t')

df %>% 
  filter(variable == 'sil_score' | variable == 'ShanE') %>%
  mutate(selection = recode(selection, "neutral" = "No selection"),
         variable = factor(variable, levels = c('sil_score', 'ShanE'))) %>%
  #mutate(category = paste(selection,rec_rate,sep = '_')) %>%
  ggplot(aes(x = selection, y = value, fill = selection)) +
  geom_boxplot(outliers = F) +
  geom_jitter(alpha = 0.3) +
  # scale_x_discrete(labels = c(
  #   "immune high" = "Immune",
  #   "neutral high" = "None",
  #   "hybrid exp_mod" = "Hybrid")) +
  # xlab(label = 'Selection') +
  facet_wrap(~variable, 
             scales = 'free_y',
             labeller = as_labeller(c("sil_score" = "Silhouette score", 
                                      "ShanE" = "Shannon Evenness")))+
  theme_bw() +
  theme(plot.title = element_text(size = 25),
        
        axis.title.y = element_text(size = 25, face = 'bold'),
        axis.text.y = element_text(size = 25, face = 'bold'),
        
        axis.title.x = element_text(size = 25, face = 'bold'),
        axis.text.x = element_text(size = 25, angle = 45, hjust = 1, face = 'bold'), 
        
        panel.spacing = unit(4, "lines"),
        legend.title = element_text(size = 25, face = 'bold'),
        legend.text = element_text(size = 25, face = 'bold'),
        strip.text = element_text(size = 25, face = 'bold'))

setwd('/Users/brandonely/Desktop/bb_strain_model_dev/evo_sim/results_tables/new_plots/')
#ggsave('cluster_results_boxplots_V2_DIS.png', width = 12, height = 8)


# dot plot in/out dists with margin density histogram
p <- df2 %>%
  select(selection, avg_in_dist, avg_out_dist) %>%
  mutate(selection = recode(selection, "neutral" = "No selection")) %>%
  ggplot(aes(x = avg_in_dist, y = avg_out_dist, color = selection)) +
  geom_point(alpha = 0.5) +
  stat_ellipse(level = 0.95, size = 1) +
  labs(x = "Mean in cluster distance", y = "Mean across cluster distance") +
  theme_bw() +
  theme(#plot.title = element_text(size = 25),
        legend.position = 'none',
        axis.title.y = element_text(size = 25, face = 'bold'),
        axis.text.y = element_text(size = 25, face = 'bold'),
        
        axis.title.x = element_text(size = 25, face = 'bold'),
        axis.text.x = element_text(size = 25, face = 'bold'))
  # theme(legend.text = element_text(size = 10),
  #       legend.position = 'bottom',
  #       legend.title = element_text(size = 0),
  #       axis.title = element_text(size = 15),
  #       axis.text = element_text(size = 15))

#ggsave('InOut_ClusterDists.png', width = 8, height = 8)

ggMarginal(p, type = "density", groupColour = F, groupFill = TRUE, alpha = 0.5)
p2 <- ggMarginal(p, type = "density", groupColour = F, groupFill = TRUE, alpha = 0.5)
ggsave('InOut_ClusterDists_v3_DIS.jpeg', plot = p2, height = 10, width = 12)

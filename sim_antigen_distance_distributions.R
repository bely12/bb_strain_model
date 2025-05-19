# sim antigen distance distribution analysis 
setwd('/Users/brandonely/Desktop/bb_strain_model_dev/')
library(tidyverse)
library(ggplot2)
library(ggpubr)
source("sim_analysis_functions.R")

### sampled distances from final population ###
### upload data
setwd('/Users/brandonely/Desktop/bb_strain_model_dev/sim_results/v5000_yrs500_results/')

cr_dists <- read.delim('immune_5000_20_500_25SET_sampled_pairwise_antigen_dists.tsv', header = F, col.names = c('run_id', 'distance')) 
cr_dists$selection <- 'immune selection' # label 
cr_dists <- sample_n(cr_dists, 1000)

hs_dists <- read.delim('hs_5000_20_500_25SET_sampled_pairwise_antigen_dists.tsv', header = F, col.names = c('run_id', 'distance')) 
hs_dists$selection <- 'adaptive selection' # label 
hs_dists <- sample_n(hs_dists, 1000)

ctrl_dists <- read.delim('neutral_5000_20_500_25SET_sampled_pairwise_antigen_dists.tsv', header = F, col.names = c('run_id', 'distance'))
ctrl_dists$selection <- 'neutral selection' # label
ctrl_dists <- sample_n(ctrl_dists, 1000)

all_dists <- bind_rows(cr_dists, hs_dists, ctrl_dists)

### plot historgrams

ggplot(all_dists, aes(x = distance, fill = selection, alpha = 0.5)) +
  geom_density(binwidth = 1, color = 'black') +
  #geom_histogram(binwidth = 1, color = 'black') +
  scale_fill_manual(values = c("immune selection" = "skyblue3", "adaptive selection" = "red3", "neutral selection" = "yellow2")) +
  facet_wrap(~ selection, nrow = 1) +
  theme_bw() +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title = "Random sampling of pairwise antigen distances", subtitle = "1000 sampled distances from n = 100 simulations", x = "Distance", y = "Frequency") +
  theme(legend.position = "none", panel.grid = element_blank())

#setwd('/Users/brandonely/Desktop/bb_strain_model_dev/sim_results/results/plots/')
#ggsave('pw_Adist_distribution_historgram.jpeg', height = 4, width = 4)


# mean pairwise antigen distance comparison 
summary_df <- all_dists %>%
  group_by(selection, run_id) %>%
  summarise(mean_distance = mean(distance, na.rm = TRUE), .groups = "drop")

ggplot(summary_df, aes(x = selection, y = mean_distance)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 1.0)) +
  geom_jitter(aes(colour = selection)) +
  scale_color_manual(values = c("immune selection" = "skyblue3", "adaptive selection" = "red3", "neutral selection" = "yellow2")) +
  labs(y = "Average of mean pairwise distances", x = '') +
  stat_compare_means(method = "t.test", label = "p.signif",
                     comparisons = list(
                       c("immune selection", "adaptive selection"),
                       c("immune selection", "neutral selection"),
                       c("adaptive selection", "neutral selection"))) +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text = element_text(size = 10), 
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 20, hjust = 0.5))

#ggsave('v5000_25set_final_sim_mean_sampled_AntDist_boxplot.jpeg', height = 5, width = 4)


# housekeeping! clears everything except functions
#rm(list = setdiff(ls(), lsf.str()))

# ks tests
ks_test_crVShs <- ks.test(cr_dists$distance, hs_dists$distance) #ks.test(cr_ks, hs_ks)
ks_test_crVSctrl <- ks.test(cr_dists$distance, ctrl_dists$distance) #ks.test(cr_ks, ctrl_ks)
ks_test_hsVSctrl <- ks.test(hs_dists$distance, ctrl_dists$distance) #ks.test(hs_ks, ctrl_ks)

# plot
ggplot(all_dists, aes(x = distance, group = selection, color = selection)) +
  scale_color_manual(values = c("immune selection" = "skyblue3", "adaptive selection" = "red3", "neutral selection" = 'yellow2')) +
  stat_ecdf(size=1) +
  labs(y = 'ECDF', x = 'Pairwise Antigen Distances') +
  theme_bw() +
  theme(legend.position ="bottom", axis.title = element_text(size = 20), axis.text = element_text(size = 20)) +
  ### PRINT KS TEST RESULT STATISTICS ONTO PLOT ###
  annotate('text', x = 0.05, y = 0.99, label = 'KS Results:', size = 4, color = 'black', fontface = 'bold') +
  # comparison 1
  annotate('text', x = 0.7, y = 0.94, 
           label = 'immune vs adaptive', size = 4, color = 'black', fontface = 'bold') +
  annotate('text', x = 0, y = 0.90, 
           label = 'stat = ', size = 4, color = 'black', fontface = 'bold') +
  annotate('text', x = 0.8, y = 0.90, 
           label = paste(round(ks_test_crVShs$statistic, 3)), size = 4, color = 'black', fontface = 'bold') +
  annotate('text', x = 0.25, y = 0.86, 
           label = 'p-value = ', size = 4, color = 'black', fontface = 'bold') +
  annotate('text', x = 1.67, y = 0.86, 
           label = format(ks_test_crVShs$p.value, scientific = TRUE, digits = 3), size = 4, color = 'black', fontface = 'bold') +
  #comparison 2
  annotate('text', x = 0.65, y = 0.8, 
           label = 'immune vs neutral', size = 4, color = 'black', fontface = 'bold') +
  annotate('text', x = 0, y = 0.76, 
           label = 'stat = ', size = 4, color = 'black', fontface = 'bold') +
  annotate('text', x = 0.8, y = 0.76, 
           label = paste(round(ks_test_crVSctrl$statistic, 3)), size = 4, color = 'black', fontface = 'bold') +
  annotate('text', x = 0.25, y = 0.72, 
           label = 'p-value = ', size = 4, color = 'black', fontface = 'bold') +
  annotate('text', x = 1.67, y = 0.72, 
           label = format(ks_test_crVSctrl$p.value, scientific = TRUE, digits = 3), size = 4, color = 'black', fontface = 'bold') +
  #comparison 3
  annotate('text', x = 0.65, y = 0.66, 
           label = 'adaptive vs neutral', size = 4, color = 'black', fontface = 'bold') +
  annotate('text', x = 0.01, y = 0.62, 
           label = 'stat = ', size = 4, color = 'black', fontface = 'bold') +
  annotate('text', x = 0.8, y = 0.62, 
           label = paste(round(ks_test_hsVSctrl$statistic, 3)), size = 4, color = 'black', fontface = 'bold') +
  annotate('text', x = 0.25, y = 0.58, 
           label = 'p-value = ', size = 4, color = 'black', fontface = 'bold') +
  annotate('text', x = 1.67, y = 0.58, 
           label = format(ks_test_hsVSctrl$p.value, scientific = TRUE, digits = 3), size = 4, color = 'black', fontface = 'bold')

#ggsave('ks_test_results.jpeg', height = 8, width = 12)


# qq plot
# cr_sampled <- sort(cr_sampled_all$V2)
# hs_sampled <- sort(hs_sampled_all$V2)
# df <- data.frame(cr = cr_sampled, hs = hs_sampled)
# ggplot(df, aes(x = cr, y = hs)) +
#   geom_point() +
#   geom_abline(aes(slope = 1, intercept = 0), linetype = 2, color = 'red') +
#   labs(title = 'Q-Q for antigen distance distributions', 
#        subtitle = 'immune cross reactivity VS host specialization',
#        x = 'immune cross reactivity',
#        y = 'host specialization') +
#   theme_bw()

#ggsave('test_set_500_allyrs_crVShs_qqplot.jpeg', height = 5, width = 8)

# base R for plotting this stuff 
# plot empirical distributions in base R
# cr_ecdf <- ecdf(cr_sampled$V2)
# hs_ecdf <- ecdf(hs_sampled$V2)
# plot(cr_ecdf, verticals=TRUE, do.points=FALSE, col="blue") 
# plot(hs_ecdf, verticals=TRUE, do.points=FALSE, col="green", add=TRUE) 

# qq plot in base R
# qqplot(cr_sampled, hs_sampled, 
#        main = "antigen distance distributions", 
#        xlab = "cross reactivity", 
#        ylab = "host specialization")
# abline(c(0,1), col = "red")

  

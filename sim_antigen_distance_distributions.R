# sim antigen distance distribution analysis 
setwd('/Users/brandonely/Desktop/bb_strain_model_dev/')
library(tidyverse)
library(ggplot2)
library(ggpubr)
source("sim_analysis_functions.R")

### sampled distances from final population ###
### upload data
setwd('/Users/brandonely/Desktop/bb_strain_model_dev/new_sim_output/set_larger_pops/')

cr_dists <- read.delim('immune_20_5000_500_sampled_pairwise_antigen_dists.tsv', header = F) 
cr_dists$selection <- 'immune selection' # label 

hs_dists <- read.delim('hs_20_5000_500_sampled_pairwise_antigen_dists.tsv', header = F) 
hs_dists$selection <- 'adaptive selection' # label 

ctrl_dists <- read.delim('ctrl_len20_vec1000_yrs500_sampled_pairwise_antigen_dists.tsv', header = F)
ctrl_dists$selection <- 'neutral selection' # label

all_dists <- bind_rows(cr_dists, hs_dists, ctrl_dists)

### plot historgrams
ggplot(all_dists, aes(x = V2, fill = selection)) +
  geom_histogram(binwidth = 1, color = 'black') +
  scale_fill_manual(values = c("immune selection" = "skyblue3", "adaptive selection" = "red3", "neutral selection" = "yellow2")) +
  facet_wrap(~ selection, nrow = 1) +
  theme_bw() +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title = "Random sampling of pairwise antigen distances", subtitle = "5000 sampled distances from n = 100 simulations", x = "Distance", y = "Frequency") +
  theme(legend.position = "none", panel.grid = element_blank())

#setwd('/Users/brandonely/Desktop/bb_strain_model_dev/new_sim_output/test_set_500yrs/plots/')
#ggsave('test_set_500_full_sim_AntDist.jpeg', height = 5, width = 8)


### Complete distances from FINAL POPULATION ###

### upload data
setwd('/Users/brandonely/Desktop/bb_strain_model_dev/new_sim_output/set_larger_pops/')
header_names <- c('variant', 'counts', 'frequency', 'run_id') # set col names

cr_final_pop_data <- read.delim('immune_5000_20_500_variant_frequencies.tsv', header = T, col.names = header_names ,colClasses = c('character','numeric','numeric','character'))
cr_final_pop_data_filtered <- cr_final_pop_data %>% filter(frequency >1.49)

hs_final_pop_data <- read.delim('hs_20_5000_500_variant_frequencies.tsv', header = T, col.names = header_names, colClasses = c('character','numeric','numeric','character'))
hs_final_pop_data_filtered <- hs_final_pop_data %>% filter(frequency >1.49)

ctrl_final_pop_data <- read.delim('neutral_5000_20_500_variant_frequencies.tsv', header = T, col.names = header_names, colClasses = c('character','numeric','numeric','character'))
ctrl_final_pop_data_filtered <- ctrl_final_pop_data %>% filter(frequency >1.49)

# check some stuff to see if filtering is reasonable
# x <- hs_final_pop_data %>%
#   filter(frequency >1.49) %>%
#   group_by(run_id) %>%
#   summarise(
#     n_variants = n()
#   )
# print(mean(x$n_variants))

# final distributions for cross reactivity
sample_catcher <- list()
index <- 1
for (run in unique(cr_final_pop_data_filtered$run_id)) {
  filtered_df <- cr_final_pop_data_filtered %>% filter(run_id == run)
  dists <- final_antigen_distances(filtered_df)
  sampled_dists <- dists[sample(nrow(dists), 100), ]
  sample_catcher[[index]] <- sampled_dists
  index <- index + 1
  print(index)
}

all_dists <- unlist(sample_catcher)
cr <- sample(all_dists, size = 1000, replace = FALSE)
#cr <- sample(all_dists, size = length(all_dists), replace = F)
cr2 <- data.frame(values = cr)
cr2$selection <- 'Immune selection'

# final distributions for host specialization
sample_catcher <- list()
index <- 1
for (run in unique(hs_final_pop_data_filtered$run_id)) {
  filtered_df <- hs_final_pop_data_filtered %>% filter(run_id == run)
  dists <- final_antigen_distances(filtered_df)
  sampled_dists <- dists[sample(nrow(dists), 100), ]
  sample_catcher[[index]] <- sampled_dists
  index <- index + 1
  print(index)
}

all_dists <- unlist(sample_catcher)
hs <- sample(all_dists, size = 1000, replace = FALSE)
#hs <- sample(all_dists, size = length(all_dists), replace = FALSE)
hs2 <- data.frame(values = hs)
hs2$selection <- 'Host specialization'

# final distributions for negative control
sample_catcher <- list()
index <- 1
for (run in unique(ctrl_final_pop_data_filtered$run_id)) {
  filtered_df <- ctrl_final_pop_data_filtered %>% filter(run_id == run)
  dists <- final_antigen_distances(filtered_df)
  sampled_dists <- dists[sample(nrow(dists), 100), ]
  sample_catcher[[index]] <- sampled_dists
  index <- index + 1
  print(index)
}

all_dists <- unlist(sample_catcher)
#ctrl <- sample(all_dists, size = 5000, replace = FALSE)
ctrl <- sample(all_dists, size = 1000, replace = FALSE)
ctrl2 <- data.frame(values = ctrl)
ctrl2$selection <- 'Neutral selection'

# combine distributions to one df
set_dists <- bind_rows(cr2, hs2, ctrl2)
colnames(set_dists) <- c('Antigen Distance', 'selection')

ggplot(set_dists, aes(x = `Antigen Distance`, fill = selection)) +
  geom_histogram(binwidth = 1, color = 'black', show.legend = F) +
  scale_fill_manual(values = c("Immune selection" = "skyblue3", "Host specialization" = "red3", "Neutral selection"= "yellow2")) +
  facet_wrap(~ selection, scales = 'free_y', nrow = 1) +
  #facet_grid(~ selection) +
  #labs(title = '5000 sampled distances, n = 100 simuations',x = 'Antigen Distance') +
  scale_y_continuous(expand = c(0, 0)) +
  #theme(legend.position = "top", panel.grid = element_blank()) +
  theme_bw() #+
  #theme(strip.text = element_text(size = 15))
  # theme(strip.text = element_blank(),
  #       plot.title = element_text(size = 20),
  #       plot.subtitle = element_text(size = 20),
  #       axis.text = element_text(size = 15),
  #       axis.title = element_text(size = 20),
  #       legend.text = element_text(size = 20),
  #       legend.position = "bottom", 
  #       legend.title = element_blank())

#setwd('/Users/brandonely/Desktop/bb_strain_model_dev/new_sim_output/test_set_500yrs/plots/')
#ggsave('v5000_set_final_sim_AntDist_filtered_freq_less1.5.jpeg', height = 5, width = 12)

# mean pairwise antigen distance comparison 
bar_df_cr <- data.frame(value = cr, selection = 'Immune selection')
bar_df_hs <- data.frame(value = hs, selection = 'Host specialization')
bar_df_ctrl <- data.frame(value = ctrl, selection = 'Neutral selection')
bar_df <- bind_rows(bar_df_cr,bar_df_hs,bar_df_ctrl)

ggplot(bar_df, aes(x = selection, y = value)) +
  geom_bar(stat = "summary", fun = mean, aes(fill = selection), width = 0.5) +
  scale_fill_manual(values = c("Immune selection" = "skyblue3", 
                               "Host specialization" = "red3", 
                               "Neutral selection" = "yellow2")) +
  scale_x_discrete(expand = expansion(mult = c(0.2, 0.2))) + #reduce space between bars
  labs(y = "Mean distance", x = '') +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  #coord_cartesian(ylim = c(0, 15)) +
  stat_compare_means(method = "t.test", label = "p.signif",
                     comparisons = list(
                       c("Immune selection", "Host specialization"),
                       c("Immune selection", "Neutral selection"),
                       c("Host specialization", "Neutral selection"))) +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text = element_text(size = 10), 
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 20, hjust = 0.5))

setwd('/Users/brandonely/Desktop/bb_strain_model_dev/new_sim_output/set_larger_pops/plots/')
ggsave('v5000_set_final_sim_mean_AntDist_bar_filtered_freq_less1.5.jpeg', height = 5, width = 4)


# housekeeping! clears everything except functions
#rm(list = setdiff(ls(), lsf.str()))

### STATISTICAL TESTING ON DISTRIBUTIONS ###
# run k-s test

#ks_test <- ks.test(cr_sampled_all$V2, hs_sampled_all$V2) # samples from any year in sim

# record ks test result to print on plot
# ks_result <- paste('k-s test', "\n", 
#                    "stat= ", round(ks_test$statistic, 3), "\n",
#                    #"pval= ", format("2.2e-16"))
#                    "pval= ", format(ks_test$p.value, scientific = TRUE, digits = 3))

# plot results for cr vs hs conditions
# ggplot(set_dists %>% filter(selection != 'control'), 
#        aes(x = values, group = selection, color = selection)) +
#   scale_color_manual(values = c("cross reactivity" = "skyblue3", "host specialization" = "red3")) +
#   stat_ecdf(size=1) +
#   theme_bw() +
#   theme(legend.position ="top") +
#   xlab("value") +
#   ylab("ECDF") +
#   annotate("text", x = 2, y = 0.8, label = ks_result, size = 4, color = "black")

# redoing this for all combinations of conditions

# sample 100 from final pops of each condition
cr_ks <- sample(cr2$values, size = 100, replace = FALSE)
hs_ks <- sample(hs2$values, size = 100, replace = FALSE)
ctrl_ks <- sample(ctrl2$values, size = 100, replace = FALSE)
# ks tests
ks_test_crVShs <- ks.test(cr_ks, hs_ks)
ks_test_crVSctrl <- ks.test(cr_ks, ctrl_ks)
ks_test_hsVSctrl <- ks.test(hs_ks, ctrl_ks)

# record ks test result to print on plot
ks_result_crVShs <- paste('CR vs HS','\n',
                          'p-value = ', format(ks_test_crVShs$p.value, scientific = TRUE, digits = 3),'\n', 
                          'stat = ', round(ks_test_crVShs$statistic, 3))
ks_test_crVShs$p.value
ks_result_crVSctrl <- paste('CR vs CTRL p-value = ', format(ks_test_crVSctrl$p.value, scientific = TRUE, digits = 3)) 
ks_result_hsVSctrl <- paste('HS vs CTRL p-value = ', format(ks_test_hsVSctrl$p.value, scientific = TRUE, digits = 3)) 

ggplot(set_dists, aes(x = `Antigen Distance`, group = selection, color = selection)) +
  scale_color_manual(values = c("Immune selection" = "skyblue3", "Host specialization" = "red3", "Neutral selection" = 'yellow2')) +
  stat_ecdf(size=1) +
  labs(y = 'ECDF', x = 'Pairwise Antigen Distances') +
  theme_bw() +
  theme(legend.position ="bottom", axis.title = element_text(size = 20), axis.text = element_text(size = 20)) +
  annotate('text', x = 0, y = 0.99, label = 'KS Results:', size = 4, color = 'black', fontface = 'bold') +
  # comparison 1
  annotate('text', x = -0.1, y = 0.94, 
           label = 'CR vs HS', size = 4, color = 'black', fontface = 'bold') +
  annotate('text', x = 0, y = 0.90, 
           label = 'stat = ', size = 4, color = 'black', fontface = 'bold') +
  annotate('text', x = 0.6, y = 0.90, 
           label = paste(round(ks_test_crVShs$statistic, 3)), size = 4, color = 'black', fontface = 'bold') +
  annotate('text', x = 0.25, y = 0.86, 
           label = 'p-value = ', size = 4, color = 'black', fontface = 'bold') +
  annotate('text', x = 1.37, y = 0.86, 
           label = format(ks_test_crVShs$p.value, scientific = TRUE, digits = 3), size = 4, color = 'black', fontface = 'bold') +
  #comparison 2
  annotate('text', x = 0.1, y = 0.8, 
           label = 'CR vs CTRL', size = 4, color = 'black', fontface = 'bold') +
  annotate('text', x = 0, y = 0.76, 
           label = 'stat = ', size = 4, color = 'black', fontface = 'bold') +
  annotate('text', x = 0.6, y = 0.76, 
           label = paste(round(ks_test_crVSctrl$statistic, 3)), size = 4, color = 'black', fontface = 'bold') +
  annotate('text', x = 0.25, y = 0.72, 
           label = 'p-value = ', size = 4, color = 'black', fontface = 'bold') +
  annotate('text', x = 1.37, y = 0.72, 
           label = format(ks_test_crVSctrl$p.value, scientific = TRUE, digits = 3), size = 4, color = 'black', fontface = 'bold') +
  #comparison 3
  annotate('text', x = 0.1, y = 0.66, 
           label = 'HS vs CTRL', size = 4, color = 'black', fontface = 'bold') +
  annotate('text', x = 0, y = 0.62, 
           label = 'stat = ', size = 4, color = 'black', fontface = 'bold') +
  annotate('text', x = 0.6, y = 0.62, 
           label = paste(round(ks_test_hsVSctrl$statistic, 3)), size = 4, color = 'black', fontface = 'bold') +
  annotate('text', x = 0.25, y = 0.58, 
           label = 'p-value = ', size = 4, color = 'black', fontface = 'bold') +
  annotate('text', x = 1.37, y = 0.58, 
           label = format(ks_test_hsVSctrl$p.value, scientific = TRUE, digits = 3), size = 4, color = 'black', fontface = 'bold')


ggsave('v5000_set_ks_test_filtered_freq_less1.5.jpeg', height = 7, width = 12)


# qq plot
cr_sampled <- sort(cr_sampled_all$V2)
hs_sampled <- sort(hs_sampled_all$V2)
df <- data.frame(cr = cr_sampled, hs = hs_sampled)
ggplot(df, aes(x = cr, y = hs)) +
  geom_point() +
  geom_abline(aes(slope = 1, intercept = 0), linetype = 2, color = 'red') +
  labs(title = 'Q-Q for antigen distance distributions', 
       subtitle = 'immune cross reactivity VS host specialization',
       x = 'immune cross reactivity',
       y = 'host specialization') +
  theme_bw()

setwd('/Users/brandonely/Desktop/bb_strain_model_dev/new_sim_output/test_set_500yrs/plots/')
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

  

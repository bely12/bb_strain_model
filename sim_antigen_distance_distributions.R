# sim antigen distance distribution analysis 
setwd('/Users/brandonely/Desktop/bb_strain_model_dev/')
library(tidyverse)
library(ggplot2)
source("sim_analysis_functions.R")

setwd('/Users/brandonely/Desktop/bb_strain_model_dev/new_sim_output/test_set_500yrs/')


### ANTIGEN DISTANCE DISTRIBUTIONS FOR ALL YEARS OF SIM ###

# upload cross reactivity data and remove failed sims
cr_ad_all <- read.delim('cr_len20_vec1000_yrs500_sampled_pairwise_antigen_dists.tsv', header = F) 
cr_data <- read.delim('cr_len20_vec1000_yrs500_sim_data.tsv', header = T)
cr_runs_to_remove <- cr_data %>% filter(active_strain_count == 0) # find failed sims
cr_ad_all <- cr_ad_all %>% filter(!V1 %in% c(unique(cr_runs_to_remove$run_tag))) # remove failed sims
cr_sampled_all <- cr_ad_all[sample(nrow(cr_ad_all), 5000), ] # sample
cr_sampled_all$selection <- 'cross reactivity' # label 

# upload host specialization data and remove failed sims
hs_ad_all <- read.delim('hs_len20_vec1000_yrs500_sampled_pairwise_antigen_dists.tsv', header = F) 
hs_data <- read.delim('hs_len20_vec1000_yrs500_sim_data.tsv', header = T)
hs_runs_to_remove <- hs_data %>% filter(active_strain_count == 0) # find failed sims
hs_ad_all <- hs_ad_all %>% filter(!V1 %in% c(unique(hs_runs_to_remove$run_tag))) # remove failed sims
hs_sampled_all <- hs_ad_all[sample(nrow(hs_ad_all), 5000), ] # sample
hs_sampled_all$selection <- 'host specialization' # label 

# upload neg control data and remove failed sims
ctrl_ad_all <- read.delim('ctrl_len20_vec1000_yrs500_sampled_pairwise_antigen_dists.tsv', header = F)
ctrl_data <- read.delim('ctrl_len20_vec1000_yrs500_sim_data.tsv', header = F, col.names = colnames(cr_data))
ctrl_runs_to_remove <- ctrl_data %>% filter(active_strain_count == 0) # find failed sims
ctrl_ad_all <- ctrl_ad_all %>% filter(!V1 %in% c(unique(ctrl_runs_to_remove$V1))) # remove failed sims
ctrl_sampled_all <- ctrl_ad_all[sample(nrow(ctrl_ad_all), 5000), ] # sample 
ctrl_sampled_all$selection <- 'control' # label 

# combine sampled antigen distances from all conditions
ad_combined_all <- bind_rows(cr_sampled_all, hs_sampled_all, ctrl_sampled_all)

# plot multi-panel histogram for ant dist distributions sampled from all years of all sims 
ggplot(ad_combined_all, aes(x = V2, fill = selection)) +
  geom_histogram(binwidth = 1, color = 'black') +
  scale_fill_manual(values = c("cross reactivity" = "skyblue3", "host specialization" = "red3", "control" = "yellow2")) +
  facet_wrap(~ selection, scales = 'free_y', nrow = 1) +
  theme_bw() +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title = "Random sampling of pairwise antigen distances", subtitle = "5000 sampled distances from n = 100 simulations", x = "Distance", y = "Frequency") +
  theme(legend.position = "none", panel.grid = element_blank())

setwd('/Users/brandonely/Desktop/bb_strain_model_dev/new_sim_output/test_set_500yrs/plots/')
#ggsave('test_set_500_full_sim_AntDist.jpeg', height = 5, width = 8)


### FINAL POPULATION ANTIGEN DISTANCE DISTRIBUTIONS ###
setwd('/Users/brandonely/Desktop/bb_strain_model_dev/new_sim_output/test_set_500yrs/')
# upload data for computing antigen distances of final populations 
cr_final_pop_data <- read.delim('cr_len20_vec1000_yrs500_variant_frequencies.tsv', header = T, colClasses = c('character','numeric','numeric','character'))
cr_final_pop_data <- cr_final_pop_data %>% filter(!run_id %in% c(unique(cr_runs_to_remove$run_tag))) # remove failed sims
hs_final_pop_data <- read.delim('hs_len20_vec1000_yrs500_variant_frequencies.tsv', header = T, colClasses = c('character','numeric','numeric','character'))
hs_final_pop_data <- hs_final_pop_data %>% filter(!run_id %in% c(unique(hs_runs_to_remove$run_tag))) # remove failed sims
ctrl_final_pop_data <- read.delim('ctrl_len20_vec1000_yrs500_variant_frequencies.tsv', header = F, col.names = colnames(cr_final_pop_data), colClasses = c('character','numeric','numeric','character'))
ctrl_final_pop_data <- ctrl_final_pop_data %>% filter(!run_id %in% c(unique(ctrl_runs_to_remove$run_tag))) # remove failed sims

# final distributions for cross reactivity
sample_catcher <- list()
index <- 1
for (run in unique(cr_final_pop_data$run_id)) {
  filtered_df <- cr_final_pop_data %>% filter(run_id == run)
  dists <- final_antigen_distances(filtered_df)
  sampled_dists <- dists[sample(nrow(dists), 100), ]
  sample_catcher[[index]] <- sampled_dists
  index <- index + 1
}

all_dists <- unlist(sample_catcher)
cr <- sample(all_dists, size = 5000, replace = FALSE)
cr2 <- data.frame(values = cr)
cr2$selection <- 'cross reactivity'

# final distributions for host specialization
sample_catcher <- list()
index <- 1
for (run in unique(hs_final_pop_data$run_id)) {
  filtered_df <- hs_final_pop_data %>% filter(run_id == run)
  dists <- final_antigen_distances(filtered_df)
  sampled_dists <- dists[sample(nrow(dists), 100), ]
  sample_catcher[[index]] <- sampled_dists
  index <- index + 1
}

all_dists <- unlist(sample_catcher)
hs <- sample(all_dists, size = 5000, replace = FALSE)
hs2 <- data.frame(values = hs)
hs2$selection <- 'host specialization'

# final distributions for negative control
sample_catcher <- list()
index <- 1
for (run in unique(ctrl_final_pop_data$run_id)) {
  filtered_df <- ctrl_final_pop_data %>% filter(run_id == run)
  dists <- final_antigen_distances(filtered_df)
  sampled_dists <- dists[sample(nrow(dists), 100), ]
  sample_catcher[[index]] <- sampled_dists
  index <- index + 1
}

all_dists <- unlist(sample_catcher)
ctrl <- sample(all_dists, size = 5000, replace = FALSE)
ctrl2 <- data.frame(values = ctrl)
ctrl2$selection <- 'control'

# combine distributions to one df
set_dists <- bind_rows(cr2, hs2, ctrl2)

ggplot(set_dists, aes(x = values, fill = selection)) +
  geom_histogram(binwidth = 1, color = 'black') +
  scale_fill_manual(values = c("cross reactivity" = "skyblue3", "host specialization" = "red3", "control"= "yellow2")) +
  facet_wrap(~ selection, scales = 'free_y', nrow = 1) +
  labs(title = '5000 sampled distances, n = 100 simuations',x = 'Antigen Distance') +
  scale_y_continuous(expand = c(0, 0)) +
  #theme(legend.position = "top", panel.grid = element_blank()) +
  theme_bw() +
  #theme(strip.text = element_text(size = 15))
  theme(strip.text = element_blank(),
        plot.title = element_text(size = 20),
        plot.subtitle = element_text(size = 20),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.position = "bottom", 
        legend.title = element_blank())

setwd('/Users/brandonely/Desktop/bb_strain_model_dev/new_sim_output/test_set_500yrs/plots/')
#ggsave('test_set_500_final_sim_AntDist_v2.jpeg', height = 5, width = 12)


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
cr_ks <- sample(cr, size = 100, replace = FALSE)
hs_ks <- sample(hs, size = 100, replace = FALSE)
ctrl_ks <- sample(ctrl, size = 100, replace = FALSE)
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

ggplot(set_dists, aes(x = values, group = selection, color = selection)) +
  scale_color_manual(values = c("cross reactivity" = "skyblue3", "host specialization" = "red3", "control" = 'yellow2')) +
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


setwd('/Users/brandonely/Desktop/bb_strain_model_dev/new_sim_output/test_set_500yrs/plots/')
ggsave('test_set_500_allyrs_crVShs_ksTest_v2.jpeg', height = 7, width = 12)


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

  

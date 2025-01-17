# sim antigen distance distribution analysis 
setwd('/Users/brandonely/Desktop/bb_strain_model_dev/')
library(tidyverse)
library(ggplot2)
source("sim_analysis_functions.R")

# upload data from sim output, format, and sample 
setwd('/Users/brandonely/Desktop/bb_strain_model_dev/new_sim_output/')
cr_ad <- read.delim('cr_len20_vec1000_yrs250_sampled_pairwise_antigen_dists.tsv', header = F) 
cr_sampled <- cr_ad[sample(nrow(cr_ad), 5000), ] # sample
cr_sampled$source <- 'cross reactivity' # label 
hs_ad <- read.delim('hs_len20_vec1000_yrs250_sampled_pairwise_antigen_dists.tsv', header = F) 
# filter to remove failed sims
hs_data <- read.delim('hs_len20_vec1000_yrs250_sim_data.tsv', header = T)
runs_to_remove <- hs_data %>% filter(active_strain_count == 0)
hs_ad <- hs_ad %>% filter(!V1 %in% c(unique(runs_to_remove$run_tag)))
# resume
hs_sampled <- hs_ad[sample(nrow(hs_ad), 5000), ] # sample
hs_sampled$source <- 'host specialization' # label 
df_combined <- bind_rows(cr_sampled, hs_sampled) # combine

# plot multi-panel histogram for ant dist distributions sampled from all years of all sims 
ggplot(df_combined, aes(x = V2, fill = source)) +
  geom_histogram(binwidth = 1) +
  scale_fill_manual(values = c("cross reactivity" = "skyblue3", "host specialization" = "red3")) +
  facet_wrap(~ source, scales = 'free_y', nrow = 1) +
  theme_bw() +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title = "Random sampling of pairwise antigen distances", subtitle = "5000 sampled distances from n = 100 simulations" ,x = "Distance", y = "Frequency") +
  theme(legend.position = "none", panel.grid = element_blank())

# ant dist distribution in final population
cr_pop_data <- read.delim('cr_len20_vec1000_yrs250_variant_frequencies.tsv', header = T, colClasses = c('character','numeric','numeric','character'))
hs_pop_data <- read.delim('hs_len20_vec1000_yrs250_variant_frequencies.tsv', header = T, colClasses = c('character','numeric','numeric','character'))

# look at final number of variants in pop, i think i have this data elsewhere though...
# result <- cr_pop_data %>%
#   group_by(run_id) %>%
#   tally() %>%
#   arrange(desc(n)) 
# mean(result$n)
# max(result$n)

cr_run_pop <- cr_pop_data %>% filter(run_id == 'run_79') # choosing a run to work with
cr_ant_dists <- final_antigen_distances(cr_run_pop) # custom function for calculating distances
cr_ant_dists$selection <- 'cross reactivity' # label
cr_sampled_dists <- cr_ant_dists[sample(nrow(cr_ant_dists), 1000), ] # sample the distances 

hs_run_pop <- hs_pop_data %>% filter(run_id == 'run_29') # choosing a run to work with
hs_ant_dists <- final_antigen_distances(hs_run_pop) # custom function for calculating distances
hs_ant_dists$selection <- 'host specialization' # label 
hs_sampled_dists <- hs_ant_dists[sample(nrow(hs_ant_dists), 1000), ] # sample the distances 

combined_samples <- bind_rows(cr_sampled_dists, hs_sampled_dists)

# plot
ggplot(combined_samples, aes(x = distance, fill = selection)) +
  geom_histogram(binwidth = 1) +
  scale_fill_manual(values = c("cross reactivity" = "skyblue3", "host specialization" = "red3")) +
  facet_wrap(~ selection, scales = 'free_y', nrow = 1) +
  labs(title = "Pairwise antigen distances in final populations", 
       subtitle = "Single simulation run for each condition, n variants in final population = 20",
       x = "Distance",
       y = "Frequency") +
  theme_bw() +
  scale_y_continuous(expand = c(0, 0)) +
  theme(legend.position = "none", panel.grid = element_blank())


# statistical testing on distributions 
# run k-s test
ks_test <- ks.test(cr_sampled$V2, hs_sampled$V2)
# record ks test result to print on plot
ks_result <- paste('k-s test', "\n", 
                   "stat= ", round(ks_test$statistic, 3), "\n", 
                   "pval= ", format(ks_test$p.value, scientific = TRUE, digits = 3))

ggplot(df_combined, aes(x = V2, group = source, color = source))+
  scale_color_manual(values = c("cross reactivity" = "skyblue3", "host specialization" = "red3")) +
  stat_ecdf(size=1) +
  theme_bw() +
  theme(legend.position ="top") +
  xlab("value") +
  ylab("ECDF") +
  annotate("text", x = 2, y = 0.8, label = ks_result, size = 4, color = "black")

# qq plot
cr_sampled <- sort(cr_sampled$V2)
hs_sampled <- sort(hs_sampled$V2)
df <- data.frame(cr = cr_sampled, hs = hs_sampled)
ggplot(df, aes(x = cr, y = hs)) +
  geom_point() +
  geom_abline(aes(slope = 1, intercept = 0), linetype = 2, color = 'red') +
  labs(title = 'Q-Q for antigen distance distributions', 
       subtitle = 'immune cross reactivity VS host specialization',
       x = 'immune cross reactivity',
       y = 'host specialization') +
  theme_bw()


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

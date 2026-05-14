### this script looks at coexistence of generalists and specialists in the final populations of simulations; to run this analysis I need: 
### 1. cluster labels from the clustering analysis
### 2. host infection frequencies 

library(tidyverse)
library(broom)
setwd('/Users/brandonely/Desktop/bb_strain_model_dev/R_analysis_scripts/')
source("sim_analysis_functions.R")

##### upload and format input data #####
setwd('/Users/brandonely/Desktop/bb_strain_model_dev/evo_sim/results_tables/cluster_analysis_results/')

# get all the cluster labels from the clustering results
labels <- read.delim('hybrid_cluster_labels.tsv', sep = '\t', header = T, colClasses = c('character','character','numeric','numeric', 'character', 'character','character', 'character','character','character','character', 'character'))

# filter labels for single parameter set
labels <- labels %>% filter(gen_fitness == 'high' & host_pop_stability == 'stable')

#host infection file names
#'hybrid_ExpMod_host_infection_freqs.tsv'
#'flux_ExpMod_host_infection_freqs.tsv'
#'hybrid_LowGenFit_host_infection_freqs.tsv'
#'hybrid_LinearGenFit_host_infection_freqs.tsv'

# host infections table
setwd('/Users/brandonely/Desktop/bb_strain_model_dev/evo_sim/results_tables/sim_output/hybrid_results/')
header_names <- c('strain',	'rodent_count',	'bird_count',	'rodent_freq',	'bird_freq',	'combined_count',	'combined_freq',	'run_id',	'selection',	'gene_type',	'rec_rate',	'mut_rate',	'gen_fitness')

hi <- read.delim('hybrid_HighGenFit_host_infection_freqs.tsv', sep = '\t', header = T, col.names = header_names, colClasses = c('character','numeric','numeric','numeric','numeric','numeric','numeric', 'character','character','character','character','character','character'))

hi$host_pop_stability <- 'stable'

# break up the binary string into immune selection and adaptive selection components;
# if not using multi gene parameter, will have to rename col of strain as immune to work with cluster_host_assoc function 
hi <- hi %>%
  mutate(
    immune = substr(strain, 1, 20),
    adaptive = substr(strain, 21, 40)
  )

##### aggregate infection counts by cluster and assign preliminary classifications #####
res <- cluster_host_assoc(hi, labels)
res2 <- cluster_host_assoc(hi, labels)
res3 <- cluster_host_assoc(hi, labels)
#res4 <- cluster_host_assoc(hi, labels)
res_all <- bind_rows(res, res2, res3)
#write.table(res_all, '_host_assoc_by_threshold_res.tsv', sep = '\t', quote = F, row.names = F)

##### binomial GLM #####
setwd('/Users/brandonely/Desktop/bb_strain_model_dev/evo_sim/results_tables/host_association_results/')
res <- read.delim('hybrid_host_assoc_by_threshold_res.tsv', sep = '\t', header = T)
res <- res %>% filter(combined_count >= 5)

# enforce bounds (avoid 0 counts and totals >= 50)
res_bounded <- res %>%
  mutate(
    bird_count   = pmin(pmax(bird_count, 1), 49),
    rodent_count = pmin(pmax(rodent_count, 1), 49),
    bird_total   = 50,   # assuming 50 hosts sampled per type
    rodent_total = 50
  )

# long format for GLM
res_long <- res_bounded %>%
  pivot_longer(
    cols = c(bird_count, rodent_count, bird_total, rodent_total),
    names_to = c("host_type", ".value"),
    names_pattern = "(bird|rodent)_(count|total)"
  ) %>%
  mutate(
    not_infected = total - count,
    host_type = factor(host_type, levels = c("bird", "rodent"))
  )

# run GLM
glm_results <- res_long %>%
  group_by(selection, gen_fitness, host_pop_stability, run_id, cluster_label) %>%
  do(
    tidy(
      glm(cbind(count, not_infected) ~ host_type, 
          data = ., 
          family = binomial)
    )
  ) %>%
  ungroup()


# classify based on coefficient and significance
glm_results <- glm_results %>%
  filter(term == "host_typerodent") %>%  # coefficient comparing rodents to birds
  mutate(
    classification = case_when( estimate > 0 & p.value < 0.05  ~ "rodent_specialist",
                                estimate < 0 & p.value < 0.05  ~ "bird_specialist",
                                TRUE ~ "generalist"),
    odds_ratio = exp(estimate)
    )


# clean output
df <- glm_results %>% select(selection, gen_fitness, host_pop_stability, run_id, cluster_label, estimate, odds_ratio, p.value, classification)

# combine all outputs to 1 df 
#df_all <- bind_rows(df, df2, df3, df4)
#rm(list = setdiff(ls(), "df_all"))
#write.table(df, 'adaptive_host_assoc_GLM_res_p05.tsv', quote = F, sep = '\t', row.names = F)

##### get community compositions #####
df_all <- read.delim('adaptive_host_assoc_GLM_res_p05.tsv', sep = '\t', header = T)

# count the number of runs for each parameter set, regardless of outcome
parameter_counts <- df_all %>%
  group_by(selection, gen_fitness, host_pop_stability) %>%
  summarise(
    n_runs = n_distinct(run_id),
    .groups = "drop"
  )

# for each run, tally the number of each type
run_counts <- df_all %>%
  group_by(selection, gen_fitness, host_pop_stability, run_id) %>%
  summarise(
    n_generalist = sum(classification == "generalist"),
    n_bird = sum(classification == "bird_specialist"),
    n_rodent = sum(classification == "rodent_specialist"),
    .groups = "drop"
  )

# grab run_id's for different categories 
gss_id <- run_counts %>%
  filter(n_generalist > 0 & n_bird > 0 & n_rodent > 0 & gen_fitness == 'high') %>%
  pull(run_id)
gss_id <- data.frame(run_id = gss_id)
#write.table(gss_id, 'hybrid_ExpMod_coexistence_runIDs.tsv',quote = F, sep = '\t', row.names = F)

# categorize community for each run
coexist_res <- run_counts %>%
  mutate(
    community = if_else(n_generalist > 0 & n_bird > 0 & n_rodent > 0, 'gss',
                if_else(n_generalist > 0 & ((n_bird > 0 & n_rodent == 0) | (n_bird == 0 & n_rodent > 0)), 'gs', 
                if_else(n_generalist > 0 & n_bird == 0 & n_rodent == 0, 'g', 
                if_else(n_generalist == 0 & n_bird > 0 & n_rodent > 0, 'ss',
                if_else(n_generalist == 0 & ((n_bird > 0 & n_rodent == 0) | (n_bird == 0 & n_rodent > 0)), 's', 'mistake')))))
  )

# tally counts for each community type
coexist_res <- coexist_res %>%
  group_by(selection, gen_fitness, host_pop_stability, community) %>%
  summarise(
    counts = n(),
    .groups = 'drop'
  )

# add col for total runs for each parameter
coexist_res <- left_join(coexist_res, parameter_counts)

# calculate proportion of runs resulting in each community type
coexist_res$probability <- coexist_res$counts/coexist_res$n_runs

# make a col for full parameter description
coexist_res$category <- paste(coexist_res$selection, 
                                      coexist_res$gen_fitness,
                                      coexist_res$host_pop_stability,
                                      sep = '_')


# drop unneccessary cols
coexist_res <- coexist_res %>% select(8,5,6,7,4)

# add communities that are missing from certain parameters, filled with 0 probability
coexist_res <- coexist_res %>%
  complete(category, community, fill = list(value = 0)) %>%
  {.[is.na(.)] <- 0; .}

# save this table for viz later
#write.table(coexist_res, 'adaptive_community_composition_probabilities_p05.tsv', quote = F, sep = '\t', row.names = F)


##### very simple community comp check #####
setwd('/Users/brandonely/Desktop/bb_strain_model_dev/evo_sim/results_tables/sim_output/adaptive_results/')

header_names <- c('strain', 'count', 'freq', 'run_id', 'v', 'w', 'x', 'y', 'z')

df <- read.delim('adaptive_ExpMod_variant_freqs_sampled.tsv', sep = '\t', header = T, col.names = header_names, colClasses = c('character','numeric','numeric', 'character', 'character','character', 'numeric','numeric', 'character'))

df <- read.delim('adaptive_HighGenFit_variant_freqs_sampled.tsv', sep = '\t', header = T, colClasses = c('character','numeric','numeric', 'character', 'character','character', 'numeric','numeric', 'character'))


df <- df %>% 
  select(1:4) %>%
  filter(freq >=1)

df <- df %>% mutate(val = str_count(pattern = '1', strain)/20)
df <- df %>% mutate(type = 
                      if_else(val <= 0.32, 'bird', 
                              if_else(val >= 0.68, 'rodent', 'generalist')))


communities <- df %>%
  group_by(run_id, type) %>%
  summarise(total = sum(count), .groups = "drop") %>%
  pivot_wider(
    names_from = type,
    values_from = total,
    values_fill = 0
  )

presence_by_run <- communities %>%
  mutate(
    bird_present = bird > 0,
    rodent_present = rodent > 0,
    generalist_present = generalist > 0
  ) %>%
  select(run_id, bird_present, rodent_present, generalist_present)

combo_summary <- presence_by_run %>%
  group_by(bird_present, rodent_present, generalist_present) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(proportion = n / sum(n))

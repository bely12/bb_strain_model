# Gupta simulation
library(tidyverse)
setwd('/Users/brandonely/Desktop/bb_strain_model_dev/gupta_sim/gupta_model_results/')

#df <- read.delim('gupta_model_results.tsv', header = T, sep = '\t', colClasses = c(genotype = 'character'))
df <- read.delim('gupta_res_frequencies.tsv', header = T, sep = '\t', colClasses = c(genotype = 'character'))

# rewrite genotypes so avoid starting with a number and getting screwy in R
df_clean <- df %>%
  mutate(genotype = paste0("g", genotype))

link_dis <- df_clean %>%
  group_by(run_id, year, cp, rec_rate) %>%
  pivot_wider(names_from = genotype, values_from = frequency, values_fill = 0) %>%
  mutate(
    A0 = g00 + g01,
    A1 = g10 + g11,
    B0 = g00 + g10,
    B1 = g01 + g11,
    D = g00 - (A0 * B0),
    r2 = ifelse(A0 * A1 * B0 * B1 > 0, (D^2) / (A0 * A1 * B0 * B1), NA_real_)) %>%
  select(run_id, year, cp, rec_rate, D, r2)


cp_colors <- c(
  high   = "red",
  medium = "orange",
  low    = "blue",
  none   = "skyblue"
)

link_dis$cp <- factor(link_dis$cp, levels = c("high", "medium", "low", "none"))

ggplot(link_dis, aes(x = year, y = r2, color = cp)) +
  #geom_line(aes(group = run_id), alpha = 0.3) + 
  stat_summary(fun = mean, geom = "line", linewidth = 1.2, aes(group = cp)) +
  labs(y = expression(r^2), x = "year", color = "cross protection") +
  scale_color_manual(values = cp_colors) +
  facet_grid(rec_rate ~ .) +
  theme_bw() +
  theme(strip.text = element_text(size = 16),
        axis.title = element_text(size = 20), 
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20))
ggsave('ld_gupta_dis.png', height = 10, width = 10)

# look at some single runs
run_samples <- df %>%
  distinct(run_id, cp, rec_rate) %>%
  group_by(cp, rec_rate) %>%
  sample_n(1) %>%
  ungroup()

df_sampled <- df %>%
  semi_join(run_samples, by = c("run_id", "cp", "rec_rate"))

# facet wrap
# ggplot(df_sampled, aes(x = year, y = frequency, color = genotype)) +
#   geom_line(aes(group = genotype), alpha = 0.8, size = 1) +
#   facet_wrap(vars(cp, rec_rate, run_id), nrow = 4) +
#   labs(title = "Genotype Frequencies Over Time by Simulation Run",
#        x = "Year", y = "Frequency", color = "Genotype") +
#   theme_bw() +
#   theme(strip.text = element_text(size = 10))

# facet grid
df_sampled$cp <- factor(df_sampled$cp, levels = c("high", "medium", "low", "none"))

ggplot(df_sampled, aes(x = year, y = frequency, color = genotype)) +
  geom_line(aes(group = genotype), alpha = 0.8, size = 1) +
  facet_grid(cp~rec_rate, 
             labeller = labeller(rec_rate = function(x) paste("recombination rate =", x))) +
  labs(x = "Year", y = "Frequency", color = "Genotype") +
  theme_bw() +
  theme(strip.text = element_text(size = 16), 
        axis.title = element_text(size = 20), 
        axis.text = element_text(size = 20),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20))

ggsave('gupta_freqs_dis.png', height = 10, width = 10)


################# GUPTA ONE2ONE MODEL #####################

df <- read.delim('gupta_one2one_res_frequencies.tsv', header = T, sep = '\t')

df <- df %>%
  mutate(specialization = if_else(substr(genotype, 2, 2) == 'X', 'rodent', if_else(substr(genotype, 2, 2) == 'Y', 'bird', 'generalist')),
         antigen = substr(genotype, 1, 1))

plot_df <- df %>% filter(run_id == 'run_5')

ggplot(plot_df, aes(x = year, y = frequency, color = antigen, linetype = specialization)) +
  geom_line(aes(group = genotype), alpha = 0.8, size = 1) +
  facet_wrap(~ antigen, ncol = 1) +
  scale_linetype_manual(values = c("rodent" = "dashed", "bird" = "dotted", "generalist" = "solid")) +
  theme_bw()


################# GUPTA plus adaptive MODEL #####################
setwd('/Users/brandonely/Desktop/bb_strain_model_dev/gupta_sim/')
df <- read.delim('gupta_plus_adprec_res_variant_frequencies.tsv', header = T, sep = '\t')

df <- df %>%
  mutate(specialization = if_else(substr(genotype, 5, 5) == 'X', 'rodent', if_else(substr(genotype, 5, 5) == 'Y', 'bird', 'generalist')),
         antigen = substr(genotype, 1, 4))

plot_df <- df %>% filter(run_id == 'run_35')

ggplot(plot_df, aes(x = year, y = frequency, color = antigen, linetype = specialization)) +
  geom_line(aes(group = genotype), alpha = 0.8, linewidth = 1) +
  facet_wrap(~ specialization, ncol = 1) +
  scale_linetype_manual(values = c("rodent" = "dashed", "bird" = "dotted", "generalist" = "solid")) +
  theme_bw()

# this gets the mean frequency for each adaptive category, separated by comp_factor and year 
stats <- df %>%
  group_by(year, specialization, comp_factor, cp, rec_rate, run_id) %>%
  summarise(total_freq = sum(frequency)) %>%
  group_by(year, specialization, comp_factor, cp, rec_rate) %>%
  summarise(mean_freq = mean(total_freq)) 

stats %>% #filter(cp == 'linear') %>%
  ggplot(aes(x = year, y = mean_freq, color = specialization)) +
  geom_line() +
  facet_grid(vars(comp_factor, rec_rate))


x <- stats %>% filter(cp == 'linear' & year == 100)
# coexist <- df %>%
#   filter(year == 100) %>%
#   group_by(specialization, comp_factor, run_id) %>%
#   summarise()



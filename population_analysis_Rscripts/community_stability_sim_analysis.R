setwd('/Users/brandonely/Desktop/bb_strain_model_dev/')
library(tidyverse)
library(circular)
library(pheatmap)
library(viridis) 
library(wesanderson)
library(stringr)

# load simulation results
setwd('/Users/brandonely/Desktop/bb_strain_model_dev/static_sim')
df <- read.delim('static_sim_bash_batch_static_sim_results.tsv', header = T)

# get circular variance
circular_variance <- function(x) {
  num_vals <- as.numeric(unlist(strsplit(gsub("\\[|\\]", "", x), ",")))
  angles <- pi * num_vals
  circular_points <- circular(angles, units = 'radians')
  return(var.circular(circular_points))
}
df$new_antigen_var <- sapply(df$antigen_vals, circular_variance)

# round antigen dispersion value to nearest 0.05
round_to_05 <- function(x) {
  rounded <- round(x / 0.05) * 0.05
  return(rounded)
}
df$antigen_dispersion <- sapply(df$new_antigen_var, round_to_05)
#df$antigen_dispersion<- round(df$new_antigen_var, 2)

# function to round hundreths place to even number
round_even_decimal <- function(x) {
  x <- x * 100 # convert value to a whole number, seems easier this way
  if (x %% 2 != 0) { 
    x <- ceiling(x / 2) * 2  # round up to nearest even number if its odd
  }
  return(x / 100) # convert back to a decimal 
}
df$spec_weight <- sapply(df$spec_weight, round_even_decimal)

#z <- df %>% filter(spec_weight == 0.1)
# cleanup
#new_df <- df %>% select(5,9,11,13)
#df <- df %>% select(13,8,11)
new_df <- df %>% select(7,12,10)

# trait value distributions
# hist(df$antigen_dispersion)
# hist(df$spec_weight)
ggplot(new_df, aes(x = df$antigen_dispersion)) +
  geom_histogram(binwidth = 0.05, color = 'black', fill = 'forestgreen')
ggplot(new_df, aes(x = df$spec_weight)) +
  geom_histogram(binwidth = 0.05, color = 'black', fill = 'purple3')

new_df %>%
  ggplot(aes(x = spec_weight, y = antigen_dispersion, color = result)) +
  geom_jitter()

# get probabilities for every pair of antigen/specialization value, plus counts for each pair
df2 <- df %>%
  group_by(spec_weight, antigen_dispersion) %>%
  summarise(stable_prob = mean(result == "stable"), counts = n()) %>%
  ungroup()

# counts and stable probability distributions
# hist(df2$counts)
# hist(df2$stable_prob)
ggplot(df2, aes(x = df2$counts)) +
  geom_histogram(binwidth = 1, color = 'black', fill = 'forestgreen')
ggplot(df2, aes(x = df2$stable_prob)) +
  geom_histogram(binwidth = 0.05, color = 'black', fill = 'purple3')

# convert to wide/matrix style for heatmap 
df_wide <- df2 %>% 
  select(1:3) %>%
  pivot_wider(names_from = spec_weight, values_from = stable_prob) %>%
  column_to_rownames(var = "antigen_dispersion")

# sort rows for heatmap
df_wide <- df_wide[rev(order(as.numeric(rownames(df_wide)))), ]

# color pallete
pal <- wes_palette("Zissou1", 100, type = "continuous")

#setwd('/Users/brandonely/Desktop/bb_strain_model_dev/static_sim/plots')
# make heatmap
pheatmap(df_wide,
         color = pal,
         cluster_rows = FALSE,  
         cluster_cols = FALSE,  
         display_numbers = FALSE,
         border_color = FALSE,
         show_rownames = T,
         show_colnames = T,
         cellwidth = 10,
         cellheight = 10)
         #filename = 'stable_communities_CircVar_heatmap.jpeg')  


# fill the grid!!!
df_wide_counts <- df2 %>% 
  select(1:2,4) %>%
  pivot_wider(names_from = spec_weight, values_from = counts) %>%
  column_to_rownames(var = "antigen_dispersion")

# sort rows for heatmap
df_wide_counts <- df_wide_counts[rev(order(as.numeric(rownames(df_wide_counts)))), ]

# color pallete
#pal2 <- colorRampPalette(c("white", "red"))(112)
pal2 <- c("red", 'orange', 'yellow2',"green3")
breaks <- c(min(df_wide_counts, na.rm = TRUE),8,16,24,max(df_wide_counts, na.rm = TRUE))
# make heatmap
pheatmap(df_wide_counts,
         color = pal2,
         breaks = breaks,
         cluster_rows = FALSE,  
         cluster_cols = FALSE,  
         display_numbers = T,
         number_format = "%.0f",
         border_color = FALSE,
         show_rownames = T,
         show_colnames = T,
         cellwidth = 20,
         cellheight = 20)  



# structured pops stability
# keep only sims that had at least one hs val at 0.5
# new_df2 <- new_df %>%
#   filter(str_detect(hs_vals, "\\b0\\.5\\b"))
# 
# # split the hs vals out of brackets into cols as numerics
# new_df3 <- new_df2 %>%
#   # remove brackets
#   mutate(hs_vals_clean = str_remove_all(hs_vals, "\\[|\\]")) %>%
#   # separate into 3 new columns
#   separate(hs_vals_clean, into = c("hs1", "hs2", "hs3"), sep = ",\\s*") %>%
#   # convert to numeric
#   mutate(across(c(hs1, hs2, hs3), as.numeric)) %>%
#   # remove brackets
#   mutate(hs_vals_clean2 = str_remove_all(final_hs, "\\[|\\]")) %>%
#   # separate into 3 new columns
#   separate(hs_vals_clean2, into = c("hsf1", "hsf2", "hsf3"), sep = ",\\s*") %>%
#   # convert to numeric
#   mutate(across(c(hsf1, hsf2, hsf3), as.numeric)) 
# 
# # make sure at least one of them is a specialists above 0.5
# new_df4 <- new_df3 %>% filter(hs1 > 0.5 | hs2 > 0.5 | hs3 > 0.5)
# 
# temp_df <- new_df4 %>% filter(hs1 == 0.5)
# x1 <- data.frame(x = temp_df$hs2, y = temp_df$hs3, result = temp_df$result, final_pop = temp_df$final_hs, antigen = temp_df$antigen_dispersion)
# temp_df <- new_df4 %>% filter(hs2 == 0.5 & hs1 != 0.5)
# x2 <- data.frame(x = temp_df$hs1, y = temp_df$hs3, result = temp_df$result, final_pop = temp_df$final_hs, antigen = temp_df$antigen_dispersion)
# temp_df <- new_df4 %>% filter(hs3 == 0.5 & hs1 != 0.5 & hs2 != 0.5)
# x3 <- data.frame(x = temp_df$hs1, y = temp_df$hs2, result = temp_df$result, final_pop = temp_df$final_hs, antigen = temp_df$antigen_dispersion)
# 
# df <- bind_rows(x1,x2,x3)
# df2 <- df %>% mutate(new_x = if_else(x < 0.5, y, x)) %>%
#   mutate(new_y = if_else(x < 0.5, x, y))
# 
# 
# ggplot(df2, aes(x=new_x, y=new_y, colour = result)) +
#   geom_jitter() +
#   xlim(0.5, 1.0) +
#   facet_wrap(~antigen)


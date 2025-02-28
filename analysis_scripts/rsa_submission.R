rm(list = ls())
filepath <- "~/nat_behav_submission/" # filepath to GitHub folder
setwd(paste(filepath, 'data/fMRI/rsa/', sep = ''))
library(ggplot2)
library(tidyverse)
library(gghalves)
library(ggdist)
library(rstatix)
library(ggpubr)
library(latex2exp)
library(dplyr)
library(scales)
library(reshape2)
library(lme4)

############# Prepare data ###############

d_rsa <- read.csv('rsa_output.csv')
head(d_rsa)
d_rsa$constant <- 1

# compute mean split according to behvaioural change between environments
mean_behav_change <- mean(d_rsa$p_accept_10_diff, na.rm = T)
d_rsa$mean_split <- ifelse(d_rsa$p_accept_10_diff > mean_behav_change, 'adaptive', 'non-adaptive')

# Create similarity matrices for later plotting
create_similarity_matrix <- function(d_rsa, d1_col, d2_col, d3_col, d4_col, env_label) {
  expand.grid(subject = unique(d_rsa$subject),
              opt_1 = c('low', 'mid', 'high'),
              opt_2 = c('low', 'mid', 'high')) %>%
    left_join(d_rsa %>% select(subject, p_accept_10_diff, mean_split, all_of(c(d1_col, d2_col, d3_col, d4_col))), by = "subject") %>%
    mutate(
      d_val = case_when(
        opt_1 == opt_2 ~ NA_real_,
        opt_1 == 'low' & opt_2 == 'mid' | opt_1 == 'mid' & opt_2 == 'low' ~ !!sym(d3_col),
        opt_1 == 'mid' & opt_2 == 'high' | opt_1 == 'high' & opt_2 == 'mid' ~ !!sym(d4_col),
        opt_1 == 'low' & opt_2 == 'high' | opt_1 == 'high' & opt_2 == 'low' ~ NA_real_,
        TRUE ~ NA_real_
      ),
      env = env_label
    ) %>%
    select(subject, opt_1, opt_2, d_val, p_accept_10_diff, env, mean_split)
}
compute_mean_distance_change <- function(df, opt1, opt2, env) {
  df %>%
    filter(behav_diff >= mean_behav_diff, opt_1 == opt1, opt_2 == opt2, env == env) %>%
    summarise(m = mean(d_val, na.rm = TRUE))
}

# create dACC matrix
dACC_rich <- create_similarity_matrix(d_rsa, "dACC_d1", "dACC_d2", "dACC_d3", "dACC_d4", "rich")
dACC_poor <- create_similarity_matrix(d_rsa, "dACC_d1", "dACC_d2", "dACC_d1", "dACC_d2", "poor")
d_rsa_dACC <- bind_rows(dACC_rich, dACC_poor)

d_rsa_dACC_summary <- d_rsa_dACC %>%
  group_by(mean_split, opt_1, opt_2, env) %>%
  summarise(mean_d_val = mean(d_val, na.rm = TRUE), .groups = "drop")

dACC_global_stats <- d_rsa_dACC_summary %>%
  summarise(
    global_min = min(mean_d_val, na.rm = TRUE),
    global_max = max(mean_d_val, na.rm = TRUE),
    global_median = median(mean_d_val, na.rm = TRUE)
  )

# create AI matrix
AI_rich <- create_similarity_matrix(d_rsa, "AI_d1", "AI_d2", "AI_d3", "AI_d4", "rich")
AI_poor <- create_similarity_matrix(d_rsa, "AI_d1", "AI_d2", "AI_d1", "AI_d2", "poor")
d_rsa_AI <- bind_rows(AI_rich, AI_poor)

d_rsa_AI_summary <- d_rsa_AI %>%
  group_by(mean_split, opt_1, opt_2, env) %>%
  summarise(mean_d_val = mean(d_val, na.rm = TRUE), .groups = "drop")

AI_global_stats <- d_rsa_AI_summary %>%
  summarise(
    global_min = min(mean_d_val, na.rm = TRUE),
    global_max = max(mean_d_val, na.rm = TRUE),
    global_median = median(mean_d_val, na.rm = TRUE)
  )

######################## Figure 3 ##########################

t.test(d_rsa$dACC_H2[d_rsa$mean_split=='adaptive'], d_rsa$dACC_H2[d_rsa$mean_split=='non-adaptive'])
t.test(d_rsa$AI_H2[d_rsa$mean_split=='adaptive'], d_rsa$dACC_H2[d_rsa$mean_split=='non-adaptive'])

# Fig. 3C-i
cor.test(d_rsa$dACC_H2, d_rsa$p_accept_10_diff)
ggplot(d_rsa, aes(x = p_accept_10_diff*100, y = dACC_H2)) + 
  geom_point(shape = 21, fill = 'red2', size = 2.5) +
  geom_smooth(method = 'lm', color = 'red4', linewidth = 0.65, alpha = 0.50, se = F, linetype = 2) + 
  scale_y_continuous(name = 'H2(d4 - d2)', labels = scales::number_format(accuracy = 0.01)) + 
  scale_x_continuous(name = TeX("$\\Delta$%acc. [Poor - Rich]")) + 
  theme_pubr() + 
  theme(
    aspect.ratio = 1, 
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 18),
    strip.text = element_text(size = 18)
  ) + 
  facet_grid(~factor(constant, labels = 'dACC')) + 
  border(color = 'black')

# Fig 3C-ii
ggplot(data = d_rsa_dACC_summary %>% filter(mean_split=='adaptive'), 
       aes(x = factor(opt_1, levels = c('high', 'mid', 'low')), 
           y = factor(opt_2, levels = c('low', 'mid', 'high')), 
           fill = mean_d_val)) +  
  geom_tile(color = 'black', size = 0.5) +  
  scale_fill_gradient2(name = 'Cosine distance', 
                       low = "#B2182B", mid = "#F7F7F7", high = "#2166AC", 
                       midpoint = dACC_global_stats$global_median, 
                       limits = c(dACC_global_stats$global_min-0.01, dACC_global_stats$global_max+0.01), 
                       breaks = round(c(dACC_global_stats$global_min-0.01, (dACC_global_stats$global_min+dACC_global_stats$global_max)/2, dACC_global_stats$global_max+0.01), 2),
                       space = "Lab") + 
  scale_x_discrete(name = '', labels = c('High', 'Mid', 'Low')) +  
  scale_y_discrete(name = '', labels = c('Low', 'Mid', 'High')) +  
  theme_pubr() +
  theme(
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    strip.text.x = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 14),
    aspect.ratio = 1,
    strip.text.y = element_blank(),
    strip.background = element_rect(colour='black', fill='grey99'),
    panel.spacing.x = unit(0.2, "lines"),
    legend.position = 'right'
  ) + 
  facet_grid(~factor(env, levels = c('poor', 'rich'), labels = c('Poor', 'Rich'))) +  
  coord_fixed()

# Plot heatmap with mean values per condition
ggplot(data = d_rsa_dACC_summary %>% filter(mean_split=='non-adaptive'), 
       aes(x = factor(opt_1, levels = c('high', 'mid', 'low')), 
           y = factor(opt_2, levels = c('low', 'mid', 'high')), 
           fill = mean_d_val)) +  
  geom_tile(color = 'black', size = 0.5) +  
  scale_fill_gradient2(name = 'Cosine distance', 
                       low = "#B2182B", mid = "#F7F7F7", high = "#2166AC", 
                       midpoint = dACC_global_stats$global_median, 
                       limits = c(dACC_global_stats$global_min-0.01, dACC_global_stats$global_max+0.01), 
                       breaks = round(c(dACC_global_stats$global_min-0.01, (dACC_global_stats$global_min+dACC_global_stats$global_max)/2, dACC_global_stats$global_max+0.01), 2),
                       space = "Lab") + 
  scale_x_discrete(name = '', labels = c('High', 'Mid', 'Low')) +  
  scale_y_discrete(name = '', labels = c('Low', 'Mid', 'High'), position = 'right') +  
  theme_pubr() +
  theme(
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    strip.text.x = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 14),
    aspect.ratio = 1,
    strip.text.y = element_blank(),
    strip.background = element_rect(colour='black', fill='grey99'),
    panel.spacing.x = unit(0.2, "lines"),
    legend.position = 'right'
  ) + 
  facet_grid(~factor(env, levels = c('poor', 'rich'), labels = c('Poor', 'Rich'))) +  
  coord_fixed() # Ensure square tiles

# Fig. 3D-i
# AI H2 vs behaviour
cor.test(d_rsa$AI_H2, d_rsa$p_accept_10_diff)
ggplot(d_rsa, aes(x = p_accept_10_diff*100, y = AI_H2)) + 
  geom_point(shape = 21, fill = 'lightblue1', size = 2.5) +
  geom_smooth(method = 'lm', color = 'blue4', linewidth = 0.65, alpha = 0.50, se = F, linetype = 2) + 
  scale_y_continuous(name = 'H2(d4 - d2)', labels = scales::number_format(accuracy = 0.01)) + 
  scale_x_continuous(name = TeX("$\\Delta$%acc. [Poor - Rich]")) + 
  theme_pubr() + 
  theme(
    aspect.ratio = 1, 
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 18),
    strip.text = element_text(size = 18)
  ) + 
  facet_grid(~factor(constant, labels = 'AI')) + 
  border(color = 'black')

ggplot(data = d_rsa_AI_summary %>% filter(mean_split=='adaptive'), 
       aes(x = factor(opt_1, levels = c('high', 'mid', 'low')), 
           y = factor(opt_2, levels = c('low', 'mid', 'high')), 
           fill = mean_d_val)) +  
  geom_tile(color = 'black', size = 0.5) +  
  scale_fill_gradient2(name = 'Cosine distance', 
                       low = "#B2182B", mid = "#F7F7F7", high = "#2166AC", 
                       midpoint = AI_global_stats$global_median, 
                       limits = c(AI_global_stats$global_min-0.01, AI_global_stats$global_max+0.01), 
                       breaks = round(c(AI_global_stats$global_min-0.01, (AI_global_stats$global_min+AI_global_stats$global_max)/2, AI_global_stats$global_max+0.01), 2),
                       space = "Lab") + 
  scale_x_discrete(name = '', labels = c('High', 'Mid', 'Low')) +  
  scale_y_discrete(name = '', labels = c('Low', 'Mid', 'High')) +  
  theme_pubr() +
  theme(
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    strip.text.x = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 14),
    aspect.ratio = 1,
    strip.text.y = element_blank(),
    strip.background = element_rect(colour='black', fill='grey99'),
    panel.spacing.x = unit(0.2, "lines"),
    legend.position = 'right'
  ) + 
  facet_grid(~factor(env, levels = c('poor', 'rich'), labels = c('Poor', 'Rich'))) +  
  coord_fixed()

# Plot heatmap with mean values per condition
ggplot(data = d_rsa_AI_summary %>% filter(mean_split=='non-adaptive'), 
       aes(x = factor(opt_1, levels = c('high', 'mid', 'low')), 
           y = factor(opt_2, levels = c('low', 'mid', 'high')), 
           fill = mean_d_val)) +  
  geom_tile(color = 'black', size = 0.5) +  
  scale_fill_gradient2(name = 'Cosine distance', 
                       low = "#B2182B", mid = "#F7F7F7", high = "#2166AC", 
                       midpoint = AI_global_stats$global_median, 
                       limits = c(AI_global_stats$global_min-0.01, AI_global_stats$global_max+0.01), 
                       breaks = round(c(AI_global_stats$global_min-0.01, (AI_global_stats$global_min+AI_global_stats$global_max)/2, AI_global_stats$global_max+0.01), 2),
                       space = "Lab") + 
  scale_x_discrete(name = '', labels = c('High', 'Mid', 'Low')) +  
  scale_y_discrete(name = '', labels = c('Low', 'Mid', 'High'), position = 'right') +  
  theme_pubr() +
  theme(
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    strip.text.x = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 14),
    aspect.ratio = 1,
    strip.text.y = element_blank(),
    strip.background = element_rect(colour='black', fill='grey99'),
    panel.spacing.x = unit(0.2, "lines"),
    legend.position = 'right'
  ) + 
  facet_grid(~factor(env, levels = c('poor', 'rich'), labels = c('Poor', 'Rich'))) +  
  coord_fixed()

#################### Supplementary Figures

# Supplementary figure S5B
cor.test(d_rsa$dACC_H1, d_rsa$p_accept_10_diff)
ggplot(d_rsa, aes(x = p_accept_10_diff*100, y = dACC_H1)) + 
  geom_point(shape = 21, fill = 'red2', size = 2.5) +
  geom_smooth(method = 'lm', color = 'red4', linewidth = 0.65, alpha = 0.50, se = F, linetype = 2) + 
  scale_y_continuous(name = 'H1(d3 - d1)', labels = scales::number_format(accuracy = 0.01)) + 
  scale_x_continuous(name = TeX("$\\Delta$%acc. [Poor - Rich]")) + 
  theme_pubr() + 
  theme(
    aspect.ratio = 1, 
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 18),
    strip.text = element_text(size = 18)
  ) + 
  facet_grid(~factor(constant, labels = 'dACC')) + 
  border(color = 'black')

cor.test(d_rsa$AI_H1, d_rsa$p_accept_10_diff)
ggplot(d_rsa, aes(x = p_accept_10_diff*100, y = AI_H1)) + 
  geom_point(shape = 21, fill = 'lightblue1', size = 2.5) +
  geom_smooth(method = 'lm', color = 'blue4', linewidth = 0.65, alpha = 0.50, se = F, linetype = 2) + 
  scale_y_continuous(name = 'H1(d3 - d1)', labels = scales::number_format(accuracy = 0.01)) + 
  scale_x_continuous(name = TeX("$\\Delta$%acc. [Poor - Rich]")) + 
  theme_pubr() + 
  theme(
    aspect.ratio = 1, 
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 18),
    strip.text = element_text(size = 18)
  ) + 
  facet_grid(~factor(constant, labels = 'AI')) + 
  border(color = 'black')

# Supplementary figure S5C
H1_long <- d_rsa %>%
  pivot_longer(cols = ends_with('_H1'), names_to = 'region', values_to = 'H1_distance_change') %>%
  select(subject, region, H1_distance_change)

ggplot(H1_long, aes(x = region, y = H1_distance_change)) + 
  stat_summary(fun = 'mean', geom = 'bar', color = 'black', fill = 'grey', width = 0.85) + 
  stat_summary(fun.data = 'mean_ci', geom = 'errorbar', color = 'black', width = 0.35) + 
  geom_point(fill = 'white', alpha = 0.50, color = 'black', shape = 21) + 
  geom_hline(yintercept = 0, linetype = 2, color = 'red') + 
  scale_x_discrete(name = 'ROI', labels = c('AI', 'dACC', 'dmPFC', 'lFPC', 'pgACC', 'sgACC')) + 
  scale_y_continuous(name = 'H1 [d3 - d1]') + 
  theme_pubr() + 
  theme(
    aspect.ratio = 0.75, 
    axis.text = element_text(size = 14),
    axis.text.x = element_text(size = 14, vjust = 0.75, angle = 30),
    axis.title = element_text(size = 18),
    strip.text = element_text(size = 18),
    legend.position = 'none'
  )

# Supplementary figure S5D

H2_long <- d_rsa %>%
  pivot_longer(cols = ends_with('_H2'), names_to = 'region', values_to = 'H2_distance_change') %>%
  select(subject, region, H2_distance_change)

ggplot(H2_long, aes(x = region, y = H2_distance_change)) + 
  stat_summary(fun = 'mean', geom = 'bar', color = 'black', fill = 'grey', width = 0.85) + 
  stat_summary(fun.data = 'mean_ci', geom = 'errorbar', color = 'black',  width = 0.35) + 
  geom_point(fill = 'white', alpha = 0.50, color = 'black', shape = 21) + 
  geom_hline(yintercept = 0, linetype = 2, color = 'red') + 
  scale_x_discrete(name = 'ROI', labels = c('AI', 'dACC', 'dmPFC', 'lFPC', 'pgACC', 'sgACC')) + 
  scale_y_continuous(name = 'H2 [d4 - d2]') + 
  theme_pubr() + 
  theme(
    aspect.ratio = 0.75, 
    axis.text = element_text(size = 14),
    axis.text.x = element_text(size = 14, vjust = 0.75, angle = 30),
    axis.title = element_text(size = 18),
    strip.text = element_text(size = 18),
    legend.position = 'none'
  )

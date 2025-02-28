rm(list = ls())
filepath <- "~/nat_behav_submission/" # filepath to GitHub folder
setwd(paste(filepath, 'data/fMRI/timecourse_outputs/', sep = ''))
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

####################### TIMECOURSE GRAPHS ######################
file_paths <- dir(pattern = '*_timecourse.csv')

# Read and process CSV files using lapply
d_list <- lapply(file_paths, function(path) {
  # Read the CSV file
  df <- read.csv(path, header = T)
  
  # Add a 'time' column
  df$timepoint <- 1:nrow(df)
  
  # Return the processed data frame
  return(df)
})

d <- bind_rows(d_list) %>%
  rename(
    m = tc_data1,
    se = tc_data2,
    region = tc_data3,
    glm = tc_data4,
    contrast = tc_data5,
  ) %>%
  mutate(
    region = recode(region,
                    "dACCsphere7" = "dACC",
                    "AIsphere7" = "AI",
                    "DRN_custom" = "DRN",
                    "HB" = "Hb",
                    "LC_Pauli" = "LC",
                    "MDB" = "MDB",
                    "PP" = "PPN",
                    "BF" = "MS"),
    timepoint = rescale(timepoint, to = c(-2, 8))
  )

###################### Figure 4 #######################
# Fig. 4C
filtered_data <- d %>%
  filter(glm %in% c('GLM_04.1', 'GLM_04.1b')
         & contrast %in% c('contr_1')
         & region %in% c('DRN'))
ggplot(filtered_data %>% filter(glm=='GLM_04.1'), aes(
  x = timepoint, y = m, ymin = m - se, ymax = m + se, fill = glm, group = glm)) +
  geom_ribbon(aes(fill = glm, group = glm), alpha = 0.65, fill = 'deeppink4') +
  geom_line(aes(group = glm), color = 'black', linewidth = 0.60) +
  geom_hline(yintercept = 0, linetype = 1, alpha = 0.25) +
  geom_vline(xintercept = 0, linetype = 1, alpha = 0.25) +
  labs(
    x = 'Time [s]',
    y = TeX("$\\beta$-coeff. [a.u.]")
  ) +
  scale_y_continuous(breaks = c(-0.05, 0,00, 0.05, 0.10), limits = c(-0.05, 0.06)) +
  scale_x_continuous(limits = c(-2,8), breaks = c(0, 4, 8)) + 
  theme_pubr() +
  facet_grid(~factor(region, labels = c('Policy-change [DRN]'))) + 
  theme(
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    strip.text.x = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.position = 'none',
    aspect.ratio = 1,
    strip.text.y = element_blank(),
    strip.background = element_rect(colour='black', fill='grey99'),
    panel.spacing.x = unit(0.2, "lines")
  )

ggplot(filtered_data %>% filter(glm=='GLM_04.1b'), aes(
  x = timepoint, y = m, ymin = m - se, ymax = m + se, fill = glm, group = glm)) +
  geom_ribbon(aes(fill = glm, group = glm), alpha = 0.65, fill = 'azure2') +
  geom_line(aes(group = glm), color = 'black', linewidth = 0.60) +
  geom_hline(yintercept = 0, linetype = 1, alpha = 0.25) +
  geom_vline(xintercept = 0, linetype = 1, alpha = 0.25) +
  labs(
    x = 'Time [s]',
    y = TeX("$\\beta$-coeff. [a.u.]")
  ) +
  scale_y_continuous(breaks = c(-0.05, 0,00, 0.05, 0.10), limits = c(-0.05, 0.06)) +
  scale_x_continuous(limits = c(-2,8), breaks = c(0, 4, 8)) + 
  theme_pubr() +
  facet_grid(~factor(region, labels = c('Action-change [DRN]'))) + 
  theme(
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    strip.text.x = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.position = 'none',
    aspect.ratio = 1,
    strip.text.y = element_blank(),
    strip.background = element_rect(colour='black', fill='grey99'),
    panel.spacing.x = unit(0.2, "lines")
  )

###################### Figure 5 #######################
filtered_data <- d %>%
  filter(glm=='GLM_04.3'
         & contrast %in% c('contr_1')
         & region == 'DRN') %>%
  mutate(region_label = 'DRN-dACC')
ggplot(filtered_data, aes(
  x = timepoint, y = m, ymin = m - se, ymax = m + se)) +
  geom_hline(yintercept = 0, linetype = 1, alpha = 0.1) +
  geom_vline(xintercept = 0, linetype = 1, alpha = 0.1) +
  geom_ribbon(aes(fill = contrast, group = contrast), alpha = 0.60, fill = 'skyblue1') +
  geom_line(aes(group = contrast), color = 'black', linewidth = 0.60) +
  labs(
    x = 'Time [s]',
    y = TeX("$\\beta$-coeff. [a.u.]")
  ) +
  scale_y_continuous(breaks = seq(-0.06, 0.03, 0.03), position = 'right') +
  scale_x_continuous(limits = c(-2,8), breaks = c(0, 4, 8)) + 
  theme_pubr() +
  facet_grid(~region_label) + 
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
    legend.position = 'none'
  ) + 
  border(color = 'black')

###################### Supplementary Figure 3 ################################

filtered_data <- d %>%
  filter(glm=='GLM_0S3a'
         & contrast %in% c('contr_1')
         & region == 'MBD') %>%
  mutate(region_label = 'MBD')
ggplot(filtered_data, aes(
  x = timepoint, y = m, ymin = m - se, ymax = m + se)) +
  geom_hline(yintercept = 0, linetype = 1, alpha = 0.1) +
  geom_vline(xintercept = 0, linetype = 1, alpha = 0.1) +
  geom_ribbon(aes(fill = contrast, group = contrast), alpha = 0.75, fill = 'darkseagreen1') +
  geom_line(aes(group = contrast), color = 'black', linewidth = 0.60) +
  labs(
    x = 'Time [s]',
    y = TeX("$\\beta$-coeff. [a.u.]")
  ) +
  #scale_y_continuous(breaks = seq(-0.06, 0.03, 0.03), position = 'right') +
  scale_x_continuous(limits = c(-2,8), breaks = c(0, 4, 8)) + 
  theme_pubr() +
  facet_grid(~region) + 
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
    legend.position = 'none'
  ) + 
  border(color = 'black')

filtered_data <- d %>%
  filter(glm %in% c('GLM_0S3b')
         & contrast %in% c('contr_1')
         & region == 'MBD') %>%
  mutate(region_label = 'MBD')
ggplot(filtered_data, aes(
  x = timepoint, y = m, ymin = m - se, ymax = m + se, group = glm, fill = glm)) +
  geom_hline(yintercept = 0, linetype = 1, alpha = 0.1) +
  geom_vline(xintercept = 0, linetype = 1, alpha = 0.1) +
  geom_ribbon(aes(fill = glm, group = glm), alpha = 0.60, fill = 'skyblue2') +
  geom_line(aes(group = glm), color = 'black', linewidth = 0.60) +
  labs(
    x = 'Time [s]',
    y = TeX("$\\beta$-coeff. [a.u.]")
  ) +
  scale_y_continuous(breaks = seq(-0.03, 0.06, 0.03)) +
  scale_x_continuous(limits = c(-2,8), breaks = c(0, 4, 8)) + 
  theme_pubr() +
  facet_grid(~region) + 
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
    legend.position = 'none'
  ) + 
  border(color = 'black')

####################### Supplementary figure 6 ################################

filtered_data <- d %>%
  filter(glm=='GLM_0S6'
         & contrast %in% c('contr_1')
         & region %in% c('MS'))
ggplot(filtered_data, aes(
  x = timepoint, y = m, ymin = m - se, ymax = m + se)) +
  geom_hline(yintercept = 0, linetype = 1, alpha = 0.1) +
  geom_vline(xintercept = 0, linetype = 1, alpha = 0.1) +
  geom_ribbon(fill = 'steelblue2', alpha = 0.60) +
  geom_line(aes(group = contrast), color = 'black', linewidth = 0.60) +
  labs(
    x = 'Time [s]',
    y = TeX("$\\beta$-coeff. [a.u.]")
  ) +
  scale_y_continuous(breaks = c(-0.05, 0,00, 0.05, 0.10)) +
  scale_x_continuous(limits = c(-2,8), breaks = c(0, 4, 8)) + 
  theme_pubr() +
  facet_grid(~region) + 
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
    legend.position = 'none'
  )

####################### PEAK GRAPHS #######################

# List CSV files matching the pattern
list_files <- dir(pattern = '*peaks.csv')

# Read and process CSV files using lapply
d_peak_list <- lapply(list_files, function(file) {
  df <- read.csv(file)
  return(df)
})

# Combine data frames into one using dplyr
d_peak <- bind_rows(d_peak_list) %>%
  pivot_longer(cols = c(1, 5:10), names_to = 'region', values_to = 'peak') %>%
  mutate(region = recode(region,
                         DRN_custom = 'DRN',
                         LC_Pauli = 'LC',
                         BF = 'MS',
                         AIsphere7 = 'AI', 
                         dACCsphere7 = 'dACC'
  ))

######################## Figure 3 #########################

# Figure 3B
filtered_data <- d_peak %>%
  filter(
    region %in% c('DRN', 'MBD', 'LC', 'MS', 'HB'),
    contrast %in% c('contr1'),
    GLM %in% c('GLM_04.1', 'GLM_04.1b')
  ) %>%
  mutate(region = factor(region, levels = c('DRN', 'MBD', 'MS', 'LC', 'HB')))

ggplot(data = filtered_data, aes(x = region, y = peak, fill = GLM, group = GLM)) + 
  stat_summary(
    fun = 'mean',
    geom = 'bar', 
    color = 'black', 
    width = 0.5, 
    position = position_dodge(width = 0.8, preserve = 'total')
  ) + 
  geom_point(
    shape = 21, 
    color = 'grey50',
    alpha = 0.25,
    position = position_dodge(width = 0.8, preserve = 'total'),
    size = 2
  ) + 
  geom_hline(yintercept = 0, linetype = 2, alpha = 1, color = 'red') + 
  scale_y_continuous(
    #breaks = seq(-0.2, 0.2, 0.2),
    labels = scales::number_format(accuracy = 0.01)
  ) + 
  scale_x_discrete(
    #name = 'Predictor',
    #labels = c('Policy', 'Action')
  ) + 
  scale_fill_manual(
    name = NULL,
    values = c('deeppink4', 'azure1'),
    labels = NULL
  ) + 
  labs(
    x = 'ROI', 
    y = TeX("$\\beta$-coeff [a.u.]")
  ) + 
  theme_pubr() + 
  theme(
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.position = 'none', 
    legend.background = element_blank(),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.direction = 'none',
    aspect.ratio = 0.3,
    strip.text.x = element_text(size = 18),
    #strip.background = element_rect(colour='black', fill='grey95'),
    strip.background = element_rect(colour='black', fill='grey99'),
    panel.spacing.x = unit(0.2, "lines")
  ) + 
  #facet_grid(~region) + 
  border(color = 'black')

# Figure 3D-i
filtered_data <- d_peak %>%
  filter(
    region=='DRN',
    contrast == 'contr1',
    GLM %in% c('GLM_04.2a', 'GLM_04.2c')
  )

ggplot(data = filtered_data, aes(x = GLM, y = peak, fill = GLM, group = GLM)) + 
  geom_hline(yintercept = 0, linetype = 2, alpha = 1, color = 'red') + 
  stat_summary(
    fun = 'mean',
    geom = 'bar', 
    color = 'black', 
    width = 0.5, 
    alpha = 0.75,
    position = position_dodge2(width = 0)
  ) +   
  geom_point(alpha = 0.50, shape = 21, color = 'grey90', size = 2.5) + 
  scale_y_continuous(name = TeX("$\\beta$-coeff [a.u.]"), labels = scales::number_format(accuracy = 0.01), limits = c(-0.30, 0.50), breaks = seq(-0.25, 0.50, 0.25)) + 
  scale_x_discrete(labels = NULL) + 
  scale_fill_manual(values = c('deeppink4', 'pink')) + 
  theme_pubr() + 
  theme(
      axis.title.y = element_text(size = 18),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text = element_text(size = 14),
      legend.position = 'none', 
      legend.background = element_blank(),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 14),
      legend.direction = 'none',
      aspect.ratio = 2.5,
      strip.text.x = element_text(size = 18),
      #strip.background = element_rect(colour='black', fill='grey95'),
      strip.background = element_rect(colour='black', fill='grey99'),
      panel.spacing.x = unit(0.2, "lines")
  ) + 
  border(color = 'black')

# Figure 3D-ii
filtered_data <- d_peak %>%
  filter(
    region=='DRN',
    contrast == 'contr1',
    GLM %in% c('GLM_04.2a', 'GLM_04.2b')
  )

ggplot(data = filtered_data, aes(x = GLM, y = peak, fill = GLM, group = GLM)) + 
  geom_hline(yintercept = 0, linetype = 2, alpha = 1, color = 'red') + 
  stat_summary(
    fun = 'mean',
    geom = 'bar', 
    color = 'black', 
    width = 0.5, 
    alpha = 0.75,
    position = position_dodge2(width = 0)
  ) +   
  geom_point(alpha = 0.50, shape = 21, color = 'grey90', size = 2.5) + 
  scale_y_continuous(name = TeX("$\\beta$-coeff [a.u.]"), labels = scales::number_format(accuracy = 0.01), limits = c(-0.30, 0.50), breaks = seq(-0.25, 0.50, 0.25)) + 
  scale_x_discrete(labels = NULL) + 
  scale_fill_manual(values = c('deeppink4', 'azure2')) + 
  theme_pubr() + 
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text = element_text(size = 14),
    legend.position = 'none', 
    legend.background = element_blank(),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.direction = 'none',
    aspect.ratio = 2.5,
    strip.text.x = element_text(size = 18),
    #strip.background = element_rect(colour='black', fill='grey95'),
    strip.background = element_rect(colour='black', fill='grey99'),
    panel.spacing.x = unit(0.2, "lines")
  ) + 
  border(color = 'black')

# Figure 3D-iii

filtered_data <- d_peak %>%
  filter(
    region=='DRN',
    contrast == 'contr1',
    GLM %in% c('GLM_04.2c', 'GLM_04.2d')
  )

ggplot(data = filtered_data, aes(x = GLM, y = peak, fill = GLM, group = GLM)) + 
  geom_hline(yintercept = 0, linetype = 2, alpha = 1, color = 'red') + 
  stat_summary(
    fun = 'mean',
    geom = 'bar', 
    color = 'black', 
    width = 0.5, 
    alpha = 0.75,
    position = position_dodge2(width = 0)
  ) +   
  geom_point(alpha = 0.50, shape = 21, color = 'grey90', size = 2.5) + 
  scale_y_continuous(name = TeX("$\\beta$-coeff [a.u.]"), labels = scales::number_format(accuracy = 0.01), limits = c(-0.30, 0.50), breaks = seq(-0.25, 0.50, 0.25)) + 
  scale_x_discrete(labels = NULL) + 
  scale_fill_manual(values = c('pink', 'azure2')) + 
  theme_pubr() + 
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text = element_text(size = 14),
    legend.position = 'none', 
    legend.background = element_blank(),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.direction = 'none',
    aspect.ratio = 2.5,
    strip.text.x = element_text(size = 18),
    #strip.background = element_rect(colour='black', fill='grey95'),
    strip.background = element_rect(colour='black', fill='grey99'),
    panel.spacing.x = unit(0.2, "lines")
  ) + 
  border(color = 'black')

# Figure 3D-iv

filtered_data <- d_peak %>%
  filter(
    region=='DRN',
    contrast == 'contr1',
    GLM %in% c('GLM_04.2a', 'GLM_04.2a_low')
  )

ggplot(data = filtered_data, aes(x = GLM, y = peak, fill = GLM, group = GLM)) + 
  geom_hline(yintercept = 0, linetype = 2, alpha = 1, color = 'red') + 
  stat_summary(
    fun = 'mean',
    geom = 'bar', 
    color = 'black', 
    width = 0.5, 
    alpha = 0.75,
    position = position_dodge2(width = 0)
  ) +   
  geom_point(alpha = 0.50, shape = 21, color = 'grey90', size = 2.5) + 
  scale_y_continuous(name = TeX("$\\beta$-coeff [a.u.]"), labels = scales::number_format(accuracy = 0.01), limits = c(-0.30, 0.50), breaks = seq(-0.25, 0.50, 0.25)) + 
  scale_x_discrete(labels = NULL) + 
  scale_fill_manual(values = c('deeppink4', 'azure2')) + 
  theme_pubr() + 
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text = element_text(size = 14),
    legend.position = 'none', 
    legend.background = element_blank(),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.direction = 'none',
    aspect.ratio = 2.5,
    strip.text.x = element_text(size = 18),
    #strip.background = element_rect(colour='black', fill='grey95'),
    strip.background = element_rect(colour='black', fill='grey99'),
    panel.spacing.x = unit(0.2, "lines")
  ) + 
  border(color = 'black')

# Fig 3D-v
filtered_data <- d_peak %>%
  filter(
    region=='DRN',
    contrast == 'contr1',
    GLM %in% c('GLM_04.2c', 'GLM_04.2c_low')
  )

ggplot(data = filtered_data, aes(x = GLM, y = peak, fill = GLM, group = GLM)) + 
  geom_hline(yintercept = 0, linetype = 2, alpha = 1, color = 'red') + 
  stat_summary(
    fun = 'mean',
    geom = 'bar', 
    color = 'black', 
    width = 0.5, 
    alpha = 0.75,
    position = position_dodge2(width = 0)
  ) +   
  geom_point(alpha = 0.50, shape = 21, color = 'grey90', size = 2.5) + 
  scale_y_continuous(name = TeX("$\\beta$-coeff [a.u.]"), labels = scales::number_format(accuracy = 0.01), limits = c(-0.30, 0.50), breaks = seq(-0.25, 0.50, 0.25)) + 
  scale_x_discrete(labels = NULL) + 
  scale_fill_manual(values = c('pink1', 'azure2')) + 
  theme_pubr() + 
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text = element_text(size = 14),
    legend.position = 'none', 
    legend.background = element_blank(),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.direction = 'none',
    aspect.ratio = 2.5,
    strip.text.x = element_text(size = 18),
    #strip.background = element_rect(colour='black', fill='grey95'),
    strip.background = element_rect(colour='black', fill='grey99'),
    panel.spacing.x = unit(0.2, "lines")
  ) + 
  border(color = 'black')

######################## Figure 4 ###################
# Fig 4E
filtered_data <- d_peak %>%
  filter(
    region %in% c('DRN', 'AI', 'dACC'),
    contrast == 'contr1',
    GLM %in% c('GLM_04.2a', 'GLM_04.2c')
  )

ggplot(data = filtered_data, aes(x = GLM, y = peak, fill = GLM, group = GLM)) + 
  geom_hline(yintercept = 0, linetype = 2, alpha = 1, color = 'red') + 
  stat_summary(
    fun = 'mean',
    geom = 'bar', 
    color = 'black', 
    width = 0.5, 
    alpha = 0.75,
    position = position_dodge2(width = 0)
  ) +   
  geom_point(alpha = 0.50, shape = 21, color = 'grey90', size = 2.5) + 
  scale_y_continuous(name = TeX("$\\beta$-coeff [a.u.]"), labels = scales::number_format(accuracy = 0.01), limits = c(-0.30, 0.60), breaks = seq(-0.25, 0.50, 0.25)) + 
  scale_x_discrete(labels = NULL) + 
  scale_fill_manual(values = c('deeppink4', 'pink')) + 
  theme_pubr() + 
  theme(
    axis.title.y = element_text(size = 18),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text = element_text(size = 14),
    legend.position = 'none', 
    legend.background = element_blank(),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.direction = 'none',
    aspect.ratio = 1.5,
    strip.text.x = element_text(size = 18),
    #strip.background = element_rect(colour='black', fill='grey95'),
    strip.background = element_rect(colour='black', fill='grey99'),
    panel.spacing.x = unit(0.2, "lines")
  ) + 
  facet_grid(~region) + 
  border(color = 'black')

########################### SUPPLEMENRARY FIGURES #######################

# Supplementary figure S3A
filtered_data <- d_peak %>%
  filter(
    region=='MBD',
    contrast == 'contr1',
    GLM=='GLM_0S3a'
  )

ggplot(data = filtered_data, aes(x = GLM, y = peak, fill = GLM, group = GLM)) + 
  geom_hline(yintercept = 0, linetype = 2, alpha = 1, color = 'red') + 
  stat_summary(
    fun = 'mean',
    geom = 'bar', 
    color = 'black', 
    width = 0.5, 
    alpha = 0.75,
    position = position_dodge2(width = 0),
    fill = 'darkseagreen1'
  ) +   
  geom_point(alpha = 0.50, shape = 21, color = 'black', size = 2.5, fill = 'darkseagreen2') + 
  scale_y_continuous(name = TeX("$\\beta$-coeff [a.u.]"), position = 'right', labels = scales::number_format(accuracy = 0.01)) + 
  scale_x_discrete(labels = NULL) + 
  labs(x = NULL, y = NULL) + 
  theme_pubr() + 
  theme(
    axis.text = element_text(size = 22),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = 'none',
    legend.background = element_blank(),
    plot.background = element_blank(),
    aspect.ratio = 2.5
  ) + 
  border(size = 1.5)

# Supplementary figure S3C

filtered_data <- d_peak %>%
  filter(
    region=='MBD',
    contrast == 'contr1',
    GLM=='GLM_0S3b'
  )

ggplot(data = filtered_data, aes(x = GLM, y = peak, fill = GLM, group = GLM)) + 
  geom_hline(yintercept = 0, linetype = 2, alpha = 1, color = 'red') + 
  stat_summary(
    fun = 'mean',
    geom = 'bar', 
    color = 'black', 
    width = 0.5, 
    alpha = 0.75,
    position = position_dodge2(width = 0),
    fill = 'skyblue1'
  ) +   
  geom_point(alpha = 0.50, shape = 21, color = 'black', size = 2.5, fill = 'skyblue2') + 
  scale_y_continuous(name = TeX("$\\beta$-coeff [a.u.]"), position = 'right', labels = scales::number_format(accuracy = 0.01), breaks = seq(-0.05, 0.15, 0.05)) + 
  scale_x_discrete(labels = NULL) + 
  labs(x = NULL, y = NULL) + 
  theme_pubr() + 
  theme(
    axis.text = element_text(size = 22),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = 'none',
    legend.background = element_blank(),
    plot.background = element_blank(),
    aspect.ratio = 2.5
  ) + 
  border(size = 1.5)

# Supplementary figure S6A
filtered_data <- d_peak %>%
  filter(
    region %in% c('DRN', 'MBD', 'MS', 'HB', 'LC'),
    contrast == 'contr1',
    GLM=='GLM_0S6'
  )

ggplot(data = filtered_data, aes(x = region, y = peak)) + 
  geom_hline(yintercept = 0, linetype = 2, alpha = 1, color = 'red') + 
  stat_summary(
    fun = 'mean',
    geom = 'bar', 
    color = 'black', 
    fill = 'steelblue1',
    width = 0.5, 
    alpha = 0.75,
    position = position_dodge2(width = 0)
  ) +   
  geom_point(alpha = 0.50, shape = 21, color = 'steelblue1', size = 2.5) + 
  scale_y_continuous(name = TeX("$\\beta$-coeff [a.u.]"), labels = scales::number_format(accuracy = 0.01)) + 
  scale_x_discrete(name = 'ROI') + 
  theme_pubr() + 
  theme(
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.position = 'none', 
    legend.background = element_blank(),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.direction = 'none',
    aspect.ratio = 0.75,
    strip.text.x = element_text(size = 18),
    #strip.background = element_rect(colour='black', fill='grey95'),
    strip.background = element_rect(colour='black', fill='grey99'),
    panel.spacing.x = unit(0.2, "lines")
  ) + 
  border(color = 'black')

# Supplementary figure S6

filtered_data <- d_peak %>%
  filter(
    region %in% c('DRN', 'MS'),
    contrast == 'contr1',
    GLM %in% c('GLM_04.1', 'GLM_0S6')
  )

ggplot(data = filtered_data, aes(x = region, y = peak, fill = GLM, color = GLM)) + 
  geom_hline(yintercept = 0, linetype = 2, alpha = 1, color = 'red') + 
  stat_summary(
    fun = 'mean',
    geom = 'bar', 
    color = 'black', 
    width = 0.5, 
    alpha = 0.75,
    position = position_dodge2(width = 0)
  ) +   
  geom_point(alpha = 0.50, shape = 21, fill = NA, size = 2.5) + 
  scale_y_continuous(name = TeX("$\\beta$-coeff [a.u.]"), labels = scales::number_format(accuracy = 0.01)) + 
  scale_x_discrete(name = 'ROI') + 
  scale_fill_manual(name = 'Contrast', values = c('steelblue1', 'purple4'), labels = c('Incongruent', 'Congruent')) + 
  scale_color_manual(name = 'Contrast', values = c('steelblue1', 'purple4'), labels = c('Incongruent', 'Congruent')) + 
  theme_pubr() + 
  theme(
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.position = 'none', 
    legend.background = element_blank(),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.direction = 'none',
    aspect.ratio = 0.75,
    strip.text.x = element_text(size = 18),
    #strip.background = element_rect(colour='black', fill='grey95'),
    strip.background = element_rect(colour='black', fill='grey99'),
    panel.spacing.x = unit(0.2, "lines")
  ) + 
  border(color = 'black') + 
  facet_grid(~factor(GLM, labels = c('Exploratory', 'Env.-driven')))

############################## STATISTICAL TESTS ###############################

# compare policy-switch and action-switch signals
policy_switch_drn <- d_peak %>% filter(GLM =='GLM_04.1', region == 'DRN', contrast== 'contr1', timelock == 'offer-timelocked') %>% select(peak)
action_switch_drn <- d_peak %>% filter(GLM =='GLM_04.1b', region == 'DRN', contrast== 'contr1', timelock == 'offer-timelocked') %>% select(peak)
t.test(policy_switch_drn, action_switch_drn)

# compare policy switch signals as a function of (i) direction and (ii) congruency
congruent_poor_drn <-  d_peak %>% filter(GLM =='GLM_04.2a', region == 'DRN', contrast== 'contr1', timelock == 'offer-timelocked') %>% select(peak)
incongruent_rich_drn <-  d_peak %>% filter(GLM =='GLM_04.2b', region == 'DRN', contrast== 'contr1', timelock == 'offer-timelocked') %>% select(peak)
congruent_rich_drn <-  d_peak %>% filter(GLM =='GLM_04.2c', region == 'DRN', contrast== 'contr1', timelock == 'offer-timelocked') %>% select(peak)
incongruent_poor_drn <-  d_peak %>% filter(GLM =='GLM_04.2d', region == 'DRN', contrast== 'contr1', timelock == 'offer-timelocked') %>% select(peak)

t.test(incongruent_rich_drn)
t.test(incongruent_poor_drn)

t.test(congruent_poor_drn, congruent_rich_drn) # are congruent changes different with respect to environment? 
t.test(congruent_poor_drn$peak, incongruent_rich_drn$peak, paired = T) # are pursue-related switches different according to congruence (i.e. as a function of environment)? 
t.test(congruent_rich_drn$peak, incongruent_poor_drn$peak, paired = T) # are pursue-related switches different according to congruence (i.e. as a function of environment)? 

# compare congruent policy switches for (i) low-value, and (ii) high-value options
congruent_poor_low_drn <- d_peak %>% filter(GLM =='GLM_04.2a_low', region == 'DRN', contrast== 'contr1', timelock == 'offer-timelocked') %>% select(peak)
congruent_rich_low_drn <- d_peak %>% filter(GLM =='GLM_04.2c_low', region == 'DRN', contrast== 'contr1', timelock == 'offer-timelocked') %>% select(peak)

t.test(congruent_poor_low_drn)
t.test(congruent_rich_low_drn)

t.test(congruent_poor_drn$peak, congruent_poor_low_drn$peak, paired = T)
t.test(congruent_rich_drn$peak, congruent_rich_low_drn$peak, paired = T)

# ANOVA comparing DRN and MS effects as a function of switch type
summary(aov(peak ~ region*GLM, data = d_peak %>% filter(region %in% c('DRN', 'MS') & contrast=='contr1' & GLM %in% c('GLM_04.1', 'GLM_0S6'))))

rm(list = ls())
filepath <- "~/nat_behav_submission/" # filepath to GitHub folder
setwd(paste(filepath, 'data/behavioural/', sep = ''))

# load packages
library(lme4)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(gghalves)
library(ggdist)
library(ggpubr)
library(latex2exp)
library(tidyverse)

################## Load data ####################
# get data
d <- read.csv('behavioural_data.csv', header = T)
###########################################################################
#§2: Reward guided decisions are modulated by the richness of the environment

GLM1.1a <- glmer(response ~ mag_z*env_bin + trial_z + (mag_z:env_bin|subject), data = d, family = 'binomial', control = glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
summary(GLM1.1a)

GLM1.1b_mid <- glmer(response ~ env_bin + trial_z + (env_bin + trial_z|subject), data = d %>% filter(mag==10), family = 'binomial', control = glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
summary(GLM1.1b_mid)

GLM1.1b_low <- glmer(response ~ env_bin + trial_z + (env_bin + trial_z|subject), data = d %>% filter(mag==5), family = 'binomial', control = glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
summary(GLM1.1b_low )

GLM1.1c_high <- glmer(response ~ env_bin + trial_z + (1|subject), data = d %>% filter(mag==50), family = 'binomial', control = glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
summary(GLM1.1c_high)

GLM1.2_mid <- glmer(response ~ env_bin*trial_n_block_z + (1|subject), data = d %>% filter(mag==10), family = 'binomial', control = glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
summary(GLM1.2_mid)

GLM1.2_low <- glmer(response ~ env_bin*trial_n_block_z + (1|subject), data = d %>% filter(mag==5), family = 'binomial', control = glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
summary(GLM1.2_low)

GLM1.2_high <- glmer(response ~ env_bin*trial_n_block_z + (1|subject), data = d %>% filter(mag==50), family = 'binomial', control = glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
summary(GLM1.2_high)

GLM1.2_mid_poor <- glmer(response ~ trial_n_block_z + (1|subject), data = d %>% filter(mag==10 & env_bin == -1), family = 'binomial', control = glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
summary(GLM1.2_mid_poor)

GLM1.2_mid_rich <- glmer(response ~ trial_n_block_z + (1|subject), data = d %>% filter(mag==10 & env_bin == 1), family = 'binomial', control = glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
summary(GLM1.2_mid_rich)

GLM1.3_mid <- glmer(response ~ mu_val_z + trial_z + (mu_val_z|subject), data = d %>% filter(mag==10), family = 'binomial', control = glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
summary(GLM1.3_mid )

GLM_1.3_low <- glmer(response ~ mu_val_z + trial_z + (mu_val_z|subject), data = d %>% filter(mag==5), family = 'binomial', control = glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
summary(GLM_1.3_low)

GLM_1.3_high <- glmer(response ~ mu_val_z + trial_z + (mu_val_z|subject), data = d %>% filter(mag==50), family = 'binomial', control = glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
summary(GLM_1.3_high)

# relationship between acceptance rates and task performance
total_rewards <- d %>%
  group_by(subject) %>%
  summarise(total_reward = sum(reward, na.rm = TRUE)) %>%
  mutate(total_reward = scale(total_reward))
acc_rates <- d %>%
  group_by(subject, mag, env) %>%
  summarise(m = mean(response, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = env, values_from = m) %>%
  pivot_wider(names_from = mag, values_from = c(poor, rich)) %>%
  left_join(total_rewards, by = "subject") %>%  # Join total_reward per subject
  ungroup()
summary(lm(total_reward ~ poor_10*rich_10, data = acc_rates))

#Figure 1B: environment-specific frequencies
plot_data <- d %>%
  filter(subject==103) %>%
  select(env_bin, trial, mu_val)
plot_data$p_low <- ifelse(plot_data$env_bin==1, 1/6+0.001,
                          ifelse(plot_data$env_bin==-1, 3/6+0.001, NA))
plot_data$p_mod <- 2/6
plot_data$p_high <- ifelse(plot_data$env_bin==-1, 1/6-0.001,
                           ifelse(plot_data$env_bin==1, 3/6-0.001, NA))
p_bar <- plot_data %>%
  pivot_longer(cols = starts_with('p_'), names_to = 'option', values_to = 'p_opt')
ggplot(p_bar, aes(x = factor(option, levels = c('p_low', 'p_mod', 'p_high')), y = p_opt, fill = option)) +
  stat_summary(fun = 'mean', geom = 'bar', color = 'black', alpha = 0.80, width = 0.65) + 
  scale_x_discrete(name = 'Rw-magnitude', labels = c('5', '10', '50')) + 
  scale_y_continuous(name = 'p(Option)', limits = c(0, 0.55)) + 
  scale_fill_manual(name = 'Rw-magnitude', values = c('#DAA520', '#CD7F32', '#C0C0C0')) + 
  theme_pubr() + 
  coord_cartesian(ylim = c(0.15, 0.55)) + 
  theme(
    aspect.ratio = 1,
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 13),
    strip.text = element_text(size = 16),
    legend.position = 'none'
  ) + 
  facet_grid(~factor(env_bin, labels = c('Poor', 'Rich'))) + 
  border(color = 'black')

# Figure 1C: example environments
example_session <- d %>% filter(subject == 106)
session_change <- which(example_session$trial_n_block==1)
session_change <- c(session_change[2], session_change[4:11])
example_session$env <- ifelse(example_session$env_bin==1, 'rich', 'poor')
example_session$env <- factor(example_session$env, levels = c('poor', 'rich'), labels = c('poor', 'rich'))

example_session$env[session_change]

plot_data <- example_session %>% filter(block_count >9 & block_count <12)
plot_data$trial <- 1:nrow(plot_data)

ggplot() +
  geom_point(data = plot_data, aes(x = trial, y = mag, fill = factor(mag)), shape = 21, size = 4) + 
  #geom_line(data = plot_data, aes(x = trial, y = mag), linewidth = 0.50, alpha = 0.25) + 
  geom_vline(xintercept = 14, color = 'darkblue', linetype = 2, size = 0.5) + 
  geom_line(data = plot_data, aes(x = trial, y = mu_val), color = 'red', alpha = 0.75, linetype = 'longdash', linewidth = 0.5) + 
  scale_x_continuous(name = 'Time [trials]', breaks = seq(0, 40, 10)) + 
  scale_y_continuous(name = 'Ave. value', breaks = seq(0, 40, 20)) + 
  scale_fill_manual(name = 'Reward-magnitude', values = c('#CD7F32', '#C0C0C0', '#DAA520')) + 
  theme_pubr() + 
  theme(    axis.title.y = element_text(size = 18),
            axis.title.x = element_text(size = 18),
            axis.text = element_text(size = 14),
            legend.text = element_text(size = 12),
            axis.text.y = element_text(size = 14),
            legend.position = 'none',
            aspect.ratio = 0.50) + 
  coord_cartesian(ylim = c(-10,60))

# Figure 1E: pursue-rate as a function of rw-magnitude and environment-type
plot_data <- d %>%
  group_by(subject, mag, env) %>%
  summarise(m = mean(response))
ggplot(plot_data, aes(x = factor(mag), fill = factor(env), color = factor(env), y = m*100)) + 
  stat_summary(fun = 'mean', geom = 'bar', color = 'black', width = 0.85, alpha = 1, position = position_dodge(width = 1)) + 
  geom_point(size = 2.5, colour = 'grey50', shape = 21, alpha = 0.75, position = position_dodge(width = 1)) +  
  scale_y_continuous(name = "Acceptance-rate [%]") + 
  scale_x_discrete(name = 'Rw-mag [points]', labels = c('5', '10', '50')) + 
  scale_fill_manual(name = 'Rw-mag [points]', values = c('lightskyblue1', 'darkblue'), labels = c('5', '10', '50')) + 
  scale_color_manual(name = 'Rw-mag [points]', values = c('lightskyblue1', 'darkblue'), labels = c('5', '10', '50')) + 
  theme_pubr() + 
  #coord_cartesian(ylim = c(-40, 40)) + 
  theme(
    text = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    legend.box.background = element_rect(colour = "black", fill = 'black'),
    legend.background = element_rect(fill = 'grey99'),
    aspect.ratio = 1,
    legend.position = 'none', 
    legend.direction="horizontal")

# Figure 1F: pursue-rate as a function of rw-magnitude, environment-type and time
d$mag_factor <- paste(d$mag, '-points', sep = '')
plot_data <- d %>% 
  filter(mag==10) %>% 
  group_by(env, trial_n_block, mag_factor) %>% 
  summarise(m = mean(response, na.rm = T)*100)
ggplot(plot_data, aes(x = trial_n_block, y = m, group = as.factor(env), fill = as.factor(env))) + 
  stat_summary(fun = 'mean', geom = 'point', size = 3, shape = 21) + 
  geom_smooth(method = 'lm', alpha = 0.35, color = 'black', linewidth = 0.5) + 
  scale_y_continuous(name = 'Accept-rate [%]', breaks = seq(0, 100, 20), limits = c(40,80)) + 
  scale_x_continuous(name = 'Time-in-env. [Trials]', limits = c(0, 15)) + 
  theme_pubr() + 
  scale_fill_manual(name = 'environment', values= c('skyblue', 'darkblue'), labels = c('Poor', 'Rich')) + 
  theme(
    text = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    legend.background = element_rect(fill = 'grey99'),
    aspect.ratio = 1,
    legend.position = 'none', 
    legend.direction="horizontal") + 
  facet_grid(.~mag_factor)

# Fig 1G: relationship between behavioural change and task performance
acc_rates$env_diff <- acc_rates$poor_10 - acc_rates$rich_10

ggplot(acc_rates, aes(x = env_diff*100, y = total_reward/max(total_reward))) +
  geom_vline(xintercept = 0, linetype = 2, color = 'blue') +
  geom_point(shape = 21, color = 'black', fill = 'azure1', size = 3) +
  geom_smooth(method = 'lm', fill = 'azure2', color = 'black', size = 0.5) +
  scale_x_continuous(name = TeX("$\\Delta$%acc. [Poor - Rich]")) +
  scale_y_continuous(name = 'Total points', breaks = c(0.8, 0.9, 1.0), limits = c(0.8, 1.0)) +
  theme_pubr() +
  theme(
    aspect.ratio = 1.3,
    text = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.title = element_blank(),
    legend.text = element_text(size = 14)
  ) +
  border(color = 'black')

###########################################################################
# §3: Option specific behavioural policies are reconciled with the richness of the environment
d$policy_switch <- ifelse(d$response==d$prev_policy, 0,
                          ifelse(d$response != d$prev_policy, 1, NA))
d$pursue_switch <- ifelse(d$response==1 & d$prev_policy==0, 1, 0)
d$reject_switch <- ifelse(d$response==0 & d$prev_policy==1, 1, 0)
d$switch_type <- ifelse(d$pursue_switch==1, 'pursue-switch', 
                        ifelse(d$reject_switch==1, 'reject-switch', 'none'))
d$switch_type_bin <- ifelse(d$pursue_switch==1, 1,
                            ifelse(d$reject_switch==1, -1, NA))

# Is there evidence of option-specific behavioural policies? 
GLM2.1 <-  glmer(response ~ mag_z + prev_policy + prev_action + trial_z + (mag_z + prev_policy + prev_action|subject), data = d, family = 'binomial', control = glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
summary(GLM2.1)

GLM2.1a <- glmer(response ~ mag_z + prev_policy + trial_z + (mag_z + prev_policy|subject), data = d, family = 'binomial', control = glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
summary(GLM2.1a)

GLM2.1b <- glmer(response ~ mag_z + prev_action + trial_z + (prev_action |subject), data = d, family = 'binomial', control = glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
summary(GLM2.1b)

# Do policy switches vary in frequency for different options? 
GLM2.3a <- glmer(policy_switch ~ factor(mag) + trial_n_block_z + (factor(mag) + trial_n_block_z|subject), data = d %>% filter(mag != 5), family = 'binomial', control = glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
summary(GLM2.3a)

GLM2.3b <- glmer(policy_switch ~ factor(mag) + trial_n_block_z + (factor(mag) + trial_n_block_z|subject), data = d %>% filter(mag != 50), family = 'binomial', control = glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
summary(GLM2.3b)

# Which type of policy switch is more likely to occur in poor (2.4a) and rich (2.4b) environments
GLM2.4a <- glmer(policy_switch ~ prev_policy + (1|subject), data = d %>% filter(mag==10 & env=='poor'), family = 'binomial', control = glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
summary(GLM2.4a)

GLM2.4b <- glmer(policy_switch ~ prev_policy + (1|subject), data = d %>% filter(mag==10 & env=='rich'), family = 'binomial', control = glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
summary(GLM2.4b)

# Which type of environment are pursue-switches (2.5a) and reject switches (2.5b) more likely to occur in? 
GLM2.5a <- glmer(pursue_switch ~ env_bin + (env_bin|subject), data = d %>% filter(mag==10 & prev_policy==0), family = 'binomial', control = glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
summary(GLM2.5a)

GLM2.5b <- glmer(reject_switch ~ env_bin + (env_bin|subject), data = d %>% filter(mag==10 & prev_policy==1), family = 'binomial', control = glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
summary(GLM2.5b)

# How does the likelihood of policy switching change according to average-value? 
GLM2.6a <- glmer(pursue_switch ~ mu_val_z + (mu_val_z|subject), data = d %>% filter(mag==10 & prev_policy==0), family = 'binomial', control = glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
summary(GLM2.5a)

GLM2.6b <- glmer(reject_switch ~ mu_val_z + (mu_val_z|subject), data = d %>% filter(mag==10 & prev_policy==1), family = 'binomial', control = glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
summary(GLM2.5b)

GLM2.6c <- glmer(pursue_switch ~ mu_val_z + trial_z + (1|subject), data = d %>% filter(mag==5 & prev_policy==0), family = 'binomial', control = glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
summary(GLM2.5c)

GLM2.6d <- glmer(reject_switch ~ mu_val_z + trial_z + (1|subject), data = d %>% filter(mag==5 & prev_policy==1), family = 'binomial', control = glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
summary(GLM2.5d)

GLM2.6e <- glm(pursue_switch ~ mu_val_z + trial_z , data = d %>% filter(mag==50 & prev_policy==0), family = 'binomial')
summary(GLM2.5e)

GLM2.6f <- glm(reject_switch ~ mu_val_z + trial_z , data = d %>% filter(mag==50 & prev_policy==1), family = 'binomial')
summary(GLM2.5f)

# Figure 3B: prev. policy and prev. action effects
policy_effect <- d %>% 
  drop_na(prev_policy) %>% 
  group_by(subject, prev_policy) %>%
  summarise(m = mean(response, na.rm = T)) %>%
  mutate(effect = 'Policy')

mean_pursue_rate <- mean(d$response)*100
ggplot(policy_effect, aes(x = factor(prev_policy), y = m*100, fill = factor(prev_policy))) + 
  geom_hline(yintercept = mean_pursue_rate, color = 'red', linetype = 2, alpha = 0.5) + 
  stat_summary(fun = 'mean', geom = 'bar', color = 'black', alpha = 0.75, position = 'dodge', width = 0.8) + 
  geom_point(size = 2.5, alpha = 0.50, shape = 21, color = 'black') + 
  scale_fill_manual(name = 'Prev. policy', values = c('grey95', 'deeppink4')) + 
  scale_x_discrete(name = 'Prev. policy', labels = c('Reject', 'Pursue')) + 
  scale_y_continuous(name = 'Pursue-rate [%]') + 
  theme_pubr() + 
  theme(
    aspect.ratio = 1.5, 
    text = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.position = 'none',
    legend.text = element_text(size = 14),
    legend.background = element_blank(), 
    strip.background = element_rect(colour='black', fill='grey99')
  ) + 
  border(color = 'black') + 
  coord_cartesian(ylim = c(0, 105)) + 
  facet_wrap(~effect)

action_effect <- d %>% 
  drop_na(prev_action) %>% 
  group_by(subject, prev_action) %>%
  summarise(m = mean(response, na.rm = T)) %>%
  mutate(effect = 'Action')

mean_pursue_rate <- mean(d$response)*100
ggplot(action_effect, aes(x = factor(prev_action), y = m*100, fill = factor(prev_action))) + 
  geom_hline(yintercept = mean_pursue_rate, color = 'red', linetype = 2, alpha = 0.5) + 
  stat_summary(fun = 'mean', geom = 'bar', color = 'black', alpha = 0.75, position = 'dodge', width = 0.8) + 
  geom_point(size = 2.5, alpha = 0.50, shape = 21, color = 'black') + 
  scale_fill_manual(name = 'Prev. action', values = c('grey95', 'skyblue3')) + 
  scale_x_discrete(name = 'Prev. action', labels = c('Reject', 'Pursue')) + 
  scale_y_continuous(name = 'Pursue-rate [%]') + 
  theme_pubr() + 
  theme(
    aspect.ratio = 1.5, 
    text = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.position = 'none',
    legend.text = element_text(size = 14),
    legend.background = element_blank(),
    strip.background = element_rect(colour='black', fill='grey99')
  ) + 
  border(color = 'black') + 
  coord_cartesian(ylim = c(0, 105)) + 
  facet_wrap(~effect)

# Figure 3C: number of policy-switches as a function of reward-magnitude
n_switch <- d %>% 
  drop_na(policy_switch) %>% 
  mutate(factor_policy_switch = factor(policy_switch)) %>%
  count(subject, mag, factor_policy_switch, .drop = FALSE) %>%
  filter(factor_policy_switch==1)

ggplot(n_switch, aes(x = factor(mag), y = n, fill = factor(mag))) + 
  stat_summary(fun = 'mean', geom = 'bar', color = 'black', alpha = 0.65) + 
  geom_point(position = position_jitter(height = 0.05, width = 0.05), size = 2.5, alpha = 0.50, shape = 21, color = 'black') + 
  scale_fill_manual(name = 'Reward-option', values = c('grey95', 'deeppink4', 'white')) + 
  scale_x_discrete(name = 'Reward-option') + 
  scale_y_continuous(name = 'Policy-switches [N]') + 
  theme_pubr() + 
  theme(
    aspect.ratio = 1, 
    text = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.position = c(0.4, 1.05),
    legend.direction = 'horizontal',
    #legend.position = 'none',
    legend.text = element_text(size = 14),
    legend.background = element_blank()
  )

# Figure 3C: frequency of switch-type as a function of environment-type

plot_data <- d %>%
  drop_na(switch_type) %>%
  filter(env == 'poor' & policy_switch==1) %>%
  group_by(subject, mag, switch_type) %>%
  summarise(n = n())

plot_data <- d %>%
  filter(env == 'poor') %>%
  group_by(subject, mag) %>%
  summarise(m_pursue = mean(pursue_switch, na.rm = T),
            m_reject = mean(reject_switch, na.rm = T)) %>%
  pivot_longer(cols = starts_with('m_'), names_to = 'switch_type', values_to = 'switch_rate')

ggplot(plot_data %>% filter(mag==10), aes(x = switch_type, y = switch_rate*100, fill = switch_type)) + 
  stat_summary(fun = 'mean', geom = 'bar', position = 'dodge', color = 'black', alpha = 0.50, width = 0.75) + 
  geom_point(shape = 21, size = 2.5) + 
  scale_fill_manual(name = 'Congruency', labels = c('Pursue', 'Reject'), values = c('deeppink4', 'grey95')) + 
  scale_x_discrete(name = 'Switch-type', labels = c('Pursue', 'Reject')) + 
  scale_y_continuous(name = 'Switch-rate[%]', limits = c(0, 40), breaks = seq(0, 30, 10)) + 
  theme_pubr() + 
  facet_wrap(~factor(mag, labels = c('Poor-env.'))) + 
  theme(
    aspect.ratio = 1.25, 
    text = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.position = 'none'
  )

plot_data <- d %>%
  filter(env == 'rich') %>%
  group_by(subject, mag) %>%
  summarise(m_pursue = mean(pursue_switch, na.rm = T),
            m_reject = mean(reject_switch, na.rm = T)) %>%
  pivot_longer(cols = starts_with('m_'), names_to = 'switch_type', values_to = 'switch_rate')

ggplot(plot_data %>% filter(mag==10), aes(x = switch_type, y = switch_rate*100, fill = switch_type)) + 
  stat_summary(fun = 'mean', geom = 'bar', position = 'dodge', color = 'black', alpha = 0.50, width = 0.75) + 
  geom_point(shape = 21, size = 2.5) + 
  scale_fill_manual(name = 'Congruency', labels = c('Pursue', 'Reject'), values = c('grey95', 'deeppink4')) + 
  scale_x_discrete(name = 'Switch-type', labels = c('Pursue', 'Reject')) + 
  scale_y_continuous(name = 'Switch-rate[%]', limits = c(0, 40), breaks = seq(0, 30, 10), position = 'right') + 
  theme_pubr() + 
  facet_wrap(~factor(mag, labels = c('Rich-env.'))) + 
  theme(
    aspect.ratio = 1.25, 
    text = element_text(size = 18),
    axis.text = element_text(size = 14),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    legend.position = 'none'
  )

# breakdown by env.
plot_data <- d %>%
  group_by(subject, env, mag) %>%
  summarise(m_pursue = mean(pursue_switch, na.rm = T),
            m_reject = mean(reject_switch, na.rm = T)) %>%
  pivot_longer(cols = starts_with('m_'), names_to = 'switch_type', values_to = 'switch_rate')

ggplot(plot_data %>% filter(mag==10 & switch_type=='m_pursue'), aes(x = env, y = switch_rate*100, fill = env)) + 
  stat_summary(fun = 'mean', geom = 'bar', position = 'dodge', color = 'black', alpha = 0.50, width = 0.75) + 
  geom_point(shape = 21, size = 2.5) + 
  scale_fill_manual(name = 'Env.', labels = c('Poor', 'Rich'), values = c('deeppink4', 'grey95')) + 
  scale_x_discrete(name = 'Environment', labels = c('Poor', 'Rich')) + 
  scale_y_continuous(name = 'Switch-rate[%]', limits = c(0, 40), breaks = seq(0, 30, 10)) + 
  theme_pubr() + 
  facet_wrap(~factor(mag, labels = c('Pursue-switch'))) + 
  theme(
    aspect.ratio = 1.25, 
    text = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.position = 'none'
  )

ggplot(plot_data %>% filter(mag==10 & switch_type=='m_reject'), aes(x = env, y = switch_rate*100, fill = env)) + 
  stat_summary(fun = 'mean', geom = 'bar', position = 'dodge', color = 'black', alpha = 0.50, width = 0.75) + 
  geom_point(shape = 21, size = 2.5) + 
  scale_fill_manual(name = 'Env.', labels = c('Poor', 'Rich'), values = c('deeppink4', 'grey95')) + 
  scale_x_discrete(name = 'Environment', labels = c('Poor', 'Rich')) + 
  scale_y_continuous(name = 'Switch-rate[%]', position = 'right', limits = c(0, 40), breaks = seq(0, 30, 10)) + 
  theme_pubr() + 
  facet_wrap(~factor(mag, labels = c('Reject-switch'))) + 
  theme(
    aspect.ratio = 1.25, 
    text = element_text(size = 18),
    axis.text = element_text(size = 14),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    legend.position = 'none'
  )

plot_data <- d %>%
  filter(env == 'rich') %>%
  group_by(subject, mag) %>%
  summarise(m_pursue = mean(pursue_switch, na.rm = T),
            m_reject = mean(reject_switch, na.rm = T)) %>%
  pivot_longer(cols = starts_with('m_'), names_to = 'switch_type', values_to = 'switch_rate')

# Fig. 3E: Switch rate as a function of average value
plot_data <- d %>%
  filter(mag==10) %>%
  drop_na(mu_val_z) %>%
  mutate(mu_val_ind = cut(mu_val_z, breaks = unique(quantile(mu_val_z, probs = seq(0, 1, length.out = 11))), labels = FALSE, include.lowest = TRUE)) %>%
  group_by(mu_val_ind) %>%
  mutate(mu_val_graph = mean(mu_val_z, na.rm = T)) %>%
  ungroup() %>%
  group_by(subject, mag, prev_policy, mu_val_graph) %>% 
  summarise(m = 100*mean(pursue_switch, na.rm = T))

ggplot(plot_data, aes(x = mu_val_graph, y = m)) + 
  geom_smooth(method = 'lm', fill = 'deeppink4', color = 'black', linewidth = 0.5) +
  stat_summary(fun.data = 'mean_se', geom = 'pointrange', fill = 'deeppink4', shape = 21, color = 'black', size = 0.65, linewidth = 0.5) + 
  theme_pubr() + 
  scale_y_continuous(name = 'Switch-rate [%]', breaks = seq(0, 100, 20)) + 
  scale_x_continuous(name = 'Ave. value [Z]') + 
  theme(
    aspect.ratio = 1.1, 
    text = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.position = 'none'
  ) + 
  coord_cartesian(ylim = c(0, 40)) + 
  facet_grid(~factor(mag, labels = 'Pursue-switch')) + 
  border('black')

plot_data <- d %>%
  filter(mag==10) %>%
  drop_na(mu_val_z) %>%
  mutate(mu_val_ind = cut(mu_val_z, breaks = unique(quantile(mu_val_z, probs = seq(0, 1, length.out = 11))), labels = FALSE, include.lowest = TRUE)) %>%
  group_by(mu_val_ind) %>%
  mutate(mu_val_graph = mean(mu_val_z, na.rm = T)) %>%
  ungroup() %>%
  group_by(subject, mag, prev_policy, mu_val_graph) %>% 
  summarise(m = 100*mean(reject_switch, na.rm = T))

ggplot(plot_data, aes(x = mu_val_graph, y = m)) + 
  geom_smooth(method = 'lm', fill = 'deeppink4', color = 'black', linewidth = 0.5) +
  stat_summary(fun.data = 'mean_se', geom = 'pointrange', fill = 'deeppink4', shape = 21, color = 'black', size = 0.65, linewidth = 0.5) + 
  theme_pubr() + 
  scale_y_continuous(name = 'Switch-rate [%]', position = 'right', breaks = seq(0, 40, 20)) + 
  scale_x_continuous(name = 'Ave. value [Z]') + 
  theme(
    aspect.ratio = 1.1, 
    text = element_text(size = 18),
    axis.text = element_text(size = 14),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    legend.position = 'none'
  ) + 
  coord_cartesian(ylim = c(0, 40)) + 
  facet_grid(~factor(mag, labels = 'Reject-switch')) + 
  border('black')

############################################################
# Supplementary Figures

# Supplementary figure S1: Additional details of task and behaviour

# fig S1B
max_trials <- d %>% 
  group_by(subject, env_bin, block_count) %>%
  summarise(max = max(trial_n_block, na.rm = T)) %>%
  filter(max > 1) %>%
  ungroup() %>%
  group_by(subject, env_bin) %>%
  summarise(m_trials = mean(max, na.rm = T))
dur_poor <- max_trials %>% filter(env_bin==-1)
dur_rich <- max_trials %>% filter(env_bin==1)
t.test(dur_poor$m_trials, dur_rich$m_trials, paired = T)

ggplot(max_trials, aes(x = factor(env_bin), fill = factor(env_bin), group = factor(env_bin), y = m_trials)) + 
  geom_hline(yintercept = mean(max_trials$m_trials), linetype = 2, color = 'red', alpha = 0.50) + 
  geom_point(shape = 21, alpha = 0.50) + 
  stat_summary(fun.data = 'mean_se', geom = 'pointrange', color = 'black', shape = 21) + 
  stat_summary(fun = 'mean', geom = 'line', group = 2, color = 'blue', linetype = 2) + 
  theme_pubr() + 
  scale_y_continuous(name = 'Max.trials [N]') + 
  scale_x_discrete(name = 'Environment-type', labels = c('Poor', 'Rich')) + 
  scale_fill_manual(name = 'Environment-type', values = c('skyblue2', 'darkblue')) + 
  theme(
    text = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    legend.background = element_rect(fill = 'grey99'),
    aspect.ratio = 1,
    legend.position = 'none', 
    legend.direction="horizontal")

# fig S1C
plot_data <- d %>%
  group_by(subject, mag, block_count, env) %>%
  summarise(m = mean(response, na.rm = T)) %>%
  group_by(subject, mag, env) %>%
  summarise(mu = mean(m, na.rm = T)) %>%
  pivot_wider(names_from = env, values_from = mu) %>%
  mutate(env_diff = poor - rich) %>%
  ungroup()
ggplot(plot_data, aes(x = factor(mag), fill = factor(mag), group = factor(mag), y = env_diff*100)) + 
  geom_hline(yintercept = 0, linetype = 2, color = 'red', alpha = 0.50) + 
  geom_point(shape = 21, alpha = 0.50) + 
  stat_summary(fun.data = 'mean_se', geom = 'pointrange', color = 'black', shape = 21) + 
  stat_summary(fun = 'mean', geom = 'line', group = 2, color = 'blue', linetype = 2) + 
  theme_pubr() + 
  scale_y_continuous(name = TeX("$\\Delta$%acc. [Poor - Rich]"), limits = c(-35, 60)) + 
  scale_x_discrete(name = 'Rw.-mag.[points]') + 
  scale_fill_manual(name = 'Rw-magnitude', values = c('azure1', 'skyblue2', 'darkblue')) + 
  theme(
    text = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    legend.background = element_rect(fill = 'grey99'),
    aspect.ratio = 1.5,
    legend.position = 'none', 
    legend.direction="horizontal")

# fig S1D
ggplot(plot_data, aes(x = mu_val_graph, y = m*100, fill = factor(mag))) + 
  stat_summary(fun.data = 'mean_se', geom = 'pointrange', size = 0.66, shape = 21, alpha = 0.65, color = 'black') + 
  geom_smooth(method = 'lm', color = 'black', linewidth = 0.5) + 
  theme_pubr() + 
  scale_fill_manual(name = 'Rw-magnitude', values = c('azure1', 'skyblue2', 'darkblue')) + 
  scale_x_continuous(name = 'Ave. value [Z]') + 
  scale_y_continuous(name = 'Pursue-rate [%]', breaks = seq(0, 100, 20)) + 
  theme(
    text = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    legend.background = element_rect(fill = 'grey99'),
    aspect.ratio = 1,
    legend.position = 'none', 
    legend.direction="horizontal") + 
  facet_wrap(.~factor(mag, labels = c('5-point', '10-point', '50-point')), scales = 'free')

# Supplementary figure S2: Policy switches for 5-point and 50-point options are 
# not related to the richness of the environment

# 5-pt option
plot_data <- d %>%
  filter(mag==5) %>%
  drop_na(mu_val_z) %>%
  mutate(mu_val_ind = cut(mu_val_z, breaks = unique(quantile(mu_val_z, probs = seq(0, 1, length.out = 11))), labels = FALSE, include.lowest = TRUE)) %>%
  group_by(mu_val_ind) %>%
  mutate(mu_val_graph = mean(mu_val_z, na.rm = T)) %>%
  ungroup() %>%
  group_by(subject, mag, prev_policy, mu_val_graph) %>% 
  summarise(m = 100*mean(pursue_switch, na.rm = T))

ggplot(plot_data, aes(x = mu_val_graph, y = m)) + 
  geom_smooth(method = 'lm', fill = 'deeppink4', color = 'black', linewidth = 0.5) +
  stat_summary(fun.data = 'mean_se', geom = 'pointrange', fill = 'deeppink4', shape = 21, color = 'black', size = 0.65, linewidth = 0.5) + 
  theme_pubr() + 
  scale_y_continuous(name = 'Switch-rate [%]', breaks = seq(0, 100, 10)) + 
  scale_x_continuous(name = 'Ave. value [Z]') + 
  theme(
    aspect.ratio = 1.1, 
    text = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.position = 'none'
  ) + 
  coord_cartesian(ylim = c(0, 30)) + 
  facet_grid(~factor(mag, labels = 'Pursue-sw. [5-pt]')) + 
  border('black')

plot_data <- d %>%
  filter(mag==5) %>%
  drop_na(mu_val_z) %>%
  mutate(mu_val_ind = cut(mu_val_z, breaks = unique(quantile(mu_val_z, probs = seq(0, 1, length.out = 11))), labels = FALSE, include.lowest = TRUE)) %>%
  group_by(mu_val_ind) %>%
  mutate(mu_val_graph = mean(mu_val_z, na.rm = T)) %>%
  ungroup() %>%
  group_by(subject, mag, prev_policy, mu_val_graph) %>% 
  summarise(m = 100*mean(reject_switch, na.rm = T))

ggplot(plot_data, aes(x = mu_val_graph, y = m)) + 
  geom_smooth(method = 'lm', fill = 'deeppink4', color = 'black', linewidth = 0.5) +
  stat_summary(fun.data = 'mean_se', geom = 'pointrange', fill = 'deeppink4', shape = 21, color = 'black', size = 0.65, linewidth = 0.5) + 
  theme_pubr() + 
  scale_y_continuous(name = 'Switch-rate [%]', position = 'right', breaks = seq(0, 40, 20)) + 
  scale_x_continuous(name = 'Ave. value [Z]') + 
  theme(
    aspect.ratio = 1.1, 
    text = element_text(size = 18),
    axis.text = element_text(size = 14),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    legend.position = 'none'
  ) + 
  coord_cartesian(ylim = c(0, 30)) + 
  facet_grid(~factor(mag, labels = 'Reject-sw. [5-pt]')) + 
  border('black')

# 50-pt option
plot_data <- d %>%
  filter(mag==50) %>%
  drop_na(mu_val_z) %>%
  mutate(mu_val_ind = cut(mu_val_z, breaks = unique(quantile(mu_val_z, probs = seq(0, 1, length.out = 11))), labels = FALSE, include.lowest = TRUE)) %>%
  group_by(mu_val_ind) %>%
  mutate(mu_val_graph = mean(mu_val_z, na.rm = T)) %>%
  ungroup() %>%
  group_by(subject, mag, prev_policy, mu_val_graph) %>% 
  summarise(m = 100*mean(pursue_switch, na.rm = T))

ggplot(plot_data, aes(x = mu_val_graph, y = m)) + 
  geom_smooth(method = 'lm', fill = 'deeppink4', color = 'black', linewidth = 0.5) +
  stat_summary(fun.data = 'mean_se', geom = 'pointrange', fill = 'deeppink4', shape = 21, color = 'black', size = 0.65, linewidth = 0.5) + 
  theme_pubr() + 
  scale_y_continuous(name = 'Switch-rate [%]', breaks = seq(0, 100, 20)) + 
  scale_x_continuous(name = 'Ave. value [Z]') + 
  theme(
    aspect.ratio = 1.1, 
    text = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.position = 'none'
  ) + 
  coord_cartesian(ylim = c(0, 20)) + 
  facet_grid(~factor(mag, labels = 'Pursue-sw. [50-pt]')) + 
  border('black')

plot_data <- d %>%
  filter(mag==50) %>%
  drop_na(mu_val_z) %>%
  mutate(mu_val_ind = cut(mu_val_z, breaks = unique(quantile(mu_val_z, probs = seq(0, 1, length.out = 11))), labels = FALSE, include.lowest = TRUE)) %>%
  group_by(mu_val_ind) %>%
  mutate(mu_val_graph = mean(mu_val_z, na.rm = T)) %>%
  ungroup() %>%
  group_by(subject, mag, prev_policy, mu_val_graph) %>% 
  summarise(m = 100*mean(reject_switch, na.rm = T))

ggplot(plot_data, aes(x = mu_val_graph, y = m)) + 
  geom_smooth(method = 'lm', fill = 'deeppink4', color = 'black', linewidth = 0.5) +
  stat_summary(fun.data = 'mean_se', geom = 'pointrange', fill = 'deeppink4', shape = 21, color = 'black', size = 0.65, linewidth = 0.5) + 
  theme_pubr() + 
  scale_y_continuous(name = 'Switch-rate [%]', position = 'right', breaks = seq(0, 40, 20)) + 
  scale_x_continuous(name = 'Ave. value [Z]') + 
  theme(
    aspect.ratio = 1.1, 
    text = element_text(size = 18),
    axis.text = element_text(size = 14),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    legend.position = 'none'
  ) + 
  coord_cartesian(ylim = c(0, 20)) + 
  facet_grid(~factor(mag, labels = 'Reject-sw. [50-pt]')) + 
  border('black')

# Supplementary figure S3: effect of value difference between current and average reward opportunity
d$mu_val_pe_z <- d$mag_z - d$mu_val_z #

plot_data <- d %>%
  drop_na(mu_val_pe_z) %>%
  mutate(mu_val_pe_ind = cut(mu_val_pe_z, breaks = unique(quantile(mu_val_pe_z, probs = seq(0, 1, length.out = 11))), labels = FALSE, include.lowest = TRUE)) %>%
  group_by(mu_val_pe_ind) %>%
  mutate(mu_val_pe_graph = mean(mu_val_pe_z, na.rm = T)) %>%
  ungroup() %>%
  group_by(subject, mag, prev_policy, mu_val_pe_graph) %>% 
  summarise(m = mean(response, na.rm = T))

ggplot(plot_data, aes(x = mu_val_pe_graph, y = m*100)) + 
  geom_smooth(method = 'lm', fill = 'skyblue1', color = 'black', linewidth = 0.5) +
  stat_summary(fun.data = 'mean_se', geom = 'pointrange', fill = 'skyblue2', shape = 21, color = 'black', size = 0.65, linewidth = 0.5) + 
  theme_pubr() + 
  scale_y_continuous(name = 'Pursue-rate [%]', breaks = seq(0, 100, 20)) + 
  scale_x_continuous(name = 'Rw-mag - Ave.value [Z]') + 
  theme(
    aspect.ratio = 1.2, 
    text = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.position = 'none'
  ) + 
  border('black')

GLM.s3 <- glmer(response ~ mu_val_pe_z + trial_z + (0 + mu_val_pe_z|subject), data = d, family = 'binomial', control = glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
summary(GLM.s3)


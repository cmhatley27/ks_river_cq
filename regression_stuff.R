# Load libraries and data ----------------------------------------------------------
library(tidyverse)
library(cluster)  
library(factoextra)
library(reshape)
library(brunnermunzel)
library(ggpubr)
source('load_and_combine.R')
source('Theme+Settings.R')


quantize <- function(x, q){
  if(q[1] == q[2]){
    return(cut(unlist(x), 
               quantile(unlist(x), c(0, q[1], 1)), 
               labels = c('low', 'high'),
               include.lowest = TRUE))
  } else{
    return(cut(unlist(x), 
               quantile(unlist(x), c(0, q, 1)), 
               labels = c('low', 'mid', 'high'),
               include.lowest = TRUE))
  }
}


# Select reduced variables ------------------------------------------------
event_summary_trim <- dplyr::select(event_summary, c(max_N, season, duration, duration_since_last, initial_Q, max_Q, delta_rat_Q, total_Q,
                                              starts_with('p_') & ends_with(c('1d_total', '7d_ante', '15d_ante', '30d_ante', '365d_ante')) & contains(c('10260015','10270104')),
                                              starts_with('d_') & contains(c('pet', 'bal')) & !ends_with(c('6','9')) & contains(c('10260015','10270104')),
                                              starts_with('r_'), reservoir_runoff_ratio))

# Regressions -------------------------------------------------------------
library(MASS)
mod_full <- lm(max_N ~ ., data = event_summary_trim)
mod_step <- stepAIC(mod_full, direction = 'backward', trace = FALSE)
mod1 <- lm(max_N~., data = mod_step$model)
summary(mod1)

mod_spring <- lm(FI_n ~ ., data = subset(event_summary_trim, season == 'Spring'))
mod_step <- stepAIC(mod_spring, direction = 'backward', trace = FALSE)
mod1 <- lm(FI_n~., data = mod_step$model)
summary(mod1)


library(MASS)
event_summary_fil <- filter(event_summary_trim) %>%
  dplyr::select(., where(is.numeric)) 
cors <- cor(event_summary_fil, use = 'pairwise.complete.obs', method = 'spearman') %>%
  melt(.) %>%
  filter(X1 == 'FI_n')

vars <- cors$X2[abs(cors$value) >= 0.1]

trimmer <- dplyr::select(event_summary, all_of(vars))
mod_full <- lm(FI_n ~ ., data = trimmer)
summary(mod_full)
mod_step <- stepAIC(mod_full, direction = 'backward', trace = FALSE)
mod1 <- lm(FI_n ~ duration_since_last + initial_Q + p_10270104_1d_total + 
             r_t_1d_delta + r_p_share + r_m_share + reservoir_runoff_ratio, data = event_summary)
summary(mod1)

mod2 <- lm(FI_n ~ (duration_since_last + initial_Q + p_10270104_1d_total + 
             r_t_1d_delta + r_p_share + r_m_share + reservoir_runoff_ratio)*reservoir_runoff_ratio*season, data = event_summary)
summary(mod2)
mod_step2 <- stepAIC(mod2, trace = FALSE)
mod_step2$anova

mod3 <- lm(FI_n ~ duration_since_last + initial_Q + p_10270104_1d_total + 
             r_t_1d_delta + r_p_share + r_m_share + reservoir_runoff_ratio + 
             season + initial_Q:reservoir_runoff_ratio + p_10270104_1d_total:reservoir_runoff_ratio + 
             r_t_1d_delta:reservoir_runoff_ratio + r_p_share:reservoir_runoff_ratio + 
             r_m_share:reservoir_runoff_ratio + initial_Q:season + p_10270104_1d_total:season + 
             r_t_1d_delta:season + r_p_share:season + r_m_share:season + 
             reservoir_runoff_ratio:season + initial_Q:reservoir_runoff_ratio:season + 
             p_10270104_1d_total:reservoir_runoff_ratio:season + r_t_1d_delta:reservoir_runoff_ratio:season + 
             r_p_share:reservoir_runoff_ratio:season + r_m_share:reservoir_runoff_ratio:season, data = event_summary)
summary(mod3)

anova(mod1, mod3)
ggplot() +
  geom_point(aes(x = mod3$model$FI_n, y = mod3$fitted.values)) +
  xlim(c(-1,1)) +
  ylim(c(-1,1)) +
  geom_abline(slope = 1)

ggplot(data = event_summary) +
  geom_point(aes(x = reservoir_runoff_ratio, y = p_10260015_30d_ante))
cor(event_summary$reservoir_runoff_ratio, event_summary$p_10270104_1d_total, use = 'pairwise.complete.obs', method = 'spearman')
# Interactions with res ---------------------------------------------------

mod1 <- lm(FI_n ~ (initial_Q)*reservoir_runoff_ratio*season, data = event_summary)
summary(mod1)

anova(mod1, mod3)
library(MASS)
mod_step <- stepAIC(mod1, direction = 'backward', trace = FALSE)
mod2 <- lm(FI_n ~ p_10260010_1d_total + initial_Q + reservoir_runoff_ratio + 
             season + p_10260010_1d_total:reservoir_runoff_ratio + initial_Q:reservoir_runoff_ratio + 
             initial_Q:season + reservoir_runoff_ratio:season + initial_Q:reservoir_runoff_ratio:season, data = mod_step$model)
summary(mod2)
mod_step$anova
mod2 <- lm(FI_n ~ (p_10270104_1d_total + initial_Q)*reservoir_runoff_ratio, data = event_summary)
summary(mod2)

anova(mod1,mod2)


mod1 <- lm(FI_n ~ p_10260015_30d_ante, data = event_summary)
summary(mod1)

mod3 <- lm(FI_n ~ p_10260015_30d_ante*reservoir_runoff_ratio + p_10260015_1d_total*reservoir_runoff_ratio, data = event_summary)
summary(mod3)

anova(mod2,mod3)
cor.test(event_summary$p_10260008_30d_ante, event_summary$p_10260008_1d_total, use = 'pairwise.complete.obs')

library(MASS)
mod <- lm(FI_n ~ (delta_rat_Q + initial_Q + duration_since_last + p_10260008_1d_total + p_10260008_30d_ante)*reservoir_runoff_ratio, data = event_summary)
summary(mod)

step <- stepAIC(mod, trace = FALSE)
step$anova

mod2 <- lm(FI_n ~ delta_rat_Q + initial_Q + duration_since_last + p_10260008_30d_ante + 
             reservoir_runoff_ratio + delta_rat_Q:reservoir_runoff_ratio + 
             p_10260008_30d_ante:reservoir_runoff_ratio, data = event_summary)
summary(mod2)




# Based on sig cors -------------------------------------------------------

library(MASS)
mod1 <- lm(FI_n ~ initial_Q + delta_rat_Q + duration_since_last + d_10260008_spei12 + p_10260008_1d_total + p_10260008_30d_ante, data = event_summary)
summary(mod1)
mod1_step <- stepAIC(mod1)
mod1_step$anova
mod1_reduce <- lm(FI_n ~ initial_Q + delta_rat_Q, data = event_summary)
summary(mod1_reduce)

mod2 <- lm(FI_n ~ (initial_Q + delta_rat_Q + duration_since_last + d_10260008_spei12 + p_10260008_1d_total + p_10260008_30d_ante)*reservoir_runoff_ratio, data = event_summary)
summary(mod2)
mod2_step <- stepAIC(mod2)
mod2_step$anova
mod2_reduce <- lm(FI_n ~ initial_Q + delta_rat_Q + duration_since_last + p_10260008_30d_ante + 
                    reservoir_runoff_ratio + delta_rat_Q:reservoir_runoff_ratio + 
                    p_10260008_30d_ante:reservoir_runoff_ratio, data = event_summary)
summary(mod2_reduce)

anova(mod1, mod2)
anova(mod1, mod1_reduce)
anova(mod2, mod2_reduce)
anova(mod1_reduce, mod2_reduce)


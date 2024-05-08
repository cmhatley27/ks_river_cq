library(tidyverse)
library(cluster)  
library(factoextra)
library(reshape)
library(brunnermunzel)
library(ggpubr)
library(sf)
library(corrplot)
library(tmap)
source('load_and_combine.R')
source('Theme+Settings.R')


# correlations ------------------------------------------------------------

event_summary_numeric <- select(event_summary, c(FI_n, HI_n, max_N, duration, total_Q, delta_Q, delta_rat_Q, initial_Q, initial_N, duration_since_last, p_10260008_30d_ante, d_10260008_spei12, p_10260008_1d_total, reservoir_runoff_ratio))
event_summary_numeric_lowres <- event_summary_numeric %>%
  filter(reservoir_runoff_ratio <= 0.55)
event_summary_numeric_highres <- event_summary_numeric %>%
  filter(reservoir_runoff_ratio > 0.55)

cors <- tibble(var = colnames(event_summary_numeric),
               cor = NA,
               p = NA,
               cor_lowres = NA,
               p_lowres = NA,
               cor_highres = NA,
               p_highres = NA)

for(variable in 1:nrow(cors)){
  variable_select <- event_summary_numeric[cors$var[variable]]
  cor_result <- cor.test(event_summary_numeric$FI_n, unlist(variable_select), method = 'spearman', use = 'pairwise.complete.obs')
  cors[variable,'cor'] <- cor_result$estimate
  cors[variable,'p'] <- cor_result$p.value
  
  variable_select <- event_summary_numeric_lowres[cors$var[variable]]
  cor_result <- cor.test(event_summary_numeric_lowres$FI_n, unlist(variable_select), method = 'spearman', use = 'pairwise.complete.obs')
  cors[variable,'cor_lowres'] <- cor_result$estimate
  cors[variable,'p_lowres'] <- cor_result$p.value
  
  variable_select <- event_summary_numeric_highres[cors$var[variable]]
  cor_result <- cor.test(event_summary_numeric_highres$FI_n, unlist(variable_select), method = 'spearman', use = 'pairwise.complete.obs')
  cors[variable,'cor_highres'] <- cor_result$estimate
  cors[variable,'p_highres'] <- cor_result$p.value
}
sig_level <- 0.05
cors_fil <- mutate(cors,
                   sig = p <= sig_level,
                   sig_lowres = p_lowres <= sig_level,
                   sig_highres = p_highres <= sig_level)# %>%
  filter(sig | sig_lowres | sig_highres)

corrplot(cor(event_summary_numeric, use = 'pairwise.complete.obs', method = 'spearman'), method = 'number', tl.cex = 0.8)
cor(event_summary_numeric) %>%
  melt(.)

cor.mtest(event_summary_numeric, use = 'pairwise.complete.obs', method = 'spearman', conf.level = 0.95)$p
cor_vars <- colnames(select(event_summary, c(starts_with('r_'), ends_with('_Q'), ends_with('ratio'), season_index, duration, duration_since_last)))

res_var <- 'reservoir_runoff_ratio'
res_name <- 'Res/Runoff'
res_q <- c(0.5,0.5)

response_vars <- c('FI_n', 'HI_n')

sig_level <- 0.05

dat <- select(event_summary, c(all_of(c(response_vars, res_var, cor_vars)))) %>%
  mutate(res_end = quantize(select(.,res_var), res_q)) %>%
  pivot_longer(all_of(cor_vars), names_to = 'cor_var', values_to = 'cor_val')

cor_calcs <- function(x, cor_var, sig, res = NULL, res_side = c('low','high')){
  if(is.null(res)){
    cor_stat <- cor.test(x,cor_var, use = 'pairwise.complete.obs')$estimate
    cor_p <- cor.test(x,cor_var, use = 'pairwise.complete.obs')$p.value
  }
  else{
    cor_stat <- cor.test(x[res == res_side],cor_var[res == res_side], use = 'pairwise.complete.obs')$estimate
    cor_p <- cor.test(x[res == res_side],cor_var[res == res_side], use = 'pairwise.complete.obs')$p.value
  }
  return(ifelse(cor_p <= sig, cor_stat, NA))
}

cors <- dat %>%
  group_by(cor_var) %>%
  summarise(across(all_of(response_vars),
                   ~cor_calcs(.x, cor_val, sig_level),
                   .names = '{.col}_stat'),
            across(all_of(response_vars),
                   ~cor_calcs(.x, cor_val, sig_level, res_end, res_side = 'low'),
                   .names = '{.col}_lowres_stat'),
            across(all_of(response_vars),
                   ~cor_calcs(.x, cor_val, sig_level, res_end, res_side = 'high'),
                   .names = '{.col}_highres_stat')) %>%
  select(c(cor_var, starts_with('FI'), starts_with('HI')))

FI_cors <- cors %>%
  select(cor_var, FI_n_lowres_stat, FI_n_highres_stat) %>%
  pivot_longer(!cor_var, names_to = 'var', values_to = 'val')
ggplot(data = FI_cors, aes(x = var, y = val, fill = var)) +
  geom_col(color = 'black') +
  geom_hline(yintercept = 0) +
  facet_wrap(vars(cor_var))



var(event_summary$FI_n, na.rm = T)
var(event_summary$FI_n[event_summary$reservoir_runoff_ratio <= 0.55], na.rm = T)
var(event_summary$FI_n[event_summary$reservoir_runoff_ratio >= 0.55], na.rm = T)



# corrplot for trimmed variables ------------------------------------------

event_summary_trim <- select(event_summary, c(FI_n, HI_n, initial_Q, max_Q, delta_Q, delta_rat_Q, total_Q, duration, duration_since_last,
                                              p_10260008_30d_ante, d_10260008_spei12, p_10260008_1d_total, reservoir_runoff_ratio))

rho <- cor(event_summary_trim, use = 'pairwise.complete.obs', method = 'spearman')
p <- cor.mtest(event_summary_trim, use = 'pairwise.complete.obs', method = 'spearman', conf.level = 0.95)$p
corrplot(rho, method = 'number', tl.cex = 0.8, p.mat = p, sig.level = 0.05)



event_summary_trim_lowres <- event_summary_trim %>%
  filter(reservoir_runoff_ratio <= 0.55)
event_summary_trim_highres <- event_summary_trim %>%
  filter(reservoir_runoff_ratio > 0.55)

cors <- tibble(var = colnames(event_summary_trim),
               cor = NA,
               p = NA,
               cor_lowres = NA,
               p_lowres = NA,
               cor_highres = NA,
               p_highres = NA)

for(variable in 1:nrow(cors)){
  variable_select <- event_summary_trim[cors$var[variable]]
  cor_result <- cor.test(event_summary_trim$HI_n, unlist(variable_select), method = 'spearman', use = 'pairwise.complete.obs')
  cors[variable,'cor'] <- cor_result$estimate
  cors[variable,'p'] <- cor_result$p.value
  
  variable_select <- event_summary_trim_lowres[cors$var[variable]]
  cor_result <- cor.test(event_summary_trim_lowres$HI_n, unlist(variable_select), method = 'spearman', use = 'pairwise.complete.obs')
  cors[variable,'cor_lowres'] <- cor_result$estimate
  cors[variable,'p_lowres'] <- cor_result$p.value
  
  variable_select <- event_summary_trim_highres[cors$var[variable]]
  cor_result <- cor.test(event_summary_trim_highres$HI_n, unlist(variable_select), method = 'spearman', use = 'pairwise.complete.obs')
  cors[variable,'cor_highres'] <- cor_result$estimate
  cors[variable,'p_highres'] <- cor_result$p.value
}
sig_level <- 0.05
cors_fil <- mutate(cors,
                   sig = p <= sig_level,
                   sig_lowres = p_lowres <= sig_level,
                   sig_highres = p_highres <= sig_level)# %>%
filter(sig | sig_lowres | sig_highres)


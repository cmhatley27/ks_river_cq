# load data ---------------------------------------------------------------
library(tidyverse)
library(ggpubr)
library(rstatix)
source('8_combine_all_event_data.R')
source('Theme+Settings.R')

plot_dat <- event_summary %>%
  mutate(cluster = factor(cluster))

cluster_means <- group_by(plot_dat, cluster) %>%
  summarise(FI = mean(FI_n),
            HI = mean(HI_n))

# calculate characteristic differences ---------------------------------------
vars <- colnames(select(event_summary, c(initial_Q, duration_since_last, 
                                         max_Q, delta_Q, delta_rat_Q, duration, 
                                         starts_with('p_') & ends_with('1d_total'),
                                         starts_with('p_') & ends_with('30d_ante'),
                                         starts_with('d_') & ends_with('spei12'), 
                                         starts_with('r_') & ends_with('_share'), reservoir_runoff_ratio)))
sig_level <- 0.05
best_k <- 7
vars_kw <- tibble(
  var = vars,
  p = NA,
  kw_sig = NA
)

for(var in 1:nrow(vars_kw)){
  kw_p <- kruskal.test(unlist(plot_dat[vars[var]]) ~ plot_dat$cluster)$p.value
  vars_kw[var,2] <- kw_p
  vars_kw[var,3] <- kw_p < sig_level
}

sig_vars_kw <- vars_kw$var[vars_kw$kw_sig == TRUE]

combos <- get_comparisons(plot_dat, cluster)
n_combos <- length(combos)

con_results_long <- tibble()
con_results_short <- tibble(
  var = sig_vars_kw,
  n_sig_con = NA,
  sig_combs = NA
)
library(conover.test)
for(var in 1:length(sig_vars_kw)){
  con_test <- conover.test(unlist(plot_dat[sig_vars_kw[var]]), plot_dat$cluster, method = 'bh')
  test_ps <- con_test$P.adjusted
  test_comb <- con_test$comparisons
  con_results_long <- rbind(con_results_long, tibble(
    characteristic = rep(sig_vars_kw[var], length(test_ps)),
    cluster1 = str_sub(test_comb,1,1),
    cluster2 = str_sub(test_comb,-1,-1),
    p.adj = test_ps,
    conover_stat = con_test$T
  ))
  
  sig_comb <- test_comb[test_ps < (sig_level/2)]
  con_results_short[var,2] <- sum(test_ps < (sig_level/2), na.rm = TRUE)
  con_results_short[var,3] <- paste(sig_comb, collapse = ', ')
  
  for(k in 1:best_k){
    k_count <- sum(str_detect(sig_comb, as.character(k)))
    con_results_short[var,3+k] <- k_count
  }
  names(con_results_short)[4:ncol(con_results_short)] <- paste0('n_',seq(1:best_k))
}

# plot scatter ------------------------------------------------------------

ggplot(plot_dat, aes(x = FI_n, y = HI_n, color = cluster)) +
  geom_hline(yintercept = 0, color = 'grey') +
  geom_vline(xintercept = 0, color = 'grey') +
  geom_point(size = 1) +
  geom_label(data = cluster_means, aes(x = FI, y = HI, label = cluster)) +
  scale_color_discrete(guide = 'none') +
  ylim(-1,1) +
  xlim(-1,1) +
  xlab('Flushing Index') +
  ylab('Hysteresis Index') +
  theme(axis.title = element_text(size = 10))
ggsave(file.path('Figures','final','clusters','clusters_scatter.png'), height = 2.625, width = 3, units = 'in')


# plot boxes --------------------------------------------------------------

box_dat <- select(plot_dat, c(cluster, initial_Q, delta_rat_Q, duration_since_last, p_10270207_1d_total, p_10270104_1d_total, r_m_share)) %>%
  pivot_longer(!cluster, names_to = 'var', values_to = 'val') %>%
  mutate(var = factor(var, levels = c('p_10270207_1d_total', 'p_10270104_1d_total', 
                                      'initial_Q', 'duration_since_last', 'delta_rat_Q', 
                                      'r_m_share')))

var_names <- list(
  'p_10270207_1d_total' = 'Event Precipitation [mm]\n(Central HUC8)',
  'p_10270104_1d_total' = 'Event Precipitation [mm]\n(Easternmost HUC8)',
  'initial_Q' = expression('Initial Q'~'['*m^3*'/s]'),
  'duration_since_last' = 'Duration Since\nLast Event [days]',
  'delta_rat_Q' = 'Delta/Initial Q [-]',
  'r_m_share' = 'Milford Lake\nOutflow-Runoff Ratio [-]'
)
var_labeller <- function(variable,value){
  return(var_names[value])
}

ggplot(data = box_dat, aes(x = cluster, fill = cluster, y = val)) +
  geom_boxplot(outliers = F, size = 0.25) +
  facet_wrap(vars(var), scales = 'free_y', labeller = var_labeller, ncol = 2) +
  scale_fill_discrete(guide = 'none') +
  xlab('Cluster') +
  scale_y_continuous(name = NULL, expand = expansion(mult = c(0.05, 0.3))) +
  theme(strip.text = element_text(size = 9))
ggsave(file.path('Figures','final','clusters','boxes.png'), height = 4.5, width = 3.5, units = 'in')

# Load libraries and data ---------------------------------------------------------------

# Libraries
library(tidyverse)
library(cluster)  
library(factoextra)

CQ_summary <- read_csv('./DataFiles/hydro_data/other_params/CQ_summary.csv') %>%
  select(FI_n, HI_n)

simple_temporal <- read_csv('./DataFiles/hydro_data/event_char/simple_temporal.csv')
simple_size <- read_csv('./DataFiles/hydro_data/event_char/simple_size.csv')

precip_inevent_gridded <- read_csv('./DataFiles/hydro_data/event_char/precip_inevent_gridded.csv') %>%
  select(contains('1d'))
precip_anteevent_gridded <- read_csv('./DataFiles/hydro_data/event_char/precip_anteevent_gridded.csv') %>%
  select(contains('30d'))
spei12 <- read_csv('./DataFiles/hydro_data/event_char/spei12.csv')

flow_ratios <- read_csv('./DataFiles/hydro_data/event_char/flow_ratios.csv')

#Combine
event_summary <- tibble(
  CQ_summary, 
  simple_temporal, 
  simple_size,
  precip_inevent_gridded,
  precip_anteevent_gridded,
  flow_ratios,
  spei12
) %>%
  mutate(wateryear = factor(wateryear),
         month = factor(month, levels = c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec')),
         season = factor(season, levels = c('Winter', 'Spring', 'Summer', 'Fall')))

rm(list = ls()[ls() != 'event_summary'])
source('Theme+Settings.R')


# set up clustering data --------------------------------------------------
hclust_data <- event_summary %>%
  filter(complete.cases(.)) %>%
  select(c(FI_n, HI_n))

dist <- daisy(hclust_data, metric = 'euclidean')

cls <- agnes(dist, method = 'ward')

# Loop through different numbers of clusters and count sig differences -----
vars <- colnames(select(event_summary, c(initial_Q, max_Q, delta_Q, 
                                         delta_rat_Q, total_Q, duration, 
                                         duration_since_last, starts_with('p_') & ends_with('1d_total'), 
                                         starts_with('p_') & ends_with('30d_ante'), 
                                         starts_with('d_') & ends_with('spei12'),  
                                         starts_with('r_') & ends_with('_share'), reservoir_runoff_ratio)))
sig_level <- 0.05
con_sigs <- tibble(
  k = seq(3,12),
  sig_n = NA
)
for(ks in seq(3,12)){
  best_k <- ks
  
  event_ids <- mutate(hclust_data, cluster = cutree(cls, k = best_k))
  cls_id <- left_join(event_summary, event_ids)$cluster
  
  event_summary_cluster <- left_join(event_summary, event_ids) %>%
    filter(!is.na(cluster)) %>%
    mutate(season = factor(season, levels = c('Winter', 'Spring', 'Summer', 'Fall')))
  
  vars_kw <- tibble(
    var = vars,
    kw_sig = NA
  )
  
  for(var in 1:nrow(vars_kw)){
    kw_p <- kruskal.test(unlist(event_summary_cluster[vars[var]]) ~ event_summary_cluster$cluster)$p.value
    vars_kw[var,2] <- kw_p < sig_level
  }
  
  vars_con <- tibble(
    var = vars_kw$var[vars_kw$kw_sig == TRUE],
    n_sig_con = NA,
    sig_combs = NA
  )
  
  library(conover.test)
  for(var in 1:nrow(vars_con)){
    con_test <- conover.test(unlist(event_summary_cluster[vars_con$var[var]]), event_summary_cluster$cluster, method = 'bh')
    test_ps <- con_test$P.adjusted
    test_comb <- con_test$comparisons
    sig_comb <- test_comb[test_ps < (sig_level/2)]
    vars_con[var,2] <- sum(test_ps < (sig_level/2), na.rm = TRUE)
    vars_con[var,3] <- paste(sig_comb, collapse = ', ')
    
    for(k in 1:best_k){
      k_count <- sum(str_detect(sig_comb, as.character(k)))
      vars_con[var,3+k] <- k_count
    }
    names(vars_con)[4:ncol(vars_con)] <- paste0('n_',seq(1:best_k))
  }
  con_sigs[ks-2,2] <- sum(vars_con$n_sig_con)
}

ggplot(data = con_sigs, aes(x = k, y = sig_n)) +
  geom_line() +
  geom_point(aes(color = factor(k)), size = 4) +
  geom_vline(xintercept = 7, linetype = 'dashed') +
  scale_x_continuous(breaks = seq(3,12))


# Examine best k ----------------------------------------------------------
vars <- colnames(select(event_summary, c(initial_Q, max_Q, delta_Q, 
                                         delta_rat_Q, total_Q, duration, 
                                         duration_since_last, starts_with('p_') & ends_with('1d_total'), 
                                         starts_with('p_') & ends_with('30d_ante'), 
                                         starts_with('d_') & ends_with('spei12'),  
                                         starts_with('r_') & ends_with('_share'), reservoir_runoff_ratio)))
sig_level <- 0.05

best_k <- 7

# Renumber clusters so they plot in order
cluster_ids <- mutate(hclust_data, cluster = cutree(cls, k = best_k)) %>%
  mutate(cluster = factor(case_when(
    cluster == 1 ~ 5,
    cluster == 2 ~ 7,
    cluster == 3 ~ 1,
    cluster == 4 ~ 2,
    cluster == 5 ~ 3,
    cluster == 6 ~ 6,
    cluster == 7 ~ 4,
  )))

event_summary_cluster <- left_join(event_summary, cluster_ids)

ggplot(data = event_summary_cluster) +
  geom_point(aes(x = FI_n, y = HI_n, color = factor(cluster))) +
  scale_color_discrete(name = 'Cluster') +
  xlab('') +
  ylab('') +
  xlim(-1,1) +
  ylim(-1,1) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0)

write_csv(select(event_summary_cluster, cluster), './DataFiles/cluster_ids.csv')


# Calculate characteristic differences between clusters -------------------

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
write_csv(con_results_long, './DataFiles/cluster_differences.csv')

library(tidyverse)
library(cluster)  
library(factoextra)

# Load data ---------------------------------------------------------------
# Basic characteristics 

CQ_summary <- read_csv('./DataFiles/hydro_data/other_params/CQ_summary.csv') %>%
  select(FI_n, HI_n)

simple_temporal <- read_csv('./DataFiles/hydro_data/event_char/simple_temporal.csv')
simple_size <- read_csv('./DataFiles/hydro_data/event_char/simple_size.csv')

N_yield <- read_csv('./DataFiles/hydro_data/event_char/N_yield.csv') %>%
  cbind(., CQ_summary)

#* Gridded Precip Predictors 
precip_inevent_gridded <- read_csv('./DataFiles/hydro_data/event_char/precip_inevent_gridded.csv') %>%
  select(contains('1d'))

precip_anteevent_gridded <- read_csv('./DataFiles/hydro_data/event_char/precip_anteevent_gridded.csv') 

precip_unintercepted_volumes <- read_csv('./DataFiles/hydro_data/event_char/precip_unintercepted_volumes.csv')

max_precip_coords <- read_csv('./DataFiles/hydro_data/event_char/max_precip_coords.csv')

avg_precip_intensity <- read_csv('./DataFiles/hydro_data/event_char/avg_precip_intensity.csv') %>%
  select(contains('1d'))

max_precip_intensity <- read_csv('./DataFiles/hydro_data/event_char/max_precip_intensity.csv')


# * Outflows Predictors 
outflows_inevent_totals <- read_csv('./DataFiles/hydro_data/event_char/outflows_inevent_totals.csv') %>%
  select(contains('1d'))

outflows_deltas <- read_csv('./DataFiles/hydro_data/event_char/outflows_deltas.csv') %>%
  select(contains('1d'))

flow_ratios <- read_csv('./DataFiles/hydro_data/event_char/flow_ratios.csv')


# Drought Predictors 
spei1 <- read_csv('./DataFiles/hydro_data/event_char/spei1.csv')

spei3 <- read_csv('./DataFiles/hydro_data/event_char/spei3.csv')

spei6 <- read_csv('./DataFiles/hydro_data/event_char/spei6.csv')

spei9 <- read_csv('./DataFiles/hydro_data/event_char/spei9.csv')

spei12 <- read_csv('./DataFiles/hydro_data/event_char/spei12.csv')

spei18 <- read_csv('./DataFiles/hydro_data/event_char/spei18.csv')

spei24 <- read_csv('./DataFiles/hydro_data/event_char/spei24.csv')

pet1 <- read_csv('./DataFiles/hydro_data/event_char/pet1.csv')

pet3 <- read_csv('./DataFiles/hydro_data/event_char/pet3.csv')

pet6 <- read_csv('./DataFiles/hydro_data/event_char/pet6.csv')

pet9 <- read_csv('./DataFiles/hydro_data/event_char/pet9.csv')

pet12 <- read_csv('./DataFiles/hydro_data/event_char/pet12.csv')

bal1 <- read_csv('./DataFiles/hydro_data/event_char/bal1.csv')

bal3 <- read_csv('./DataFiles/hydro_data/event_char/bal3.csv')

bal6 <- read_csv('./DataFiles/hydro_data/event_char/bal6.csv')

bal9 <- read_csv('./DataFiles/hydro_data/event_char/bal9.csv')

bal12 <- read_csv('./DataFiles/hydro_data/event_char/bal12.csv')

droughts <- tibble(spei1, spei3, spei6, spei9, spei12, spei18, spei24, 
                   pet1, pet3, pet6, pet9, pet12,
                   bal1, bal3, bal6, bal9, bal12)


# * Correlated Predictors 

correlated_predictors <- unlist(as.vector(read_csv('./DataFiles/hydro_data/event_char/correlated_predictors.csv'))) 

#* Combine

event_summary <- tibble(
  CQ_summary, 
  simple_temporal, 
  simple_size,
  precip_inevent_gridded,
  precip_anteevent_gridded,
  precip_unintercepted_volumes,
  max_precip_coords,
  avg_precip_intensity,
  max_precip_intensity,
  outflows_inevent_totals,
  outflows_deltas,
  flow_ratios,
  droughts
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

fviz_nbclust(hclust_data, hcut, method = 'wss', k.max = 10)


# Set k manually ----------------------------------------------------------
k <- 4
event_ids <- mutate(hclust_data, cluster = cutree(cls, k = k))
cls_id <- left_join(event_summary, event_ids)$cluster
#write_csv(as_tibble(cls_id), 'C:/School/SAFE KAW/Data/DataFiles/cluster_data/cluster_assignment.csv')

#Add cluster ID back in
event_summary_cluster <- left_join(event_summary, event_ids) %>%
  filter(!is.na(cluster)) %>%
  mutate(season = factor(season, levels = c('Winter', 'Spring', 'Summer', 'Fall')))
cluster_summaries <- event_summary_cluster %>%
  group_by(cluster) %>%
  summarise(across(where(is.numeric), median))

#FI HI scatter
ggplot(data = event_summary_cluster) +
  geom_point(aes(x = FI_n, y = HI_n, color = factor(cluster))) +
  scale_color_discrete(name = 'Cluster') +
  xlab('') +
  ylab('') +
  xlim(-1,1) +
  ylim(-1,1) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0)


# Loop ks -----------------------------------------------------------------
vars <- colnames(select(event_summary, c(initial_Q, max_Q, delta_Q, 
                                delta_rat_Q, total_Q, duration, 
                                duration_since_last, starts_with('p_') & ends_with('1d_total'),
                                starts_with('p_') & ends_with('30d_ante'), 
                                starts_with('p_') & ends_with('90d_ante'),
                                starts_with('p_') & ends_with('15d_ante'),
                                starts_with('d_') & ends_with('spei12'), 
                                starts_with('r_') & !contains('_c_'), reservoir_runoff_ratio)))

# vars <- colnames(select(event_summary, where(is.numeric) & !c(FI_n,HI_n)))

sig_level <- 0.05
k_summary <- tibble(
  k_n = seq(2, 12),
  diff_vars = NA,
  kw_count = NA
)

for(k in k_summary$k_n){
  event_ids <- mutate(hclust_data, cluster = cutree(cls, k = k))
  cls_id <- left_join(event_summary, event_ids)$cluster

  event_summary_cluster <- left_join(event_summary, event_ids) %>%
    filter(!is.na(cluster)) %>%
    mutate(season = factor(season, levels = c('Winter', 'Spring', 'Summer', 'Fall')))
  
  ps <- 0
  kw_n <- 0
  for(var in 1:length(vars)) {
    test_ps <- pairwise.wilcox.test(unlist(event_summary_cluster[vars[var]]), event_summary_cluster$cluster, p.adjust.method = 'none')$p.value
    ps <- ps + sum(test_ps <= sig_level, na.rm = TRUE)
    kw_p <- kruskal.test(unlist(event_summary_cluster[vars[var]]) ~ event_summary_cluster$cluster)$p.value
    kw_sig <- kw_p < sig_level
    kw_n <- kw_n + kw_sig
  }
  
  k_summary[k-1,2] <- ps
  k_summary[k-1,3] <- kw_n
}

k_summary <- k_summary %>%
  mutate(ncomb = choose(k_n,2)*length(vars),
         pct = diff_vars/ncomb,
         kw_pct = kw_count/length(vars))

ggplot(data = k_summary, aes(x = k_n, y = kw_count)) +
  scale_x_continuous(breaks = seq(1,12)) +
  geom_line()


# Loop conovers -----------------------------------------------------------
vars <- colnames(select(event_summary, c(initial_Q, max_Q, delta_Q, 
                                         delta_rat_Q, total_Q, duration, 
                                         duration_since_last, starts_with('p_') & ends_with('1d_total'),
                                         starts_with('p_') & ends_with('30d_ante'), 
                                         starts_with('p_') & ends_with('90d_ante'),
                                         starts_with('p_') & ends_with('15d_ante'),
                                         starts_with('d_') & ends_with('spei12'), 
                                         starts_with('r_') & !contains('_c_'), reservoir_runoff_ratio)))
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
  geom_vline(xintercept = 7, linetype = 'dashed') +
  scale_x_continuous(breaks = seq(3,12))


# Examine best k ----------------------------------------------------------
vars <- colnames(select(event_summary, c(initial_Q, max_Q, delta_Q, 
                                         delta_rat_Q, total_Q, duration, 
                                         duration_since_last, starts_with('p_') & ends_with('1d_total'),
                                         starts_with('p_') & ends_with('30d_ante'), 
                                         starts_with('p_') & ends_with('90d_ante'),
                                         starts_with('p_') & ends_with('15d_ante'),
                                         starts_with('d_') & ends_with('spei12'), 
                                         starts_with('r_') & !contains('_c_'), reservoir_runoff_ratio)))
sig_level <- 0.05

best_k <- 7

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


ggplot(data = event_summary_cluster, aes(x = factor(month), fill = factor(season))) +
  geom_bar(color = 'black') +
  scale_fill_manual(values = c('#8DA9C4','#60992D','#C33C54','#F4B860')) +
  facet_wrap(vars(factor(cluster)), nrow = 2)

test_var <- 'reservoir_runoff_ratio'
conover.test(event_summary_cluster$reservoir_runoff_ratio, event_summary_cluster$cluster, method = 'bh')
ggplot(data = event_summary_cluster, aes(x = factor(cluster), fill = factor(cluster), y = r_m_share)) +
  geom_boxplot() #+
  ylim(c(0,2))
  

month_counts <- event_summary_cluster %>% count(month, cluster) %>%
  group_by(month) %>%
  mutate(month_sum = sum(n),
         prop = n/month_sum)
ggplot(data = month_counts) +
  geom_col(aes(x = month, y = n, fill = factor(cluster), group = factor(cluster)), color = 'black')

ggplot(data = event_summary_cluster) +
  geom_point(aes(x = FI_n, y = HI_n, color = factor(cluster))) +
  scale_color_discrete(name = 'Cluster') +
  xlab('') +
  ylab('') +
  xlim(-1,1) +
  ylim(-1,1) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0)
#ggsave('C:/School/SAFE KAW/Data/Figures/maps/cluster_summaries/scatterk9.tiff', height = 6, width = 7, units = 'in', compression = 'lzw')


# Seasonal diffs ----------------------------------------------------------

library(reporttools)

pairwise.fisher.test(event_summary_cluster$month, event_summary_cluster$cluster, p.adjust.method = 'fdr')

conover.test(event_summary_cluster$HI_n, event_summary_cluster$season, method = 'bh')

counts <- count(event_summary_cluster, event_summary_cluster$month)$n
counts_clus12 <- count(filter(event_summary_cluster, cluster %in% c(1,2)), month)$n

chisq.test(counts_clus12, p = counts, rescale.p = TRUE)

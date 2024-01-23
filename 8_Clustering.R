library(tidyverse)
library(cluster)  
library(factoextra)



# Load data ---------------------------------------------------------------


# Basic characteristics --------------------------------------------------

CQ_summary <- read_csv('./DataFiles/hydro_data/event_char/CQ_summary.csv') %>%
  select(FI, HI = HI_mean)
  
simple_temporal <- read_csv('./DataFiles/hydro_data/event_char/simple_temporal.csv')
simple_size <- read_csv('./DataFiles/hydro_data/event_char/simple_size.csv')

N_yield <- read_csv('./DataFiles/hydro_data/event_char/N_yield.csv') %>%
  cbind(., CQ_summary)

#* Gridded Precip Predictors -----------------------------------------------
precip_inevent_gridded <- read_csv('./DataFiles/hydro_data/event_char/precip_inevent_gridded.csv') %>%
  select(contains('1d')) %>%
  select(contains(c('10260008', '10260010', '10260015', '102701')))

precip_anteevent_gridded <- read_csv('./DataFiles/hydro_data/event_char/precip_anteevent_gridded.csv') %>%
  select(contains(c('10260008', '10260010', '10260015', '102701')))

precip_unintercepted_volumes <- read_csv('./DataFiles/hydro_data/event_char/precip_unintercepted_volumes.csv')

max_precip_coords <- read_csv('./DataFiles/hydro_data/event_char/max_precip_coords.csv')

avg_precip_intensity <- read_csv('./DataFiles/hydro_data/event_char/avg_precip_intensity.csv') %>%
  select(contains('1d')) %>%
  select(contains(c('10260008', '10260010', '10260015', '102701')))

max_precip_intensity <- read_csv('./DataFiles/hydro_data/event_char/max_precip_intensity.csv') %>%
  select(contains(c('10260008', '10260010', '10260015', '102701')))


# * Outflows Predictors ---------------------------------------------------
outflows_inevent_totals <- read_csv('./DataFiles/hydro_data/event_char/outflows_inevent_totals.csv') %>%
  select(contains('1d'))

outflows_deltas <- read_csv('./DataFiles/hydro_data/event_char/outflows_deltas.csv') %>%
  select(contains('1d'))

flow_ratios <- read_csv('./DataFiles/hydro_data/event_char/flow_ratios.csv')


# Drought Predictors ------------------------------------------------------
spei1 <- read_csv('./DataFiles/hydro_data/event_char/spei1.csv') %>%
  select(contains(c('10260015', '10270104')))

spei3 <- read_csv('./DataFiles/hydro_data/event_char/spei3.csv') %>%
  select(contains(c('10260015', '10270104')))

spei6 <- read_csv('./DataFiles/hydro_data/event_char/spei6.csv') %>%
  select(contains(c('10260015', '10270104')))

spei9 <- read_csv('./DataFiles/hydro_data/event_char/spei9.csv') %>%
  select(contains(c('10260015', '10270104')))

spei12 <- read_csv('./DataFiles/hydro_data/event_char/spei12.csv') %>%
  select(contains(c('10260015', '10270104')))

spei18 <- read_csv('./DataFiles/hydro_data/event_char/spei18.csv') %>%
  select(contains(c('10260015', '10270104')))

spei24 <- read_csv('./DataFiles/hydro_data/event_char/spei24.csv') %>%
  select(contains(c('10260015', '10270104')))

# * Correlated Predictors -------------------------------------------------

correlated_predictors <- unlist(as.vector(read_csv('./DataFiles/hydro_data/event_char/correlated_predictors.csv'))) 

#* Combine ----------------------------------------------

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
  spei1,
  spei3,
  spei6,
  spei9,
  spei12,
  spei18,
  spei24
) %>%
  mutate(wateryear = factor(wateryear),
         month = factor(month),
         season = factor(season)) %>%
  select(!contains(correlated_predictors[correlated_predictors != 'reservoir_runoff_ratio']))

rm(CQ_summary, 
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
   spei1,
   spei3,
   spei6,
   spei9,
   spei12,
   spei18,
   spei24,
   correlated_predictors)

# Helpful functions -------------------------------------------------------

factor_mode <- function(x) {
  names(sort(table(x), decreasing = TRUE))[1]
}


# Cluster by CQ Indices --------------------------------------------------

# * create clusters -------------------------------------------------------

hclust_data <- event_summary %>%
  filter(complete.cases(.)) %>%
  select(FI, HI)

dist <- daisy(hclust_data, metric = 'euclidean')

cls <- agnes(dist, method = 'ward')

#fviz_nbclust(hclust_data, hcut, method = 'wss', k.max = 20)

# * for loop for calculating number of significant characteristics and average cluster size for each  number of clusters
# k_summary <- tibble(
#   k = seq(2, 25),
#   diff_vars = NA,
#   avg_n = NA,
#   sd_n = NA,
# )
# 
# for(k in k_summary$k){
k <- 8
# cls_id <- kmeans(hclust_data, centers = k, nstart = 10)$cluster
cls_id <- cutree(cls, k = k)
#write_csv(as_tibble(cls_id), 'C:/School/SAFE KAW/Data/DataFiles/cluster_data/cluster_assignment.csv')

#Add cluster ID back in
event_summary_cluster <- event_summary %>%
  filter(complete.cases(.)) %>%
  mutate(cluster = cls_id,
         season = factor(season, levels = c('Winter', 'Spring', 'Summer', 'Fall')))

#summary statistics of variables for each cluster
cluster_summary <- event_summary_cluster %>%
  group_by(cluster) %>%
  summarise(across(where(is.numeric), list(mean = mean)),
            across(where(is.factor), factor_mode))


# #ALTERNATIVE: MANUAL CLUSTERING
# #QUADRANTS + MIDDLE
# k <- 5
# event_summary_cluster <- event_summary_trim %>%
#   filter(complete.cases(.)) %>%
#   mutate(cluster = if_else(FI < 0 & HI >= 0, 1,
#                            if_else(FI >= 0 & HI >= 0, 2,
#                                    if_else(FI < 0 & HI < 0, 3, 4)))) %>%
#   mutate(cluster = if_else(FI^2 + HI^2 <= 0.2^2, 5, cluster),
#          season = factor(season, levels = c('Winter', 'Spring', 'Summer', 'Fall')))
# 
#9 WAYS
# k <- 9
# cut <- 0.2
# event_summary_cluster <- event_summary_trim %>%
#   filter(complete.cases(.)) %>%
#   mutate(cluster = if_else(FI < cut & HI > cut, 1,
#                            if_else(abs(FI) < cut & HI > cut, 2,
#                                    if_else(FI > cut & HI > cut, 3,
#                                            if_else(FI < cut & abs(HI) < cut, 4,
#                                                    if_else(FI > cut & abs(HI) < cut, 6,
#                                                            if_else(FI < cut & HI < cut, 7,
#                                                                    if_else(abs(FI) < cut & HI < cut, 8,
#                                                                            if_else(FI > cut & HI < cut, 9, NA))))))))) %>%
#   mutate(cluster = if_else(abs(FI) < cut & HI > cut, 2, cluster),
#          cluster = if_else(abs(FI) < cut & HI < cut, 8, cluster),
#          cluster = if_else(abs(FI) < cut & abs(HI) < cut, 5, cluster),
#          season = factor(season, levels = c('Winter', 'Spring', 'Summer', 'Fall')))



#FI HI scatter
ggplot(data = event_summary_cluster) +
  geom_point(aes(x = FI, y = HI, color = factor(cluster))) +
  scale_color_discrete(name = 'Cluster') +
  xlab('') +
  ylab('') +
  xlim(-1,1) +
  ylim(-1,1) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0)


# * k-s test for sig different vars ---------------------------------------
event_cdfs <- event_summary_cluster %>%
  select(!c(FI, HI, where(is.factor))) %>%
  pivot_longer(!cluster, names_to = 'var', values_to = 'value') %>%
  group_by(cluster, var) %>%
  arrange((value)) %>%
  mutate(rank = seq(1:length(value)),
         P = rank/(length(value) + 1)*100) %>%
  ungroup()

combs <- combn(k, 2)

cluster_ks <- tibble(
  var = sort(unique(event_cdfs$var))
)

for(comb in seq(ncol(combs))){
comb_ks <- event_cdfs %>%
  group_by(var) %>%
  summarize(ks = ks.test(value[cluster == combs[1,comb]], value[cluster == combs[2,comb]])$statistic,
            sig = ks.test(value[cluster == combs[1,comb]], value[cluster == combs[2,comb]])$p.value)

cluster_ks[,comb+1] <- comb_ks[,3]
colnames(cluster_ks)[comb+1] <- paste0('sig_',combs[1,comb],'_',combs[2,comb])
}

#All significant pairwise differences
diff_vars <- cluster_ks %>%
  pivot_longer(!var, names_to = 'combo', values_to = 'sig') %>%
  filter(sig <= 0.05)

same_vars <- cluster_ks %>%
  pivot_longer(!var, names_to = 'combo', values_to = 'sig') %>%
  filter(sig >= 0.95)

#Number of significant differencess per var
diff_vars_obs <- diff_vars %>%
  count(var) %>%
  arrange(desc(n))

same_vars_obs <- same_vars %>%
  count(var) %>%
  arrange(desc(n))


#cluster size check
clst_obs <- event_summary_cluster %>%
  count(cluster)

# k_summary$diff_vars[k-1] <- nrow(diff_vars)
# k_summary$avg_n[k-1] <- mean(clst_obs$n)
# k_summary$sd_n[k-1] <- sqrt(var(clst_obs$n))
# 
# print(paste0('k = ', k, '/', max(k_summary$k), ' finished'))
# }
# write_csv(k_summary, 'C:/School/SAFE KAW/Data/DataFiles/cluster_data/k_summary_ks/k_summary_kmeans.csv')
# write_csv(diff_vars_obs, 'C:/School/SAFE KAW/Data/DataFiles/cluster_data/k_summary_ks/sig_vars_kmeans.csv')
# 
# ggplot(data = k_summary) +
#    geom_point(aes(x = k, y = diff_vars/k))


#Rank vars by number of pairwise differences per cluster
diff_vars_clust_summary_var <- tibble()
diff_vars_clust_summary_n <- tibble()
for(cluster in seq(k)) {
cluster_char <- as.character(cluster)
diff_vars_clust <- cluster_ks %>%
  select(var, c(contains(paste0('_',cluster_char,'_')), ends_with(paste0('_',cluster_char)))) %>%
  pivot_longer(-var, names_to = 'combo', values_to = 'sig') %>%
  filter(sig <= 0.05) %>%
  count(var) %>%
  arrange(desc(n))
diff_vars_clust_summary_var[1:nrow(diff_vars_clust),cluster] <- diff_vars_clust$var
diff_vars_clust_summary_n[1:nrow(diff_vars_clust),cluster] <- diff_vars_clust$n
}

#box plots of significantly different vars, remove outliers
diff_cutoff <- 2
event_summary_cluster %>%
  select(all_of(diff_vars_obs$var[diff_vars_obs$n >= diff_cutoff]), cluster) %>%
  pivot_longer(!cluster, names_to = 'var', values_to = 'value') %>%
  group_by(var) %>%
  filter(value <= IQR(value)+quantile(value, 0.75) &
           value >= quantile(value, 0.25)-IQR(value)) %>%
  ggplot() +
  geom_boxplot(aes(y = value, x = cluster, fill = factor(cluster))) +
  scale_fill_discrete(guide = 'none') +
  facet_wrap(vars(var), scales = 'free_y')

#CDFs of significantly different vars
ggplot(data = subset(event_cdfs, var %in% unique(diff_vars$var))) +
  geom_line(aes(x = value, y = P, group = cluster, color = factor(cluster))) +
  facet_wrap(vars(var), scales = 'free_x')

ggplot(data = subset(event_cdfs, var %in% unique(diff_vars$var))) +
  geom_line(aes(x = value, y = P, group = cluster, color = factor(cluster))) +
  facet_wrap(vars(var), scales = 'free_x')

#Distributions of season across each cluster
ggplot(data = event_summary_cluster) +
  geom_bar(aes(x = season, fill = season), color = 'black') +
  scale_fill_manual(values = c('#8DA9C4','#60992D','#C33C54','#F4B860'), guide = 'none') +
  xlab('') +
  facet_wrap(vars(cluster))
ggplot(data = event_summary_cluster) +
  geom_bar(aes(x = season, fill = season), color = 'black') +
  scale_fill_manual(values = c('#8DA9C4','#60992D','#C33C54','#F4B860'), guide = 'none')

ggplot(data = event_summary_cluster) +
  geom_bar(aes(x = wateryear, fill = wateryear), color = 'black') +
  xlab('') +
  facet_wrap(vars(cluster))
ggplot(data = event_summary_cluster) +
  geom_bar(aes(x = wateryear, fill = wateryear), color = 'black')



# * Using other tests for  sig different vars ----------------------------------------
library(brunnermunzel)

# dist <- daisy(hclust_data, metric = 'gower')
# cls <- agnes(dist, method = "ward")
# 
# # * for loop for calculating number of significant characteristics and average cluster size for each  number of clusters
# k_summary <- tibble(
#   k = seq(2, 25),
#   diff_vars = NA,
#   avg_n = NA,
#   sd_n = NA,
# )
# 
# for(k in k_summary$k){
# k <- 8
# # cls_id <- kmeans(hclust_data, centers = k, nstart = 10)$cluster
# cls_id <- cutree(cls, k = k)
# 
# event_summary_cluster <- event_summary %>%
#   filter(complete.cases(.)) %>%
#   mutate(cluster = cls_id,
#          season = factor(season, levels = c('Winter', 'Spring', 'Summer', 'Fall')))

event_cdfs <- event_summary_cluster %>%
  select(!c(FI, HI, where(is.factor))) %>%
  pivot_longer(!cluster, names_to = 'var', values_to = 'value') %>%
  group_by(cluster, var) %>%
  arrange((value)) %>%
  mutate(rank = seq(1:length(value)),
         P = rank/(length(value) + 1)*100) %>%
  ungroup()

combs <- combn(k, 2)

cluster_bm <- tibble(
  var = sort(unique(event_cdfs$var))
)

for(comb in seq(ncol(combs))){
  comb_bm <- event_cdfs %>%
    group_by(var) %>%
    summarize(bm = brunnermunzel.test(value[cluster == combs[1,comb]], value[cluster == combs[2,comb]])$statistic,
              sig = brunnermunzel.test(value[cluster == combs[1,comb]], value[cluster == combs[2,comb]])$p.value)
  
  cluster_bm[,comb+1] <- comb_bm[,3]
  colnames(cluster_bm)[comb+1] <- paste0('sig_',combs[1,comb],'_',combs[2,comb])
}

#All significant pairwise differences
diff_vars <- cluster_bm %>%
  pivot_longer(!var, names_to = 'combo', values_to = 'sig') %>%
  filter(sig <= 0.05)

clst_obs <- event_summary_cluster %>%
  count(cluster)

# k_summary$diff_vars[k-1] <- nrow(diff_vars)
# k_summary$avg_n[k-1] <- mean(clst_obs$n)
# k_summary$sd_n[k-1] <- sqrt(var(clst_obs$n))
# 
# print(paste0('k = ', k, '/', max(k_summary$k), ' finished'))
# }
# write_csv(k_summary, 'C:/School/SAFE KAW/Data/DataFiles/cluster_data/k_summary_bm/k_summary_kmeans.csv')
# 
# ggplot(data = k_summary) +
#    geom_point(aes(x = k, y = diff_vars/k))


#Number of significant differencess per var
diff_vars_obs <- diff_vars %>%
  count(var) %>%
  arrange(desc(n))

diff_vars_clust_summary_var <- tibble()
diff_vars_clust_summary_n <- tibble()
for(cluster in seq(k)) {
  cluster_char <- as.character(cluster)
  diff_vars_clust <- cluster_ks %>%
    select(var, c(contains(paste0('_',cluster_char,'_')), ends_with(paste0('_',cluster_char)))) %>%
    pivot_longer(-var, names_to = 'combo', values_to = 'sig') %>%
    filter(sig <= 0.05) %>%
    count(var) %>%
    arrange(desc(n))
  diff_vars_clust_summary_var[1:nrow(diff_vars_clust),cluster] <- diff_vars_clust$var
  diff_vars_clust_summary_n[1:nrow(diff_vars_clust),cluster] <- diff_vars_clust$n
}





# Look at precips more closely --------------------------------------------
huc8_centroids <- read_csv('C:/School/SAFE KAW/Data/DataFiles/shapefiles/created/huc8_centroids.csv')
precip_order_x <- paste0('p_', as.character(huc8_centroids$huc8), '_30d_ante')
precip_order_y <- paste0('p_', as.character(arrange(huc8_centroids, y)$huc8), '_1d_total')


#box plots of just precip
event_summary_cluster %>%
  select(starts_with('d_'), ends_with('_ante'), cluster) %>%
  pivot_longer(!cluster, names_to = 'var', values_to = 'value') %>%
  group_by(var) %>%
  filter(value <= IQR(value)+quantile(value, 0.75) &
           value >= quantile(value, 0.25)-IQR(value)) %>%
  ggplot() +
  geom_boxplot(aes(y = value, x = cluster, fill = factor(cluster))) +
  scale_fill_discrete(guide = 'none') +
  facet_wrap(vars(var), scales = 'free_y')

#box plots of precip by var instead of by cluster
event_summary_cluster %>%
  select(ends_with('30d_ante'), cluster) %>%
  pivot_longer(!cluster, names_to = 'var', values_to = 'value') %>%
  group_by(var) %>%
  filter(value <= IQR(value)+quantile(value, 0.75) &
           value >= quantile(value, 0.25)-IQR(value)) %>%
  ggplot() +
  geom_boxplot(aes(y = value, x = factor(var, levels = precip_order_x), fill = factor(var, levels = precip_order_x))) +
  scale_fill_viridis_d(guide = 'none') +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  facet_wrap(vars(cluster), scales = 'fixed')

#Compare precips between vars to see if each HUC8 is basically giving the same stuff
event_cdfs <- event_summary_cluster %>%
  select(!c(FI, HI, where(is.factor))) %>%
  pivot_longer(!cluster, names_to = 'var', values_to = 'value') %>%
  group_by(cluster, var) %>%
  arrange((value)) %>%
  mutate(rank = seq(1:length(value)),
         P = rank/(length(value) + 1)*100) %>%
  ungroup()

precip_cdfs <- event_cdfs %>%
  filter(str_detect(var, '1d_total'))

combs_precip <- combn(unique(precip_cdfs$var), 2)

precip_ks <- tibble(
  cluster = seq_len(k)
)

for(comb in seq(ncol(combs_precip))){
  comb_ks <- precip_cdfs %>%
    group_by(cluster) %>%
    summarize(ks = ks.test(value[var == combs_precip[1,comb]], value[var == combs_precip[2,comb]])$statistic,
              sig = ks.test(value[var == combs_precip[1,comb]], value[var == combs_precip[2,comb]])$p.value)
  
  precip_ks[,comb+1] <- comb_ks[,3]
  colnames(precip_ks)[comb+1] <- paste0('sig_',combs_precip[1,comb],'_',combs_precip[2,comb])
}

diff_precips <- precip_ks %>%
  pivot_longer(!cluster, names_to = 'combo', values_to = 'sig') %>%
  filter(sig <= 0.05)

same_precips <- precip_ks %>%
  pivot_longer(!cluster, names_to = 'combo', values_to = 'sig') %>%
  filter(sig >= 0.95)

same_precips_summary_clust <- tibble()
same_precips_summary_n <- tibble()
for(variable in unique(precip_cdfs$var)) {
  same_precips_clust <- precip_ks %>%
    select(cluster, contains(variable)) %>%
    pivot_longer(-cluster, names_to = 'combo', values_to = 'sig') %>%
    filter(sig >= 0.95) %>%
    count(cluster) %>%
    arrange(desc(n))
  same_precips_summary_clust[1:nrow(same_precips_clust),variable] <- same_precips_clust$cluster
  same_precips_summary_n[1:nrow(same_precips_clust),variable] <- same_precips_clust$n
}

# Compare each cluster's spatial precip distribution to the combined set to see which is diff

all_clusters_dist <- event_summary_cluster %>%
  select(ends_with('1d_total'))

spatial_dist_ks <- tibble(
  cluster = seq(k)
)
for(cluster_sel in seq(k)){
  cluster_dist <- event_summary_cluster %>%
    filter(cluster == cluster_sel) %>%
    select(ends_with('1d_total'))
  for(huc8 in seq(ncol(all_clusters_dist))){
    spatial_dist_ks[cluster_sel,huc8+1] <- ks.test(cluster_dist[,huc8], all_clusters_dist[,huc8])$p.value
  }
}
colnames(spatial_dist_ks)[2:ncol(spatial_dist_ks)] <- colnames(all_clusters_dist)

diff_spatial_precip <- spatial_dist_ks %>%
  pivot_longer(!cluster, names_to = 'huc8', values_to = 'sig') %>%
  filter(sig <= 0.05)

event_summary_cluster %>%
  select(ends_with('1d_total'), cluster) %>%
  pivot_longer(!cluster, names_to = 'var', values_to = 'value') %>%
  group_by(var) %>%
  filter(value <= IQR(value)+quantile(value, 0.75) &
           value >= quantile(value, 0.25)-IQR(value)) %>%
  filter(cluster == 8) %>%
  ggplot() +
  geom_boxplot(aes(y = value, x = factor(var, levels = precip_order_x), fill = factor(var, levels = precip_order_x))) +
  scale_fill_viridis_d(guide = 'none') +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  facet_wrap(vars(cluster), scales = 'fixed')

event_summary_cluster %>%
  select(ends_with('1d_total'), cluster) %>%
  pivot_longer(!cluster, names_to = 'var', values_to = 'value') %>%
  group_by(var) %>%
  filter(value <= IQR(value)+quantile(value, 0.75) &
           value >= quantile(value, 0.25)-IQR(value)) %>%
  ggplot() +
  geom_boxplot(aes(y = value, x = factor(var, levels = precip_order_x), fill = factor(var, levels = precip_order_x))) +
  scale_fill_viridis_d(guide = 'none') +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

# * summary plots of clustered vars ---------------------------------------

#Distributions of season across each cluster
ggplot(data = event_summary_cluster) +
  geom_bar(aes(x = season, fill = season), color = 'black') +
  scale_fill_manual(values = c('#8DA9C4','#60992D','#C33C54','#F4B860'), guide = 'none') +
  xlab('') +
  facet_wrap(vars(cluster))
ggplot(data = event_summary_cluster) +
  geom_bar(aes(x = season, fill = season), color = 'black') +
  scale_fill_manual(values = c('#8DA9C4','#60992D','#C33C54','#F4B860'), guide = 'none')


#box plots by cluster
event_summary_cluster %>%
  select(!c(FI, HI, where(is.factor))) %>%
  pivot_longer(!cluster, names_to = 'var', values_to = 'value') %>%
  ggplot() +
  geom_boxplot(aes(y = value, group = cluster, fill = factor(cluster))) +
  scale_fill_discrete(guide = 'none') +
  facet_wrap(vars(var), scales = 'free_y')

#remove outliers
event_summary_cluster %>%
  select(!c(FI, HI, where(is.factor))) %>%
  pivot_longer(!cluster, names_to = 'var', values_to = 'value') %>%
  group_by(var) %>%
  filter(value <= IQR(value)+quantile(value, 0.75) &
           value >= quantile(value, 0.25)-IQR(value)) %>%
  ggplot() +
  geom_boxplot(aes(y = value, x = cluster, fill = factor(cluster))) +
  scale_fill_discrete(guide = 'none') +
  facet_wrap(vars(var), scales = 'free_y')

#Single variable
event_summary_cluster %>%
  select(!c(FI, HI, where(is.factor))) %>%
  select(cluster, max_precip_x) %>%
  pivot_longer(!cluster, names_to = 'var', values_to = 'value') %>%
  group_by(var) %>%
  filter(value <= IQR(value)+quantile(value, 0.75) &
           value >= quantile(value, 0.25)-IQR(value)) %>%
  ggplot() +
  geom_boxplot(aes(y = value, x = cluster, fill = factor(cluster))) +
  scale_fill_discrete(guide = 'none') +
  scale_x_continuous(breaks = seq(6))

#cdfs by cluster
event_cdfs <- event_summary_cluster %>%
  select(!c(FI, HI, where(is.factor))) %>%
  pivot_longer(!cluster, names_to = 'var', values_to = 'value') %>%
  group_by(cluster, var) %>%
  arrange((value)) %>%
  mutate(rank = seq(1:length(value)),
         P = rank/(length(value) + 1)*100) %>%
  ungroup()

ggplot(data = event_cdfs) +
  geom_line(aes(x = value, y = P, group = cluster, color = factor(cluster))) +
  facet_wrap(vars(var), scales = 'free_x')


# Cluster by Event Characteristics -------------------------------------------------

hclust_data <- event_summary %>%
  filter(complete.cases(.)) %>%
  select(d_10260015_spei12, reservoir_precip_ratio) %>%
  filter(reservoir_precip_ratio < 2)

dist <- daisy(hclust_data, metric = 'gower')

cls <- agnes(dist, method = "ward")

fviz_nbclust((hclust_data %>% select(where(is.numeric))), hcut, method = 'wss')

k <- 4
cls_id <- cutree(cls, k = k)

#Add cluster ID back in
event_summary_cluster <- event_summary %>%
  filter(complete.cases(.)) %>%
  mutate(cluster = cls_id)

cluster_summary <- event_summary_cluster %>%
  group_by(cluster) %>%
  summarise(across(where(is.numeric), list(mean = mean)),
            across(where(is.factor), factor_mode))

ggplot(data = event_summary_cluster) +
  geom_point(aes(x = FI, y = HI, color = factor(cluster))) +
  geom_point(data = cluster_summary, aes(x = FI_mean, y = HI_mean, color = factor(cluster)), size = 3, shape = 15) +
  xlim(-1,1) +
  ylim(-1,1) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0)

ggplot(data = event_summary_cluster) +
  geom_point(aes(x = FI, y = HI, color = factor(cluster))) +
  geom_point(data = cluster_summary, aes(x = FI_mean, y = HI_mean, color = factor(cluster)), size = 3, shape = 15) +
  xlim(-1,1) +
  ylim(-1,1) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  facet_wrap(vars(cluster))

ggplot(data = event_summary_cluster) +
  geom_point(aes(x = FI, y = HI, color = factor(season))) +
  xlim(-1,1) +
  ylim(-1,1) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  facet_wrap(vars(season))

event_cdfs <- event_summary_cluster %>%
  select(cluster, FI, HI) %>%
  pivot_longer(!cluster, names_to = 'var', values_to = 'value') %>%
  group_by(cluster, var) %>%
  arrange((value)) %>%
  mutate(rank = seq(1:length(value)),
         P = rank/(length(value) + 1)*100) %>%
  ungroup()

combs <- combn(k, 2)

cluster_ks <- tibble(
  var = sort(unique(event_cdfs$var))
)

for(comb in seq(ncol(combs))){
  comb_ks <- event_cdfs %>%
    group_by(var) %>%
    summarize(ks = ks.test(value[cluster == combs[1,comb]], value[cluster == combs[2,comb]])$statistic,
              sig = ks.test(value[cluster == combs[1,comb]], value[cluster == combs[2,comb]])$p.value)
  
  cluster_ks[,comb+1] <- comb_ks[,3]
  colnames(cluster_ks)[comb+1] <- paste0('sig_',combs[1,comb],'_',combs[2,comb])
}

#All significant pairwise differences
diff_vars <- cluster_ks %>%
  pivot_longer(!var, names_to = 'combo', values_to = 'sig') %>%
  filter(sig <= 0.05)

same_vars <- cluster_ks %>%
  pivot_longer(!var, names_to = 'combo', values_to = 'sig') %>%
  filter(sig >= 0.95)

#Number of significant differencess per var
diff_vars_obs <- diff_vars %>%
  count(var) %>%
  arrange(desc(n))

same_vars_obs <- same_vars %>%
  count(var) %>%
  arrange(desc(n))

event_summary_cluster %>%
  select(FI, HI, cluster) %>%
  pivot_longer(!cluster, names_to = 'var', values_to = 'value') %>%
  group_by(var) %>%
  ggplot() +
  geom_boxplot(aes(y = value, x = cluster, fill = factor(cluster))) +
  scale_fill_discrete(guide = 'none') +
  facet_wrap(vars(var), scales = 'free_y')




# cluster combinations test -----------------------------------------------

k_combs <- tibble(
    k = seq(2, 25),
    diff_vars = NA
  )
chars <- 20

for(k_hyp in k_combs$k){
  k_hyp <- k_hyp
  combos <- combn(k_hyp,2)
  k_combs$diff_vars[k_hyp-1] <- ncol(combos)*chars
}

k_combs <- k_combs %>%
  mutate(slope = diff_vars/k,
         formula = factorial(k)/(factorial(2)*factorial(k-2)))

ggplot(data = k_combs) +
  geom_point(aes(x = k, y = diff_vars/k))

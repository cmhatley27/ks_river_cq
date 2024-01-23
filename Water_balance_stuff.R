
# Load libraries ----------------------------------------------------------
library(tidyverse)
library(cluster)  
library(factoextra)
library(reshape)
library(brunnermunzel)
library(ggpubr)

# Load data ---------------------------------------------------------------
# Basic characteristics 

CQ_summary <- read_csv('./DataFiles/hydro_data/event_char/CQ_summary.csv') %>%
  select(FI, HI = HI_mean)

simple_temporal <- read_csv('./DataFiles/hydro_data/event_char/simple_temporal.csv')
simple_size <- read_csv('./DataFiles/hydro_data/event_char/simple_size.csv')

N_yield <- read_csv('./DataFiles/hydro_data/event_char/N_yield.csv') %>%
  cbind(., CQ_summary)

#* Gridded Precip Predictors 
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


# * Outflows Predictors 
outflows_inevent_totals <- read_csv('./DataFiles/hydro_data/event_char/outflows_inevent_totals.csv') %>%
  select(contains('1d'))

outflows_deltas <- read_csv('./DataFiles/hydro_data/event_char/outflows_deltas.csv') %>%
  select(contains('1d'))

flow_ratios <- read_csv('./DataFiles/hydro_data/event_char/flow_ratios.csv')


# Drought Predictors 
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
         season = factor(season))


hclust_data <- event_summary %>%
  filter(complete.cases(.)) %>%
  select(FI, HI)

dist <- daisy(hclust_data, metric = 'euclidean')

cls <- agnes(dist, method = "ward")

k <- 8

cls_id <- cutree(cls, k = k)

#Add cluster ID back in and add N yield data
event_summary_cluster <- event_summary %>%
  filter(complete.cases(.)) %>%
  mutate(cluster = cls_id,
         season = factor(season, levels = c('Winter', 'Spring', 'Summer', 'Fall'))) %>%
  left_join(., filter(N_yield, !is.na(initial_N) & !is.na(HI))) %>%
  mutate(loadrate = N_load/duration,
         log_loadrate = log10(loadrate),
         initial_loadrate = initial_Q*initial_N*60*60*24/1000,
         log_initial_loadrate = log10(initial_loadrate),
         delta_N = (max_N - initial_N) - (initial_N - min_N),
         range_N = max_N - min_N,
         new_HI = if_else(FI < 0, HI*-1, HI))


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

# Calculate water balance data -----------------------------------------------------

#EKSRB area in m^2
library(sf)
huc8_boundaries <- st_read('C:/School/SAFE KAW/Data/DataFiles/shapefiles/created/huc8_boundaries.shp') 
huc8_areas <- tibble(HUC8 = huc8_boundaries$HUC8,
                     area_sqkm = huc8_boundaries$AREASQKM) %>%
  filter(HUC8 %in% c('10260008', '10260010', '10260015', '10270101', '10270102', '10270104'))
area <- sum(huc8_areas$area_sqkm)*1000000

#water balance in mm
water_bal <- event_summary_cluster %>%
  #select(FI, HI, cluster, total_Q, p_unint_vol, r_sum_1d_total, reservoir_runoff_ratio, reservoir_precip_ratio, delta_N, loadrate) %>%
  mutate(p_depth = p_unint_vol/area*1000,
         bal = (r_sum_1d_total + p_unint_vol - total_Q)/area*1000,
         bal_onlyprecip = (p_unint_vol - total_Q)/area*1000,
         bal_overprecip = bal/(p_unint_vol/area*1000),
         bal_overinputs = bal/((p_unint_vol+r_sum_1d_total)/area*1000),
         runoff_ratio = event_summary_cluster$runoff_precip_ratio,
         spei24_w = event_summary_cluster$d_10260015_spei24,
         spei24_e = event_summary_cluster$d_10270104_spei24,
         spei24_mean = (spei24_w + spei24_e)/2,
         spei18_w = event_summary_cluster$d_10260015_spei18,
         spei18_e = event_summary_cluster$d_10270104_spei18,
         spei18_mean = (spei18_w + spei18_e)/2,
         spei12_w = event_summary_cluster$d_10260015_spei12,
         spei12_e = event_summary_cluster$d_10270104_spei12,
         spei12_mean = (spei12_w + spei12_e)/2,
         spei9_w = event_summary_cluster$d_10260015_spei9,
         spei9_e = event_summary_cluster$d_10270104_spei9,
         spei9_mean = (spei9_w + spei9_e)/2,
         spei6_w = event_summary_cluster$d_10260015_spei6,
         spei6_e = event_summary_cluster$d_10270104_spei6,
         spei6_mean = (spei6_w + spei6_e)/2,
         spei3_w = event_summary_cluster$d_10260015_spei3,
         spei3_e = event_summary_cluster$d_10270104_spei3,
         spei3_mean = (spei3_w + spei3_e)/2,
         spei1_w = event_summary_cluster$d_10260015_spei1,
         spei1_e = event_summary_cluster$d_10270104_spei1,
         spei1_mean = (spei1_w + spei1_e)/2,
         log_loadrate = event_summary_cluster$log_loadrate,
         ap_15d_median = unlist((select(event_summary_cluster, ends_with('15d_ante')) %>%
                                   mutate(event = seq(1:nrow(event_summary_cluster))) %>%
                                   pivot_longer(!event, names_to = 'huc8', values_to = 'ap') %>%
                                   group_by(event) %>%
                                   summarise(median = median(ap)))[,2]),
         ap_30d_median = unlist((select(event_summary_cluster, ends_with('30d_ante')) %>%
                                   mutate(event = seq(1:nrow(event_summary_cluster))) %>%
                                   pivot_longer(!event, names_to = 'huc8', values_to = 'ap') %>%
                                   group_by(event) %>%
                                   summarise(median = median(ap)))[,2]),
         ap_60d_median = unlist((select(event_summary_cluster, ends_with('60d_ante')) %>%
                                   mutate(event = seq(1:nrow(event_summary_cluster))) %>%
                                   pivot_longer(!event, names_to = 'huc8', values_to = 'ap') %>%
                                   group_by(event) %>%
                                   summarise(median = median(ap)))[,2]),
         ap_90d_median = unlist((select(event_summary_cluster, ends_with('90d_ante')) %>%
                                   mutate(event = seq(1:nrow(event_summary_cluster))) %>%
                                   pivot_longer(!event, names_to = 'huc8', values_to = 'ap') %>%
                                   group_by(event) %>%
                                   summarise(median = median(ap)))[,2]),
         ap_180d_median = unlist((select(event_summary_cluster, ends_with('180d_ante')) %>%
                                   mutate(event = seq(1:nrow(event_summary_cluster))) %>%
                                   pivot_longer(!event, names_to = 'huc8', values_to = 'ap') %>%
                                   group_by(event) %>%
                                   summarise(median = median(ap)))[,2]),
         ap_365d_median = unlist((select(event_summary_cluster, ends_with('365d_ante')) %>%
                                   mutate(event = seq(1:nrow(event_summary_cluster))) %>%
                                   pivot_longer(!event, names_to = 'huc8', values_to = 'ap') %>%
                                   group_by(event) %>%
                                   summarise(median = median(ap)))[,2]))
           

# Plot event data broken down by Climate and Reservoir info ---------------
climate_var <- 'spei12_e'
climate_name <- climate_var
climate_q <- c(0.25,0.75)

res_var <- 'reservoir_runoff_ratio'
res_name <- 'Res/Runoff'
res_q <- c(0.25,0.75)

response_vars <- c('FI', 'HI', 'log_loadrate')

climate_vals <- unlist(water_bal[climate_var])
climate_quant <- quantile(climate_vals, climate_q)
res_vals <- unlist(water_bal[res_var])
res_quant <- quantile(res_vals, res_q)


plot_dat <- water_bal %>%
  mutate(climate_end = if_else(climate_vals <= climate_quant[1], 'low', 'mid'),
         climate_end = if_else(climate_vals >= climate_quant[2], 'high', climate_end),
         climate_end = factor(climate_end, levels = c('low', 'mid', 'high')),
         res_end = if_else(res_vals <= res_quant[1], 'low', 'mid'),
         res_end = if_else(res_vals >= res_quant[2], 'high', res_end),
         res_end = factor(res_end, levels = c('low', 'mid', 'high'))) %>%
  pivot_longer(all_of(response_vars), names_to = 'var', values_to = 'val') %>%
  select(var, val, climate_end, res_end)

#Plot only climate boxes
if(climate_q[1] != 0.5) {
  ggplot(data = plot_dat, aes(x = climate_end, y = val)) +
    geom_boxplot(aes(fill = climate_end)) +
    scale_fill_manual(limits = c('low', 'mid', 'high'), values = c('coral2', '#E9ECEE', 'skyblue'), na.value = NA,
                      labels = c(paste0('< ',res_q[1]*100, 'th'), paste0(res_q[1]*100,'th - ',res_q[2]*100,'th'), paste0('> ',res_q[2]*100, 'th')),
                      name = paste0(climate_name, ' Percentile')) +
    scale_x_discrete(limits = c('low', 'mid', 'high'), labels = c(paste0('< ',climate_q[1]*100, 'th'), paste0(climate_q[1]*100,'th - ',climate_q[2]*100,'th'), paste0('> ',climate_q[2]*100, 'th')),
                     name = paste(climate_name, ' Percentile')) +
    facet_wrap(vars(var), scales = 'free_y') +
    ylab(NULL)
} else {
  ggplot(data = plot_dat, aes(x = climate_end, y = val)) +
    geom_boxplot(aes(fill = climate_end)) +
    scale_fill_manual(limits = c('low', 'mid', 'high'), values = c('coral2', '#E9ECEE', 'skyblue'), na.value = NA,
                      labels = c(paste0('< ',res_q[1]*100, 'th'), paste0(res_q[1]*100,'th - ',res_q[2]*100,'th'), paste0('> ',res_q[2]*100, 'th')),
                      name = paste0(climate_name, ' Percentile')) +
    scale_x_discrete(limits = c('low', 'high'), labels = c(paste0('< ',climate_q[1]*100, 'th'), paste0('> ',climate_q[2]*100, 'th')),
                     name = paste0(climate_name, ' Percentile')) +
    facet_wrap(vars(var), scales = 'free_y')
}

#Plot reservoir boxes
if(climate_q[1] != 0.5) {
  ggplot(data = plot_dat, aes(x = climate_end, y = val)) +
    geom_boxplot(aes(fill = climate_end, alpha = res_end)) +
    scale_fill_manual(limits = c('low', 'mid', 'high'), values = c('coral2', '#E9ECEE', 'skyblue'), na.value = NA,
                      labels = c(paste0('< ',res_q[1]*100, 'th'), paste0(res_q[1]*100,'th - ',res_q[2]*100,'th'), paste0('> ',res_q[2]*100, 'th')),
                      name = paste0(climate_name, ' Percentile')) +
    scale_x_discrete(limits = c('low', 'mid', 'high'), labels = c(paste0('< ',climate_q[1]*100, 'th'), paste0(climate_q[1]*100,'th - ',climate_q[2]*100,'th'), paste0('> ',climate_q[2]*100, 'th')),
                     name = paste(climate_name, ' Percentile')) +
    facet_wrap(vars(var), scales = 'free_y') +
    ylab(NULL)
  } else {
  ggplot(data = plot_dat, aes(x = climate_end, y = val)) +
    geom_boxplot(aes(fill = climate_end, alpha = res_end)) +
    scale_fill_manual(limits = c('low', 'mid', 'high'), values = c('coral2', '#E9ECEE', 'skyblue'), na.value = NA,
                      labels = c(paste0('< ',res_q[1]*100, 'th'), paste0(res_q[1]*100,'th - ',res_q[2]*100,'th'), paste0('> ',res_q[2]*100, 'th')),
                      name = paste0(climate_name, ' Percentile')) +
    scale_x_discrete(limits = c('low', 'high'), labels = c(paste0('< ',climate_q[1]*100, 'th'), paste0('> ',climate_q[2]*100, 'th')),
                     name = paste0(climate_name, ' Percentile')) +
    facet_wrap(vars(var), scales = 'free_y')
  }

#Calculate significant differences
climate_combs <- combn(unique(plot_dat$climate_end), 2)
climate_sig <- tibble(
  response_var = rep(response_vars, each = ncol(climate_combs)),
  climate1 = rep(climate_combs[1,], length.out = length(response_var)),
  climate2 = rep(climate_combs[2,], length.out = length(response_var)),
  p = NA,
  climate_var = climate_var
)

for(combo in 1:nrow(climate_sig)){
  dat_fil <- filter(plot_dat, var == unlist(climate_sig[combo, 'response_var']))
  
  climate_sig[combo, 'p'] <- round(brunnermunzel.test(filter(dat_fil, climate_end == unlist(climate_sig[combo,'climate1']))$val,
                                                      filter(dat_fil, climate_end == unlist(climate_sig[combo,'climate2']))$val)$p.value,
                                   digits = 3)
}


res_combs <- combn(unique(plot_dat$res_end), 2)
res_sig <- tibble(
  response_var = rep(response_vars, each = ncol(res_combs)*length(unique(plot_dat$climate_end))),
  climate = rep(rep(unique(plot_dat$climate_end), each = ncol(res_combs)), length.out = length(response_var)),
  res1 = rep(res_combs[1,], length.out = length(response_var)),
  res2 = rep(res_combs[2,], length.out = length(response_var)),
  p = NA,
  climate_var = climate_var
)

for(combo in 1:nrow(res_sig)){
  dat_fil <- filter(plot_dat, var == unlist(res_sig[combo, 'response_var']) & climate_end == unlist(res_sig[combo, 'climate']))
  
  res_sig[combo, 'p'] <- round(brunnermunzel.test(filter(dat_fil, res_end == unlist(res_sig[combo,'res1']))$val,
                                                  filter(dat_fil, res_end == unlist(res_sig[combo,'res2']))$val)$p.value,
                               digits = 3)
}

#Print distribution summary
print(climate_quant)
print(res_quant)
print(table(select(plot_dat, climate_end, res_end))/length(response_vars))


# Event Characteristics -------------------------------------------------------
climate_var <- 'p_10260015_90d_ante'
climate_name <- '12 Month SPEI (West & East Mean)'
climate_q <- c(0.25,0.75)

res_var <- 'reservoir_runoff_ratio'
res_name <- 'Res/Runoff'
res_q <- c(0.25,0.75)

response_vars <- colnames(event_summary_cluster %>%
  select(p_10270104_1d_total, duration_since_last, duration, initial_Q, delta_rat_Q, HI, 
         r_t_1d_delta, r_m_1d_delta, r_t_share, r_m_share, initial_N, p_10260015_1d_total, r_c_1d_delta, r_c_share,
         r_p_share, r_p_1d_delta, reservoir_runoff_ratio,
         FI, new_HI))

climate_vals <- unlist(water_bal[climate_var])
climate_quant <- quantile(climate_vals, climate_q)
res_vals <- unlist(water_bal[res_var])
res_quant <- quantile(res_vals, res_q)


plot_dat <- water_bal %>%
  mutate(climate_end = if_else(climate_vals <= climate_quant[1], 'low', 'mid'),
         climate_end = if_else(climate_vals >= climate_quant[2], 'high', climate_end),
         climate_end = factor(climate_end, levels = c('low', 'mid', 'high')),
         res_end = if_else(res_vals <= res_quant[1], 'low', 'mid'),
         res_end = if_else(res_vals >= res_quant[2], 'high', res_end),
         res_end = factor(res_end, levels = c('low', 'mid', 'high'))) %>%
  pivot_longer(all_of(response_vars), names_to = 'var', values_to = 'val')


if(climate_q[1] != 0.5) {
  ggplot(data = plot_dat, aes(x = climate_end, y = val)) +
    geom_boxplot(aes(fill = climate_end)) +
    scale_fill_manual(limits = c('low', 'mid', 'high'), values = c('coral2', '#E9ECEE', 'skyblue'), na.value = NA,
                      labels = c(paste0('< ',res_q[1]*100, 'th'), paste0(res_q[1]*100,'th - ',res_q[2]*100,'th'), paste0('> ',res_q[2]*100, 'th')),
                      name = paste0(climate_name, ' Percentile')) +
    scale_x_discrete(limits = c('low', 'mid', 'high'), labels = c(paste0('< ',climate_q[1]*100, 'th'), paste0(climate_q[1]*100,'th - ',climate_q[2]*100,'th'), paste0('> ',climate_q[2]*100, 'th')),
                     name = paste(climate_name, ' Percentile')) +
    facet_wrap(vars(var), scales = 'free_y') +
    ylab(NULL)
} else {
  ggplot(data = plot_dat, aes(x = climate_end, y = val)) +
    geom_boxplot(aes(fill = climate_end)) +
    scale_fill_manual(limits = c('low', 'mid', 'high'), values = c('coral2', '#E9ECEE', 'skyblue'), na.value = NA,
                      labels = c(paste0('< ',res_q[1]*100, 'th'), paste0(res_q[1]*100,'th - ',res_q[2]*100,'th'), paste0('> ',res_q[2]*100, 'th')),
                      name = paste0(res_name, ' Percentile')) +
    scale_x_discrete(limits = c('low', 'high'), labels = c(paste0('< ',climate_q[1]*100, 'th'), paste0('> ',climate_q[2]*100, 'th')),
                     name = paste0(climate_name, ' Percentile')) +
    facet_wrap(vars(var), scales = 'free_y')
}


if(climate_q[1] != 0.5) {
  ggplot(data = plot_dat, aes(x = climate_end, y = val)) +
    geom_boxplot(aes(fill = climate_end, alpha = res_end)) +
    scale_fill_manual(limits = c('low', 'mid', 'high'), values = c('coral2', '#E9ECEE', 'skyblue'), na.value = NA,
                      labels = c(paste0('< ',res_q[1]*100, 'th'), paste0(res_q[1]*100,'th - ',res_q[2]*100,'th'), paste0('> ',res_q[2]*100, 'th')),
                      name = paste0(res_name, ' Percentile')) +
    scale_x_discrete(limits = c('low', 'mid', 'high'), labels = c(paste0('< ',climate_q[1]*100, 'th'), paste0(climate_q[1]*100,'th - ',climate_q[2]*100,'th'), paste0('> ',climate_q[2]*100, 'th')),
                     name = paste(climate_name, ' Percentile')) +
    facet_wrap(vars(var), scales = 'free_y') +
    ylab(NULL)
} else {
  ggplot(data = plot_dat, aes(x = climate_end, y = val)) +
    geom_boxplot(aes(fill = climate_end, alpha = res_end)) +
    scale_fill_manual(limits = c('low', 'mid', 'high'), values = c('coral2', '#E9ECEE', 'skyblue'), na.value = NA,
                      labels = c(paste0('< ',res_q[1]*100, 'th'), paste0(res_q[1]*100,'th - ',res_q[2]*100,'th'), paste0('> ',res_q[2]*100, 'th')),
                      name = paste0(res_name, ' Percentile')) +
    scale_x_discrete(limits = c('low', 'high'), labels = c(paste0('< ',climate_q[1]*100, 'th'), paste0('> ',climate_q[2]*100, 'th')),
                     name = paste0(climate_name, ' Percentile')) +
    facet_wrap(vars(var), scales = 'free_y')
}

table(select(plot_dat, climate_end, res_end))/length(response_vars)
table(select(plot_dat, climate_end, res_end, season))/length(response_vars)


# FI HI scatter -----------------------------------------------------------
climate_var <- 'spei12_mean'
climate_name <- '12 Month SPEI (West & East Mean)'
climate_q <- c(0.25,0.75)

res_var <- 'reservoir_runoff_ratio'
res_name <- 'Res/Runoff'
res_q <- c(0.25,0.75)

climate_vals <- unlist(water_bal[climate_var])
climate_quant <- quantile(climate_vals, climate_q)
res_vals <- unlist(water_bal[res_var])
res_quant <- quantile(res_vals, res_q)


plot_dat <- water_bal %>%
  mutate(climate_end = if_else(climate_vals <= climate_quant[1], 'dry', 'avg'),
         climate_end = if_else(climate_vals >= climate_quant[2], 'wet', climate_end),
         climate_end = factor(climate_end, levels = c('dry', 'avg', 'wet')),
         res_end = if_else(res_vals <= res_quant[1], 'precip', 'mix'),
         res_end = if_else(res_vals >= res_quant[2], 'reservoir', res_end),
         res_end = factor(res_end, levels = c('precip', 'mix', 'reservoir')),
         grouping = factor(paste0(climate_end,'-',res_end),
                           levels = c('dry-reservoir', 'avg-reservoir', 'wet-reservoir',
                                      'dry-mix', 'avg-mix', 'wet-mix',
                                      'dry-precip', 'avg-precip', 'wet-precip')))

group_means <- group_by(plot_dat, grouping) %>%
  summarise(FI_mean = mean(FI),
            HI_mean = mean(HI))

ggplot(data = plot_dat, aes(x = FI, y = HI, color = grouping)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point() +
  geom_point(data = group_means, aes(x = FI_mean, y = HI_mean, color = grouping), shape = 15, size = 3) +
  facet_wrap(vars(grouping))


# Distributions within categories -----------------------------------------
climate_var <- 'spei12_w'
climate_name <- 'SPEI West'
climate_q <- c(0.25,0.75)

res_var <- 'reservoir_runoff_ratio'
res_name <- 'Res/Runoff'
res_q <- c(0.25,0.75)

climate_vals <- unlist(water_bal[climate_var])
climate_quant <- quantile(climate_vals, climate_q)
res_vals <- unlist(water_bal[res_var])
res_quant <- quantile(res_vals, res_q)

plot_dat <- water_bal %>%
  mutate(climate_end = if_else(climate_vals <= climate_quant[1], 'low', 'mid'),
         climate_end = if_else(climate_vals >= climate_quant[2], 'high', climate_end),
         climate_end = factor(climate_end, levels = c('low', 'mid', 'high')),
         res_end = if_else(res_vals <= res_quant[1], 'low', 'mid'),
         res_end = if_else(res_vals >= res_quant[2], 'high', res_end),
         res_end = factor(res_end, levels = c('low', 'mid', 'high')))

table(select(plot_dat, season, climate_end))

westplot <- ggplot(data = plot_dat, aes(x = unlist(plot_dat[res_var]), fill = climate_end)) +
  geom_histogram(binwidth = 0.1) +
  geom_vline(xintercept = res_quant) +
  scale_fill_manual(limits = c('low', 'mid', 'high'), values = c('coral2', 'grey75', 'skyblue'), na.value = NA,
                    labels = c(paste0('< ',res_q[1]*100, 'th'), paste0(res_q[1]*100,'th - ',res_q[2]*100,'th'), paste0('> ',res_q[2]*100, 'th')),
                    name = paste0(climate_name, ' Percentile')) +
  scale_x_continuous(name = 'Reservoir/Runoff Ratio') +
  ggtitle(climate_name) +
  facet_wrap(vars(factor(climate_end, labels = c('Dry', 'Avg', 'Wet'))))

ggarrange(meanplot, westplot, eastplot, nrow = 3, align = 'v')



# Significant Differences for each climate category -----------------------
climate_indicators <- colnames(select(water_bal, starts_with('spei'), ends_with('_ante'), starts_with('ap'), initial_Q))

climate_summary <- tibble()
res_summary <- tibble()

for(indicator in 1:length(climate_indicators)){
  climate_var <- climate_indicators[indicator]
  climate_q <- c(0.25,0.75)
  
  res_var <- 'reservoir_runoff_ratio'
  res_q <- c(0.25,0.75)
  
  response_vars <- c('FI', 'HI', 'log_loadrate')
  
  climate_vals <- unlist(water_bal[climate_var])
  climate_quant <- quantile(climate_vals, climate_q)
  res_vals <- unlist(water_bal[res_var])
  res_quant <- quantile(res_vals, res_q)
  

  plot_dat <- water_bal %>%
    mutate(climate_end = if_else(climate_vals <= climate_quant[1], 'low', 'mid'),
           climate_end = if_else(climate_vals >= climate_quant[2], 'high', climate_end),
           climate_end = factor(climate_end, levels = c('low', 'mid', 'high')),
           res_end = if_else(res_vals <= res_quant[1], 'low', 'mid'),
           res_end = if_else(res_vals >= res_quant[2], 'high', res_end),
           res_end = factor(res_end, levels = c('low', 'mid', 'high'))) %>%
    pivot_longer(all_of(response_vars), names_to = 'var', values_to = 'val') %>%
    select(var, val, climate_end, res_end)
  
  #Calculate significant differences
  climate_combs <- combn(unique(plot_dat$climate_end), 2)
  climate_sig <- tibble(
    response_var = rep(response_vars, each = ncol(climate_combs)),
    climate1 = rep(climate_combs[1,], length.out = length(response_var)),
    climate2 = rep(climate_combs[2,], length.out = length(response_var)),
    p = NA,
    stat = NA,
    climate_var = climate_var
  )
  
  for(combo in 1:nrow(climate_sig)){
    dat_fil <- filter(plot_dat, var == unlist(climate_sig[combo, 'response_var']))
    test <- brunnermunzel.test(filter(dat_fil, climate_end == unlist(climate_sig[combo,'climate1']))$val,
                               filter(dat_fil, climate_end == unlist(climate_sig[combo,'climate2']))$val)
    
    climate_sig[combo, 'p'] <- round(test$p.value, digits = 3)
    climate_sig[combo, 'stat'] <- round(test$statistic, digits = 3)
  }
  
  climate_summary <- rbind(climate_summary, climate_sig)
  
  
  res_combs <- combn(unique(plot_dat$res_end), 2)
  res_sig <- tibble(
    response_var = rep(response_vars, each = ncol(res_combs)*length(unique(plot_dat$climate_end))),
    climate = rep(rep(unique(plot_dat$climate_end), each = ncol(res_combs)), length.out = length(response_var)),
    res1 = rep(res_combs[1,], length.out = length(response_var)),
    res2 = rep(res_combs[2,], length.out = length(response_var)),
    p = NA,
    stat = NA,
    climate_var = climate_var
  )
  
  for(combo in 1:nrow(res_sig)){
    dat_fil <- filter(plot_dat, var == unlist(res_sig[combo, 'response_var']) & climate_end == unlist(res_sig[combo, 'climate']))
    test <- brunnermunzel.test(filter(dat_fil, res_end == unlist(res_sig[combo,'res1']))$val,
                               filter(dat_fil, res_end == unlist(res_sig[combo,'res2']))$val)
    res_sig[combo, 'p'] <- round(test$p.value, digits = 3)
    res_sig[combo, 'stat'] <- round(test$statistic, digits = 3)
  }
  
  res_summary <- rbind(res_summary, res_sig)
}

filter(res_summary, p <= 0.05) %>%
  group_by(climate_var) %>%
  summarise(n = n()) %>%
  arrange(desc(n)) %>%
  mutate(combos = ((length(unique(climate_quant))+1)*(length(unique(res_quant))+1)*length(response_vars)),
         pct = n/combos)




# Box Plots for ALL climate variables -------------------------------------
climate_vars <- colnames(select(water_bal, starts_with('spei'), ends_with('_ante'), starts_with('ap'), initial_Q))
climate_q <- c(0.25,0.75)

res_var <- 'reservoir_runoff_ratio'
res_name <- 'Res/Runoff'
res_q <- c(0.25,0.75)

response_var <- 'FI'

plot_dat <- select(water_bal, all_of(c(response_var, res_var, climate_vars))) %>%
  dplyr::rename(response_var = response_var) %>%
  mutate(res_end = quantize(select(.,res_var), res_q)) %>%
  mutate(across(all_of(climate_vars), ~quantize(.x, climate_q))) %>%
  pivot_longer(c(all_of(climate_vars)), names_to = 'climate_var', values_to = 'climate_end')

ggplot(data = subset(plot_dat, climate_end != 'mid'), aes(x = climate_end, y = response_var)) +
  geom_boxplot(aes(fill = climate_end)) +
  facet_wrap(vars(climate_var)) +
  ylab(response_var)

ggplot(data = subset(plot_dat, climate_end != 'mid' & res_end != 'mid'), aes(x = climate_end, y = response_var)) +
  geom_boxplot(aes(fill = climate_end, alpha = res_end)) +
  facet_wrap(vars(climate_var)) +
  ylab(response_var)


# Median response variable value across ALL climate variables --------------
climate_vars <- colnames(select(water_bal, starts_with('spei') & !ends_with('_mean'), 
                                ends_with('_ante') & contains('10270104'), ends_with('_ante') & contains('10260015'),
                                initial_Q))

climate_q <- c(0.25,0.75)

res_var <- 'reservoir_runoff_ratio'
res_name <- 'Res/Runoff'
res_q <- c(0.25,0.75)

response_var <- 'log_loadrate'

plot_dat <- select(water_bal, all_of(c(response_var, res_var, climate_vars))) %>%
  dplyr::rename(response_var = response_var) %>%
  mutate(res_end = quantize(select(.,res_var), res_q)) %>%
  mutate(across(all_of(climate_vars), ~quantize(.x, climate_q))) %>%
  pivot_longer(c(all_of(climate_vars)), names_to = 'climate_var', values_to = 'climate_end')


ggplot(data = subset(plot_dat, climate_end != 'mid'), aes(x = climate_end, y = response_var)) +
  geom_boxplot(aes(fill = climate_end))

ggplot(data = subset(plot_dat, climate_end != 'mid' & res_end != 'mid'), aes(x = climate_end, y = response_var)) +
  geom_boxplot(aes(fill = climate_end, alpha = res_end)) +
  ylab(response_var)


climate_medians <- group_by(plot_dat, climate_end) %>%
  summarise(median = median(response_var, na.rm = TRUE))
res_medians <- group_by(plot_dat, climate_end, res_end) %>%
  summarise(median = median(response_var, na.rm = TRUE))

ggplot(data = subset(climate_medians,
                     climate_end != 'mid')) +
  geom_col(aes(x = climate_end, y = median, fill = climate_end), position = position_dodge())

ggplot(data = subset(res_medians,
                     climate_end != 'mid' & res_end != 'mid')) +
  geom_col(aes(x = climate_end, y = median, fill = climate_end, alpha = res_end), position = position_dodge()) +
  scale_alpha_manual(values = c(0.6,1))

#B-M tests probably not legit bc pooling events from multiple climate indicators
#means that the groupings are not independent - e.g. an event in the 'dry' group
#according to SPEI 12 will most likely also be put in the 'dry' group by 12-month AP.
#Essentially inflating the number of samples in each group to reduce the p-value.
brunnermunzel.test(plot_dat$response_var[plot_dat$climate_end == 'low'], plot_dat$response_var[plot_dat$climate_end == 'high'])$p.value
brunnermunzel.test(plot_dat$response_var[plot_dat$climate_end == 'low' & plot_dat$res_end == 'low'], 
                   plot_dat$response_var[plot_dat$climate_end == 'low' & plot_dat$res_end == 'high'])$p.value
brunnermunzel.test(plot_dat$response_var[plot_dat$climate_end == 'high' & plot_dat$res_end == 'low'], 
                   plot_dat$response_var[plot_dat$climate_end == 'high' & plot_dat$res_end == 'high'])$p.value

# Compare co-occurence of climate categories ------------------------------
ggplot(data = water_bal) +
  geom_jitter(aes(x = spei3_w, y = spei3_e, color = cut(reservoir_runoff_ratio, breaks = quantile(reservoir_runoff_ratio, c(0,0.25,0.75,1)))), width = 0.1) +
  geom_abline(slope = 1) +
  geom_hline(yintercept = quantile(water_bal$spei3_w, c(0.25,0.75))) +
  geom_vline(xintercept = quantile(water_bal$spei3_e, c(0.25,0.75))) +
  scale_color_discrete(name = 'res_cat')

table(
water_bal %>%
  mutate(res_cat = cut(water_bal$reservoir_runoff_ratio, breaks = quantile(water_bal$reservoir_runoff_ratio, c(0,0.25,0.75,1))),
         clim_cat = cut(water_bal$spei18_w, breaks = quantile(water_bal$spei18_w, c(0,0.25,0.75,1)))) %>%
  # filter(spei18_w < quantile(spei18_w, 0.75),
  #        spei1_w > quantile(spei1_w, 0.75)) %>%
  select(res_cat, clim_cat, wateryear)
)

view(
  water_bal %>%
    mutate(res_cat = cut(water_bal$reservoir_runoff_ratio, breaks = quantile(water_bal$reservoir_runoff_ratio, c(0,0.25,0.75,1)),
                         labels = c('precip','mix','res')),
           clim_cat = cut(water_bal$spei18_w, breaks = quantile(water_bal$spei18_w, c(0,0.25,0.75,1)),
                          labels = c('dry', 'avg', 'wet'))) %>%
    # filter(spei18_w < quantile(spei18_w, 0.75),
    #        spei1_w > quantile(spei1_w, 0.75)) %>%
    select(res_cat, clim_cat, wateryear, month, season, initial_Q, total_Q)
)

water_bal %>%
  mutate(res_cat = cut(water_bal$reservoir_runoff_ratio, breaks = quantile(water_bal$reservoir_runoff_ratio, c(0,0.25,0.75,1)),
                       labels = c('precip','mix','res')),
         clim_cat = cut(water_bal$spei1_w, breaks = quantile(water_bal$spei1_w, c(0,0.25,0.75,1)),
                        labels = c('dry', 'avg', 'wet'))) %>%
  # filter(spei18_w < quantile(spei18_w, 0.75),
  #        spei1_w > quantile(spei1_w, 0.75)) %>%
  select(res_cat, clim_cat, wateryear, month, season, initial_Q, total_Q, delta_Q, delta_rat_Q, duration) %>%
  group_by(clim_cat, res_cat) %>%
  summarise(across(c(ends_with('_Q'), duration), median))

select(water_bal, spei3_w, spei3_e) %>% cor(.)
summary(lm(spei3_e ~ spei3_w, data = water_bal))
           
select(water_bal, ends_with('_ante')) %>%
  select(contains('10270104'), contains('10260015')) %>%
  mutate(event = seq(1:nrow(water_bal))) %>%
  pivot_longer(!event, names_to = 'var', values_to = 'ap') %>%
  ggplot() +
  geom_histogram(aes(x = ap)) +
  facet_wrap(vars(var))




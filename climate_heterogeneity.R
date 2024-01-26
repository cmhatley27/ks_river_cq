
# Load libraries ----------------------------------------------------------
library(tidyverse)
library(cluster)  
library(factoextra)
library(reshape)
library(brunnermunzel)
library(ggpubr)
library(sf)
library(tmap)

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


# Correlations between climate and responses ------------------------------
fi_cors <- cor(select(event_summary_cluster, FI, HI, starts_with('d_'), ends_with('_ante')), use = 'pairwise.complete.obs') %>%
  melt(.) %>%
  filter(X1 == 'FI')

hi_cors <- cor(select(event_summary_cluster, FI, HI, starts_with('d_'), ends_with('_ante')), use = 'pairwise.complete.obs') %>%
  melt(.) %>%
  filter(X1 == 'HI')


# BM differences between climate/reservoir groupings --------------------------------------------
climate_vars <- colnames(select(event_summary_cluster, starts_with('d_'), ends_with('_ante')))

climate_q <- c(0.25,0.75)

res_var <- 'reservoir_runoff_ratio'
res_name <- 'Res/Runoff'
res_q <- c(0.25,0.75)

response_vars <- c('FI', 'HI', 'log_loadrate')

dat <- select(event_summary_cluster, all_of(c(response_vars, res_var, climate_vars))) %>%
  #dplyr::rename(response_var = response_var) %>%
  mutate(res_end = quantize(select(.,res_var), res_q)) %>%
  mutate(across(all_of(climate_vars), ~quantize(.x, climate_q))) %>%
  pivot_longer(c(all_of(climate_vars)), names_to = 'climate_var', values_to = 'climate_end')

bm <- dat %>%
  group_by(climate_var) %>%
  summarise(fi_climate_p = brunnermunzel.test(FI[climate_end == 'low'], FI[climate_end == 'high'])$p.value,
            fi_climate_stat = brunnermunzel.test(FI[climate_end == 'low'], FI[climate_end == 'high'])$statistic,
            hi_climate_p = brunnermunzel.test(HI[climate_end == 'low'], HI[climate_end == 'high'])$p.value,
            hi_climate_stat = brunnermunzel.test(HI[climate_end == 'low'], HI[climate_end == 'high'])$statistic,
            load_climate_p = brunnermunzel.test(log_loadrate[climate_end == 'low'], log_loadrate[climate_end == 'high'])$p.value,
            load_climate_stat = brunnermunzel.test(log_loadrate[climate_end == 'low'], log_loadrate[climate_end == 'high'])$statistic) %>%
  mutate(HUC8 = str_sub(climate_var, 3, 10),
         var = str_sub(climate_var, 12, -1))


# Map statistical significance in BM differences --------------------------

huc8s <- st_read('./DataFiles/shapefiles/created/huc8_boundaries.shp')
streams <- st_read('./DataFiles/shapefiles/created/rivers_ksrb.shp')
climate_var_order <- c('15d_ante', '30d_ante', '60d_ante', '90d_ante', '180d_ante', '365d_ante',
                       'spei1', 'spei3', 'spei6', 'spei9', 'spei12', 'spei18', 'spei24')
sig_level <- 0.05

huc8_dat <- list()
for(element in 1:(length(response_vars)*2)){
  huc8_dat[[element]] <- pivot_wider(bm, id_cols = HUC8, names_from = var, values_from = colnames(bm)[element+1]) %>%
                                   left_join(huc8s, .) %>%
                                   pivot_longer(!colnames(huc8s), names_to = 'var', values_to = 'val') %>%
                                   mutate(var = factor(var, levels = climate_var_order))
  if(str_detect(colnames(bm)[element+1], '_p')){
    huc8_dat[[element]] <- mutate(huc8_dat[[element]], sig = ifelse(val <= sig_level, TRUE, FALSE))
  }
}
names(huc8_dat) <- colnames(bm)[2:7]

tm_shape(huc8_dat$fi_climate_p) +
  tm_polygons('sig', palette = c('grey85', 'coral2')) +
  tm_facets('var')

tm_shape(huc8_dat$hi_climate_p) +
  tm_polygons('sig', palette = c('grey85', 'coral2')) +
  tm_facets('var')

tm_shape(huc8_dat$load_climate_p) +
  tm_polygons('sig', palette = c('grey85', 'coral2')) +
  tm_facets('var')

tm_shape(huc8_dat$fi_climate_stat) +
  tm_polygons('val') +
  tm_facets('var')

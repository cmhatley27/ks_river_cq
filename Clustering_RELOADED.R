
# Load libraries ----------------------------------------------------------
library(tidyverse)
library(cluster)  
library(factoextra)
library(brunnermunzel)
library(ggnewscale)

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

#EKSRB area in m^2
library(sf)
huc8_boundaries <- st_read('C:/School/SAFE KAW/Data/DataFiles/shapefiles/created/huc8_boundaries.shp') 
huc8_areas <- tibble(HUC8 = huc8_boundaries$HUC8,
                     area_sqkm = huc8_boundaries$AREASQKM) %>%
  filter(HUC8 %in% c('10260008', '10260010', '10260015', '10270101', '10270102', '10270104'))
area <- sum(huc8_areas$area_sqkm)*1000000

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
         season = factor(season),
         p_depth = p_unint_vol/area*1000,
         bal = (r_sum_1d_total + p_unint_vol - total_Q)/area*1000,
         bal_onlyprecip = (p_unint_vol - total_Q)/area*1000,
         bal_overprecip = bal/(p_unint_vol/area*1000),
         bal_overinputs = bal/((p_unint_vol+r_sum_1d_total)/area*1000)) %>%
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
   correlated_predictors,
   huc8_areas,
   huc8_boundaries)



# Cluster --------------------------------------------------------------------

hclust_data <- event_summary %>%
  filter(complete.cases(.)) %>%
  select(d_10270104_spei12, reservoir_runoff_ratio) %>%
  mutate(across(everything(), ~(.-min(.))/(max(.)-min(.))))

dist <- daisy(hclust_data, metric = 'euclidean')

cls <- agnes(dist, method = 'ward')

fviz_nbclust(hclust_data, hcut, method = 'wss', k.max = 20)

k <- 6

cls_id <- cutree(cls, k = k)

#Add cluster ID back in and add N yield data
event_summary_cluster <- event_summary %>%
  filter(complete.cases(.)) %>%
  mutate(cluster = factor(cls_id),
         season = factor(season, levels = c('Winter', 'Spring', 'Summer', 'Fall'))) %>%
  left_join(., filter(N_yield, !is.na(initial_N) & !is.na(HI))) %>%
  mutate(loadrate = N_load/duration,
         log_loadrate = log10(loadrate),
         initial_loadrate = initial_Q*initial_N*60*60*24/1000,
         log_initial_loadrate = log10(initial_loadrate),
         delta_N = (max_N - initial_N) - (initial_N - min_N))

factor_mode <- function(x) {
  names(sort(table(x), decreasing = TRUE))[1]
}
cluster_summary <- event_summary_cluster %>%
  group_by(cluster) %>%
  summarise(across(where(is.numeric), mean),
            across(where(is.factor), factor_mode)) 

ggplot(data = event_summary_cluster) +
  geom_point(aes(x = d_10270104_spei12, y = reservoir_runoff_ratio, color = cluster))

event_summary_cluster %>%
  select(!c(wateryear, month, season)) %>%
  pivot_longer(!c(d_10270104_spei12, reservoir_runoff_ratio, cluster), names_to = 'var', values_to = 'val') %>%
  filter(!is.na(val)) %>%
  filter(val <= IQR(val)+quantile(val, 0.75) &
           val >= quantile(val, 0.25)-IQR(val)) %>%
  ggplot() +
  geom_boxplot(aes(x = cluster, y = val, fill = cluster)) +
  facet_wrap(vars(var), scales = 'free_y')
  
library(tidyverse)

# Load data ---------------------------------------------------------------
# Basic characteristics 

CQ_summary <- read_csv('./DataFiles/hydro_data/other_params/CQ_summary.csv')

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

# * Cluster Assignment

clusters <- read_csv('C:/School/SAFE KAW/Data/DataFiles/cluster_data/cluster_assignment.csv')

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
  droughts,
  clusters
) %>%
  dplyr::rename(cluster = value) %>%
  mutate(wateryear = factor(wateryear),
         month = factor(month),
         season = factor(season, levels = c('Winter', 'Spring', 'Summer', 'Fall'))) %>%
  left_join(., filter(N_yield, !is.na(initial_N) & !is.na(HI_n))) %>%
  mutate(loadrate = N_load/duration,
         log_loadrate = log10(loadrate),
         initial_loadrate = initial_Q*initial_N*60*60*24/1000,
         log_initial_loadrate = log10(initial_loadrate),
         delta_N = (max_N - initial_N) - (initial_N - min_N),
         range_N = max_N - min_N)
  
#clustering data
# hclust_data <- event_summary %>%
#   select(FI_n, HI_n) %>%
#   filter(complete.cases(.))
# 
# dist <- daisy(hclust_data, metric = 'euclidean')
# 
# cls <- agnes(dist, method = "ward")
# 
# k <- 8
# kmean <- kmeans(hclust_data, centers = k)[1]
# event_ids <- mutate(hclust_data, cluster = cutree(cls, k = k))
# event_ids <- mutate(hclust_data, cluster = kmeans(hclust_data, centers = k)$cluster)
# event_summary_cluster <- left_join(event_summary, event_ids)
# ggplot(data = event_summary_cluster, aes(x = FI_n, y = HI_n, color = factor(cluster))) +
#   geom_point()
# 
rm(list = ls()[ls() != 'event_summary'])

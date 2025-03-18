library(tidyverse)

# Load data ---------------------------------------------------------------
# Basic characteristics 

CQ_summary <- read_csv('./DataFiles/CQ_summary.csv')

simple_temporal <- read_csv('./DataFiles/hydro_data/event_char/simple_temporal.csv')
simple_size <- read_csv('./DataFiles/hydro_data/event_char/simple_size.csv')

N_yield <- read_csv('./DataFiles/hydro_data/event_char/N_yield.csv') %>%
  cbind(., CQ_summary)

precip_inevent_huc8 <- read_csv('./DataFiles/hydro_data/event_char/precip_inevent_huc8.csv') %>%
  select(contains('1d'))
precip_inevent_global <- read_csv('./DataFiles/hydro_data/event_char/precip_inevent_global.csv')

precip_anteevent_huc8 <- read_csv('./DataFiles/hydro_data/event_char/precip_anteevent_huc8.csv') %>%
  select(contains('30d'))
precip_anteevent_global <- read_csv('./DataFiles/hydro_data/event_char/precip_anteevent_global.csv') 

flow_ratios <- read_csv('./DataFiles/hydro_data/event_char/flow_ratios.csv')

spei12 <- read_csv('./DataFiles/hydro_data/event_char/spei12.csv')
spei12_global <-  data.frame(d_global_spei12 = rowMeans(spei12))

clusters <- read_csv('./DataFiles/cluster_ids.csv')

# Combine -----------------------------------------------------------------
event_summary <- tibble(
  CQ_summary, 
  simple_temporal, 
  simple_size,
  precip_inevent_huc8,
  precip_inevent_global,
  precip_anteevent_huc8,
  precip_anteevent_global,
  spei12,
  spei12_global,
  flow_ratios,
  clusters
) %>%
  mutate(wateryear = factor(wateryear),
         month = factor(month),
         season = factor(season, levels = c('Winter', 'Spring', 'Summer', 'Fall'))) %>%
  filter(complete.cases(.))

# write_csv(event_summary, './DataFiles/event_summary.csv')

rm(list = ls()[ls() != 'event_summary'])
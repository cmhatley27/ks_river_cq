library(tidyverse)
library(nnet)

# Load data ---------------------------------------------------------------


# * Basic characteristics

CQ_summary <- read_csv('./DataFiles/hydro_data/event_char/CQ_summary.csv') %>%
  select(FI, HI = HI_mean)

simple_temporal <- read_csv('./DataFiles/hydro_data/event_char/simple_temporal.csv')
simple_size <- read_csv('./DataFiles/hydro_data/event_char/simple_size.csv')

N_yield <- read_csv('./DataFiles/hydro_data/event_char/N_yield.csv') %>%
  cbind(., CQ_summary)

# * Gridded Precip Predictors
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


# * Drought Predictors 
spei1 <- read_csv('./DataFiles/hydro_data/event_char/spei1.csv') %>%
  select(contains(c('10260015', '10270104')))

spei3 <- read_csv('./DataFiles/hydro_data/event_char/spei3.csv') %>%
  select(contains(c('10260015', '10270104')))

spei6 <- read_csv('./DataFiles/hydro_data/event_char/spei6.csv') %>%
  select(contains(c('10260015', '10270104')))

spei12 <- read_csv('./DataFiles/hydro_data/event_char/spei12.csv') %>%
  select(contains(c('10260015', '10270104')))

# * Correlated Predictors 

correlated_predictors <- unlist(as.vector(read_csv('./DataFiles/hydro_data/event_char/correlated_predictors.csv'))) 

# * Combine 

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
  spei12
) %>%
  mutate(wateryear = factor(wateryear),
         month = factor(month),
         season = factor(season)) %>%
  select(!contains(correlated_predictors))

# * create clusters 

hclust_data <- event_summary %>%
  filter(complete.cases(.)) %>%
  select(FI, HI)

dist <- daisy(hclust_data, metric = 'euclidean')

cls <- agnes(dist, method = "ward")

k <- 8

cls_id <- cutree(cls, k = k)

# * Add cluster ID back in
event_summary_cluster <- event_summary %>%
  filter(complete.cases(.)) %>%
  mutate(cluster = cls_id)

rm(list = ls()[ls() != 'event_summary_cluster'])

#Set Reference categories
event_summary_cluster$cluster <- relevel(factor(event_summary_cluster$cluster), ref = 8)
# event_summary_cluster$wateryear <- relevel(factor(event_summary_cluster$wateryear), ref = '2017')
event_summary_cluster$season <- relevel(event_summary_cluster$season, ref = 'Spring')



#Remove irrelevants
event_summary_cluster <- event_summary_cluster %>%
  select(!c(yday, month, wateryear))


sig_preds <- c('season', 'delta_Q', 'delta_rat_Q', 'duration', 'duration_since_last', 
               'initial_N', 'initial_Q', 'p_10260015_1d_total', 'p_10260015_maxint',
               'Eastern Precip Total' = 'p_10270104_1d_total', 'd_10270104_spei3', 'r_c_1d_delta', 'r_m_share',
               'r_t_1d_delta', 'reservoir_precip_ratio')

mod_frame <- select(event_summary_cluster, cluster, sig_preds) %>%
  mutate(initial_Q_log = log(ifelse(initial_Q <= 0, 0.01, initial_Q)),
         r_t_1d_delta_log = log(ifelse(r_t_1d_delta <= 0, 0.01, r_t_1d_delta))) %>%
  select(!c(initial_Q, r_t_1d_delta))

mod_frame$cluster <- relevel(factor(mod_frame$cluster, levels = seq(1:8)), ref = 8)
mod <- multinom(cluster ~ 1 + . + 
                  delta_rat_Q*season + delta_rat_Q*duration + p_10260015_1d_total*reservoir_precip_ratio,
                data = mod_frame)

mod_coef <- (exp(summary(mod)$coefficients))

mod_std <- summary(mod)$standard.errors
mod_wald <- abs(summary(mod)$coefficients/mod_std)
mod_wald_p <- 2*(1 - pnorm(mod_wald))
mod_sig <- ifelse(mod_coef > 0, 1, -1) * ifelse(abs(mod_wald_p) > 0.05, 0, 1)




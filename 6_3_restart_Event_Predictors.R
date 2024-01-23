library(tidyverse)
library(reshape)
library(corrplot)

# Load all relevant data --------------------------------------------------

all_hydro_15min <- with_tz(read_csv('./DataFiles/hydro_data/all_hydro_15min.csv'), 'America/Chicago')
all_hydro_daily <- read_csv('./DataFiles/hydro_data/all_hydro_daily.csv')
events <- with_tz(read_csv('DataFiles/hydro_data/event_delineations/BFLOW_events_bfi_threshold_adjusted.csv'), 'America/Chicago') %>%
  select(event_number, start_dateTime, end_dateTime)
CQ_summary <- with_tz(read_csv('./DataFiles/hydro_data/CQ_summary.csv'), 'America/Chicago') %>%
  select(FI, HI = HI_mean)
event_data <- with_tz(read_csv('./DataFiles/hydro_data/event_data/event_data.csv'), 'America/Chicago')
rise_limbs <- with_tz(read_csv('./DataFiles/hydro_data/event_data/rise_limbs.csv'), 'America/Chicago')
fall_limbs <- with_tz(read_csv('./DataFiles/hydro_data/event_data/fall_limbs.csv'), 'America/Chicago')

#*Gridded Precip Predictors -----------------------------------------------
precip_inevent_gridded <- read_csv('./DataFiles/hydro_data/event_char/precip_inevent_gridded.csv') 
precip_anteevent_gridded <- read_csv('./DataFiles/hydro_data/event_char/precip_anteevent_gridded.csv')
max_precip_coords <- read_csv('./DataFiles/hydro_data/event_char/max_precip_coords.csv')
avg_precip_intensity <- read_csv('./DataFiles/hydro_data/event_char/avg_precip_intensity.csv')
precip_unintercepted_volumes <- read_csv('./DataFiles/hydro_data/event_char/precip_unintercepted_volumes.csv')
max_precip_intensity <- read_csv('./DataFiles/hydro_data/event_char/max_precip_intensity.csv')


# Times/seasonality -------------------------------------------------------
simple_temporal <- events %>%
  mutate(
    #Water Year
    wateryear = factor(ifelse(as.integer(month(start_dateTime)) >= 10, as.integer(year(start_dateTime)) + 1, year(start_dateTime))),
    #Month
    month = factor(month(start_dateTime, label = TRUE)),
    #Season (Winter = Dec - Feb, Spring = Mar - May, etc.)
    season = factor(month(floor_date(start_dateTime, unit = 'season')), levels = c(12,3,6,9), labels = c('Winter', 'Spring', 'Summer', 'Fall')),
    #Julian day
    yday = yday(start_dateTime),
    #Event duration (in days)
    duration = time_length(as.duration(interval(start = start_dateTime, end = end_dateTime)))/(60*60*24),
    #Duration since previous event (in days)
    duration_since_last = time_length(as.duration(interval(start = lag(end_dateTime), end = start_dateTime)))/(60*60*24)
  ) %>%
  select(-c(event_number, start_dateTime, end_dateTime))

write_csv(simple_temporal, './DataFiles/hydro_data/event_char/simple_temporal.csv')


# Simple size/loads/initial conditions ---------------------------------------------------------
simple_size <- event_data %>%
  group_by(event) %>%
  summarize(
    #Initial discharge
    initial_Q = riverQ_linterp[1],
    #Maximum discharge
    max_Q = max(riverQ_linterp, na.rm = TRUE),
    #Delta discharge
    delta_Q = max_Q - initial_Q,
    #Delta percent discharge
    delta_rat_Q = delta_Q/initial_Q,
    #Total discharge
    total_Q = sum(riverQ_linterp*60*15),
    #Initial N
    initial_N = N_linterp_smooth[1]) %>%
  select(-event)

write_csv(simple_size, './DataFiles/hydro_data/event_char/simple_size.csv')

#N Yields
N_yield <- event_data %>%
  group_by(event) %>%
  summarize(
    initial_N = N_linterp_smooth[1],
    #load in kg/15min, summed across event
    N_load = sum(N_linterp_smooth*(riverQ_linterp*1000)*60*15/1000/1000),
    max_N = max(N, na.rm = TRUE),
    min_N = min(N, na.rm = TRUE)) %>%
  select(-event)
write_csv(N_yield, './DataFiles/hydro_data/event_char/N_yield.csv')

# Total outflows from each reservoir within events ------------------------

outflows <- all_hydro_15min %>%
  select(dateTime, c = clintonQ_linterp, m = milfordQ_linterp, p = perryQ_linterp, t = tuttleQ_linterp, event)

#Function for adding up 15min flows between event start minus a buffer to event end
event_flow <- function(dateTime_start, dateTime_end, buffer_length, dateTimes, flow_series){
  df <- data.frame(dateTime = dateTimes, cum_flows = cumsum(flow_series))
  flow_total <- df$cum_flows[df$dateTime == dateTime_end] - df$cum_flows[df$dateTime == dateTime_start - days(buffer_length) - minutes(15)]
  return(ifelse(date(dateTime_start) < ymd('2013-10-01'), NA, flow_total))
}

#**Loops within-event flow function over every gage and buffer length combo
within_event_flows <- function(flows_frame, buffer_lengths, event_frame = event_data){
  #extract dates
  dateTimes <- flows_frame$dateTime
  #extract just precip series
  flows <- flows_frame %>%
    select(!c(dateTime, event))
  #initialize list of output vectors
  window_totals <- list()
  #Loop through each combo of gage+buffer length and calculate cumulative precip
  for(gage in 1:ncol(flows)){
    for(buffer in 1:length(buffer_lengths)){
      summary_frame <- event_frame %>%
        group_by(event) %>%
        summarize(dateTime_start = dateTime[1],
                  dateTime_end = max(dateTime),
                  flow_in_window = event_flow(dateTime_start = dateTime_start, dateTime_end = dateTime_end, buffer_length = buffer_lengths[buffer], dateTimes = flows_frame$dateTime, flow_series = flows[[gage]])
        )
      window_totals[[paste0('r_',colnames(flows[gage]), '_', buffer_lengths[buffer], 'd_total')]] <- summary_frame$flow_in_window
    }
  }
  return(window_totals)
}

outflows_inevent_totals <- (data.frame(within_event_flows(flows_frame = outflows, buffer_lengths = c(0,1,2))) * 60 * 15) %>%
  mutate(r_sum_0d_total = r_c_0d_total + r_p_0d_total + r_m_0d_total + r_t_0d_total,
         r_sum_1d_total = r_c_1d_total + r_p_1d_total + r_m_1d_total + r_t_1d_total,
         r_sum_2d_total = r_c_2d_total + r_p_2d_total + r_m_2d_total + r_t_2d_total)

write_csv(outflows_inevent_totals, './DataFiles/hydro_data/event_char/outflows_inevent_totals.csv')


# Ratios and shares of total outflows -------------------------------------

flow_ratios <- tibble(
  runoff_precip_ratio = ifelse(event_size$total_Q/precip_unintercepted_volumes$p_unint_vol > 10, 10, event_size$total_Q/precip_unintercepted_volumes$p_unint_vol),
  reservoir_precip_ratio = ifelse(outflows_inevent_totals$r_sum_1d_total/precip_unintercepted_volumes$p_unint_vol > 10, 10, outflows_inevent_totals$r_sum_1d_total/precip_unintercepted_volumes$p_unint_vol),
  reservoir_runoff_ratio = outflows_inevent_totals$r_sum_1d_total/event_size$total_Q,
  r_c_share = outflows_inevent_totals$r_c_1d_total/outflows_inevent_totals$r_sum_1d_total,
  r_p_share = outflows_inevent_totals$r_p_1d_total/outflows_inevent_totals$r_sum_1d_total,
  r_m_share = outflows_inevent_totals$r_m_1d_total/outflows_inevent_totals$r_sum_1d_total,
  r_t_share = outflows_inevent_totals$r_t_1d_total/outflows_inevent_totals$r_sum_1d_total)

write_csv(flow_ratios, './DataFiles/hydro_data/event_char/flow_ratios.csv')

# Outflow deltas ----------------------------------------------------------

delta_flows <- function(flows_frame, buffer_lengths){
  #extract just flow series
  flows <- flows_frame %>%
    select(!event)
  max_flows <- flows_frame %>%
    group_by(event) %>%
    summarize(across(!dateTime, max)) %>%
    filter(!is.na(event))
  deltas <- list()
  ratios <- list()
  for(buffer in 1:length(buffer_lengths)){
    start_dateTimes <- flows_frame %>%
      group_by(event) %>%
      summarize(dateTime = dateTime[1] - days(buffer_lengths[buffer]) + minutes(90)) %>%
      filter(!is.na(event))
    start_flows <- start_dateTimes %>%
      left_join(flows) %>%
      select(-dateTime) %>%
      filter(!is.na(event))
    delta_frame <- max_flows - start_flows
    colnames(delta_frame) <- paste0('r_', colnames(delta_frame), '_', buffer_lengths[buffer], 'd_delta')
    ratio_frame <- delta_frame / start_flows
    colnames(ratio_frame) <- paste0(colnames(delta_frame), '_rat')
    deltas[[buffer]] <- delta_frame
    ratios[[buffer]] <- ratio_frame
  }
  results <- list(deltas, ratios)
  return(results)
}

outflows_deltas <- data.frame(delta_flows(flows_frame = outflows, buffer_lengths = c(0,1,2))) %>%
  select(!contains('event'))
write_csv(outflows_deltas, './DataFiles/hydro_data/event_char/outflows_deltas.csv')



# Check for colinearity -------------------------------------------------

#Load all and combine
#Basics

CQ_summary <- read_csv('./DataFiles/hydro_data/event_char/CQ_summary.csv') %>%
  select(FI, HI = HI_mean)

simple_temporal <- read_csv('./DataFiles/hydro_data/event_char/simple_temporal.csv')
simple_size <- read_csv('./DataFiles/hydro_data/event_char/simple_size.csv')

#Gridded precip
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

#Outflows
outflows_inevent_totals <- read_csv('./DataFiles/hydro_data/event_char/outflows_inevent_totals.csv') %>%
  select(contains('1d'))

outflows_deltas <- read_csv('./DataFiles/hydro_data/event_char/outflows_deltas.csv') %>%
  select(contains('1d'))

flow_ratios <- read_csv('./DataFiles/hydro_data/event_char/flow_ratios.csv')

#SPEIs
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

#Combine
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

#p cors
p_cors <- event_summary %>%
  select(starts_with('p_')) %>%
  cor(.) %>%
  melt(.) %>%
  filter(X1 != X2)
view(filter(p_cors, X1 == 'p_10270104_180d_ante'))

corrplot(cor(select(event_summary, starts_with('p_') & ends_with('_total'))),
         type = 'lower', method = 'number')
#p totals - 10270104 & 10260008/10/15

corrplot(cor(select(event_summary, starts_with('p_') & ends_with('_180d_ante'))),
         type = 'lower', method = 'number')
#p antes - 15d 10270104 & 10260008/10/15
# 30d 10270104 & 10260015
# 60d 10270104 only
# 90d 10270104 only
#180d 10270104 only
#365d 10270104 & 10260008/10/15

#p maxint - 10270104 & 10260008/10/15 + 10270101
#p avgint - same as above

#also remove p_unint_vol

#r cors
r_cors <- event_summary %>%
  select(starts_with('r_')) %>%
  cor(.) %>%
  melt(.) %>%
  filter(X1 != X2)

corrplot(cor(select(event_summary, starts_with('r_'))),
         type = 'lower', method = 'number')

#remove r_sum_1d_total

cor(event_summary$reservoir_precip_ratio, event_summary$runoff_precip_ratio)

#remove runoff precip ratio

#spei cors
spei_comb <- cbind(spei1, spei3, spei6, spei9, spei12, spei18, spei24)
spei_cors <- melt(cor(select(spei_comb, contains('10270104'), contains('10260015')))) %>%
  filter(X1 != X2)

#remove spei 9 and 18s, western 24 and 6

correlated_predictors <- as_data_frame(c(colnames(select(event_summary, contains(c('10260008', '10260010', '10270101', '10270102', '10270103')))),
                     colnames(select(event_summary, c(p_10260015_30d_ante, p_10260015_60d_ante, p_10260015_90d_ante, p_10260015_180d_ante))),
                     colnames(select(event_summary, p_unint_vol)),
                     colnames(select(event_summary, c(r_sum_1d_total, runoff_precip_ratio))),
                     'd_10260015_spei9', 'd_10260015_spei18', 'd_10260015_spei6', 'd_10260015_spei24', 'd_10270104_spei9', 'd_10270104_spei18'))

write_csv(correlated_predictors, './DataFiles/hydro_data/event_char/correlated_predictors.csv')













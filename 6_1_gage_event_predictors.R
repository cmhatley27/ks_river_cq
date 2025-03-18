library(tidyverse)

# Load all relevant data --------------------------------------------------

all_hydro_15min <- with_tz(read_csv('./DataFiles/hydro_data/all_hydro_15min.csv'), 'America/Chicago')
events <- with_tz(read_csv('DataFiles/hydro_data/event_delineations/adjusted_events.csv'), 'America/Chicago') %>%
  select(event_number, start_dateTime, end_dateTime)
event_data <- with_tz(read_csv('./DataFiles/hydro_data/event_data/event_data.csv'), 'America/Chicago')


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
  select(!event_number)

# write_csv(simple_temporal, './DataFiles/hydro_data/event_char/simple_temporal.csv')


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

# write_csv(simple_size, './DataFiles/hydro_data/event_char/simple_size.csv')

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
# write_csv(N_yield, './DataFiles/hydro_data/event_char/N_yield.csv')

# Total outflows from each reservoir within events ------------------------

outflows <- all_hydro_15min %>%
  select(dateTime, event, 
         c = clintonQ_linterp, m = milfordQ_linterp, p = perryQ_linterp, t = tuttleQ_linterp, 
         k = kanopolisQ_linterp, wi = wilsonQ_linterp, wa = wacondaQ_linterp)

#Function for adding up 15min flows between event start minus a buffer to event end
event_flow <- function(dateTime_start, dateTime_end, buffer_length, dateTimes, flow_series){
  df <- data.frame(dateTime = dateTimes, cum_flows = cumsum(flow_series))
  flow_total <- df$cum_flows[df$dateTime == dateTime_end - days(buffer_length) - minutes(15)] - df$cum_flows[df$dateTime == dateTime_start - days(buffer_length) - minutes(15)]
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

outflows_inevent_totals <- (data.frame(within_event_flows(flows_frame = outflows, buffer_lengths = c(1))) * 60 * 15) %>%
  mutate(r_sum_1d_total = r_c_1d_total + r_p_1d_total + r_m_1d_total + r_t_1d_total + r_k_1d_total + r_wi_1d_total + r_wa_1d_total)

write_csv(outflows_inevent_totals, './DataFiles/hydro_data/event_char/outflows_inevent_totals.csv')


# Ratios and shares of total outflows -------------------------------------

flow_ratios <- tibble(
  reservoir_runoff_ratio = outflows_inevent_totals$r_sum_1d_total/simple_size$total_Q,
  r_c_share = outflows_inevent_totals$r_c_1d_total/simple_size$total_Q,
  r_p_share = outflows_inevent_totals$r_p_1d_total/simple_size$total_Q,
  r_m_share = outflows_inevent_totals$r_m_1d_total/simple_size$total_Q,
  r_t_share = outflows_inevent_totals$r_t_1d_total/simple_size$total_Q,
  r_k_share = outflows_inevent_totals$r_k_1d_total/simple_size$total_Q,
  r_wi_share = outflows_inevent_totals$r_wi_1d_total/simple_size$total_Q,
  r_wa_share = outflows_inevent_totals$r_wa_1d_total/simple_size$total_Q)

write_csv(flow_ratios, './DataFiles/hydro_data/event_char/flow_ratios.csv')
library(tidyverse)
library(sf)
library(terra)


# Load in data ------------------------------------------------------------

all_hydro_15min <- with_tz(read_csv('./DataFiles/hydro_data/all_hydro_15min.csv'), 'America/Chicago')
events <- with_tz(read_csv('DataFiles/hydro_data/event_delineations/BFLOW_events_bfi_threshold_adjusted.csv'), 'America/Chicago') %>%
  select(event_number, start_dateTime, end_dateTime)
event_data <- with_tz(read_csv('./DataFiles/hydro_data/event_data/event_data.csv'), 'America/Chicago')

precip_huc8 <- read_csv('./DataFiles/hydro_data/precip/precip_huc8.csv')
precip_global_mean <- read_csv('./DataFiles/hydro_data/precip/precip_global_mean.csv')

spei12_huc8 <- read_csv('./DataFiles/hydro_data/drought/spei12_huc8.csv')

# Within-event total precip ---------------------------------------------

#Function for adding up precip in a window from event start minus a buffer to event max
event_precip <- function(date_start, date_max, buffer_length, dates, precip_series){
  df <- data.frame(date = dates, precip = precip_series)
  df <- df %>%
    mutate(cumul_precip = cumsum(precip))
  precip_total <- df$cumul_precip[df$date == date_max] - max(df$cumul_precip[df$date == date_start - buffer_length - 1], 0)
  return(ifelse(date_start %within% interval(dates[1], dates[1] + buffer_length), NA, precip_total))
}

#Loops within-event precip function over every gage+buffer length combo
within_event_precip <- function(precip_frame, buffer_lengths, event_frame = event_data){
  #extract dates
  dates <- precip_frame$date
  #extract just precip series
  precips <- precip_frame %>%
    select(!date)
  #initialize list of output vectors
  window_totals <- list()
  #Loop through each combo of gage+buffer length and calculate cumulative precip
  for(gage in 1:ncol(precips)){
    for(buffer in 1:length(buffer_lengths)){
      summary_frame <- event_frame %>%
        group_by(event) %>%
        summarize(date_start = date(dateTime[1]),
                  date_max = date(dateTime[which.max(riverQ_linterp)]),
                  precip_in_window = event_precip(date_start = date_start, date_max = date_max, buffer_length = buffer_lengths[buffer], dates = precip_frame$date, precip_series = precips[[gage]])
        )
      window_totals[[paste('p_', colnames(precips[gage]), '_', buffer_lengths[buffer], 'd_total', sep = "")]] <- summary_frame$precip_in_window
    }
    print(paste0('huc8 ', gage, '/', ncol(precips), ' completed'))
  }
  return(window_totals)
}

#by huc8
precip_inevent_huc8 <- data.frame(within_event_precip(precip_frame = precip_huc8, buffer_lengths = c(1)))
write_csv(precip_inevent_huc8, './DataFiles/hydro_data/event_char/precip_inevent_huc8.csv')


# global
precip_inevent_global <- data.frame(within_event_precip(precip_frame = precip_global_mean, buffer_lengths = 1)) %>%
  rename(p_global_1d_total = p_precip_1d_total)
write_csv(precip_inevent_global, './DataFiles/hydro_data/event_char/precip_inevent_global.csv')


# Ante-event total precip by huc8 -----------------------------------------

#Function for adding up precip in a window from event start
cumul_precip <- function(date_start, window_length, dates, precip_series){
  df <- data.frame(date = dates, precip = precip_series)
  df <- df %>%
    mutate(cumul_precip = cumsum(precip))
  precip_total <- df$cumul_precip[df$date == date_start] - max(df$cumul_precip[df$date == date_start - window_length - 1], 0)
  return(ifelse(date_start %within% interval(dates[1], dates[1] + window_length), NA, precip_total))
}

#**Loops ante-event precip function over every gage+window length combo
before_event_precip <- function(precip_frame, window_lengths, event_frame = event_data){
  #extract dates
  dates <- precip_frame$date
  #extract just precip series
  precips <- precip_frame %>%
    select(!date)
  #initialize list of output vectors
  window_totals <- list()
  #Loop through each combo of gage+window length and calculate cumulative precip
  for(gage in 1:ncol(precips)){
    for(window in 1:length(window_lengths)){
      summary_frame <- event_frame %>%
        group_by(event) %>%
        summarize(date_start = date(dateTime[1]),
                  date_max = date(dateTime[which.max(riverQ_linterp)]),
                  precip_in_window = cumul_precip(date_start = date_start, window_length = window_lengths[window], dates = precip_frame$date, precip_series = precips[[gage]])
        )
      window_totals[[paste('p_', colnames(precips[gage]), '_', window_lengths[window], 'd_ante', sep = "")]] <- summary_frame$precip_in_window
    }
    print(paste0('huc8 ', gage, '/', ncol(precips), ' completed'))
  }
  return(window_totals)
}

#by HUC 8
precip_anteevent_huc8 <- data.frame(before_event_precip(precip_frame = precip_huc8, window_lengths = 30))
write_csv(precip_anteevent_huc8, './DataFiles/hydro_data/event_char/precip_anteevent_huc8.csv')

#global
precip_anteevent_global <- data.frame(before_event_precip(precip_frame = precip_global_mean, window_lengths = 30)) %>%
  rename(p_global_30d_ante = p_precip_30d_ante)
write_csv(precip_anteevent_global, './DataFiles/hydro_data/event_char/precip_anteevent_global.csv')


# 1-Year SPEI for each event by huc8 ------------------------------------------
spei12 <- as_tibble(matrix(nrow = nrow(events), ncol = ncol(spei12_huc8) - 1))

for(event_num in seq(events$event_number)) {
  event_filter <- subset(event_data, event == event_num)
  start_date <- floor_date(event_filter$dateTime[1], 'month')
  if(day(event_filter$dateTime[1]) < 15){start_date <- start_date - months(1)}
  
  spei12_select <- spei12_huc8 %>%
    filter(date == start_date)
  
  spei12[event_num,] <- spei12_select[1,-1]
}
colnames(spei12) <- paste0('d_', colnames(spei12_huc8)[-1], '_spei12')

write_csv(spei12, './DataFiles/hydro_data/event_char/spei12.csv')


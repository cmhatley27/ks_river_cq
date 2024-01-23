library(tidyverse)
library(sf)
library(terra)
library(reshape)


# Load in data ------------------------------------------------------------

all_hydro_15min <- with_tz(read_csv('./DataFiles/hydro_data/all_hydro_15min.csv'), 'America/Chicago')
events <- with_tz(read_csv('DataFiles/hydro_data/event_delineations/BFLOW_events_bfi_threshold_adjusted.csv'), 'America/Chicago') %>%
  select(event_number, start_dateTime, end_dateTime)
event_data <- with_tz(read_csv('./DataFiles/hydro_data/event_data/event_data.csv'), 'America/Chicago')
event_summary_trim <- with_tz(read_csv('./DataFiles/hydro_data/event_summary_trim.csv'), 'America/Chicago') %>%
  mutate(wateryear = factor(wateryear),
         month = factor(month),
         season = factor(season))
event_summary <- with_tz(read_csv('./DataFiles/hydro_data/event_summary.csv'), 'America/Chicago')

precip_huc8 <- read_csv('./DataFiles/hydro_data/precip/precip_huc8.csv')

spei1_huc8 <- read_csv('./DataFiles/hydro_data/drought/spei1_huc8.csv')
spei3_huc8 <- read_csv('./DataFiles/hydro_data/drought/spei3_huc8.csv')
spei6_huc8 <- read_csv('./DataFiles/hydro_data/drought/spei6_huc8.csv')
spei9_huc8 <- read_csv('./DataFiles/hydro_data/drought/spei9_huc8.csv')
spei12_huc8 <- read_csv('./DataFiles/hydro_data/drought/spei12_huc8.csv')
spei18_huc8 <- read_csv('./DataFiles/hydro_data/drought/spei18_huc8.csv')
spei24_huc8 <- read_csv('./DataFiles/hydro_data/drought/spei24_huc8.csv')

# Within-event total precip by HUC8 ---------------------------------------------

#Function for adding up precip in a window from event start minus a buffer to event max
event_precip <- function(date_start, date_max, buffer_length, dates, precip_series){
  df <- data.frame(date = dates, precip = precip_series)
  df <- df %>%
    mutate(cumul_precip = cumsum(precip))
  precip_total <- df$cumul_precip[df$date == date_max] - max(df$cumul_precip[df$date == date_start - buffer_length - 1], 0)
  return(ifelse(date_start %within% interval(dates[1], dates[1] + buffer_length), NA, precip_total))
}

#**Loops within-event precip function over every gage+buffer length combo
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

precip_inevent_gridded <- data.frame(within_event_precip(precip_frame = precip_huc8, buffer_lengths = c(0,1,2)))
write_csv(precip_inevent_gridded, './DataFiles/hydro_data/event_char/precip_inevent_gridded.csv')

##Compare results with station data
precip_comp <- precip_inevent_gridded %>%
  select(contains('1d')) %>%
  cbind(., select(event_summary, contains('precip') & contains('1day_buffer')))

ggplot(data = precip_comp) +
  geom_point(aes(x = total_event_lawrence_precip_1day_buffer, y = p_10270104_1d_total)) +
  geom_abline(slope = 1)

precip_corr <- cor(precip_comp)
library(corrplot)
corrplot(precip_corr, method = 'pie', tl.cex = 0.7, tl.col = 'black', type = 'lower', diag = FALSE, col.lim = c(0, 1))


##Find coordinates of HUC8 with max precip over course of event for 1d buffer
huc8_centroids <- read_csv('C:/School/SAFE KAW/Data/DataFiles/shapefiles/created/huc8_centroids.csv')

max_precip_coords <- select(precip_inevent_gridded, contains('1d_total')) %>%
  mutate(event = seq(1,nrow(precip_inevent_gridded))) %>%
  pivot_longer(!event, names_to = 'huc8', values_to = 'total_precip') %>%
  group_by(event) %>%
  summarise(max_huc8_1d = huc8[which.max(total_precip)]) %>%
  select(max_huc8_1d) %>%
  mutate(max_huc8_1d = as.double(str_sub(max_huc8_1d, 3, 10))) %>%
  left_join(huc8_centroids, by = join_by(max_huc8_1d == huc8)) %>%
  select(p_max_x = x, p_max_y = y)

write_csv(max_precip_coords, './DataFiles/hydro_data/event_char/max_precip_coords.csv')


# Avg and max precip intensity per huc8 during each event ----------------

#avg intensity
event_durations <- vector(length = nrow(events))
for(event_num in seq(events$event_number)) {
  event_filter <- subset(event_data, event == event_num)
  start_date <- (event_filter$dateTime[1]) - days(1)
  max_date <- (event_filter$dateTime[which.max(event_filter$riverQ_linterp)])
  event_durations[event_num] <- time_length(as.duration(interval(start = start_date, end = max_date)))/(60*60*24)
}

avg_precip_intensity <- precip_inevent_gridded/event_durations

colnames(avg_precip_intensity) <- paste0(str_sub(colnames(precip_inevent_gridded), 1, 14), 'avgint')

write_csv(avg_precip_intensity, './DataFiles/hydro_data/event_char/avg_precip_intensity.csv')


#max intensity

max_precip_intensity <- as_tibble(matrix(nrow = nrow(events), ncol = ncol(precip_huc8) - 1))

for(event_num in seq(events$event_number)) {
  event_filter <- subset(event_data, event == event_num)
  start_date <- date(event_filter$dateTime[1]) - 1
  max_date <- date(event_filter$dateTime[which.max(event_filter$riverQ_linterp)])
  precip_filter <- precip_huc8 %>%
    filter(date %within% interval(start = start_date, end = max_date))
  max_intensity <- precip_filter %>%
    pivot_longer(!date, names_to = 'huc8', values_to = 'precip') %>%
    group_by(huc8) %>%
    summarize(max_intensity = max(precip)) %>%
    pivot_wider(names_from = huc8, values_from = max_intensity)
  
  max_precip_intensity[event_num,] <- max_intensity[1,]
}
colnames(max_precip_intensity) <- paste0('p_', colnames(max_intensity), '_maxint')

write_csv(max_precip_intensity, './DataFiles/hydro_data/event_char/max_precip_intensity.csv')


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

precip_anteevent_gridded <- data.frame(before_event_precip(precip_frame = precip_huc8, window_lengths = c(15,30,60,90,180,365)))
write_csv(precip_anteevent_gridded, './DataFiles/hydro_data/event_char/precip_anteevent_gridded.csv')


# Within-event precip VOLUME across basin ---------------------------------

#Multiplies 1d buffer in-event precip depths in each UNINTERCEPTED huc8 by that
#huc8's area, and adds all together. Gives volume of precip that directly fell
#into river during event

huc8_boundaries <- st_read('C:/School/SAFE KAW/Data/DataFiles/shapefiles/created/huc8_boundaries.shp') 
huc8_areas <- tibble(HUC8 = huc8_boundaries$HUC8,
                     area_sqkm = huc8_boundaries$AREASQKM)

precip_inevent_unint <- read_csv('./DataFiles/hydro_data/event_char/precip_inevent_gridded.csv') %>%
  select(contains('1d')) %>%
  select(contains(c('10260008', '10260010', '10260015', '10270101', '10270102', '10270104')))

colnames(precip_inevent_unint) <- str_sub(colnames(precip_inevent_unint), 3, 10)

precip_unintercepted_volumes <- precip_inevent_unint %>%
  mutate(event = seq(nrow(precip_inevent_unint))) %>%
  pivot_longer(!event, names_to = 'HUC8', values_to = 'precip_depth') %>%
  left_join(huc8_areas) %>%
  mutate(precip_volume = (precip_depth/1000)*(area_sqkm*1000000)) %>%
  group_by(event) %>%
  summarize(p_unint_vol = sum(precip_volume)) %>%
  select(p_unint_vol)
  
write_csv(precip_unintercepted_volumes, './DataFiles/hydro_data/event_char/precip_unintercepted_volumes.csv')


# Precip distributions per huc8 -------------------------------------------
#### Shows that daily spatial means basically equal daily spatial medians,
#### i.e. precip across each huc8 is roughly normally distributed

huc8 <- 43
huc8_ids <- huc8_boundaries$HUC8

local_boundary <- huc8_boundaries %>%
  filter(HUC8 == huc8_ids[huc8])

local_precip <- mask(crop(precip_ksrb, local_boundary), vect(local_boundary)) %>%
  r2df(., 'precip')

rain_dates <- precip_huc8 %>%
  select(date, precip = huc8_ids[huc8]) %>%
  filter(precip >= 1)

dist_stats <- tibble(date = rain_dates$date, mean = NA, median = NA)

for(rain_date in 1:nrow(rain_dates)){
  
  date_precip <- local_precip %>%
    filter(date == date[rain_date])
  
  # ggplot(data = date_precip) +
  #   geom_histogram(aes(x = precip)) +
  #   geom_vline(xintercept = mean(date_precip$precip), color = 'blue') +
  #   geom_vline(xintercept = median(date_precip$precip), color = 'red')
  
  dist_stats$mean[rain_date] <- mean(date_precip$precip)
  dist_stats$median[rain_date] <- median(date_precip$precip)
  
  print(paste0('date ', rain_date, '/', nrow(rain_dates), ' completed'))
}

ggplot(data = dist_stats) +
  geom_point(aes(x = mean, y = median))




# Drought indices for each event ------------------------------------------

#SPEI1

spei1 <- as_tibble(matrix(nrow = nrow(events), ncol = ncol(spei1_huc8) - 1))

for(event_num in seq(events$event_number)) {
  event_filter <- subset(event_data, event == event_num)
  start_date <- floor_date(event_filter$dateTime[1], 'month')
  if(day(event_filter$dateTime[1]) < 15){start_date <- start_date - months(1)}

  spei1_select <- spei1_huc8 %>%
    filter(date == start_date)

  spei1[event_num,] <- spei1_select[1,-1]
}
colnames(spei1) <- paste0('d_', colnames(spei1_huc8)[-1], '_spei1')

write_csv(spei1, './DataFiles/hydro_data/event_char/spei1.csv')

#SPEI3

spei3 <- as_tibble(matrix(nrow = nrow(events), ncol = ncol(spei3_huc8) - 1))

for(event_num in seq(events$event_number)) {
  event_filter <- subset(event_data, event == event_num)
  start_date <- floor_date(event_filter$dateTime[1], 'month')
  if(day(event_filter$dateTime[1]) < 15){start_date <- start_date - months(1)}
  
  spei3_select <- spei3_huc8 %>%
    filter(date == start_date)
  
  spei3[event_num,] <- spei3_select[1,-1]
}
colnames(spei3) <- paste0('d_', colnames(spei3_huc8)[-1], '_spei3')

write_csv(spei3, './DataFiles/hydro_data/event_char/spei3.csv')

#SPEI6

spei6 <- as_tibble(matrix(nrow = nrow(events), ncol = ncol(spei6_huc8) - 1))

for(event_num in seq(events$event_number)) {
  event_filter <- subset(event_data, event == event_num)
  start_date <- floor_date(event_filter$dateTime[1], 'month')
  if(day(event_filter$dateTime[1]) < 15){start_date <- start_date - months(1)}
  
  spei6_select <- spei6_huc8 %>%
    filter(date == start_date)
  
  spei6[event_num,] <- spei6_select[1,-1]
}
colnames(spei6) <- paste0('d_', colnames(spei6_huc8)[-1], '_spei6')

write_csv(spei6, './DataFiles/hydro_data/event_char/spei6.csv')

#SPEI9

spei9 <- as_tibble(matrix(nrow = nrow(events), ncol = ncol(spei9_huc8) - 1))

for(event_num in seq(events$event_number)) {
  event_filter <- subset(event_data, event == event_num)
  start_date <- floor_date(event_filter$dateTime[1], 'month')
  if(day(event_filter$dateTime[1]) < 15){start_date <- start_date - months(1)}
  
  spei9_select <- spei9_huc8 %>%
    filter(date == start_date)
  
  spei9[event_num,] <- spei9_select[1,-1]
}
colnames(spei9) <- paste0('d_', colnames(spei9_huc8)[-1], '_spei9')

write_csv(spei9, './DataFiles/hydro_data/event_char/spei9.csv')

#SPEI12

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

#SPEI18

spei18 <- as_tibble(matrix(nrow = nrow(events), ncol = ncol(spei18_huc8) - 1))

for(event_num in seq(events$event_number)) {
  event_filter <- subset(event_data, event == event_num)
  start_date <- floor_date(event_filter$dateTime[1], 'month')
  if(day(event_filter$dateTime[1]) < 15){start_date <- start_date - months(1)}
  
  spei18_select <- spei18_huc8 %>%
    filter(date == start_date)
  
  spei18[event_num,] <- spei18_select[1,-1]
}
colnames(spei18) <- paste0('d_', colnames(spei18_huc8)[-1], '_spei18')

write_csv(spei18, './DataFiles/hydro_data/event_char/spei18.csv')

#SPEI24

spei24 <- as_tibble(matrix(nrow = nrow(events), ncol = ncol(spei24_huc8) - 1))

for(event_num in seq(events$event_number)) {
  event_filter <- subset(event_data, event == event_num)
  start_date <- floor_date(event_filter$dateTime[1], 'month')
  if(day(event_filter$dateTime[1]) < 15){start_date <- start_date - months(1)}
  
  spei24_select <- spei24_huc8 %>%
    filter(date == start_date)
  
  spei24[event_num,] <- spei24_select[1,-1]
}
colnames(spei24) <- paste0('d_', colnames(spei24_huc8)[-1], '_spei24')

write_csv(spei24, './DataFiles/hydro_data/event_char/spei24.csv')

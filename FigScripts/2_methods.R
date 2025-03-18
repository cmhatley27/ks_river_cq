library(tidyverse)
source('Theme+Settings.R')
source('8_combine_all_event_data.R')

all_hydro_15min <- with_tz(read_csv('./DataFiles/hydro_data/all_hydro_15min.csv'), 'America/Chicago')
all_hydro_daily <- read_csv('./DataFiles/hydro_data/all_hydro_daily.csv')
events <- with_tz(read_csv('DataFiles/hydro_data/event_delineations/BFLOW_events_bfi_threshold_adjusted.csv'), 'America/Chicago') %>%
  select(event_number, start_dateTime, end_dateTime)


# Plot long time series ---------------------------------------------------
selected_events <- c(20,139)
#2 events = c(10,139)
#4 events = c(20,139,167,255)

ggplot(data = all_hydro_15min) +
  #highligh all events
  geom_rect(data = subset(events, event_number %in% event_summary$event), aes(xmin = start_dateTime, xmax = end_dateTime),
            ymin = -200, ymax = 1.05*max(all_hydro_15min$riverQ_linterp), fill = 'grey90') +
  #highlight selected events
  geom_rect(data = subset(events, event_number %in% selected_events[1:2]), aes(xmin = start_dateTime, xmax = end_dateTime),
            ymin = -200, ymax = 1.05*max(all_hydro_15min$riverQ_linterp), fill = 'gold', color = 'gold') +
  #label numbers for each event
  # geom_text(data = events, aes(x = start_dateTime, label = row_number(start_dateTime)), y = 600) +
  #N series
  geom_line(aes(x = dateTime, y = N*100), color = 'red', linewidth = 0.25) +
  #Q series
  geom_line(aes(x = dateTime, y = riverQ_linterp), linewidth = 0.25) +
  #scales
  scale_y_continuous(limits = c(0, 2000), name = expression(Q~'['*m^3*'/s]'),
                     sec.axis = sec_axis(transform = ~. / 100, name = expression(NO[3]^'-'-N~'[mg/L]'))) +
  scale_x_datetime(date_breaks = '1 year', date_labels = '%Y', expand = c(0,0), name = NULL)
  

ggsave(file.path('Figures','final','methods','long_time_series.png'), height = 1.5, width = 6.5, units = 'in', dpi = 500)


# Event time series -------------------------------------------------------
event_data <- read_csv(file.path('DataFiles','hydro_data','event_data','event_data.csv'))
event_summary$reservoir_runoff_ratio[event_summary$event %in% selected_events]

selected_events <- c(20,139)
#2 events = c(10,139)
#4 events = c(20,139,167,255)

event_sel <- 2
event_dat <- filter(event_data, event == selected_events[event_sel])
n_scale <- c(50,100,700,100)
n_shift <- c(0,0,0,1.5)
date_breaks <- c('1 day', '2 weeks', '3 days', '1 day')

ggplot(event_dat, aes(x = dateTime)) +
  {if(event_sel %in% c(2)) geom_line(aes(y = reservoir_sum_linterp), linetype = 'dashed', color = 'grey50')} +
  geom_line(aes(y = (N_linterp_smooth + n_shift[event_sel])*n_scale[event_sel]), color = 'red') +
  geom_line(aes(y = riverQ_linterp)) +
  scale_y_continuous(name = expression(Q~'['*m^3*'/s]'),
                     sec.axis = sec_axis(transform = ~./n_scale[event_sel] - n_shift[event_sel], name = expression(NO[3]^'-'-N~'[mg/L]'))) +
  scale_x_datetime(name = NULL, date_labels = '%d %b', date_breaks = date_breaks[event_sel]) +
  #add margin to left to line up with long time series if needed
  if(event_sel %in% c(1,2,4)){
    theme(plot.margin = unit(c(1,0,0,2.9),'mm'))
  } else{
    theme(plot.margin = unit(c(1,0,0,1),'mm'))
  } 
  
ggsave(file.path('Figures','final','methods',paste0('event_ts_',event_sel,'.png')), height = 1.5, width = 3, units = 'in', dpi = 500)

  
# Event CQ loops ----------------------------------------------------------
rise <- read_csv(file.path('DataFiles','hydro_data','event_data','rise_limbs.csv')) %>%
  select(!dateTime) %>%
  dplyr::rename(N_norm = rise_N_norm) %>%
  mutate(timestep = Q_norm)
fall <- read_csv(file.path('DataFiles','hydro_data','event_data','fall_limbs.csv')) %>%
  select(!dateTime) %>%
  dplyr::rename(N_norm = fall_N_norm) %>%
  mutate(timestep = 2-Q_norm)
rise_fall <- rbind(rise, fall)

selected_events <- c(20,139,167,255)

event_sel <- 2
event_dat <- filter(rise_fall, event == selected_events[event_sel]) %>%
  dplyr::arrange(timestep)

ggplot(data = event_dat, aes(x = Q_norm, y = N_norm, color = timestep)) +
  geom_path() +
  scale_color_gradient2(low = 'black',mid = 'blue', high = 'red', midpoint = 1,guide = 'none') +
  scale_x_continuous(name = NULL, breaks = c(0,1)) +
  scale_y_continuous(name = expression(NO[3]^'-'-N~'(normalized)'), breaks = c(0,1),
                     expand = c(0.05,0.05), position = 'left') +
  theme(plot.margin = unit(c(1,1,0,1),'mm'))


ggsave(file.path('Figures','final','methods',paste0('event_loop_',event_sel,'.png')), height = 1.5, width = 2.5, units = 'in', dpi = 500)

event_summary$FI_n[event_summary$event %in% selected_events]
event_summary$HI_n[event_summary$event %in% selected_events]
event_summary$reservoir_runoff_ratio[event_summary$event %in% selected_events]

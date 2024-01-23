library(tidyverse)
library(lubridate)
library(dataRetrieval)
library(hydroEvents)
library(rnoaa)
library(ggpubr)
library(FluMoDL)

# Create folders and load data ---------------------------------------------------------------

if(!dir.exists('DataFiles/hydro_data/event_delineations')) {
  dir.create('DataFiles/hydro_data/event_delineations')
}

all_hydro_15min <- with_tz(read_csv('./DataFiles/hydro_data/all_hydro_15min.csv'), tz = 'America/Chicago')
all_hydro_daily <- read_csv('./DataFiles/hydro_data/all_hydro_daily.csv')
DeSoto_hourly <- with_tz(read_csv('./DataFiles/hydro_data/DeSoto_gage/DeSoto_hourly.csv'), 'America/Chicago')



# Prepare Q time series for delineation -----------------------------------

## Selected linterp time series because delineation can't take any NAs
## Added smoothed time series to 15min data


eventdelin_15min <- all_hydro_15min %>%
  select(dateTime, riverQ_linterp) %>%
  mutate(riverQ_linterp_smooth = (lag(riverQ_linterp)+2*riverQ_linterp+lead(riverQ_linterp))/4)

eventdelin_hourly <- DeSoto_hourly %>%
  select(dateHour, riverQ_linterp = Q_linterp)

eventdelin_daily <- all_hydro_daily %>%
  select(date, riverQ_linterp)


# Calculate baseflow and quickflow using Lyne and Hollick filter ----------

alpha <- 0.925
passes <- 3
reflected <- 30

eventdelin_15min <- eventdelin_15min %>%
  mutate(bf = baseflowB(eventdelin_15min$riverQ_linterp, alpha = alpha, passes = passes, r = reflected)[["bf"]], # 
         qf = riverQ_linterp - bf,
         bfi = bf/riverQ_linterp)

eventdelin_hourly <- eventdelin_hourly %>%
  mutate(bf = baseflowB(eventdelin_hourly$riverQ_linterp, alpha = alpha, passes = passes, r = reflected)[["bf"]],
         qf = riverQ_linterp - bf,
         bfi = bf/riverQ_linterp)

eventdelin_daily <- eventdelin_daily %>%
  mutate(bf = baseflowB(eventdelin_daily$riverQ_linterp, alpha = alpha, passes = passes, r = reflected)[["bf"]],
         qf = riverQ_linterp - bf,
         bfi = bf/riverQ_linterp)


# Event Delineation - BFI Method ------------------------------------------

#*** 15 minute frequency ---------------------------------------------------------

bfi_threshold <- 0.85
mindiff <- 1

eventdelin_15min.bfi <- eventBaseflow(eventdelin_15min$riverQ_linterp, 
                                      BFI_Th = bfi_threshold, min.diff = mindiff) 
eventdelin_15min.bfi <- eventdelin_15min.bfi %>%
  mutate(event_start = as_datetime(eventdelin_15min$dateTime[eventdelin_15min.bfi$srt]),
         event_end = as_datetime(eventdelin_15min$dateTime[eventdelin_15min.bfi$end]),
         event_max_date = as_datetime(eventdelin_15min$dateTime[eventdelin_15min.bfi$which.max])) %>%
  select(event_start, event_end, event_max_date, event_max_Q = max)
write_csv(eventdelin_15min.bfi, './DataFiles/hydro_data/event_delineations/bfi_15min.csv')

startyear <- 2019
endyear <- 2019

ggplot(data = all_hydro_15min) +
  #Event limits and labels
  geom_rect(data = eventdelin_15min.bfi, aes(xmin = event_start, xmax = event_end), 
            ymin = -200, ymax = 1.05*max(eventdelin_15min$riverQ_linterp), fill = 'grey90', color = 'grey50', linetype = 'dashed') +
  geom_text(data = eventdelin_15min.bfi, aes(x = event_start, label = row_number(event_start)), y = 0.9*max(eventdelin_15min$riverQ_linterp)) +
  #River Q and baseflow
  geom_line(aes(x = dateTime, y = riverQ_linterp)) +
  geom_line(data = eventdelin_15min, aes(x = dateTime, y = bf), linetype = 'dashed') +
  theme_classic() +
  labs(title = ifelse(startyear == endyear, 
                      str_c('BFI Method, ', 'Water Year ', startyear), 
                      str_c('BFI Method, ', 'Water Years ', startyear, '-', endyear)),
       x = NULL) +
  scale_x_datetime(limits = as_datetime(c(str_c(startyear-1, '-10-01'), str_c(endyear, '-09-30'))))


#*** Hourly frequency --------------------------------------------------------

bfi_threshold <- 0.85
mindiff <- 1

eventdelin_hourly.bfi <- eventBaseflow(eventdelin_hourly$riverQ_linterp, 
                                       BFI_Th = bfi_threshold, min.diff = mindiff)
eventdelin_hourly.bfi <- eventdelin_hourly.bfi %>%
  mutate(event_start = as_datetime(eventdelin_hourly$dateHour[eventdelin_hourly.bfi$srt]),
         event_end = as_datetime(eventdelin_hourly$dateHour[eventdelin_hourly.bfi$end]),
         event_max_date = as_datetime(eventdelin_hourly$dateHour[eventdelin_hourly.bfi$which.max])) %>%
  select(event_start, event_end, event_max_date, event_max_Q = max)
write_csv(eventdelin_hourly.bfi, './DataFiles/hydro_data/event_delineations/bfi_hourly.csv')

startyear <- 2019
endyear <- 2019

ggplot() +
  geom_rect(data = eventdelin_hourly.bfi, aes(xmin = event_start, xmax = event_end), 
            ymin = -200, ymax = 1.05*max(eventdelin_hourly$riverQ_linterp), fill = 'grey90', color = 'grey50', linetype = 'dashed') +
  geom_text(data = eventdelin_hourly.bfi, aes(x = event_start, label = row_number(event_start)), y = 0.9*max(eventdelin_hourly$riverQ_linterp)) +
  geom_line(data = eventdelin_hourly, aes(x = dateHour, y = riverQ_linterp)) +
  theme_classic() +
  labs(title = ifelse(startyear == endyear, 
                      str_c('BFI Method, ', 'Water Year ', startyear), 
                      str_c('BFI Method, ', 'Water Years ', startyear, '-', endyear)),
       x = NULL) +
  scale_x_datetime(limits = as_datetime(c(str_c(startyear-1, '-10-01'), str_c(endyear, '-09-30'))))


#*** Daily frequency ---------------------------------------------------------

bfi_threshold <- 0.85
mindiff <- 1

eventdelin_daily.bfi <- eventBaseflow(eventdelin_daily$riverQ_linterp, 
                                      BFI_Th = bfi_threshold, min.diff = mindiff)
eventdelin_daily.bfi <- eventdelin_daily.bfi %>%
  mutate(event_start = date(eventdelin_daily$date[eventdelin_daily.bfi$srt]),
         event_end = date(eventdelin_daily$date[eventdelin_daily.bfi$end]),
         event_max_date = date(eventdelin_daily$date[eventdelin_daily.bfi$which.max])) %>%
  select(event_start, event_end, event_max_date, event_max_Q = max)
write_csv(eventdelin_daily.bfi, './DataFiles/hydro_data/event_delineations/bfi_daily.csv')

startyear <- 2019
endyear <- 2019

ggplot(data = all_hydro_daily) +
  #Event Highlights & Labels
  geom_rect(data = eventdelin_daily.bfi, aes(xmin = event_start, xmax = event_end), 
            ymin = -200, ymax = 1.05*max(eventdelin_daily$riverQ_linterp), fill = 'grey95', color = 'grey50', linetype = 'dashed') +
  geom_text(data = eventdelin_daily.bfi, aes(x = event_start + 1, label = row_number(event_start)), y = 0.9*max(eventdelin_daily$riverQ_linterp)) +
  #Precip Bars
  geom_col(aes(x = date, y = lawrence_precip*10), fill = 'lightblue', alpha = 1) +
  #River Q and baseflow
  geom_line(aes(x = date, y = riverQ)) +
  geom_line(data = eventdelin_daily, aes(x = date, y = bf), color = 'grey20', linetype = 'dashed') +
  ##Reservoir Outflows from USACE dam info
  #geom_line(data = clinton_dam_daily, aes(x = date, y = outflow), color = 'red') +
  #geom_line(data = milford_dam_daily, aes(x = date, y = outflow), color = 'gold') +
  #geom_line(data = perry_dam_daily, aes(x = date, y = outflow), color = 'blue') +
  #geom_line(data = tuttle_dam_daily, aes(x = date, y = outflow), color = 'purple') +
  ##Reservoir Outflows from USGS gages downstream of dams
  geom_line(aes(x = date, y = clintonQ), color = 'red') +
  geom_line(aes(x = date, y = milfordQ), color = 'gold') +
  geom_line(aes(x = date, y = perryQ), color = 'blue') +
  geom_line(aes(x = date, y = tuttleQ), color = 'purple') +
  #Theme and Axes
  theme_classic() +
  labs(title = ifelse(startyear == endyear, 
                      str_c('BFI Method, ', 'Water Year ', startyear), 
                      str_c('BFI Method, ', 'Water Years ', startyear, '-', endyear)),
       x = NULL) +
  scale_x_date(limits = as_date(c(str_c(startyear-1, '-10-01'), str_c(endyear, '-09-30'))),
               date_breaks = str_c(endyear-startyear+1,' months'),
               date_labels = '%b %y')


# Event Delineation - Minima Method ---------------------------------------

#*** 15 minute frequency -----------------------------------------------------

deltay <- 40
deltax <- 1
threshold <- 15

eventdelin_15min.min <- eventMinima(eventdelin_15min$qf, 
                                  delta.y = deltay, delta.x = deltax*96, threshold = threshold) 
eventdelin_15min.min <- eventdelin_15min.min %>%
  mutate(event_start = as_datetime(eventdelin_15min$dateTime[eventdelin_15min.min$srt]),
         event_end = as_datetime(eventdelin_15min$dateTime[eventdelin_15min.min$end]),
         event_max_date = as_datetime(eventdelin_15min$dateTime[eventdelin_15min.min$which.max])) %>%
  select(event_start, event_end, event_max_date, event_max_Q = max)
write_csv(eventdelin_15min.min, './DataFiles/hydro_data/event_delineations/min_15min.csv')

startyear <- 2019
endyear <- 2019

ggplot() +
  geom_rect(data = eventdelin_15min.min, aes(xmin = event_start, xmax = event_end), 
            ymin = -200, ymax = 1.05*max(eventdelin_15min$riverQ_linterp), fill = 'grey90', color = 'grey50', linetype = 'dashed') +
  geom_text(data = eventdelin_15min.min, aes(x = event_start, label = row_number(event_start)), y = 0.9*max(eventdelin_15min$riverQ_linterp)) +
  geom_line(data = eventdelin_15min, aes(x = dateTime, y = riverQ_linterp)) +
  theme_classic() +
  labs(title = ifelse(startyear == endyear, 
                      str_c('Local Minima, ', 'Water Year ', startyear), 
                      str_c('Local Minima, ', 'Water Years ', startyear, '-', endyear)),
       x = NULL) +
  scale_x_datetime(limits = as_datetime(c(str_c(startyear-1, '-10-01'), str_c(endyear, '-09-30'))))


#*** Hourly frequency --------------------------------------------------------

deltay <- 5
deltax <- 1
threshold <- 15

eventdelin_hourly.min <- eventMinima(eventdelin_hourly$qf, 
                                  delta.y = deltay, delta.x = deltax*24, threshold = threshold)
eventdelin_hourly.min <- eventdelin_hourly.min %>%
  mutate(event_start = as_datetime(eventdelin_hourly$dateHour[eventdelin_hourly.min$srt]),
         event_end = as_datetime(eventdelin_hourly$dateHour[eventdelin_hourly.min$end]),
         event_max_date = as_datetime(eventdelin_hourly$dateHour[eventdelin_hourly.min$which.max])) %>%
  select(event_start, event_end, event_max_date, event_max_Q = max)
write_csv(eventdelin_hourly.min, './DataFiles/hydro_data/event_delineations/min_hourly.csv')

startyear <- 2019
endyear <- 2019

ggplot() +
  geom_rect(data = eventdelin_hourly.min, aes(xmin = event_start, xmax = event_end), 
            ymin = -200, ymax = 1.05*max(eventdelin_hourly$riverQ_linterp), fill = 'grey90', color = 'grey50', linetype = 'dashed') +
  geom_text(data = eventdelin_hourly.min, aes(x = event_start, label = row_number(event_start)), y = 0.9*max(eventdelin_hourly$riverQ_linterp)) +
  geom_line(data = eventdelin_hourly, aes(x = dateHour, y = riverQ_linterp)) +
  theme_classic() +
  labs(title = ifelse(startyear == endyear, 
                      str_c('Local Minima, ', 'Water Year ', startyear), 
                      str_c('Local Minima, ', 'Water Years ', startyear, '-', endyear)),
       x = NULL) +
  scale_x_datetime(limits = as_datetime(c(str_c(startyear-1, '-10-01'), str_c(endyear, '-09-30'))))


#*** Daily frequency ---------------------------------------------------------

deltay <- 5
deltax <- 1
threshold <- 15

eventdelin_daily.min <- eventMinima(eventdelin_daily$qf, 
                                  delta.y = deltay, delta.x = deltax, threshold = threshold)
eventdelin_daily.min <- eventdelin_daily.min %>%
  mutate(event_start = date(eventdelin_daily$date[eventdelin_daily.min$srt]),
         event_end = date(eventdelin_daily$date[eventdelin_daily.min$end]),
         event_max_date = date(eventdelin_daily$date[eventdelin_daily.min$which.max])) %>%
  select(event_start, event_end, event_max_date, event_max_Q = max)
write_csv(eventdelin_daily.min, './DataFiles/hydro_data/event_delineations/min_daily.csv')

startyear <- 2019
endyear <- 2019

ggplot() +
  geom_rect(data = eventdelin_daily.min, aes(xmin = event_start, xmax = event_end), 
            ymin = -200, ymax = 1.05*max(eventdelin_daily$riverQ_linterp), fill = 'grey90', color = 'grey50', linetype = 'dashed') +
  geom_text(data = eventdelin_daily.min, aes(x = event_start + 1, label = row_number(event_start)), y = 0.9*max(eventdelin_daily$riverQ_linterp)) +
  geom_line(data = eventdelin_daily, aes(x = date, y = riverQ_linterp)) +
  theme_classic() +
  labs(title = ifelse(startyear == endyear, 
                      str_c('Local Minima, ', 'Water Year ', startyear), 
                      str_c('Local Minima, ', 'Water Years ', startyear, '-', endyear)),
       x = NULL) +
  scale_x_date(limits = as_date(c(str_c(startyear-1, '-10-01'), str_c(endyear, '-09-30'))))



# Event Delineation - Maxima Method ---------------------------------------

#*** 15 minute frequency -----------------------------------------------------

deltay <- 20
deltax <- 3
threshold <- 2

eventdelin_15min.max <- eventMaxima(eventdelin_15min$qf, 
                                  delta.y = deltay, delta.x = deltax*96, threshold = threshold) 
eventdelin_15min.max <- eventdelin_15min.max %>%
  mutate(event_start = as_datetime(eventdelin_15min$dateTime[eventdelin_15min.max$srt]),
         event_end = as_datetime(eventdelin_15min$dateTime[eventdelin_15min.max$end]),
         event_max_date = as_datetime(eventdelin_15min$dateTime[eventdelin_15min.max$which.max])) %>%
  select(event_start, event_end, event_max_date, event_max_Q = max)
write_csv(eventdelin_15min.max, './DataFiles/hydro_data/event_delineations/max_15min.csv')

startyear <- 2019
endyear <- 2019

ggplot() +
  geom_rect(data = eventdelin_15min.max, aes(xmin = event_start, xmax = event_end), 
            ymin = -200, ymax = 1.05*max(eventdelin_15min$riverQ_linterp), fill = 'grey90', color = 'grey50', linetype = 'dashed') +
  geom_text(data = eventdelin_15min.max, aes(x = event_start, label = row_number(event_start)), y = 0.9*max(eventdelin_15min$riverQ_linterp)) +
  geom_line(data = eventdelin_15min, aes(x = dateTime, y = riverQ_linterp)) +
  theme_classic() +
  xlim(as_datetime(c(str_c(startyear-1, '-10-01'), str_c(endyear, '-09-30'))))


#*** Hourly frequency -----------------------------------------------------

deltay <- 1
deltax <- 3
threshold <- 10

eventdelin_hourly.max <- eventMaxima(eventdelin_hourly$qf, 
                                     delta.y = deltay, delta.x = deltax*24, threshold = threshold)
eventdelin_hourly.max <- eventdelin_hourly.max %>%
  mutate(event_start = as_datetime(eventdelin_hourly$dateHour[eventdelin_hourly.max$srt]),
         event_end = as_datetime(eventdelin_hourly$dateHour[eventdelin_hourly.max$end]),
         event_max_date = as_datetime(eventdelin_hourly$dateHour[eventdelin_hourly.max$which.max])) %>%
  select(event_start, event_end, event_max_date, event_max_Q = max)
write_csv(eventdelin_hourly.max, './DataFiles/hydro_data/event_delineations/max_hourly.csv')

startyear <- 2018
endyear <- 2018

ggplot() +
  geom_rect(data = eventdelin_hourly.max, aes(xmin = event_start, xmax = event_end), 
            ymin = -200, ymax = 1.05*max(eventdelin_hourly$riverQ_linterp), fill = 'grey90', color = 'grey50', linetype = 'dashed') +
  geom_text(data = eventdelin_hourly.max, aes(x = event_start, label = row_number(event_start)), y = 0.9*max(eventdelin_hourly$riverQ_linterp)) +
  geom_line(data = eventdelin_hourly, aes(x = dateHour, y = riverQ_linterp)) +
  theme_classic() +
  xlim(as_datetime(c(str_c(startyear-1, '-10-01'), str_c(endyear, '-09-30'))))


#*** Daily frequency -----------------------------------------------------

deltay <- 40
deltax <- 3
threshold <- -1

eventdelin_daily.max <- eventMaxima(eventdelin_daily$qf, 
                                     delta.y = deltay, delta.x = deltax, threshold = threshold)
eventdelin_daily.max <- eventdelin_daily.max %>%
  mutate(event_start = date(eventdelin_daily$date[eventdelin_daily.max$srt]),
         event_end = date(eventdelin_daily$date[eventdelin_daily.max$end]),
         event_max_date = date(eventdelin_daily$date[eventdelin_daily.max$which.max])) %>%
  select(event_start, event_end, event_max_date, event_max_Q = max)
write_csv(eventdelin_daily.max, './DataFiles/hydro_data/event_delineations/max_daily.csv')

startyear <- 2021
endyear <- 2021

ggplot() +
  geom_rect(data = eventdelin_daily.max, aes(xmin = event_start, xmax = event_end), 
            ymin = -200, ymax = 1.05*max(eventdelin_daily$riverQ_linterp), fill = 'grey90', color = 'grey50', linetype = 'dashed') +
  geom_text(data = eventdelin_daily.max, aes(x = event_start + 1, label = row_number(event_start)), y = 0.9*max(eventdelin_daily$riverQ_linterp)) +
  geom_line(data = eventdelin_daily, aes(x = date, y = riverQ_linterp)) +
  theme_classic() +
  labs(title = ifelse(startyear == endyear, 
                      str_c('Local Maxima, ', 'Water Year ', startyear), 
                      str_c('Local Maxima, ', 'Water Years ', startyear, '-', endyear)),
       x = NULL) +
  scale_x_date(limits = as_date(c(str_c(startyear-1, '-10-01'), str_c(endyear, '-09-30'))))


# Event Delineation - Peaks over Threshold Method -------------------------

#*** 15 minute frequency -----------------------------------------------------

threshold <- 50
mindiff <- 1

eventdelin_15min.pot <- eventPOT(eventdelin_15min$qf, 
                                  threshold = threshold, min.diff = mindiff) 
eventdelin_15min.pot <- eventdelin_15min.pot %>%
  mutate(event_start = as_datetime(eventdelin_15min$dateTime[eventdelin_15min.pot$srt]),
         event_end = as_datetime(eventdelin_15min$dateTime[eventdelin_15min.pot$end]),
         event_max_date = as_datetime(eventdelin_15min$dateTime[eventdelin_15min.pot$which.max])) %>%
  select(event_start, event_end, event_max_date, event_max_Q = max)
write_csv(eventdelin_15min.pot, './DataFiles/hydro_data/event_delineations/pot_15min.csv')

startyear <- 2019
endyear <- 2019

ggplot() +
  geom_rect(data = eventdelin_15min.pot, aes(xmin = event_start, xmax = event_end), 
            ymin = -200, ymax = 1.05*max(eventdelin_15min$riverQ_linterp), fill = 'grey90', color = 'grey50', linetype = 'dashed') +
  geom_text(data = eventdelin_15min.pot, aes(x = event_start, label = row_number(event_start)), y = 0.9*max(eventdelin_15min$riverQ_linterp)) +
  geom_line(data = eventdelin_15min, aes(x = dateTime, y = riverQ_linterp)) +
  theme_classic() +
  labs(title = ifelse(startyear == endyear, 
                      str_c('Peaks over Threshold, ', 'Water Year ', startyear), 
                      str_c('Peaks over Threshold, ', 'Water Years ', startyear, '-', endyear)),
       x = NULL) +
  scale_x_datetime(limits = as_datetime(c(str_c(startyear-1, '-10-01'), str_c(endyear, '-09-30'))))


#*** Hourly frequency -----------------------------------------------------

threshold <- 50
mindiff <- 1

eventdelin_hourly.pot <- eventPOT(eventdelin_hourly$qf, 
                                  threshold = threshold, min.diff = mindiff)
eventdelin_hourly.pot <- eventdelin_hourly.pot %>%
  mutate(event_start = as_datetime(eventdelin_hourly$dateHour[eventdelin_hourly.pot$srt]),
         event_end = as_datetime(eventdelin_hourly$dateHour[eventdelin_hourly.pot$end]),
         event_max_date = as_datetime(eventdelin_hourly$dateHour[eventdelin_hourly.pot$which.max])) %>%
  select(event_start, event_end, event_max_date, event_max_Q = max)
write_csv(eventdelin_hourly.pot, './DataFiles/hydro_data/event_delineations/pot_hourly.csv')

startyear <- 2019
endyear <- 2019

ggplot() +
  geom_rect(data = eventdelin_hourly.pot, aes(xmin = event_start, xmax = event_end), 
            ymin = -200, ymax = 1.05*max(eventdelin_hourly$riverQ_linterp), fill = 'grey90', color = 'grey50', linetype = 'dashed') +
  geom_text(data = eventdelin_hourly.pot, aes(x = event_start, label = row_number(event_start)), y = 0.9*max(eventdelin_hourly$riverQ_linterp)) +
  geom_line(data = eventdelin_hourly, aes(x = dateHour, y = riverQ_linterp)) +
  theme_classic() +
  labs(title = ifelse(startyear == endyear, 
                      str_c('Peaks over Threshold, ', 'Water Year ', startyear), 
                      str_c('Peaks over Threshold, ', 'Water Years ', startyear, '-', endyear)),
       x = NULL) +
  scale_x_datetime(limits = as_datetime(c(str_c(startyear-1, '-10-01'), str_c(endyear, '-09-30'))))


#*** Daily frequency -----------------------------------------------------

threshold <- 50
mindiff <- 1

eventdelin_daily.pot <- eventPOT(eventdelin_daily$qf, 
                                  threshold = threshold, min.diff = mindiff)
eventdelin_daily.pot <- eventdelin_daily.pot %>%
  mutate(event_start = date(eventdelin_daily$date[eventdelin_daily.pot$srt]),
         event_end = date(eventdelin_daily$date[eventdelin_daily.pot$end]),
         event_max_date = date(eventdelin_daily$date[eventdelin_daily.pot$which.max])) %>%
  select(event_start, event_end, event_max_date, event_max_Q = max)
write_csv(eventdelin_daily.pot, './DataFiles/hydro_data/event_delineations/pot_daily.csv')

startyear <- 2019
endyear <- 2019

ggplot() +
  geom_rect(data = eventdelin_daily.pot, aes(xmin = event_start, xmax = event_end), 
            ymin = -200, ymax = 1.05*max(eventdelin_daily$riverQ_linterp), fill = 'grey90', color = 'grey50', linetype = 'dashed') +
  geom_text(data = eventdelin_daily.pot, aes(x = event_start + 1, label = row_number(event_start)), y = 0.9*max(eventdelin_daily$riverQ_linterp)) +
  geom_line(data = eventdelin_daily, aes(x = date, y = riverQ_linterp)) +
  theme_classic() +
  labs(title = ifelse(startyear == endyear, 
                      str_c('Peaks over Threshold, ', 'Water Year ', startyear), 
                      str_c('Peaks over Threshold, ', 'Water Years ', startyear, '-', endyear)),
       x = NULL) +
  scale_x_date(limits = as_date(c(str_c(startyear-1, '-10-01'), str_c(endyear, '-09-30'))))
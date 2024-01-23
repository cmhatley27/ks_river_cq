library(tidyverse)
library(lubridate)
library(dataRetrieval)
library(hydroEvents)
library(rnoaa)
library(ggpubr)
library(FluMoDL)
library(zoo)
library(magrittr)
library(quantreg)


all_hydro_15min <- with_tz(read_csv('./DataFiles/hydro_data/all_hydro_15min.csv'), tz = 'America/Chicago')
all_hydro_daily <- read_csv('./DataFiles/hydro_data/all_hydro_daily.csv')



k <- baseflow_RecessionConstant(all_hydro_15min$riverQ_linterp)
BFImax <- 0.8 #baseflow_BFImax(all_hydro_15min$riverQ_linterp, k)
alpha <- 0.99
passes <- 12

baseflow_15min <- data.frame(
  dateTime = all_hydro_15min$dateTime,
  riverQ_linterp = all_hydro_15min$riverQ_linterp,
  #HYSEP_fixed = baseflow_HYSEP(all_hydro_15min$riverQ_linterp, 59756, 'fixed'),
  #HYSEP_sliding = baseflow_HYSEP(all_hydro_15min$riverQ_linterp, 59756, 'sliding'),
  #HYSEP_local = baseflow_HYSEP(all_hydro_15min$riverQ_linterp, 59756, 'local'),
  #UKIH = baseflow_UKIH(all_hydro_15min$riverQ_linterp, 'B'),
  #BFLOW = baseflow_BFLOW(all_hydro_15min$riverQ_linterp, alpha, passes),
  Eckhardt = baseflow_Eckhardt(all_hydro_15min$riverQ_linterp, BFImax, k)) %>%
  mutate(Eckhardt_bfi = Eckhardt/riverQ_linterp,
         Eckhardt_qf = riverQ_linterp - Eckhardt)


# BFI threshold method ----------------------------------------------------



ggplot(data = baseflow_15min) +
  geom_histogram(aes(x = Eckhardt_bfi), binwidth = 0.01) +
  geom_vline(xintercept = bfi_threshold, color = 'red')
summary(baseflow_15min$Eckhardt_bfi)

bfi_threshold <- 0.8

Eckhardt_events_bfi_threshold <- data.frame(
  start_dateTime = baseflow_15min$dateTime[baseflow_15min$Eckhardt_bfi< bfi_threshold & lag(baseflow_15min$Eckhardt_bfi) >= bfi_threshold],
  end_dateTime = baseflow_15min$dateTime[baseflow_15min$Eckhardt_bfi< bfi_threshold & lead(baseflow_15min$Eckhardt_bfi) >= bfi_threshold],
  event_length = rle(baseflow_15min$Eckhardt_bfi< bfi_threshold)$lengths[rle(baseflow_15min$Eckhardt_bfi< bfi_threshold)$values]
) %>%
  filter(event_length > 3)
Eckhardt_events$start_dateTime[1] <- baseflow_15min$dateTime[1]

startyear <- 2019
endyear <- 2019
ggplot(data = all_hydro_15min) +
  #Event limits and labels
  geom_rect(data = Eckhardt_events_bfi_threshold, aes(xmin = start_dateTime, xmax = end_dateTime), 
            ymin = -200, ymax = 1.05*max(baseflow_15min$riverQ_linterp), fill = 'grey90', color = 'grey50', linetype = 'dashed') +
  geom_text(data = Eckhardt_events_bfi_threshold, aes(x = start_dateTime, label = row_number(start_dateTime)), y = 0.9*max(baseflow_15min$riverQ_linterp)) +
  #Precip
  geom_col(aes(x = dateTime, y = belvue_precip*100), color = 'blue', alpha = 0.5) +
  #River Q
  geom_line(aes(x = dateTime, y = riverQ_linterp)) +
  ##Reservoir Outflows
  #geom_line(aes(x = dateTime, y = clintonQ_linterp), color = 'red') +
  #geom_line(aes(x = dateTime, y = milfordQ_linterp), color = 'gold') +
  #geom_line(aes(x = dateTime, y = perryQ_linterp), color = 'lightblue') +
  #geom_line(aes(x = dateTime, y = tuttleQ_linterp), color = 'purple') +
  geom_line(aes(x = dateTime, y = reservoir_sum_linterp), color = 'darkgreen') +
  ##Baseflows
  #geom_line(aes(x = dateTime, y = HYSEP_fixed), linetype = 'dashed', color = 'red') +
  #geom_line(aes(x = dateTime, y = HYSEP_sliding), linetype = 'dashed', color = 'blue') +
  #geom_line(aes(x = dateTime, y = HYSEP_local), linetype = 'dashed', color = 'purple') +
  #geom_line(aes(x = dateTime, y = UKIH), linetype = 'dashed', color = 'orange') +
  #geom_line(aes(x = dateTime, y = BFLOW), linetype = 'dashed', color = 'red') +
  geom_line(data = baseflow_15min, aes(x = dateTime, y = Eckhardt), linetype = 'dashed', color = 'black') +
  ##BFI target line
  #geom_line(data = baseflow_15min,aes(x = dateTime, y = riverQ_linterp*bfi_threshold), linetype = 'dashed', color = 'red') +
  #Theme and scale
  theme_classic() +
  scale_x_datetime(limits = as_datetime(c(str_c(startyear-1, '-10-01'), str_c(endyear, '-09-30'))),
                   date_breaks = '2 weeks', date_labels = '%d %b')


# Minima method -----------------------------------------------------------
deltay <- 30
deltax <- 96
threshold <- 15

Eckhardt_events_minima <- eventMinima(baseflow_15min$Eckhardt_qf, deltay, deltax, threshold)

Eckhardt_events_minima <- Eckhardt_events_minima %>%
  mutate(start_dateTime = as_datetime(baseflow_15min$dateTime[Eckhardt_events_minima$srt]),
         end_dateTime = as_datetime(baseflow_15min$dateTime[Eckhardt_events_minima$end]),
         event_max_date = as_datetime(baseflow_15min$dateTime[Eckhardt_events_minima$which.max])) %>%
  select(start_dateTime, end_dateTime, event_max_date, event_max_Q = max)

startyear <- 2019
endyear <- 2019
ggplot(data = baseflow_15min) +
  #Event limits and labels
  geom_rect(data = Eckhardt_events_minima, aes(xmin = start_dateTime, xmax = end_dateTime), 
            ymin = -200, ymax = 1.05*max(baseflow_15min$riverQ_linterp), fill = 'grey90', color = 'grey50', linetype = 'dashed') +
  geom_text(data = Eckhardt_events_minima, aes(x = start_dateTime, label = row_number(start_dateTime)), y = 0.9*max(baseflow_15min$riverQ_linterp)) +
  #River Q
  geom_line(aes(x = dateTime, y = riverQ_linterp)) +
  ##Baseflows
  #geom_line(aes(x = dateTime, y = HYSEP_fixed), linetype = 'dashed', color = 'red') +
  #geom_line(aes(x = dateTime, y = HYSEP_sliding), linetype = 'dashed', color = 'blue') +
  #geom_line(aes(x = dateTime, y = HYSEP_local), linetype = 'dashed', color = 'purple') +
  #geom_line(aes(x = dateTime, y = UKIH), linetype = 'dashed', color = 'orange') +
  #geom_line(aes(x = dateTime, y = BFLOW), linetype = 'dashed', color = 'red') +
  geom_line(aes(x = dateTime, y = Eckhardt), linetype = 'dashed', color = 'black') +
  #Theme and scale
  theme_classic() +
  scale_x_datetime(limits = as_datetime(c(str_c(startyear-1, '-10-01'), str_c(endyear, '-09-30'))))
  
  


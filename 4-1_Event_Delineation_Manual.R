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


# Attempt 1 - use daily bfi interpolated to 15min -------------------------

alpha <- 0.925
passes <- 3
reflected <- 30

daily_bf_sep <- all_hydro_daily %>%
  mutate(bf = baseflowB(all_hydro_daily$riverQ_linterp, alpha = alpha, passes = passes, r = reflected)[["bf"]],
         qf = riverQ_linterp - bf,
         bfi = bf/riverQ_linterp) %>%
  select(dateTime = date, bf, bfi)

daily_to_15min <- all_hydro_15min %>%
  left_join(daily_bf_sep)
daily_to_15min$bfi[nrow(daily_to_15min)] <- daily_to_15min$bfi[nrow(daily_to_15min) - 95]
daily_to_15min$bf[nrow(daily_to_15min)] <- daily_to_15min$bf[nrow(daily_to_15min) - 95]
daily_to_15min <- daily_to_15min %>%
  mutate(bfi = linterp(bfi, max_allow = NULL),
         bf = linterp(bf, max_allow = NULL)) %>%
  mutate(bf_adjust = ifelse(bf > riverQ_linterp, riverQ_linterp, bf),
         bfi_calc = bf_adjust/riverQ_linterp)
startyear <- 2019
endyear <- 2019

ggplot(data = daily_to_15min) +
  geom_point(aes(x = bfi, y = bfi_calc))

ggplot(data = daily_to_15min) +
  geom_line(aes(x = dateTime, y = riverQ_linterp)) +
  geom_line(aes(x = dateTime, y = bf), linetype = 'dashed') +
  theme_classic() +
  scale_x_datetime(limits = as_datetime(c(str_c(startyear-1, '-10-01'), str_c(endyear, '-09-30'))))

bfi_threshold <- 0.85

event_delin_bfi <- data.frame(
  start_dateTime = daily_to_15min$dateTime[daily_to_15min$bfi < 0.85 & lag(daily_to_15min$bfi) >= 0.85],
  end_dateTime = daily_to_15min$dateTime[daily_to_15min$bfi < 0.85 & lead(daily_to_15min$bfi) >= 0.85],
  event_length = rle(daily_to_15min$bfi < 0.85)$lengths[rle(daily_to_15min$bfi < 0.85)$values]
)
event_delin_bfi$start_dateTime[1] <- daily_to_15min$dateTime[1]

startyear <- 2019
endyear <- 2019

ggplot(data = daily_to_15min) +
  #Event limits and labels
  geom_rect(data = event_delin_bfi, aes(xmin = start_dateTime, xmax = end_dateTime), 
            ymin = -200, ymax = 1.05*max(daily_to_15min$riverQ_linterp), fill = 'grey90', color = 'grey50', linetype = 'dashed') +
  geom_text(data = event_delin_bfi, aes(x = start_dateTime, label = row_number(start_dateTime)), y = 0.9*max(daily_to_15min$riverQ_linterp)) +
  #River Q and baseflow
  geom_line(aes(x = dateTime, y = riverQ_linterp)) +
  geom_line(aes(x = dateTime, y = bf), linetype = 'dashed') +
  theme_classic() +
  labs(title = ifelse(startyear == endyear, 
                      str_c('BFI Method, ', 'Water Year ', startyear), 
                      str_c('BFI Method, ', 'Water Years ', startyear, '-', endyear)),
       x = NULL) +
  scale_x_datetime(limits = as_datetime(c(str_c(startyear-1, '-10-01'), str_c(endyear, '-09-30'))))



# Using bfi_calc ----------------------------------------------------------

bfi_threshold <- 0.85

event_delin_bfi_calc<- data.frame(
  start_dateTime = daily_to_15min$dateTime[daily_to_15min$bfi_calc< bfi_threshold & lag(daily_to_15min$bfi_calc) >= bfi_threshold],
  end_dateTime = daily_to_15min$dateTime[daily_to_15min$bfi_calc< bfi_threshold & lead(daily_to_15min$bfi_calc) >= bfi_threshold],
  event_length = rle(daily_to_15min$bfi_calc< bfi_threshold)$lengths[rle(daily_to_15min$bfi_calc< bfi_threshold)$values]
)
event_delin_bfi_calc$start_dateTime[1] <- daily_to_15min$dateTime[1]

startyear <- 2019
endyear <- 2019

ggplot(data = daily_to_15min) +
  #Event limits and labels
  geom_rect(data = event_delin_bfi_calc, aes(xmin = start_dateTime, xmax = end_dateTime), 
            ymin = -200, ymax = 1.05*max(daily_to_15min$riverQ_linterp), fill = 'grey90', color = 'grey50', linetype = 'dashed') +
  geom_text(data = event_delin_bfi_calc, aes(x = start_dateTime, label = row_number(start_dateTime)), y = 0.9*max(daily_to_15min$riverQ_linterp)) +
  #River Q and baseflow
  geom_line(aes(x = dateTime, y = riverQ_linterp)) +
  geom_line(aes(x = dateTime, y = bf_adjust), linetype = 'dashed') +
  theme_classic() +
  labs(title = ifelse(startyear == endyear, 
                      str_c('BFI Method, ', 'Water Year ', startyear), 
                      str_c('BFI Method, ', 'Water Years ', startyear, '-', endyear)),
       x = NULL) +
  scale_x_datetime(limits = as_datetime(c(str_c(startyear-1, '-10-01'), str_c(endyear, '-09-30'))))



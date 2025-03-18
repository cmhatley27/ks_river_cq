library(tidyverse)
library(dataRetrieval)
library(rnoaa)
library(FluMoDL)
library(jsonlite)
library(curl)

# Master date(time) vectors ------------------------------------------------------

study_dates <- data.frame(date = seq(ymd('2013-10-01'), ymd('2022-09-30'), by = 'days'))
study_dateTimes <- data.frame(dateTime = seq.POSIXt(as.POSIXlt('2013-10-01 00:00:00'), as.POSIXlt('2022-09-30 23:45:00'), by = '15 min'))

# Download and prepare DeSoto gage data -----------------------------------

#download at 15min frequency
DeSoto_raw <- readNWISuv("06892350", c("00060", "99133"), "2013-10-01", "2022-09-30", tz = 'America/Chicago')

#clean up and rename columns, interpolate between missing Q and N
DeSoto_15min <- left_join(x = study_dateTimes, y = with_tz(DeSoto_raw, 'America/Chicago'), by = 'dateTime') %>%
  mutate(Q = X_00060_00000/35.314666212661,
         N = X_99133_00000,
         wateryear = if_else(month(dateTime) >= 10, year(dateTime) + 1, year(dateTime)),
         yday = yday(dateTime)) %>%
  mutate(Q_linterp = linterp(Q, max_allow = NULL),
         N_linterp = linterp(N, max_allow = NULL)) %>%
  select(dateTime, wateryear, yday, Q, Q_linterp, N, N_linterp)
write_csv(DeSoto_15min, './DataFiles/hydro_data/DeSoto_gage/DeSoto_15min.csv')

#aggregate to daily
DeSoto_daily <- DeSoto_15min %>%
  group_by(date = date(dateTime)) %>%
  summarise(Q = mean(Q, na.rm=TRUE),
            Q_linterp = mean(Q_linterp),
            N = mean(N, na.rm = TRUE)) %>%
  mutate(wateryear = if_else(month(date) >= 10, year(date) + 1, year(date)),
         yday = yday(date),
         Q = replace(Q, is.nan(Q), NA),
         N = replace(N, is.nan(N), NA)) %>%
  select(date, wateryear, yday, Q, Q_linterp, N)
write_csv(DeSoto_daily, './DataFiles/hydro_data/DeSoto_gage/DeSoto_daily.csv')


# Download and prepare reservoir outflow data -----------------------------

# Clinton
clinton_gage_raw <- readNWISuv("06891500", "00060", "2013-10-01", "2022-09-30", tz = 'America/Chicago')

clinton_gage_15min <- left_join(x = study_dateTimes, y = clinton_gage_raw, by = 'dateTime') %>%
  mutate(Q = X_00060_00000/35.314666212661,
         wateryear = if_else(month(dateTime) >= 10, year(dateTime) + 1, year(dateTime)),
         yday = yday(dateTime)) %>%
  mutate(Q_linterp = linterp(Q, max_allow = NULL)) %>%
  select(dateTime, wateryear, yday, Q, Q_linterp)
write_csv(clinton_gage_15min, './DataFiles/hydro_data/reservoir_outflows/clinton_gage_15min.csv')

clinton_gage_daily <- clinton_gage_15min %>%
  group_by(date = date(dateTime)) %>%
  summarise(Q = mean(Q, na.rm=TRUE),
            Q_linterp = mean(Q_linterp)) %>%
  mutate(wateryear = if_else(month(date) >= 10, year(date) + 1, year(date)),
         yday = yday(date),
         Q = replace(Q, is.nan(Q), NA)) %>%
  select(date, wateryear, yday, Q, Q_linterp)
write_csv(clinton_gage_daily, './DataFiles/hydro_data/reservoir_outflows/clinton_gage_daily.csv')


# Milford
milford_gage_raw <- readNWISuv("06857100", "00060", "2013-10-01", "2022-09-30", tz = 'America/Chicago')

milford_gage_15min <- left_join(x = study_dateTimes, y = milford_gage_raw, by = 'dateTime') %>%
  mutate(Q = X_00060_00000/35.314666212661,
         wateryear = if_else(month(dateTime) >= 10, year(dateTime) + 1, year(dateTime)),
         yday = yday(dateTime)) %>%
  mutate(Q_linterp = linterp(Q, max_allow = NULL)) %>%
  select(dateTime, wateryear, yday, Q, Q_linterp)
write_csv(milford_gage_15min, './DataFiles/hydro_data/reservoir_outflows/milford_gage_15min.csv')

milford_gage_daily <- milford_gage_15min %>%
  group_by(date = date(dateTime)) %>%
  summarise(Q = mean(Q, na.rm=TRUE),
            Q_linterp = mean(Q_linterp)) %>%
  mutate(wateryear = if_else(month(date) >= 10, year(date) + 1, year(date)),
         yday = yday(date),
         Q = replace(Q, is.nan(Q), NA)) %>%
  select(date, wateryear, yday, Q, Q_linterp)
write_csv(milford_gage_daily, './DataFiles/hydro_data/reservoir_outflows/milford_gage_daily.csv')


# Perry

perry_gage_raw <- readNWISuv("06890900", "00060", "2013-10-01", "2022-09-30", tz = 'America/Chicago')

perry_gage_15min <- left_join(x = study_dateTimes, y = perry_gage_raw, by = 'dateTime') %>%
  mutate(Q = X_00060_00000/35.314666212661,
         wateryear = if_else(month(dateTime) >= 10, year(dateTime) + 1, year(dateTime)),
         yday = yday(dateTime)) %>%
  mutate(Q_linterp = linterp(Q, max_allow = NULL)) %>%
  select(dateTime, wateryear, yday, Q, Q_linterp)
write_csv(perry_gage_15min, './DataFiles/hydro_data/reservoir_outflows/perry_gage_15min.csv')

perry_gage_daily <- perry_gage_15min %>%
  group_by(date = date(dateTime)) %>%
  summarise(Q = mean(Q, na.rm=TRUE),
            Q_linterp = mean(Q_linterp)) %>%
  mutate(wateryear = if_else(month(date) >= 10, year(date) + 1, year(date)),
         yday = yday(date),
         Q = replace(Q, is.nan(Q), NA)) %>%
  select(date, wateryear, yday, Q, Q_linterp)
write_csv(perry_gage_daily, './DataFiles/hydro_data/reservoir_outflows/perry_gage_daily.csv')



# Tuttle

tuttle_gage_raw <- readNWISuv("06887000", "00060", "2013-10-01", "2022-09-30", tz = 'America/Chicago')

tuttle_gage_15min <- left_join(x = study_dateTimes, y = tuttle_gage_raw, by = 'dateTime') %>%
  mutate(Q = X_00060_00000/35.314666212661,
         wateryear = if_else(month(dateTime) >= 10, year(dateTime) + 1, year(dateTime)),
         yday = yday(dateTime)) %>%
  mutate(Q_linterp = linterp(Q, max_allow = NULL)) %>%
  select(dateTime, wateryear, yday, Q, Q_linterp)
write_csv(tuttle_gage_15min, './DataFiles/hydro_data/reservoir_outflows/tuttle_gage_15min.csv')

tuttle_gage_daily <- tuttle_gage_15min %>%
  group_by(date = date(dateTime)) %>%
  summarise(Q = mean(Q, na.rm=TRUE),
            Q_linterp = mean(Q_linterp)) %>%
  mutate(wateryear = if_else(month(date) >= 10, year(date) + 1, year(date)),
         yday = yday(date),
         Q = replace(Q, is.nan(Q), NA)) %>%
  select(date, wateryear, yday, Q, Q_linterp)
write_csv(tuttle_gage_daily, './DataFiles/hydro_data/reservoir_outflows/tuttle_gage_daily.csv')


# Kanopolis

kanopolis_gage_raw <- readNWISuv("06865500", "00060", "2013-10-01", "2022-09-30", tz = 'America/Chicago')

kanopolis_gage_15min <- left_join(x = study_dateTimes, y = kanopolis_gage_raw, by = 'dateTime') %>%
  mutate(Q = X_00060_00000/35.314666212661,
         wateryear = if_else(month(dateTime) >= 10, year(dateTime) + 1, year(dateTime)),
         yday = yday(dateTime)) %>%
  mutate(Q_linterp = linterp(Q, max_allow = NULL)) %>%
  select(dateTime, wateryear, yday, Q, Q_linterp)
write_csv(kanopolis_gage_15min, './DataFiles/hydro_data/reservoir_outflows/kanopolis_gage_15min.csv')

kanopolis_gage_daily <- kanopolis_gage_15min %>%
  group_by(date = date(dateTime)) %>%
  summarise(Q = mean(Q, na.rm=TRUE),
            Q_linterp = mean(Q_linterp)) %>%
  mutate(wateryear = if_else(month(date) >= 10, year(date) + 1, year(date)),
         yday = yday(date),
         Q = replace(Q, is.nan(Q), NA)) %>%
  select(date, wateryear, yday, Q, Q_linterp)
write_csv(kanopolis_gage_daily, './DataFiles/hydro_data/reservoir_outflows/kanopolis_gage_daily.csv')


# Wilson

wilson_gage_raw <- readNWISuv("06868200", "00060", "2013-10-01", "2022-09-30", tz = 'America/Chicago')

wilson_gage_15min <- left_join(x = study_dateTimes, y = wilson_gage_raw, by = 'dateTime') %>%
  mutate(Q = X_00060_00000/35.314666212661,
         wateryear = if_else(month(dateTime) >= 10, year(dateTime) + 1, year(dateTime)),
         yday = yday(dateTime)) %>%
  mutate(Q_linterp = linterp(Q, max_allow = NULL)) %>%
  select(dateTime, wateryear, yday, Q, Q_linterp)
write_csv(wilson_gage_15min, './DataFiles/hydro_data/reservoir_outflows/wilson_gage_15min.csv')

wilson_gage_daily <- wilson_gage_15min %>%
  group_by(date = date(dateTime)) %>%
  summarise(Q = mean(Q, na.rm=TRUE),
            Q_linterp = mean(Q_linterp)) %>%
  mutate(wateryear = if_else(month(date) >= 10, year(date) + 1, year(date)),
         yday = yday(date),
         Q = replace(Q, is.nan(Q), NA)) %>%
  select(date, wateryear, yday, Q, Q_linterp)
write_csv(wilson_gage_daily, './DataFiles/hydro_data/reservoir_outflows/wilson_gage_daily.csv')


# Waconda

waconda_gage_raw <- readNWISuv("06875900", "00060", "2013-10-01", "2022-09-30", tz = 'America/Chicago')

waconda_gage_15min <- left_join(x = study_dateTimes, y = waconda_gage_raw, by = 'dateTime') %>%
  mutate(Q = X_00060_00000/35.314666212661,
         wateryear = if_else(month(dateTime) >= 10, year(dateTime) + 1, year(dateTime)),
         yday = yday(dateTime)) %>%
  mutate(Q_linterp = linterp(Q, max_allow = NULL)) %>%
  select(dateTime, wateryear, yday, Q, Q_linterp)
write_csv(waconda_gage_15min, './DataFiles/hydro_data/reservoir_outflows/waconda_gage_15min.csv')

waconda_gage_daily <- waconda_gage_15min %>%
  group_by(date = date(dateTime)) %>%
  summarise(Q = mean(Q, na.rm=TRUE),
            Q_linterp = mean(Q_linterp)) %>%
  mutate(wateryear = if_else(month(date) >= 10, year(date) + 1, year(date)),
         yday = yday(date),
         Q = replace(Q, is.nan(Q), NA)) %>%
  select(date, wateryear, yday, Q, Q_linterp)
write_csv(waconda_gage_daily, './DataFiles/hydro_data/reservoir_outflows/waconda_gage_daily.csv')


# Aggregate all hydro data -----------------------------------------------

DeSoto_15min <- read_csv('./DataFiles/hydro_data/DeSoto_gage/DeSoto_15min.csv')
clinton_gage_15min <- read_csv('./DataFiles/hydro_data/reservoir_outflows/clinton_gage_15min.csv')
milford_gage_15min <- read_csv('./DataFiles/hydro_data/reservoir_outflows/milford_gage_15min.csv')
perry_gage_15min <- read_csv('./DataFiles/hydro_data/reservoir_outflows/perry_gage_15min.csv')
tuttle_gage_15min <- read_csv('./DataFiles/hydro_data/reservoir_outflows/tuttle_gage_15min.csv')
kanopolis_gage_15min <- read_csv('./DataFiles/hydro_data/reservoir_outflows/kanopolis_gage_15min.csv')
wilson_gage_15min <- read_csv('./DataFiles/hydro_data/reservoir_outflows/wilson_gage_15min.csv')
waconda_gage_15min <- read_csv('./DataFiles/hydro_data/reservoir_outflows/waconda_gage_15min.csv')

#function for flagging measurement points that occur within a data gap longer
#than some specified threshold (here we use 12 15-minute measurements = 3 hours)
gap_flagger <- function(x, max_gap = 12){
  gap_start <- which(is.na(x) & !is.na(lag(x)))
  gap_end <- which(is.na(x) & !is.na(lead(x)))
  gap_length <- gap_end - gap_start
  
  gap_start_trimmed <- gap_start[gap_length > max_gap]
  gap_end_trimmed <- gap_end[gap_length > max_gap]
  
  in_gap_indices <- double()
  for(gap in 1:length(gap_start_trimmed)){
    in_gap_indices <- c(in_gap_indices, gap_start_trimmed[gap]:gap_end_trimmed[gap])
  }
  x_in_gap <- rep(0,length(x))
  x_in_gap[in_gap_indices] <- 1
  return(x_in_gap)
}

all_hydro_15min <- DeSoto_15min %>%
  rename(riverQ = Q, riverQ_linterp = Q_linterp) %>%
  left_join(clinton_gage_15min) %>% rename(clintonQ = Q, clintonQ_linterp = Q_linterp) %>%
  left_join(milford_gage_15min) %>% rename(milfordQ = Q, milfordQ_linterp = Q_linterp) %>% 
  left_join(perry_gage_15min) %>% rename(perryQ = Q, perryQ_linterp = Q_linterp) %>% 
  left_join(tuttle_gage_15min) %>% rename(tuttleQ = Q, tuttleQ_linterp = Q_linterp) %>%
  left_join(kanopolis_gage_15min) %>% rename(kanopolisQ = Q, kanopolisQ_linterp = Q_linterp) %>%
  left_join(waconda_gage_15min) %>% rename(wacondaQ = Q, wacondaQ_linterp = Q_linterp) %>%
  left_join(wilson_gage_15min) %>% rename(wilsonQ = Q, wilsonQ_linterp = Q_linterp) %>%
  mutate(reservoir_sum = select(., clintonQ, milfordQ, perryQ, tuttleQ, kanopolisQ, wacondaQ, wilsonQ) %>% rowSums(na.rm = TRUE),
         reservoir_sum_linterp = select(., clintonQ_linterp, milfordQ_linterp, perryQ_linterp, tuttleQ_linterp, kanopolisQ_linterp, wacondaQ_linterp, wilsonQ_linterp) %>% rowSums(na.rm = TRUE)) %>%
  mutate(across(c(ends_with('Q'),N), ~gap_flagger(.x), .names = '{.col}_ingap'))
write_csv(all_hydro_15min, './DataFiles/hydro_data/all_hydro_15min.csv')


DeSoto_daily <- read_csv('./DataFiles/hydro_data/DeSoto_gage/DeSoto_daily.csv')
clinton_gage_daily <- read_csv('./DataFiles/hydro_data/reservoir_outflows/clinton_gage_daily.csv')
milford_gage_daily <- read_csv('./DataFiles/hydro_data/reservoir_outflows/milford_gage_daily.csv')
perry_gage_daily <- read_csv('./DataFiles/hydro_data/reservoir_outflows/perry_gage_daily.csv')
tuttle_gage_daily <- read_csv('./DataFiles/hydro_data/reservoir_outflows/tuttle_gage_daily.csv')
kanopolis_gage_daily <- read_csv('./DataFiles/hydro_data/reservoir_outflows/kanopolis_gage_daily.csv')
wilson_gage_daily <- read_csv('./DataFiles/hydro_data/reservoir_outflows/wilson_gage_daily.csv')
waconda_gage_daily <- read_csv('./DataFiles/hydro_data/reservoir_outflows/waconda_gage_daily.csv')

all_hydro_daily <- DeSoto_daily %>%
  dplyr::rename(riverQ = Q, riverQ_linterp = Q_linterp) %>%
  left_join(clinton_gage_daily) %>% dplyr::rename(clintonQ = Q, clintonQ_linterp = Q_linterp) %>%
  left_join(milford_gage_daily) %>% dplyr::rename(milfordQ = Q, milfordQ_linterp = Q_linterp) %>% 
  left_join(perry_gage_daily) %>% dplyr::rename(perryQ = Q, perryQ_linterp = Q_linterp) %>% 
  left_join(tuttle_gage_daily) %>% dplyr::rename(tuttleQ = Q, tuttleQ_linterp = Q_linterp) %>%
  left_join(kanopolis_gage_daily) %>% dplyr::rename(kanopolisQ = Q, kanopolisQ_linterp = Q_linterp) %>%
  left_join(waconda_gage_daily) %>% dplyr::rename(wacondaQ = Q, wacondaQ_linterp = Q_linterp) %>%
  left_join(wilson_gage_daily) %>% dplyr::rename(wilsonQ = Q, wilsonQ_linterp = Q_linterp) %>%
  mutate(reservoir_sum = select(., clintonQ, milfordQ, perryQ, tuttleQ, kanopolisQ, wacondaQ, wilsonQ) %>% rowSums(na.rm = TRUE),
         reservoir_sum_linterp = select(., clintonQ_linterp, milfordQ_linterp, perryQ_linterp, tuttleQ_linterp, kanopolisQ_linterp, wacondaQ_linterp, wilsonQ_linterp) %>% rowSums(na.rm = TRUE))
write_csv(all_hydro_daily, './DataFiles/hydro_data/all_hydro_daily.csv')


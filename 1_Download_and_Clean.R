library(tidyverse)
library(dataRetrieval)
library(rnoaa)
library(FluMoDL)
library(jsonlite)
library(curl)

# Create Data Folders -----------------------------------------------------

if(!dir.exists('DataFiles')) {
  dir.create('DataFiles')
}

if(!dir.exists('DataFiles/hydro_data')) {
  dir.create('DataFiles/hydro_data')
}

if(!dir.exists('DataFiles/hydro_data/DeSoto_gage')) {
  dir.create('DataFiles/hydro_data/DeSoto_gage')
}

if(!dir.exists('DataFiles/hydro_data/lawrence_gage')) {
  dir.create('DataFiles/hydro_data/lawrence_gage')
}

if(!dir.exists('DataFiles/hydro_data/precip')) {
  dir.create('DataFiles/hydro_data/precip')
}

if(!dir.exists('DataFiles/hydro_data/reservoir_outflows')) {
  dir.create('DataFiles/hydro_data/reservoir_outflows')
}

if(!dir.exists('Figures')) {
  dir.create('Figures')
}


# Master date(time) vectors ------------------------------------------------------

study_dates <- data.frame(date = seq(ymd('2013-10-01'), ymd('2022-09-30'), by = 'days'))
study_dateTimes <- data.frame(dateTime = seq.POSIXt(as.POSIXlt('2013-10-01 00:00:00'), as.POSIXlt('2022-09-30 23:45:00'), by = '15 min'))
study_dateTimes_extrayear <- data.frame(dateTime = seq.POSIXt(as.POSIXlt('2012-10-01 00:00:00'), as.POSIXlt('2022-09-30 23:45:00'), by = '15 min'))


# Download and prepare DeSoto gage data -----------------------------------

DeSoto_raw <- readNWISuv("06892350", c("00060", "99133"), "2013-10-01", "2022-09-30", tz = 'America/Chicago')
DeSoto_raw_extrayear <- readNWISuv("06892350", c("00060", "99133"), "2012-10-01", "2022-09-30", tz = 'America/Chicago')
write_csv(DeSoto_raw, './DataFiles/hydro_data/DeSoto_gage/DeSoto_raw.csv')
write_csv(DeSoto_raw_extrayear, './DataFiles/hydro_data/DeSoto_gage/DeSoto_raw_extrayear.csv')

DeSoto_raw <- with_tz(read_csv('./DataFiles/hydro_data/DeSoto_gage/DeSoto_raw.csv'), 'America/Chicago')
##All data - 15 minute frequency
DeSoto_15min <- left_join(x = study_dateTimes, y = with_tz(DeSoto_raw, 'America/Chicago'), by = 'dateTime') %>%
  mutate(Q = X_00060_00000/35.314666212661,
         N = X_99133_00000,
         wateryear = if_else(month(dateTime) >= 10, year(dateTime) + 1, year(dateTime)),
         yday = yday(dateTime)) %>%
  mutate(Q_linterp = linterp(Q, max_allow = NULL),
         N_linterp = linterp(N, max_allow = NULL)) %>%
  select(dateTime, wateryear, yday, Q, Q_linterp, N, N_linterp)
write_csv(DeSoto_15min, './DataFiles/hydro_data/DeSoto_gage/DeSoto_15min.csv')


##Hourly averages
DeSoto_hourly <- DeSoto_15min %>%
  group_by(dateHour = floor_date(dateTime, unit = 'hours')) %>%
  summarise(Q = mean(Q, na.rm=TRUE),
            Q_linterp = mean(Q_linterp),
            N = mean(N, na.rm = TRUE)) %>%
  mutate(wateryear = if_else(month(dateHour) >= 10, year(dateHour) + 1, year(dateHour)),
         yday = yday(dateHour),
         Q = replace(Q, is.nan(Q), NA),
         N = replace(N, is.nan(N), NA)) %>%
  select(dateHour, wateryear, yday, Q, Q_linterp, N)
write_csv(DeSoto_hourly, './DataFiles/hydro_data/DeSoto_gage/DeSoto_hourly.csv')

##Daily averages
#Pull in extra year of data for historical precip data
DeSoto_15min_extrayear <- left_join(x = study_dateTimes_extrayear, y = with_tz(DeSoto_raw_extrayear, 'America/Chicago'), by = 'dateTime') %>%
  mutate(Q = X_00060_00000/35.314666212661,
         N = X_99133_00000,
         wateryear = if_else(month(dateTime) >= 10, year(dateTime) + 1, year(dateTime)),
         yday = yday(dateTime)) %>%
  mutate(N = ifelse(dateTime == '2012-10-01 00:00:00', 0, N)) %>%
  mutate(Q_linterp = linterp(Q, max_allow = NULL),
         N_linterp = linterp(N, max_allow = NULL)) %>%
  select(dateTime, wateryear, yday, Q, Q_linterp, N, N_linterp)

DeSoto_daily <- DeSoto_15min_extrayear %>%
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



# Download and prepare Lawrence gage data ---------------------------------

lawrence_gage_raw <- readNWISuv("06891080", "00060", "2013-10-01", "2022-09-30", tz = 'America/Chicago')

lawrence_gage_15min <- left_join(x = study_dateTimes, y = lawrence_gage_raw, by = 'dateTime') %>%
  mutate(Q = X_00060_00000/35.314666212661,
         wateryear = if_else(month(dateTime) >= 10, year(dateTime) + 1, year(dateTime)),
         yday = yday(dateTime)) %>%
  mutate(Q_linterp = linterp(Q, max_allow = NULL)) %>%
  select(dateTime, wateryear, yday, Q, Q_linterp)
write_csv(lawrence_gage_15min, './DataFiles/hydro_data/lawrence_gage/lawrence_gage_15min.csv')


# Download and prepare precip data ----------------------------------------

noaa_stations <- ghcnd_stations()
ks_stations <- meteo_distance(noaa_stations, 39.0473, -95.6752, radius = 200) %>%
  filter(element == 'PRCP', last_year == 2022, first_year <= 2013)


#*** Lawrence - daily --------------------------------------------------------

lawrence_precip_raw <- meteo_pull_monitors("USW00003997", date_min = "2012-10-01", date_max = '2022-09-30', var = 'all')
lawrence_precip_daily <- left_join(x = study_dates, y = lawrence_precip_raw, by = 'date') %>%
  mutate(precip = prcp/10,
         precip = replace_na(precip, 0),
         wateryear = if_else(month(date) >= 10, year(date) + 1, year(date)),
         yday = yday(date)) %>%
  select(date, wateryear, yday, precip)
write_csv(lawrence_precip_daily, './DataFiles/hydro_data/precip/lawrence_precip_daily.csv')


#*** Topeka - daily --------------------------------------------------------

topeka_precip_raw <- meteo_pull_monitors("USW00013996", date_min = "2012-10-01", date_max = '2022-09-30', var = 'all')
topeka_precip_daily <- left_join(x = study_dates, y = topeka_precip_raw, by = 'date') %>%
  mutate(precip = prcp/10,
         precip = replace_na(precip, 0),
         wateryear = if_else(month(date) >= 10, year(date) + 1, year(date)),
         yday = yday(date)) %>%
  select(date, wateryear, yday, precip)
write_csv(topeka_precip_daily, './DataFiles/hydro_data/precip/topeka_precip_daily.csv')


#*** Manhattan - daily -------------------------------------------------------

manhattan_precip_raw <- meteo_pull_monitors("USW00003936", date_min = "2012-10-01", date_max = '2022-09-30', var = 'all')
manhattan_precip_daily <- left_join(x = study_dates, y = manhattan_precip_raw, by = 'date') %>%
  mutate(precip = prcp/10,
         precip = replace_na(precip, 0),
         wateryear = if_else(month(date) >= 10, year(date) + 1, year(date)),
         yday = yday(date)) %>%
  select(date, wateryear, yday, precip)
write_csv(manhattan_precip_daily, './DataFiles/hydro_data/precip/manhattan_precip_daily.csv')


#*** Salina - daily -------------------------------------------------------

salina_precip_raw <- meteo_pull_monitors("USW00003919", date_min = "2012-10-01", date_max = '2022-09-30', var = 'all')
salina_precip_daily <- left_join(x = study_dates, y = salina_precip_raw, by = 'date') %>%
  mutate(precip = prcp/10,
         precip = replace_na(precip, 0),
         wateryear = if_else(month(date) >= 10, year(date) + 1, year(date)),
         yday = yday(date)) %>%
  select(date, wateryear, yday, precip)
write_csv(salina_precip_daily, './DataFiles/hydro_data/precip/salina_precip_daily.csv')


#*** Belvue - 15 min ---------------------------------------------------------

belvue_precip_raw <- readNWISuv('06888350', '00045', "2012-10-01", "2022-09-30", tz = 'America/Chicago')
belvue_precip_15min <- left_join(x = study_dateTimes_extrayear, y = belvue_precip_raw, by = 'dateTime') %>%
  mutate(precip = X_00045_00000*25.4,
         precip = replace_na(precip, 0),
         wateryear = if_else(month(dateTime) >= 10, year(dateTime) + 1, year(dateTime)),
         yday = yday(dateTime)) %>%
  select(dateTime, wateryear, yday, precip)
write_csv(belvue_precip_15min, './DataFiles/hydro_data/precip/belvue_precip_15min.csv')

belvue_precip_daily <- belvue_precip_15min %>%
  group_by(date = date(dateTime)) %>%
  summarise(precip = sum(precip)) %>%
  mutate(wateryear = if_else(month(date) >= 10, year(date) + 1, year(date)),
         yday = yday(date))
write_csv(belvue_precip_daily, './DataFiles/hydro_data/precip/belvue_precip_daily.csv')


# Download and prepare reservoir outflow data -----------------------------

#*** USACE dam release data --------------------------------------------------

clinton_dam_data_raw <- fromJSON('https://water.usace.army.mil/a2w/CWMS_CRREL.cwms_data_api.get_report_json?p_location_id=5488030&p_parameter_type=Flow%3AStor%3APrecip%3AStage%3AElev&p_last=10&p_last_unit=years&p_unit_system=EN&p_format=JSON')
clinton_dam_outflow_raw <- clinton_dam_data_raw[['Outflow']][[3]]
clinton_dam_daily <- clinton_dam_outflow_raw %>%
  mutate(datetime = dmy_hms(time),
         date = as_date(datetime),
         outflow = value*0.028316847) %>%
  select(date, outflow) %>%
  filter(date < as_date('2022-10-01'))
write_csv(clinton_dam_daily, './DataFiles/hydro_data/reservoir_outflows/clinton_dam_daily.csv')

tuttle_dam_data_raw <- fromJSON('https://water.usace.army.mil/a2w/CWMS_CRREL.cwms_data_api.get_report_json?p_location_id=5722030&p_parameter_type=Flow%3AStor%3APrecip%3AStage%3AElev&p_last=10&p_last_unit=years&p_unit_system=EN&p_format=JSON')
tuttle_dam_outflow_raw <- tuttle_dam_data_raw[['Outflow']][[3]]
tuttle_dam_daily <- tuttle_dam_outflow_raw %>%
  mutate(datetime = dmy_hms(time),
         date = as_date(datetime),
         outflow = value*0.028316847) %>%
  select(date, outflow) %>%
  filter(date < as_date('2022-10-01'))
write_csv(tuttle_dam_daily, './DataFiles/hydro_data/reservoir_outflows/tuttle_dam_daily.csv')

perry_dam_data_raw <- fromJSON('https://water.usace.army.mil/a2w/CWMS_CRREL.cwms_data_api.get_report_json?p_location_id=5651030&p_parameter_type=Flow%3AStor%3APrecip%3AStage%3AElev&p_last=10&p_last_unit=years&p_unit_system=EN&p_format=JSON')
perry_dam_outflow_raw <- perry_dam_data_raw[['Outflow']][[3]]
perry_dam_daily <- perry_dam_outflow_raw %>%
  mutate(datetime = dmy_hms(time),
         date = as_date(datetime),
         outflow = value*0.028316847) %>%
  select(date, outflow) %>%
  filter(date < as_date('2022-10-01'))
write_csv(perry_dam_daily, './DataFiles/hydro_data/reservoir_outflows/perry_dam_daily.csv')

milford_dam_data_raw <- fromJSON('https://water.usace.army.mil/a2w/CWMS_CRREL.cwms_data_api.get_report_json?p_location_id=5604030&p_parameter_type=Flow%3AStor%3APrecip%3AStage%3AElev&p_last=10&p_last_unit=years&p_unit_system=EN&p_format=JSON')
milford_dam_outflow_raw <- milford_dam_data_raw[['Outflow']][[3]]
milford_dam_daily <- milford_dam_outflow_raw %>%
  mutate(datetime = dmy_hms(time),
         date = as_date(datetime),
         outflow = value*0.028316847) %>%
  select(date, outflow) %>%
  filter(date < as_date('2022-10-01'))
write_csv(milford_dam_daily, './DataFiles/hydro_data/reservoir_outflows/milford_dam_daily.csv')


#*** USGS downstream gage data ---------------------------------------------------------------

#*** *** Clinton -----------------------------------------------------------------

clinton_gage_raw <- readNWISuv("06891500", "00060", "2013-10-01", "2022-09-30", tz = 'America/Chicago')

clinton_gage_15min <- left_join(x = study_dateTimes, y = clinton_gage_raw, by = 'dateTime') %>%
  mutate(Q = X_00060_00000/35.314666212661,
         wateryear = if_else(month(dateTime) >= 10, year(dateTime) + 1, year(dateTime)),
         yday = yday(dateTime)) %>%
  mutate(Q_linterp = linterp(Q, max_allow = NULL)) %>%
  select(dateTime, wateryear, yday, Q, Q_linterp)
write_csv(clinton_gage_15min, './DataFiles/hydro_data/reservoir_outflows/clinton_gage_15min.csv')


clinton_gage_hourly <- clinton_gage_15min %>%
  group_by(dateHour = floor_date(dateTime, unit = 'hours')) %>%
  summarise(Q = mean(Q, na.rm=TRUE),
            Q_linterp = mean(Q_linterp)) %>%
  mutate(wateryear = if_else(month(dateHour) >= 10, year(dateHour) + 1, year(dateHour)),
         yday = yday(dateHour),
         Q = replace(Q, is.nan(Q), NA)) %>%
  select(dateHour, wateryear, yday, Q, Q_linterp)
write_csv(clinton_gage_hourly, './DataFiles/hydro_data/reservoir_outflows/clinton_gage_hourly.csv')


clinton_gage_daily <- clinton_gage_15min %>%
  group_by(date = date(dateTime)) %>%
  summarise(Q = mean(Q, na.rm=TRUE),
            Q_linterp = mean(Q_linterp)) %>%
  mutate(wateryear = if_else(month(date) >= 10, year(date) + 1, year(date)),
         yday = yday(date),
         Q = replace(Q, is.nan(Q), NA)) %>%
  select(date, wateryear, yday, Q, Q_linterp)
write_csv(clinton_gage_daily, './DataFiles/hydro_data/reservoir_outflows/clinton_gage_daily.csv')


#*** *** Milford -----------------------------------------------------------------

milford_gage_raw <- readNWISuv("06857100", "00060", "2013-10-01", "2022-09-30", tz = 'America/Chicago')

milford_gage_15min <- left_join(x = study_dateTimes, y = milford_gage_raw, by = 'dateTime') %>%
  mutate(Q = X_00060_00000/35.314666212661,
         wateryear = if_else(month(dateTime) >= 10, year(dateTime) + 1, year(dateTime)),
         yday = yday(dateTime)) %>%
  mutate(Q_linterp = linterp(Q, max_allow = NULL)) %>%
  select(dateTime, wateryear, yday, Q, Q_linterp)
write_csv(milford_gage_15min, './DataFiles/hydro_data/reservoir_outflows/milford_gage_15min.csv')


milford_gage_hourly <- milford_gage_15min %>%
  group_by(dateHour = floor_date(dateTime, unit = 'hours')) %>%
  summarise(Q = mean(Q, na.rm=TRUE),
            Q_linterp = mean(Q_linterp)) %>%
  mutate(wateryear = if_else(month(dateHour) >= 10, year(dateHour) + 1, year(dateHour)),
         yday = yday(dateHour),
         Q = replace(Q, is.nan(Q), NA)) %>%
  select(dateHour, wateryear, yday, Q, Q_linterp)
write_csv(milford_gage_hourly, './DataFiles/hydro_data/reservoir_outflows/milford_gage_hourly.csv')


milford_gage_daily <- milford_gage_15min %>%
  group_by(date = date(dateTime)) %>%
  summarise(Q = mean(Q, na.rm=TRUE),
            Q_linterp = mean(Q_linterp)) %>%
  mutate(wateryear = if_else(month(date) >= 10, year(date) + 1, year(date)),
         yday = yday(date),
         Q = replace(Q, is.nan(Q), NA)) %>%
  select(date, wateryear, yday, Q, Q_linterp)
write_csv(milford_gage_daily, './DataFiles/hydro_data/reservoir_outflows/milford_gage_daily.csv')


#*** *** Perry -----------------------------------------------------------------

perry_gage_raw <- readNWISuv("06890900", "00060", "2013-10-01", "2022-09-30", tz = 'America/Chicago')

perry_gage_15min <- left_join(x = study_dateTimes, y = perry_gage_raw, by = 'dateTime') %>%
  mutate(Q = X_00060_00000/35.314666212661,
         wateryear = if_else(month(dateTime) >= 10, year(dateTime) + 1, year(dateTime)),
         yday = yday(dateTime)) %>%
  mutate(Q_linterp = linterp(Q, max_allow = NULL)) %>%
  select(dateTime, wateryear, yday, Q, Q_linterp)
write_csv(perry_gage_15min, './DataFiles/hydro_data/reservoir_outflows/perry_gage_15min.csv')


perry_gage_hourly <- perry_gage_15min %>%
  group_by(dateHour = floor_date(dateTime, unit = 'hours')) %>%
  summarise(Q = mean(Q, na.rm=TRUE),
            Q_linterp = mean(Q_linterp)) %>%
  mutate(wateryear = if_else(month(dateHour) >= 10, year(dateHour) + 1, year(dateHour)),
         yday = yday(dateHour),
         Q = replace(Q, is.nan(Q), NA)) %>%
  select(dateHour, wateryear, yday, Q, Q_linterp)
write_csv(perry_gage_hourly, './DataFiles/hydro_data/reservoir_outflows/perry_gage_hourly.csv')


perry_gage_daily <- perry_gage_15min %>%
  group_by(date = date(dateTime)) %>%
  summarise(Q = mean(Q, na.rm=TRUE),
            Q_linterp = mean(Q_linterp)) %>%
  mutate(wateryear = if_else(month(date) >= 10, year(date) + 1, year(date)),
         yday = yday(date),
         Q = replace(Q, is.nan(Q), NA)) %>%
  select(date, wateryear, yday, Q, Q_linterp)
write_csv(perry_gage_daily, './DataFiles/hydro_data/reservoir_outflows/perry_gage_daily.csv')


#*** *** Tuttle -----------------------------------------------------------------

tuttle_gage_raw <- readNWISuv("06887000", "00060", "2013-10-01", "2022-09-30", tz = 'America/Chicago')

tuttle_gage_15min <- left_join(x = study_dateTimes, y = tuttle_gage_raw, by = 'dateTime') %>%
  mutate(Q = X_00060_00000/35.314666212661,
         wateryear = if_else(month(dateTime) >= 10, year(dateTime) + 1, year(dateTime)),
         yday = yday(dateTime)) %>%
  mutate(Q_linterp = linterp(Q, max_allow = NULL)) %>%
  select(dateTime, wateryear, yday, Q, Q_linterp)
write_csv(tuttle_gage_15min, './DataFiles/hydro_data/reservoir_outflows/tuttle_gage_15min.csv')


tuttle_gage_hourly <- tuttle_gage_15min %>%
  group_by(dateHour = floor_date(dateTime, unit = 'hours')) %>%
  summarise(Q = mean(Q, na.rm=TRUE),
            Q_linterp = mean(Q_linterp)) %>%
  mutate(wateryear = if_else(month(dateHour) >= 10, year(dateHour) + 1, year(dateHour)),
         yday = yday(dateHour),
         Q = replace(Q, is.nan(Q), NA)) %>%
  select(dateHour, wateryear, yday, Q, Q_linterp)
write_csv(tuttle_gage_hourly, './DataFiles/hydro_data/reservoir_outflows/tuttle_gage_hourly.csv')


tuttle_gage_daily <- tuttle_gage_15min %>%
  group_by(date = date(dateTime)) %>%
  summarise(Q = mean(Q, na.rm=TRUE),
            Q_linterp = mean(Q_linterp)) %>%
  mutate(wateryear = if_else(month(date) >= 10, year(date) + 1, year(date)),
         yday = yday(date),
         Q = replace(Q, is.nan(Q), NA)) %>%
  select(date, wateryear, yday, Q, Q_linterp)
write_csv(tuttle_gage_daily, './DataFiles/hydro_data/reservoir_outflows/tuttle_gage_daily.csv')


#*** *** Kanopolis ---------------------------------------------------------------

kanopolis_gage_raw <- readNWISuv("06865500", "00060", "2013-10-01", "2022-09-30", tz = 'America/Chicago')

kanopolis_gage_15min <- left_join(x = study_dateTimes, y = kanopolis_gage_raw, by = 'dateTime') %>%
  mutate(Q = X_00060_00000/35.314666212661,
         wateryear = if_else(month(dateTime) >= 10, year(dateTime) + 1, year(dateTime)),
         yday = yday(dateTime)) %>%
  mutate(Q_linterp = linterp(Q, max_allow = NULL)) %>%
  select(dateTime, wateryear, yday, Q, Q_linterp)
write_csv(kanopolis_gage_15min, './DataFiles/hydro_data/reservoir_outflows/kanopolis_gage_15min.csv')


kanopolis_gage_hourly <- kanopolis_gage_15min %>%
  group_by(dateHour = floor_date(dateTime, unit = 'hours')) %>%
  summarise(Q = mean(Q, na.rm=TRUE),
            Q_linterp = mean(Q_linterp)) %>%
  mutate(wateryear = if_else(month(dateHour) >= 10, year(dateHour) + 1, year(dateHour)),
         yday = yday(dateHour),
         Q = replace(Q, is.nan(Q), NA)) %>%
  select(dateHour, wateryear, yday, Q, Q_linterp)
write_csv(kanopolis_gage_hourly, './DataFiles/hydro_data/reservoir_outflows/kanopolis_gage_hourly.csv')


kanopolis_gage_daily <- kanopolis_gage_15min %>%
  group_by(date = date(dateTime)) %>%
  summarise(Q = mean(Q, na.rm=TRUE),
            Q_linterp = mean(Q_linterp)) %>%
  mutate(wateryear = if_else(month(date) >= 10, year(date) + 1, year(date)),
         yday = yday(date),
         Q = replace(Q, is.nan(Q), NA)) %>%
  select(date, wateryear, yday, Q, Q_linterp)
write_csv(kanopolis_gage_daily, './DataFiles/hydro_data/reservoir_outflows/kanopolis_gage_daily.csv')


#*** *** Wilson ---------------------------------------------------------------

wilson_gage_raw <- readNWISuv("06868200", "00060", "2013-10-01", "2022-09-30", tz = 'America/Chicago')

wilson_gage_15min <- left_join(x = study_dateTimes, y = wilson_gage_raw, by = 'dateTime') %>%
  mutate(Q = X_00060_00000/35.314666212661,
         wateryear = if_else(month(dateTime) >= 10, year(dateTime) + 1, year(dateTime)),
         yday = yday(dateTime)) %>%
  mutate(Q_linterp = linterp(Q, max_allow = NULL)) %>%
  select(dateTime, wateryear, yday, Q, Q_linterp)
write_csv(wilson_gage_15min, './DataFiles/hydro_data/reservoir_outflows/wilson_gage_15min.csv')


wilson_gage_hourly <- wilson_gage_15min %>%
  group_by(dateHour = floor_date(dateTime, unit = 'hours')) %>%
  summarise(Q = mean(Q, na.rm=TRUE),
            Q_linterp = mean(Q_linterp)) %>%
  mutate(wateryear = if_else(month(dateHour) >= 10, year(dateHour) + 1, year(dateHour)),
         yday = yday(dateHour),
         Q = replace(Q, is.nan(Q), NA)) %>%
  select(dateHour, wateryear, yday, Q, Q_linterp)
write_csv(wilson_gage_hourly, './DataFiles/hydro_data/reservoir_outflows/wilson_gage_hourly.csv')


wilson_gage_daily <- wilson_gage_15min %>%
  group_by(date = date(dateTime)) %>%
  summarise(Q = mean(Q, na.rm=TRUE),
            Q_linterp = mean(Q_linterp)) %>%
  mutate(wateryear = if_else(month(date) >= 10, year(date) + 1, year(date)),
         yday = yday(date),
         Q = replace(Q, is.nan(Q), NA)) %>%
  select(date, wateryear, yday, Q, Q_linterp)
write_csv(wilson_gage_daily, './DataFiles/hydro_data/reservoir_outflows/wilson_gage_daily.csv')


#*** *** Waconda ---------------------------------------------------------------

waconda_gage_raw <- readNWISuv("06875900", "00060", "2013-10-01", "2022-09-30", tz = 'America/Chicago')

waconda_gage_15min <- left_join(x = study_dateTimes, y = waconda_gage_raw, by = 'dateTime') %>%
  mutate(Q = X_00060_00000/35.314666212661,
         wateryear = if_else(month(dateTime) >= 10, year(dateTime) + 1, year(dateTime)),
         yday = yday(dateTime)) %>%
  mutate(Q_linterp = linterp(Q, max_allow = NULL)) %>%
  select(dateTime, wateryear, yday, Q, Q_linterp)
write_csv(waconda_gage_15min, './DataFiles/hydro_data/reservoir_outflows/waconda_gage_15min.csv')


waconda_gage_hourly <- waconda_gage_15min %>%
  group_by(dateHour = floor_date(dateTime, unit = 'hours')) %>%
  summarise(Q = mean(Q, na.rm=TRUE),
            Q_linterp = mean(Q_linterp)) %>%
  mutate(wateryear = if_else(month(dateHour) >= 10, year(dateHour) + 1, year(dateHour)),
         yday = yday(dateHour),
         Q = replace(Q, is.nan(Q), NA)) %>%
  select(dateHour, wateryear, yday, Q, Q_linterp)
write_csv(waconda_gage_hourly, './DataFiles/hydro_data/reservoir_outflows/waconda_gage_hourly.csv')


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
lawrence_gage_15min <- read_csv('./DataFiles/hydro_data/lawrence_gage/lawrence_gage_15min.csv')
clinton_gage_15min <- read_csv('./DataFiles/hydro_data/reservoir_outflows/clinton_gage_15min.csv')
milford_gage_15min <- read_csv('./DataFiles/hydro_data/reservoir_outflows/milford_gage_15min.csv')
perry_gage_15min <- read_csv('./DataFiles/hydro_data/reservoir_outflows/perry_gage_15min.csv')
tuttle_gage_15min <- read_csv('./DataFiles/hydro_data/reservoir_outflows/tuttle_gage_15min.csv')
kanopolis_gage_15min <- read_csv('./DataFiles/hydro_data/reservoir_outflows/kanopolis_gage_15min.csv')
wilson_gage_15min <- read_csv('./DataFiles/hydro_data/reservoir_outflows/wilson_gage_15min.csv')
waconda_gage_15min <- read_csv('./DataFiles/hydro_data/reservoir_outflows/waconda_gage_15min.csv')

belvue_precip_15min <- read_csv('./DataFiles/hydro_data/precip/belvue_precip_15min.csv')

DeSoto_daily <- read_csv('./DataFiles/hydro_data/DeSoto_gage/DeSoto_daily.csv')
clinton_gage_daily <- read_csv('./DataFiles/hydro_data/reservoir_outflows/clinton_gage_daily.csv')
milford_gage_daily <- read_csv('./DataFiles/hydro_data/reservoir_outflows/milford_gage_daily.csv')
perry_gage_daily <- read_csv('./DataFiles/hydro_data/reservoir_outflows/perry_gage_daily.csv')
tuttle_gage_daily <- read_csv('./DataFiles/hydro_data/reservoir_outflows/tuttle_gage_daily.csv')
kanopolis_gage_daily <- read_csv('./DataFiles/hydro_data/reservoir_outflows/kanopolis_gage_daily.csv')
wilson_gage_daily <- read_csv('./DataFiles/hydro_data/reservoir_outflows/wilson_gage_daily.csv')
waconda_gage_daily <- read_csv('./DataFiles/hydro_data/reservoir_outflows/waconda_gage_daily.csv')

lawrence_precip_daily <- read_csv('./DataFiles/hydro_data/precip/lawrence_precip_daily.csv')
topeka_precip_daily <- read_csv('./DataFiles/hydro_data/precip/topeka_precip_daily.csv')
manhattan_precip_daily <- read_csv('./DataFiles/hydro_data/precip/manhattan_precip_daily.csv')
salina_precip_daily <- read_csv('./DataFiles/hydro_data/precip/salina_precip_daily.csv')
belvue_precip_daily <- read_csv('./DataFiles/hydro_data/precip/belvue_precip_daily.csv')

all_hydro_15min <- DeSoto_15min %>%
  rename(riverQ = Q, riverQ_linterp = Q_linterp) %>%
  left_join(lawrence_gage_15min) %>% rename(lawrenceQ = Q, lawrenceQ_linterp = Q_linterp) %>%
  left_join(clinton_gage_15min) %>% rename(clintonQ = Q, clintonQ_linterp = Q_linterp) %>%
  left_join(milford_gage_15min) %>% rename(milfordQ = Q, milfordQ_linterp = Q_linterp) %>% 
  left_join(perry_gage_15min) %>% rename(perryQ = Q, perryQ_linterp = Q_linterp) %>% 
  left_join(tuttle_gage_15min) %>% rename(tuttleQ = Q, tuttleQ_linterp = Q_linterp) %>%
  left_join(kanopolis_gage_15min) %>% rename(kanopolisQ = Q, kanopolisQ_linterp = Q_linterp) %>%
  left_join(waconda_gage_15min) %>% rename(wacondaQ = Q, wacondaQ_linterp = Q_linterp) %>%
  left_join(wilson_gage_15min) %>% rename(wilsonQ = Q, wilsonQ_linterp = Q_linterp) %>%
  left_join(belvue_precip_15min) %>% rename(belvue_precip = precip) %>%
  mutate(reservoir_sum = select(., clintonQ, milfordQ, perryQ, tuttleQ, kanopolisQ, wacondaQ, wilsonQ) %>% rowSums(na.rm = TRUE),
         reservoir_sum_linterp = select(., clintonQ_linterp, milfordQ_linterp, perryQ_linterp, tuttleQ_linterp, kanopolisQ_linterp, wacondaQ_linterp, wilsonQ_linterp) %>% rowSums(na.rm = TRUE))
write_csv(all_hydro_15min, './DataFiles/hydro_data/all_hydro_15min.csv')

all_hydro_daily <- DeSoto_daily %>%
  dplyr::rename(riverQ = Q, riverQ_linterp = Q_linterp) %>%
  left_join(clinton_gage_daily) %>% dplyr::rename(clintonQ = Q, clintonQ_linterp = Q_linterp) %>%
  left_join(milford_gage_daily) %>% dplyr::rename(milfordQ = Q, milfordQ_linterp = Q_linterp) %>% 
  left_join(perry_gage_daily) %>% dplyr::rename(perryQ = Q, perryQ_linterp = Q_linterp) %>% 
  left_join(tuttle_gage_daily) %>% dplyr::rename(tuttleQ = Q, tuttleQ_linterp = Q_linterp) %>%
  left_join(kanopolis_gage_daily) %>% dplyr::rename(kanopolisQ = Q, kanopolisQ_linterp = Q_linterp) %>%
  left_join(waconda_gage_daily) %>% dplyr::rename(wacondaQ = Q, wacondaQ_linterp = Q_linterp) %>%
  left_join(wilson_gage_daily) %>% dplyr::rename(wilsonQ = Q, wilsonQ_linterp = Q_linterp) %>%
  left_join(lawrence_precip_daily) %>% dplyr::rename(lawrence_precip = precip) %>%
  left_join(topeka_precip_daily) %>% dplyr::rename(topeka_precip = precip) %>%
  left_join(manhattan_precip_daily) %>% dplyr::rename(manhattan_precip = precip) %>%
  left_join(salina_precip_daily) %>% dplyr::rename(salina_precip = precip) %>%
  left_join(belvue_precip_daily) %>% dplyr::rename(belvue_precip = precip) %>%
  mutate(reservoir_sum = select(., clintonQ, milfordQ, perryQ, tuttleQ, kanopolisQ, wacondaQ, wilsonQ) %>% rowSums(na.rm = TRUE),
         reservoir_sum_linterp = select(., clintonQ_linterp, milfordQ_linterp, perryQ_linterp, tuttleQ_linterp, kanopolisQ_linterp, wacondaQ_linterp, wilsonQ_linterp) %>% rowSums(na.rm = TRUE))
write_csv(all_hydro_daily, './DataFiles/hydro_data/all_hydro_daily.csv')


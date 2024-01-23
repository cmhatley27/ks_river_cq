library(tidyverse)
library(dataRetrieval)

study_dateTimes <- data.frame(dateTime = seq.POSIXt(as.POSIXlt('2013-10-01 00:00:00'), as.POSIXlt('2022-09-30 23:45:00'), by = '15 min'))
all_hydro_15min <- with_tz(read_csv('./DataFiles/hydro_data/all_hydro_15min.csv'), 'America/Chicago')

ftriley_raw <- readNWISuv("06879100", c("00060", "99133"), "2013-10-01", "2022-09-30", tz = 'America/Chicago')

ftriley_15min <- left_join(x = study_dateTimes, y = with_tz(ftriley_raw, 'America/Chicago'), by = 'dateTime') %>%
  mutate(Q = X_00060_00000/35.314666212661,
         wateryear = if_else(month(dateTime) >= 10, year(dateTime) + 1, year(dateTime)),
         yday = yday(dateTime)) %>%
  select(dateTime, wateryear, yday, Q) %>%
  mutate(milfordQ = all_hydro_15min$milfordQ,
         DeSotoQ = all_hydro_15min$riverQ_linterp) %>%
  mutate(WestKSQ = Q - milfordQ)

ggplot(data = ftriley_15min) +
  geom_line(aes(x = dateTime, y = WestKSQ)) +
  geom_line(aes(x = dateTime, y = DeSotoQ), color = 'blue')


diff_ccf <- ccf(ftriley_15min$DeSotoQ, ftriley_15min$WestKSQ, lag.max = 672, na.action = na.pass)
(which.max(diff_ccf$acf) - 672)

ggplot(data = ftriley_15min) +
  geom_histogram(aes(x = WestKSQ/DeSotoQ)) +
  xlim(c(0,1))
  


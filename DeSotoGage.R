library(tidyverse)
library(lubridate)
library(dataRetrieval)
library(hydroEvents)
library(rnoaa)
library(ggpubr)

###Search for KS N data###

pcodes <- readNWISpCode("all")

NQcodes <- pcodes %>%
  filter(str_detect(parameter_nm, regex("nitrate", ignore_case = TRUE)) | 
           str_detect(parameter_nm, regex("discharge", ignore_case = TRUE)))

KS_N <- whatNWISsites(stateCd = "KS", parameterCd = c("99133", "99136", "99137"))

#####Download and clean river gage data#####

DeS.raw <- readNWISuv("06892350", c("00060", "99133"), "2013-10-01", "2022-09-30")

DeS.rename <- renameNWISColumns(DeS.raw) %>%
  rename(NO3 = X_99133_Inst) %>%
  mutate(Q_m3s = Flow_Inst/35.314666212661) 

DeS.clean <- DeS.rename %>%
  select(dateTime, Q_m3s, NO3)

DeS.all <- DeS.clean %>%
  mutate(logQ = log10(Q_m3s), logN = log10(NO3)) %>%
  mutate(date = date(dateTime),
         year = factor(year(dateTime)),
         month = month(dateTime),
         day = day(dateTime),
         hour = floor_date(dateTime, unit = "hour"),
         season.yearly = floor_date(dateTime, unit ="season"),
         season = factor(month(season.yearly)),
         yday = yday(dateTime))

##Subsetting by time variables
DeS_byDate <- DeS.all %>%
  group_by(date) %>%
  summarise(avgQ = mean(Q_m3s),
            avgC = mean(NO3)) %>%
  mutate(logQ = log10(avgQ),
         logC = log10(avgC),
         year = factor(year(date)),
         yday = yday(date),
         season = (month(floor_date(date, unit = 'season')))) %>%
  mutate(season = replace(season, season == 9, 'Fall'),
         season = replace(season, season == 12, 'Winter'),
         season = replace(season, season == 3, 'Spring'),
         season = replace(season, season == 6, 'Summer'),
         season = factor(season, levels = c('Winter','Spring','Summer','Fall')))
?reorder
DeS_byHour <- DeS.all %>%
  group_by(hour) %>%
  summarise(avgQ = mean(Q_m3s, na.rm=TRUE),
            avgC = mean(NO3, na.rm = TRUE)) %>%
  mutate(logQ = log10(avgQ),
         logC = log10(avgC))

DeS_byyday <- DeS.all %>%
  group_by(yday) %>%
  summarise(avgQ = mean(Q_m3s, na.rm=TRUE),
            avgC = mean(NO3, na.rm = TRUE),
            sdC = sd(NO3, na.rm = TRUE),
            minC = min(NO3, na.rm = TRUE),
            maxC = max(NO3, na.rm = TRUE)) %>%
  mutate(logQ = log10(avgQ),
         logC = log10(avgC),
         season = ifelse(yday >= 335, "Winter",
                         ifelse(yday >= 244, "Fall",
                                ifelse(yday >= 152, "Summer",
                                       ifelse(yday >= 60, "Spring", "Winter")))),
         season = factor(season, levels = c('Winter','Spring','Summer','Fall')))

DeS_bySeason <- DeS.all %>%
  group_by(season) %>%
  summarise(avgC = mean(NO3, na.rm = TRUE),
            avgQ = mean(Q_m3s, na.rm = TRUE))

## Subset byDate into over/under 2.0 logQ ##
DeS_byDate.lowQ <- DeS_byDate %>% filter(logQ < 2)
lm(logC~logQ, data = DeS_byDate.lowQ)

DeS_byDate.highQ <- DeS_byDate %>% filter(logQ >= 2)
lm(logC~logQ, data = DeS_byDate.highQ)

###Calculate windowed CQ slopes###

windowCQslopes <- function(data, c, q, winlength = 30) {
  
  datalength <- nrow(data) - (winlength - 1)
  CQ.slopes <- vector("double", length = datalength)
  
  for (i in 1:datalength) { 
    window.subset <- data[i:(i+(winlength - 1)),]
    try(window.lm <- lm(c[i:(i+(winlength - 1))]~q[i:(i+(winlength - 1))], data = data, na.action = na.omit))
    window.slope <- unname(coef(window.lm)[2])
    
    CQ.slopes[[i]] <- window.slope
  }
  CQ.slopes <- c(rep(NA, winlength/2), CQ.slopes, rep(NA, (winlength/2) - 1))
 
  return(CQ.slopes)
}

windowCQnas <- function(data, c, q, winlength = 30) {
  
  datalength <- nrow(data) - (winlength - 1)
  numNAs <- vector( "double", length = datalength)
  hasNAc <- is.na(c)
  hasNAq <- is.na(q)
  hasNA <- hasNAc | hasNAq
  
  for (i in 1:datalength) { 
    window.subset <- data[i:(i+(winlength - 1)),]
    numNAs[[i]] <- sum(hasNA[i:(i+(winlength - 1))])
  }
  numNAs <- c(rep(NA, winlength/2), numNAs, rep(NA, (winlength/2) - 1))
  
  return(numNAs)
}

CQslopes.30 <- windowCQslopes(DeS_byDate, DeS_byDate$logC, DeS_byDate$logQ)
CQnas.30 <- windowCQnas(DeS_byDate, DeS_byDate$logC, DeS_byDate$logQ)

CQnas.30
CQslopes.30

CQslopes.90 <- windowCQslopes(DeS_byDate, DeS_byDate$logC, DeS_byDate$logQ, winlength = 90)
CQnas.90 <- windowCQnas(DeS_byDate, DeS_byDate$logC, DeS_byDate$logQ, winlength = 90)

CQslopes.90                    
CQnas.90

DeS_byDate <- DeS_byDate %>%
  mutate(CQslope = CQslopes.30,
         numNAs = CQnas.30,
         CQslope.trim = ifelse(numNAs >= 15, NA, CQslope))

outliers <- DeS_byDate %>%
  filter(abs(CQslope) >= 3)

ggplot(data = DeS_byDate) +
  geom_line(aes(x = yday, y = CQslope, color = year)) +
  coord_cartesian(ylim = c(-2,2))

ggplot(data = DeS_byDate) +
  geom_line(aes(x = yday, y = CQslope.trim, color = year)) +
  coord_cartesian(ylim = c(-2,2)) +
  geom_hline(yintercept = 0.5083, linetype = 'dashed') +
  geom_hline(yintercept = 0)

ggplot(data = DeS_byDate) +
  geom_histogram(aes(x = numNAs))


####Download Precip####
precip <- meteo_pull_monitors("USW00003997", date_min = "2013-01-01", var = "PRCP")
precip <- precip %>%
  mutate(prcp = prcp/10,
         year = year(date),
         month = month(date),
         day = day(date),
         season.yearly = floor_date(date, "season"),
         season = factor(month(season.yearly)))

precip_bySeason <- precip %>%
  group_by(season.yearly) %>%
  summarise(totP = sum(prcp, na.rm = TRUE)) %>%
  mutate(season = factor(month(season.yearly))) %>%
  ungroup() %>% group_by(season) %>%
  summarise(avgP = mean(totP, na.rm = TRUE)) %>%
  mutate(yday = c(60, 152, 244, 335))

####Delineate events using 2021 data####

###Daily Q###
dailyq_2021 <- DeS.all %>%
  filter(year == 2021) %>%
  group_by(date) %>%
  summarise(daily_q = mean(Q_m3s, na.rm = TRUE))

#BF calcs
dailyq_2021$bf <- baseflowB(dailyq_2021$daily_q, alpha = 0.9)[["bf"]]
dailyq_2021$HYSEP_bf_fixed <- baseflow_HYSEP(Q = dailyq_2021$daily_q, area_mi2 = 59756, method = "fixed")
dailyq_2021$HYSEP_bf_sliding <- baseflow_HYSEP(Q = dailyq_2021$daily_q, area_mi2 = 59756, method = "sliding")
dailyq_2021$HYSEP_bf_local <- baseflow_HYSEP(Q = dailyq_2021$daily_q, area_mi2 = 59756, method = "local")

dailyq_2021$qf <- dailyq_2021$daily_q - dailyq_2021$bf

bfis <- c("digital filter" = sum(dailyq_2021$bf)/sum(dailyq_2021$daily_q),
          "HYSEP fixed" = sum(dailyq_2021$HYSEP_bf_fixed)/sum(dailyq_2021$daily_q),
          "HYSEP sliding" = sum(dailyq_2021$HYSEP_bf_sliding)/sum(dailyq_2021$daily_q),
          "HYSEP local" = sum(dailyq_2021$HYSEP_bf_local)/sum(dailyq_2021$daily_q))

ggplot(data = dailyq_2021, aes(x = date)) +
  geom_line(aes(y = daily_q)) +
  geom_line(aes(y = bf), color = "grey") +
  geom_line(aes(y = HYSEP_bf_fixed), color = "red") +
  geom_line(aes(y = HYSEP_bf_sliding), color = "green") +
  geom_line(aes(y = HYSEP_bf_local), color = "blue")

ggplot(data = dailyq_2021, aes(x = date)) +
  geom_line(aes(y = qf)) +
  geom_hline(yintercept = 50)

ggplot(data = dailyq_2021, aes(x = date)) +
  geom_line(aes(y = bf/daily_q))

q_events.min <- eventMinima(dailyq_2021$qf, delta.y = 40, delta.x = 4, threshold = 20)
plotEvents(data = dailyq_2021$qf, dates = dailyq_2021$date, events = q_events.min, 
           type = "bound", colline = "black", colbound = "grey")

q_events.max <- eventMaxima(dailyq_2021$qf, delta.y = 50, delta.x = 1, threshold = 50)
plotEvents(data = dailyq_2021$daily_q, dates = dailyq_2021$date, events = q_events.max, 
           type = "bound", colline = "black", colbound = "grey")

q_events.bf <- eventBaseflow(dailyq_2021$daily_q, BFI_Th = 0.75)
plotEvents(data = dailyq_2021$daily_q, dates = dailyq_2021$date, events = q_events.bf, 
           type = "bound", colline = "black", colbound = "grey")

q_events.pot <- eventPOT(dailyq_2021$qf, threshold = 10, min.diff = 7)
plotEvents(data = dailyq_2021$daily_q, dates = dailyq_2021$date, events = q_events.pot, 
           type = "bound", colline = "black", colbound = "grey")

###Hourly Q###
hourlyq_2021 <- DeS %>%
  filter(year == 2021) %>%
  group_by(hour) %>%
  summarise(hourly_q = mean(Q_m3s, na.rm = TRUE)) %>%
  filter(!is.na(hourly_q)) %>%
  mutate(bf = baseflowB(hourly_q, alpha = 0.925, passes = 12)[["bf"]]) %>%
  mutate(qf = hourly_q - bf)

ggplot(data = hourlyq_2021, aes(x = hour)) +
  geom_line(aes(y = hourly_q)) +
  geom_line(aes(y = bf), color = "grey")

q_events.min <- eventMinima(hourlyq_2021$qf, delta.y = 100, delta.x = 24*7, threshold = 50)
plotEvents(data = hourlyq_2021$hourly_q, dates = hourlyq_2021$hour, events = q_events.min, 
           type = "bound", colline = "black", colbound = "grey")

q_events.max <- eventMaxima(hourlyq_2021$qf, delta.y = 100, delta.x = 24, threshold = 50)
plotEvents(data = hourlyq_2021$hourly_q, dates = hourlyq_2021$hour, events = q_events.max, 
           type = "bound", colline = "black", colbound = "grey")

q_events.bf <- eventBaseflow(hourlyq_2021$hourly_q, BFI_Th = 0.9, min.diff = 24*7)
plotEvents(data = hourlyq_2021$hourly_q, dates = hourlyq_2021$hour, events = q_events.bf, 
           type = "bound", colline = "black", colbound = "grey")

q_events.pot <- eventPOT(hourlyq_2021$qf, threshold = 0, min.diff = 24)
plotEvents(data = hourlyq_2021$hourly_q, dates = hourlyq_2021$hour, events = q_events.pot, 
           type = "bound", colline = "black", colbound = "grey")


dailyp_2021 <- precip %>%
  filter(year == 2021)

ggplot(data = dailyp_2021, aes(x = date, y = prcp)) +
  geom_line()



#####Misc plots
#by Date CQ plots
ggplot(data = DeS_byDate) +
  geom_point(aes(x = avgQ, y = avgC))

### log - log plot of all dates ###
ggplot(data = DeS_byDate) +
  geom_point(aes(x = logQ, y = logC, color = season)) +
  geom_smooth(aes(x = logQ, y = logC), method = 'lm', color = 'black', se = FALSE, linetype = 'dashed') +
  annotate('text', x = 2.8, y = -0.25, label = "Slope = 0.538", color = 'black') +
  annotate('text', x = 2.8, y = -0.35, label = "CVc/CVq = 0.472", color = 'black') +
  scale_color_manual(values = c('#8DA9C4','#60992D','#C33C54','#F4B860'))

 
logQdist <- ggplot(data = DeS_byDate) +
  geom_violin(aes(x = logQ, y = ""), orientation = 'y') +
  coord_cartesian(xlim = c(1, 3.5)) 
ggarrange(logQdist, logCQall, ncol = 1, align = 'v', heights = c(1,2.5))
?ggarrange
?geom_violin
CVratio <- function(c, q) {
  cvc <- sd(c, na.rm = TRUE)/mean(c, na.rm = TRUE)
  cvq <- sd(q, na.rm = TRUE)/mean(q, na.rm = TRUE)
  cvc/cvq
}
CVratio(DeS_byDate$avgC, DeS_byDate$avgQ)
lm(logC~logQ, data = DeS_byDate)


?geom_smooth
#by Hour CQ plots
ggplot(data = DeS_byHour) +
  geom_point(aes(x = avgQ, y = avgC))
ggplot(data = DeS_byHour) +
  geom_point(aes(x = logQ, y = logC)) +
  geom_smooth(aes(x = logQ, y = logC))

#by yday CQ plots
ggplot(data = DeS_byyday) +
  geom_point(aes(x = avgQ, y = avgC, color = season)) +
  geom_vline(xintercept = 312.5, color = 'blue') +
  scale_color_manual(values = c('#8DA9C4','#60992D','#C33C54','#F4B860'))


ggplot(data = DeS_byyday) +
  geom_point(aes(x = logQ, y = logC, color = season)) +
  geom_vline(xintercept = log10(312.5), color = 'blue') +
  scale_color_manual(values = c('#8DA9C4','#60992D','#C33C54','#F4B860'))


### N by yday time series ###
ggplot(data = DeS_byyday) +
  geom_line(aes(x = yday, y = avgC)) +
  geom_line(aes(x = yday, y = avgQ/1000), color = 'blue') +
  scale_y_continuous(sec.axis = sec_axis(~ . *1000))
?sec_axis
####Time series N
#full observation period
ggplot(data = DeS.all) +
  geom_line(aes(x = dateTime, y = NO3))
#grouped by year
ggplot(data = DeS_byDate) +
  geom_line(aes(x = date, y = avgC, color = year), alpha = 1) +
  #coord_cartesian(xlim = date(c("2016-10-01","2016-12-31"))) +
  geom_line(aes(x = date, y = logQ))
#by yday
ggplot(data = DeS_byyday) +
  geom_line(aes(x = yday, y = (avgC), color = "C")) +
  geom_line(aes(x = yday, y = (maxC), color = "minmax", alpha = 0.3)) +
  geom_line(aes(x = yday, y = (minC), color = "minmax", alpha = 0.3)) +
  geom_line(aes(x = yday, y = log10(avgQ), color = "logQ")) +
  scale_color_manual(values = c("minmax" = "red", "logQ" = "blue", "C" = "black")) +
  scale_alpha_identity()

ggplot(data = DeS_byyday) +
  geom_smooth(aes(x = yday, y = avgC)) +
  geom_smooth(aes(x = yday, y = avgC + sdC, color = "red")) +
  geom_smooth(aes(x = yday, y = avgC - sdC, color = "red"))



###Seasonal Data###
seasons <- precip_bySeason %>%
  left_join(DeS_bySeason, by = 'season') %>%
  mutate(season = factor(c("Spring", "Summer", "Fall", "Winter"), 
                         levels = c("Winter", "Spring", "Summer", "Fall")))

ggplot(data = seasons) +
  geom_col(aes(x = season, y = avgP, fill = avgC)) +
  scale_fill_gradient(low = 'darkgrey', high = '#C33C54', limits = c(0.75, max(seasons$avgC)))

ggplot(data = seasons) +
  geom_col(aes(x = season, y = avgQ, fill = avgC)) +
  scale_fill_gradient(low = 'darkgrey', high = '#C33C54', limits = c(0.75, max(seasons$avgC)))


median(DeS_byDate$avgC, na.rm = TRUE) - IQR(DeS_byDate$avgC, na.rm = TRUE)/2
?geom_col
ggplot(data = DeS) +
  geom_point(aes(x = logQ, y = logC))

ggplot(data = DeS) +
  geom_point(aes(x = Q_m3s, y = NO3NO2, alpha = 0.05)) +
  geom_smooth(aes(x = Q_m3s, y = NO3NO2))

ggplot(data = DeS) +
  geom_histogram(aes(x = NO3NO2), bins = 200)

ggplot(data = DeS) +
  geom_histogram(aes(x = logC), bins = 50)

ggplot(data = DeS) +
  geom_histogram(aes(x = Q_m3s), bins = 50)

ggplot(data = DeS) +
  geom_histogram(aes(x = logQ))




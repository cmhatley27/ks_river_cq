library(tidyverse)
library(lubridate)
library(dataRetrieval)
library(hydroEvents)
library(rnoaa)
library(ggpubr)
library(FluMoDL)


#####Load in gage and precip data#####
DeS.all <- read_csv('./DataFiles/DeS_all.csv')
precip <- read_csv('./DataFiles/precip.csv')


#####Subset river gage data into different time groupings####
DeS_byDate <- DeS.all %>%
  group_by(date) %>%
  summarise(totQ = sum(Q_m3s*60*15),
            avgQ = mean(Q_m3s, na.rm=TRUE),
            avgC = mean(NO3),
            loadC = avgC/1000*totQ) %>%
  mutate(logQ = log10(avgQ),
         logC = log10(avgC),
         year = year(date),
         wateryear = if_else(month(date) >= 10, year + 1, year),
         yday = yday(date),
         season = (month(floor_date(date, unit = 'season')))) %>%
  mutate(season = replace(season, season == 9, 'Fall'),
         season = replace(season, season == 12, 'Winter'),
         season = replace(season, season == 3, 'Spring'),
         season = replace(season, season == 6, 'Summer'),
         season = factor(season, levels = c('Winter','Spring','Summer','Fall')))
write_csv(DeS_byDate, './DataFiles/DeS_byDate.csv')


?floor_date
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
            avgQ = mean(Q_m3s, na.rm = TRUE)) %>%
  mutate(season = factor(c("Spring", "Summer", "Fall", "Winter"), 
                         levels = c("Winter", "Spring", "Summer", "Fall")))


YearlyCQSlope <- DeS_byDate %>%
  group_by(wateryear) %>%
  summarise(slope = coef(lm(logC~logQ,na.action = na.omit))[2])

DeS_byYear <- DeS.all %>%
  group_by(wateryear) %>%
  summarise(avgC = mean(NO3, na.rm = TRUE),
            avgQ = mean(Q_m3s, na.rm = TRUE),
            logC = log10(avgC),
            logQ = log10(avgQ),
            totQ = sum(Q_m3s*60*15, na.rm = TRUE),
            medianQ = median(Q_m3s, na.rm = TRUE),
            iqrQ = IQR(Q_m3s, na.rm = TRUE)) %>%
  filter(wateryear < 2023) %>%
  left_join(YearlyCQSlope)


####Calculate Flow Duration Curves####

###From Daily Averages###
FDC <- DeS_byDate %>%
  select(date, avgQ) %>%
  filter(!is.na(avgQ)) %>%
  arrange(desc(avgQ)) %>%
  mutate(rank = seq(1:length(avgQ)),
         P = rank/(length(avgQ)+1)*100,
         percentile = 100 - P)

DeS_byDate <- DeS_byDate %>%
  left_join(FDC) %>%
  mutate(percentile.bin = cut(percentile, breaks = c(0, 10, 40, 60, 90, 100), labels = c('0 - 10', '10 - 40', '40 - 60', '60 - 90', '90 - 100')))

###Calculate % of annual loads by flow percentile###
loadbypercentile <- DeS_byDate %>%
  group_by(wateryear, percentile.bin) %>%
  summarise(load.percentile = sum(loadC, na.rm = TRUE)) %>%
  filter(!is.na(percentile.bin)) %>%
  group_by(wateryear) %>%
  mutate(load.frac = load.percentile/sum(load.percentile))
  

####Calculate windowed CQ slopes####

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

windowCQcvratio <- function(data, c, q, winlength = 30) {
  
  datalength <- nrow(data) - (winlength - 1)
  CVc <- vector("double", length = datalength)
  CVq <- vector("double", length = datalength)
  
  for (i in 1:datalength) { 
    window.subset <- data[i:(i+(winlength - 1)),]
    CVc[[i]] <- sd(c[i:(i+(winlength - 1))], na.rm = TRUE)/mean(c[i:(i+(winlength - 1))], na.rm = TRUE)
    CVq[[i]] <- sd(q[i:(i+(winlength - 1))], na.rm = TRUE)/mean(q[i:(i+(winlength - 1))], na.rm = TRUE)
  }
  CVc <- c(rep(NA, winlength/2), CVc, rep(NA, (winlength/2) - 1))
  CVq <- c(rep(NA, winlength/2), CVq, rep(NA, (winlength/2) - 1))
  CVratio <- CVc/CVq
  
  return(CVratio)
}

windowCQpearson <- function(data, c, q, winlength = 30) {
  
  datalength <- nrow(data) - (winlength - 1)
  pearson <- vector("double", length = datalength)
  
  for (i in 1:datalength) { 
    window.subset <- data[i:(i+(winlength - 1)),]
    pearson[[i]] <- as.double(try(unname(cor.test(c[i:(i+(winlength - 1))], q[i:(i+(winlength - 1))],
                                                  method = 'pearson', alternative = c('two.sided'))[["estimate"]])))
  }
  pearson <- c(rep(NA, winlength/2), pearson, rep(NA, (winlength/2) - 1))
  
  return(pearson)
}

CQslopes.30 <- windowCQslopes(DeS_byDate, DeS_byDate$logC, DeS_byDate$logQ)
CQnas.30 <- windowCQnas(DeS_byDate, DeS_byDate$logC, DeS_byDate$logQ)
CQcvratio.30 <- windowCQcvratio(DeS_byDate, DeS_byDate$avgC, DeS_byDate$avgQ)
CQpearson.30 <- windowCQpearson(DeS_byDate, DeS_byDate$avgC, DeS_byDate$avgQ)

DeS_byDate <- DeS_byDate %>%
  mutate(CQslope = CQslopes.30,
         numNAs = CQnas.30,
         cvratio = CQcvratio.30,
         pearson = CQpearson.30,
         significance = ifelse(is.na(cvratio), "NA", ifelse(cvratio >= 1, "Sig", "NotSig")),
         CQslope.trim = ifelse(numNAs >= 15, NA, CQslope),
         cvratio.trim = ifelse(numNAs >= 15, NA, cvratio),
         pearson.trim = ifelse(numNAs >= 15, NA, pearson),
         significance.trim = ifelse(is.na(cvratio.trim), 0, ifelse(cvratio.trim >= 1, 2, 1)))

DeS_byYear <- DeS_byDate %>%
  group_by(year) %>%
  summarise(meanC = mean(avgC, na.rm = TRUE),
            meanQ = mean(avgQ, na.rm = TRUE),
            meanlogC = mean(logC, na.rm = TRUE),
            meanlogQ = mean(logQ, na.rm = TRUE),
            medianQ = median(avgQ, na.rm = TRUE),
            medianlogQ = median(logQ, na.rm = TRUE),
            iqrQ = IQR(avgQ, na.rm = TRUE),
            meanCQslope.trim = mean(CQslope.trim, na.rm = TRUE),
            medianCQslope.trim = median(CQslope.trim, na.rm = TRUE)) %>%
  mutate(year = as.factor(year))

DeS_bySeason <- DeS_byDate %>%
  group_by(season) %>%
  summarise(meanC = mean(avgC, na.rm = TRUE),
            meanQ = mean(avgQ, na.rm = TRUE),
            meanlogC = mean(logC, na.rm = TRUE),
            meanlogQ = mean(logQ, na.rm = TRUE),
            medianC = median(avgC, na.rm = TRUE),
            medianlogC = median(logC, na.rm = TRUE),
            medianQ = median(avgQ, na.rm = TRUE),
            medianlogQ = median(logQ, na.rm = TRUE),
            iqrQ = IQR(avgQ, na.rm = TRUE),
            meanCQslope.trim = mean(CQslope.trim, na.rm = TRUE),
            medianCQslope.trim = median(CQslope.trim, na.rm = TRUE))

outliers <- DeS_byDate %>%
  filter(abs(CQslope) >= 3)

#Count number of observations in dataset that have CVratio > 1, < 1, and NA
sig_cvs <- DeS_byDate %>%
  group_by(significance) %>%
  summarise(n = n())
sig_cvs.trim <- DeS_byDate %>%
  group_by(significance.trim) %>%
  summarise(n = n())




#####Subset and calculate for precip#####

precip_bySeason <- precip %>%
  group_by(season.yearly) %>%
  summarise(totP = sum(prcp, na.rm = TRUE)) %>%
  mutate(season = factor(month(season.yearly))) %>%
  ungroup() %>% group_by(season) %>%
  summarise(avgP = mean(totP, na.rm = TRUE)) %>%
  mutate(yday = c(60, 152, 244, 335)) %>%
  mutate(season = factor(c("Spring", "Summer", "Fall", "Winter"), 
                         levels = c("Winter", "Spring", "Summer", "Fall")))


####Create Seasonal Summary####
seasons <- precip_bySeason %>%
  left_join(DeS_bySeason, by = 'season') %>%
  mutate(season = factor(c("Spring", "Summer", "Fall", "Winter"), 
                         levels = c("Winter", "Spring", "Summer", "Fall")))

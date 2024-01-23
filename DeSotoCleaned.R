library(tidyverse)
library(lubridate)
library(dataRetrieval)
library(hydroEvents)
library(rnoaa)
library(ggpubr)


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
write_csv(DeS.all, 'DeS_all.csv')


#####Subset into different time groupings####

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
            avgQ = mean(Q_m3s, na.rm = TRUE)) %>%
  mutate(season = factor(c("Spring", "Summer", "Fall", "Winter"), 
                         levels = c("Winter", "Spring", "Summer", "Fall")))

DeS_byYear <- DeS.all %>%
  group_by(year) %>%
  summarise(avgC = mean(NO3, na.rm = TRUE),
            avgQ = mean(Q_m3s, na.rm = TRUE),
            logC = log10(avgC),
            logQ = log10(avgQ),
            medianQ = median(Q_m3s, na.rm = TRUE),
            iqrQ = IQR(Q_m3s, na.rm = TRUE))

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

####Download and clean Precip and subset into groups####

precip <- meteo_pull_monitors("USW00003997", date_min = "2013-01-01", var = c("PRCP", "TAVG"))
precip <- precip %>%
  mutate(prcp = prcp/10,
         year = year(date),
         month = month(date),
         day = day(date),
         season.yearly = floor_date(date, "season"),
         season = factor(month(season.yearly)))
write_csv(precip, './DataFiles/precip.csv')

precip_bySeason <- precip %>%
  group_by(season.yearly) %>%
  summarise(totP = sum(prcp, na.rm = TRUE)) %>%
  mutate(season = factor(month(season.yearly))) %>%
  ungroup() %>% group_by(season) %>%
  summarise(avgP = mean(totP, na.rm = TRUE)) %>%
  mutate(yday = c(60, 152, 244, 335)) %>%
  mutate(season = factor(c("Spring", "Summer", "Fall", "Winter"), 
                         levels = c("Winter", "Spring", "Summer", "Fall")))

seasons <- precip_bySeason %>%
  left_join(DeS_bySeason, by = 'season') %>%
  mutate(season = factor(c("Spring", "Summer", "Fall", "Winter"), 
                         levels = c("Winter", "Spring", "Summer", "Fall")))


############
# Plotting #
############

#### LogC v LogQ####


### Log - Log by Date ###
ggplot(data = DeS_byDate) +
  geom_point(aes(x = logQ, y = logC, color = season)) +
  geom_smooth(aes(x = logQ, y = logC), method = 'lm', color = 'black', se = FALSE, linetype = 'dashed') +
  annotate('text', x = 3.1, y = -0.2, label = "Slope = 0.54\nCVc/CVq = 0.47", color = 'black') +
  scale_color_manual(values = c('#8DA9C4','#60992D','#C33C54','#F4B860'), name = "Season") +
  coord_cartesian(xlim = c(min(DeS_byDate$logQ, na.rm = TRUE),max(DeS_byDate$logQ, na.rm = TRUE))) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  xlab("Log Discharge (m^3/s)") +
  ylab("Log Concentration of Nitrate (mg/L)")
ggsave("cqloglogscatter.tiff", device = "tiff", height = 6, width = 8, units = "in", compression = "lzw", dpi = 700)

#Inset
ggplot(data = DeS_byDate) +
  geom_smooth(aes(x = logQ, y = logC), method = 'lm', color = 'black', se = FALSE, linetype = 'dashed') +
  geom_smooth(data = DeS_byDate[DeS_byDate$season == 'Winter',], aes(x = logQ, y = logC), method = 'lm', color = '#8DA9C4', se = FALSE, linetype = 'solid', fullrange = TRUE, alpha = 0.5) +
  geom_smooth(data = DeS_byDate[DeS_byDate$season == 'Spring',], aes(x = logQ, y = logC), method = 'lm', color = '#60992D', se = FALSE, linetype = 'solid', fullrange = TRUE, alpha = 0.5) +
  geom_smooth(data = DeS_byDate[DeS_byDate$season == 'Summer',], aes(x = logQ, y = logC), method = 'lm', color = '#C33C54', se = FALSE, linetype = 'solid', fullrange = TRUE, alpha = 0.5) +
  geom_smooth(data = DeS_byDate[DeS_byDate$season == 'Fall',], aes(x = logQ, y = logC), method = 'lm', color = '#F4B860', se = FALSE, linetype = 'solid', fullrange = TRUE, alpha = 0.5) +
  scale_color_manual(values = c('#8DA9C4','#60992D','#C33C54','#F4B860'), name = "Season") +
  annotate('text', x = 1.9, y = 0.5, label = "Seasonal\nRegressions", color = 'black', size = 3) +
  coord_cartesian(xlim = c(min(DeS_byDate$logQ, na.rm = TRUE),max(DeS_byDate$logQ, na.rm = TRUE))) +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) 
ggsave("cqloglogscatter_inset.tiff", device = "tiff", height = 2, width = 2, units = "in", compression = "lzw", dpi = 700)

### Log - Log by yday ###
ggplot(data = DeS_byyday) +
  geom_point(aes(x = logQ, y = logC, color = season)) +
  geom_vline(xintercept = log10(312.5), color = 'blue') +
  scale_color_manual(values = c('#8DA9C4','#60992D','#C33C54','#F4B860'))


#### 30 Day Windows ####

### Time Series colored by year ###
#All points
ggplot(data = DeS_byDate) +
  geom_line(aes(x = yday, y = CQslope, color = year)) +
  coord_cartesian(ylim = c(-2,2))

#Trimmed of windows with > 15 NAs
ggplot(data = DeS_byDate) +
  geom_line(aes(x = yday, y = CQslope.trim, color = year), size = 0.75) +
  coord_cartesian(ylim = c(-2,2)) +
  geom_hline(yintercept = 0.5083, linetype = 'dashed') +
  geom_hline(yintercept = 0)


### Slope v CVratio ###
#All points
ggplot(data = DeS_byDate) +
  geom_point(aes(x = cvratio, y = CQslope, color = year))

#Trimmed of windows with > 15 NAs
ggplot(data = DeS_byDate) +
  geom_point(aes(x = cvratio.trim, y = CQslope.trim, color = year), size = 1) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 1)

#Trimmed, colored by season instead of year
ggplot(data = DeS_byDate) +
  geom_rect(aes(xmin = 0, xmax = 1, ymin = -10, ymax = 10), fill = 'lightgrey', alpha = 0.05) +
  geom_point(aes(x = cvratio.trim, y = CQslope.trim, color = season), size = 1, alpha = 1) +
  scale_color_manual(values = c('#8DA9C4','#60992D','#C33C54','#F4B860'), name = "Season") +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 1, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'solid') +
  scale_x_continuous(limits = c(0, 8), expand = c(0,0)) +
  scale_y_continuous(breaks = seq(-8, 8, by = 2)) +
  coord_cartesian(ylim = c(-6,6)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  xlab("CVc/CVq") +
  ylab("30-Day CQ Slope")
ggsave("cq_slope_cvratio_combined.tiff", device = "tiff", height = 6, width = 8, units = "in", compression = "lzw", dpi = 700)

scatter.winter <- ggplot(data = DeS_byDate) +
  geom_rect(aes(xmin = 0, xmax = 1, ymin = -10, ymax = 10), fill = 'lightgrey', alpha = 0.05) +
  geom_point(data = DeS_byDate[DeS_byDate$season != 'Winter',],aes(x = cvratio.trim, y = CQslope.trim), color = 'darkgrey', size = 1, alpha = 1) +
  geom_point(data = DeS_byDate[DeS_byDate$season == 'Winter',], aes(x = cvratio.trim, y = CQslope.trim), 
             color = 'black', fill = '#8DA9C4', size = 1, alpha = 1, shape = 21) +
  #scale_color_manual(values = c('#8DA9C4','#60992D','#C33C54','#F4B860'), name = "Season") +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 1, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'solid') +
  scale_x_continuous(limits = c(0, 8), expand = c(0,0)) +
  scale_y_continuous(breaks = seq(-8, 8, by = 2)) +
  coord_cartesian(ylim = c(-6,6)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  xlab("") +
  ylab("30-Day CQ Slope") +
  annotate('text', x = 4, y = 5.5, size = 5, label = 'Winter', color = '#8DA9C4')
scatter.winter
scatter.spring <- ggplot(data = DeS_byDate) +
  geom_rect(aes(xmin = 0, xmax = 1, ymin = -10, ymax = 10), fill = 'lightgrey', alpha = 0.05) +
  geom_point(data = DeS_byDate[DeS_byDate$season != 'Spring',],aes(x = cvratio.trim, y = CQslope.trim), color = 'darkgrey', size = 1, alpha = 1) +
  geom_point(data = DeS_byDate[DeS_byDate$season == 'Spring',], aes(x = cvratio.trim, y = CQslope.trim), 
             fill = '#60992D', color = 'black', shape = 21, size = 1, alpha = 1) +
  #scale_color_manual(values = c('#8DA9C4','#60992D','#C33C54','#F4B860'), name = "Season") +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 1, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'solid') +
  scale_x_continuous(limits = c(0, 8), expand = c(0,0)) +
  scale_y_continuous(breaks = seq(-8, 8, by = 2)) +
  coord_cartesian(ylim = c(-6,6)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  xlab("") +
  ylab("")+
  annotate('text', x = 4, y = 5.5, size = 5, label = 'Spring', color = '#60992D')
scatter.summer <- ggplot(data = DeS_byDate) +
  geom_rect(aes(xmin = 0, xmax = 1, ymin = -10, ymax = 10), fill = 'lightgrey', alpha = 0.05) +
  geom_point(data = DeS_byDate[DeS_byDate$season != 'Summer',],aes(x = cvratio.trim, y = CQslope.trim), color = 'darkgrey', size = 1, alpha = 1) +
  geom_point(data = DeS_byDate[DeS_byDate$season == 'Summer',], aes(x = cvratio.trim, y = CQslope.trim), 
             fill = '#C33C54', color = 'black', shape = 21, size = 1, alpha = 1) +
  #scale_color_manual(values = c('#8DA9C4','#60992D','#C33C54','#F4B860'), name = "Season") +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 1, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'solid') +
  scale_x_continuous(limits = c(0, 8), expand = c(0,0)) +
  scale_y_continuous(breaks = seq(-8, 8, by = 2)) +
  coord_cartesian(ylim = c(-6,6)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  xlab("CVc/CVq") +
  ylab("30-Day CQ Slope")+
  annotate('text', x = 4, y = 5.5, size = 5, label = 'Summer', color = '#C33C54')
scatter.fall <- ggplot(data = DeS_byDate) +
  geom_rect(aes(xmin = 0, xmax = 1, ymin = -10, ymax = 10), fill = 'lightgrey', alpha = 0.05) +
  geom_point(data = DeS_byDate[DeS_byDate$season != 'Fall',],aes(x = cvratio.trim, y = CQslope.trim), color = 'darkgrey', size = 1, alpha = 1) +
  geom_point(data = DeS_byDate[DeS_byDate$season == 'Fall',], aes(x = cvratio.trim, y = CQslope.trim), 
             fill = '#F4B860', color = 'black', shape = 21, size = 1, alpha = 1) +
  #scale_color_manual(values = c('#8DA9C4','#60992D','#C33C54','#F4B860'), name = "Season") +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 1, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'solid') +
  scale_x_continuous(limits = c(0, 8), expand = c(0,0)) +
  scale_y_continuous(breaks = seq(-8, 8, by = 2)) +
  coord_cartesian(ylim = c(-6,6)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  xlab("CVc/CVq") +
  ylab("")+
  annotate('text', x = 4, y = 5.5, size = 5, label = 'Fall', color = '#F4B860')
ggarrange(scatter.winter, scatter.spring, scatter.summer, scatter.fall, nrow = 2, ncol = 2, align = 'hv')
ggsave("cq_slope_cvratio.tiff", device = "tiff", height = 6, width = 8, units = "in", compression = "lzw", dpi = 700)

### Slope v Pearson R ###
#All points
ggplot(data = DeS_byDate) +
  geom_point(aes(x = pearson, y = CQslope, color = year))

#Trimmed of windows with > 15 NAs
ggplot(data = DeS_byDate) +
  geom_point(aes(x = pearson.trim, y = CQslope.trim, color = year), size = 1) +
  coord_cartesian(xlim = c(-1, 1), ylim = c(-5,5)) 
  

#Trimmed, colored by season instead of year
ggplot(data = DeS_byDate) +
  geom_point(aes(x = pearson.trim, y = CQslope.trim, color = season), size = 1)


#### Time Series ####

### C and Q Parallel ###
#By yday
ggplot(data = DeS_byyday) +
  geom_line(aes(x = yday, y = avgC)) +
  geom_line(aes(x = yday, y = avgQ/1000), color = 'blue') +
  scale_y_continuous(sec.axis = sec_axis(~ . *1000))

#By date
ggplot(data = DeS_byDate) +
  geom_line(aes(x = date, y = avgC/max(avgC, na.rm=TRUE)), color = 'darkgreen') +
  geom_line(aes(x = date, y = avgQ/max(avgQ, na.rm = TRUE)*5), color = 'blue') +
  geom_col(data = precip, aes(x = date, y = prcp/max(prcp, na.rm = TRUE))) +
  #scale_y_continuous(sec.axis = sec_axis(~ . *1000)) +
  coord_cartesian(xlim = date(c("2016-10-01", "2017-02-01")), ylim = c(0,1)) +
  scale_x_date(date_labels = "%j")


#### Yearly ####

### Distributions of Q and CQ for each year ###
ggplot(data = DeS_byDate) +
  geom_violin(aes(x = year, y = logQ), fill = 'lightblue') +
  geom_violin(aes(x = year, y = CQslope.trim), fill = 'pink') +
  geom_line(data = DeS_byYear, aes(x = year, y = medianlogQ, group = 1), color = 'darkblue', linetype = 'dashed') +
  geom_point(data = DeS_byYear, aes(x = year, y = medianlogQ), color = 'darkblue') +
  geom_line(data = DeS_byYear, aes(x = year, y = medianCQslope.trim, group = 1), color = 'red', linetype = 'dashed') +
  geom_point(data = DeS_byYear, aes(x = year, y = medianCQslope.trim), color = 'red') +
  coord_cartesian(ylim = c(-3, 5)) +
  theme_bw() +
  theme(panel.grid = element_blank())

ggplot(data = DeS_byYear, aes(x = year, group = 1)) +
  geom_line(aes(y = medianlogQ), color = 'lightblue') +
  geom_line(aes(y = medianCQslope.trim), color = 'pink')

ggplot(data = DeS_byYear) +
  geom_point(aes(x = medianlogQ, y = medianCQslope.trim))


#### Seasonal ####

## Seasonal Bars - Q and C ##

ggplot(data = seasons) +
  geom_col(aes(x = season, y = meanQ*35.32, fill = meanC)) +
  scale_fill_gradient(low = '#c3c3c3', high = '#b288c0', limits = c(0.7, max(seasons$meanC)), name = "Average\nNitrate (mg/L)") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_cartesian(expand = TRUE, ylim = c(18*35.32, 375*35.32)) +
  geom_text(aes(label = str_c(round(meanC, digits = 2)," mg/L"), x = season, y = meanQ*35.32/2), size = 3) +
  scale_y_continuous(breaks = c(24, 5000, 10000), labels = c(0, 5000, 10000)) +
  xlab("") +
  ylab("Average Discharge (cfs)")
ggsave("cqseasonal_bar.tiff", device = "tiff", height = 3.5, width = 5, units = "in", compression = "lzw", dpi = 700)

## Q and CQ slope distributions by Season ##

ggplot(data = DeS_byDate) +
  geom_violin(aes(x = season, y = logQ), fill = 'lightblue') +
  geom_violin(aes(x = season, y = CQslope.trim, fill = season), alpha = 0.6) +
  scale_fill_manual(values = c('#8DA9C4','#60992D','#C33C54','#F4B860'), guide = 'none') +
  geom_hline(yintercept = 0, linetype = 'solid') +
  geom_hline(yintercept = 0.54, linetype = 'dotted') +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_y_continuous(limits = c(-3, 5), breaks = seq(-3, 5), 
                     sec.axis = sec_axis(~., name = "30-Day CQ Slope", breaks = seq(-3, 5))) +
  ylab("Log Discharge (m^3/s)") +
  xlab("")
ggsave("cqseasonal_violin.tiff", device = "tiff", height = 3.5, width = 5, units = "in", compression = "lzw", dpi = 700)


## Precip by Season ##
ggplot(data = precip_bySeason) +
  geom_col(aes(x = season, y = avgP), fill = 'lightblue') +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_cartesian(expand = TRUE, ylim = c(18, 340)) +
  labs(title = "Historical Data") +
  xlab("Season") +
  ylab("Average Precipitation (mm)")






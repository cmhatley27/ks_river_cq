library(tidyverse)
library(lubridate)
library(dataRetrieval)
library(hydroEvents)
library(rnoaa)
library(ggpubr)
library(ggforce)

####Create figures subfolder####
if(!dir.exists('Figures')) {
  dir.create('Figures')
}


############
# Plotting #
############


#### LogC v LogQ####


### Log - Log by Date ####
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
ggsave("./Figures/cqloglogscatter.tiff", device = "tiff", height = 6, width = 8, units = "in", compression = "lzw", dpi = 700)

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
ggsave("./Figures/cqloglogscatter_inset.tiff", device = "tiff", height = 2, width = 2, units = "in", compression = "lzw", dpi = 700)

### Log - Log by yday ###
ggplot(data = DeS_byyday) +
  geom_point(aes(x = logQ, y = logC, color = season)) +
  geom_vline(xintercept = log10(312.5), color = 'blue') +
  scale_color_manual(values = c('#8DA9C4','#60992D','#C33C54','#F4B860'))


#### 30 Day Windows ####

#* Time Series colored by year ####
#All points #
ggplot(data = DeS_byDate) +
  geom_line(aes(x = yday, y = CQslope, color = factor(wateryear))) +
  coord_cartesian(ylim = c(-2,2))

#Trimmed of windows with > 15 NAs #
ggplot(data = DeS_byDate) +
  geom_line(aes(x = yday, y = CQslope.trim, color = factor(wateryear)), size = 0.75) +
  coord_cartesian(ylim = c(-2,2)) +
  geom_hline(yintercept = 0.5083, linetype = 'dashed') +
  geom_hline(yintercept = 0)


#* Slope v CVratio ####
#All points #
ggplot(data = DeS_byDate) +
  geom_point(aes(x = cvratio, y = CQslope, color = factor(wateryear)))

#**Trimmed of windows with > 15 NAs ####
ggplot(data = DeS_byDate) +
  geom_point(aes(x = cvratio.trim, y = CQslope.trim, color = factor(wateryear)), size = 1) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 1)

#***Trimmed, colored by season instead of year ####
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
ggsave("./Figures/cq_slope_cvratio_combined.tiff", device = "tiff", height = 6, width = 8, units = "in", compression = "lzw", dpi = 700)

#**Paneled by season ####
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
ggsave("./Figures/cq_slope_cvratio.tiff", device = "tiff", height = 6, width = 8, units = "in", compression = "lzw", dpi = 700)

### Slope v Pearson R ####
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

#NAs in Q, every point
ggplot(data = DeSoto_15min) +
  geom_line(aes(x = dateTime, y = N), alpha = 0.3) +
  geom_line(aes(x = dateTime, y = Q/100), color = 'blue') +
  geom_point(data = DeSoto_15min[is.na(DeSoto_15min$Q),], aes(x = dateTime, y = 0), color = 'red', size = 0.5) +
  xlim(as_datetime(c('2013-10-01', '2018-09-30')))
#NAs in Q, hourly average
ggplot(data = DeS_byHour) +
  geom_line(aes(x = hour, y = avgC), alpha = 0.3) +
  geom_line(aes(x = hour, y = avgQ/100), color = 'blue') +
  geom_point(data = DeS_byHour[is.na(DeS_byHour$avgQ),], aes(x = hour, y = 0), color = 'red', size = 0.5) +
  xlim(as_datetime(c('2013-10-01', '2022-09-30')))
#NAs in Q, daily average
ggplot(data = DeS_byDate) +
  geom_line(aes(x = date, y = avgC), alpha = 0.3) +
  geom_line(aes(x = date, y = avgQ/100), color = 'blue') +
  geom_point(data = DeS_byDate[is.na(DeS_byDate$avgQ),], aes(x = date, y = 0), color = 'red', size = 0.5) +
  xlim(as_date(c('2017-10-01', '2022-09-30')))
  
#NAs in N, every point
ggplot(data = DeS.all) +
  geom_line(aes(x = dateTime, y = NO3)) +
  geom_line(aes(x = dateTime, y = Q_m3s/100), color = 'blue', alpha = 0.3) +
  geom_point(data = DeS.all[is.na(DeS.all$NO3),], aes(x = dateTime, y = 0), color = 'red', size = 0.5) +
  xlim(as_datetime(c('2018-10-01', '2022-09-30')))

#### Yearly ####

### Distributions of Q and CQ for each year ###
ggplot(data = DeS_byDate) +
  geom_violin(aes(x = factor(wateryear), y = logQ), fill = 'lightblue') +
  geom_violin(aes(x = factor(wateryear), y = CQslope.trim), fill = 'pink') +
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


###Discharge Statistics vs CQ slope###
ggplot(data = DeS_byYear) +
  geom_point(aes(x = totQ, y = slope, color = factor(wateryear)))

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
ggsave("./Figures/cqseasonal_bar.tiff", device = "tiff", height = 3.5, width = 5, units = "in", compression = "lzw", dpi = 700)

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
ggsave("./Figures/cqseasonal_violin.tiff", device = "tiff", height = 3.5, width = 5, units = "in", compression = "lzw", dpi = 700)

## Q and CV ratio distributions by Season ##
ggplot(data = DeS_byDate) +
  geom_violin(aes(x = season, y = logQ), fill = 'lightblue') +
  geom_violin(aes(x = season, y = cvratio.trim, fill = season), alpha = 0.6) +
  scale_fill_manual(values = c('#8DA9C4','#60992D','#C33C54','#F4B860'), guide = 'none') +
  geom_hline(yintercept = 1, linetype = 'solid') +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_y_continuous(limits = c(0, 8), breaks = seq(-3, 8), expand = c(0,0), 
                     sec.axis = sec_axis(~., name = "30-Day CVc/CVq", breaks = seq(-3, 5))) +
  ylab("Log Discharge (m^3/s)") +
  xlab("")
ggsave("./Figures/cq_cv_seasonal_violin.tiff", device = "tiff", height = 3.5, width = 5, units = "in", compression = "lzw", dpi = 700)
## Precip by Season ##
ggplot(data = precip_bySeason) +
  geom_col(aes(x = season, y = avgP), fill = 'lightblue') +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_cartesian(expand = TRUE, ylim = c(18, 340)) +
  labs(title = "Historical Data") +
  xlab("Season") +
  ylab("Average Precipitation (mm)")


#### Flow Percentile ####

##FDC##
ggplot(data = FDC) +
  geom_line(aes(x = P, y = log10(avgQ), group = 1))

##CQ vs CV##
#Trimmed CQ slopes
ggplot(data = DeS_byDate[!is.na(DeS_byDate$percentile.bin),]) +
  geom_rect(aes(xmin = 0, xmax = 1, ymin = -10, ymax = 10), fill = 'lightgrey', alpha = 0.05) +
  geom_point(aes(x = cvratio.trim, y = CQslope.trim, color = percentile.bin), size = 1, alpha = 1) +
  scale_color_manual(values = c('indianred3', 'indianred1', 'grey', 'skyblue', 'deepskyblue3')) +
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
?scale_color_discrete
##Distribution of CQ Slopes with Seasonality
ggplot(data = DeS_byDate[!is.na(DeS_byDate$percentile.bin),]) +
  geom_jitter(aes(x = percentile.bin, y = CQslope.trim, color = season), alpha = 1, width = 0.2, height = 0, size = 0.75) +
  geom_violin(aes(x = percentile.bin, y = CQslope.trim), alpha = 0, size = 1) +
  scale_color_manual(values = c('#8DA9C4','#60992D','#C33C54','#F4B860'), guide = 'none') +
  geom_hline(yintercept = 0, linetype = 'solid') +
  geom_hline(yintercept = 0.54, linetype = 'dotted') +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_cartesian(ylim = c(-6, 6)) +
  xlab("Flow Percentile") +
  ylab("CQ slope for 30-day windows")
ggsave("./Figures/CQ_slope_byPercentile_violin.tiff", device = "tiff", height = 3.5, width = 5, units = "in", compression = "lzw", dpi = 700)

##Distribution of CV ratios with Seasonality
ggplot(data = DeS_byDate[!is.na(DeS_byDate$percentile.bin),]) +
  geom_jitter(aes(x = percentile.bin, y = cvratio.trim, color = season), alpha = 1, width = 0.2, height = 0, size = 0.75) +
  geom_violin(aes(x = percentile.bin, y = cvratio.trim), alpha = 0, size = 1) +
  scale_color_manual(values = c('#8DA9C4','#60992D','#C33C54','#F4B860'), guide = 'none') +
  geom_hline(yintercept = 1, linetype = 'solid') +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_cartesian(ylim = c(0, 8)) +
  scale_y_continuous(expand = c(0,0)) +
  xlab("Flow Percentile") +
  ylab("CV ratio fro 30-day windows")
ggsave("./Figures/CV_ratio_byPercentile_violin.tiff", device = "tiff", height = 3.5, width = 5, units = "in", compression = "lzw", dpi = 700)


##Load by percentile and year
ggplot(data = loadbypercentile) +
  geom_col(aes(x = factor(wateryear), y = load.frac, fill = percentile.bin), color = 'black', position = position_fill(reverse = TRUE))+
  scale_fill_manual(values = c('#282935','#455063','#607c94','#79acc4','#94def2'), guide = guide_legend(reverse = TRUE), name = "Flow Percentile") +
  xlab("") +
  ylab("Fraction of Annual Nitrate Load")
ggsave("./Figures/N_load_byPercentile_bars.tiff", device = "tiff", height = 6, width = 8, units = "in", compression = "lzw", dpi = 700)



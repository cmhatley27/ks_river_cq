
# Load common libraries and data ------------------------------------------
library(tidyverse)
library(ggpubr)
source('load_and_combine.R')
source('Theme+Settings.R')

all_hydro_15min <- with_tz(read_csv('./DataFiles/hydro_data/all_hydro_15min.csv'), 'America/Chicago')
all_hydro_daily <- read_csv('./DataFiles/hydro_data/all_hydro_daily.csv')

# Log CQ ------------------------------------------------------------------

#Create data
cq_date <- select(all_hydro_daily, date, wateryear, riverQ, N) %>%
  mutate(logQ = log10(riverQ), logC = log10(N),
         wateryear = factor(wateryear)) %>%
  filter(wateryear != '2013' & wateryear != '2015') 

ggplot(data = cq_date) +
  geom_point(aes(x = logQ, y = logC, color = wateryear), alpha = 0.18, size = 1.5, shape = 19) +
  geom_smooth(aes(x = logQ, y = logC, group = wateryear, color = wateryear), method = 'lm', se = FALSE, fullrange = TRUE, linewidth = 1.1) +
  scale_color_viridis_d(name = 'Water Year') +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  xlab("Log Discharge (m^3/s)") +
  ylab("Log Concentration of Nitrate (mg/L)")
ggsave('./Figures/LogCQ_byyear.tiff', device = "tiff", height = 6, width = 8, units = "in", compression = "lzw", dpi = 700)

summary(lm(logC ~ logQ, data = cq_date))

ggplot(data = cq_date) +
  geom_point(aes(x = logQ, y = logC, color = wateryear), alpha = 0.6, size = 1.5, shape = 19) +
  geom_smooth(aes(x = logQ, y = logC, group = wateryear, color = wateryear), method = 'lm', se = TRUE, fullrange = TRUE, linewidth = 1.1) +
  scale_color_viridis_d(name = 'Water Year', end = 1) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  xlab("Log Discharge (m^3/s)") +
  ylab("Log Concentration of Nitrate (mg/L)") +
  facet_wrap(vars(wateryear))
ggsave('./Figures/LogCQ_byyear_grid.tiff', device = "tiff", height = 9, width = 12, units = "in", compression = "lzw", dpi = 700)

#Do year regressions and compare to total water/N fluxes from next section [NEED TO RUN NEXT SECTION FOR THIS]
cq_date_reg <- cq_date %>%
  group_by(wateryear) %>%
  do(b = unlist(as.numeric(coef(lm(logC ~ logQ, data = .))[2])),
     a = unlist(as.numeric(coef(lm(logC ~ logQ, data = .))[1]))) %>%
  mutate(slope = unlist(b),
         intercept = unlist(a)) %>%
  select(wateryear, slope, intercept)

cumul_flux_yearly <- flux_data %>%
  mutate(wateryear = factor(wateryear)) %>%
  filter(wateryear != '2015') %>%
  group_by(wateryear) %>%
  summarize(total_Q = tail(Q_cumul, 1),
            total_N = tail(N_cumul, 1)) %>%
  left_join(cq_date_reg)

ggplot(data = cumul_flux_yearly) +
  geom_point(aes(x = total_Q, y = slope))



# Cumulative fluxes --------------------------------------------------

##Create data
flux_data <- all_hydro_15min %>%
  select(dateTime, wateryear, Q = riverQ_linterp, N, res_sum = reservoir_sum_linterp) %>%
  mutate(dateTime_yearadj = ymd_hms(paste0(if_else(month(dateTime) > 9, '2019', '2020'),'-',month(dateTime),'-',day(dateTime),' ',hour(dateTime),':',minute(dateTime),':',second(dateTime))),
         Q_load = Q*60*15,
         N_load = N*Q*1000*60*15/1000000,
         res_load = res_sum*60*15) %>%
  mutate(N_load = if_else(is.na(N_load), 0, N_load)) %>%
  group_by(wateryear) %>%
  mutate(Q_cumul = cumsum(Q_load),
         N_cumul = cumsum(N_load),
         res_cumul = cumsum(res_load))


##Plots of each series
Qplot <- ggplot(data = flux_data, aes(x = dateTime_yearadj)) +
  geom_line(aes(y = Q_cumul, group = wateryear, color = as_factor(wateryear))) +
  scale_x_datetime(date_breaks = '1 month', date_labels = '%b') +
  scale_color_discrete(name = 'Water Year') +
  theme_classic() +
  xlab('') +
  ylab('Cumulative Q (m3)')
Qplot

Nplot <- ggplot(data = flux_data, aes(x = dateTime_yearadj)) +
  geom_line(aes(y = N_cumul, group = wateryear, color = as_factor(wateryear))) +
  scale_x_datetime(date_breaks = '1 month', date_labels = '%b') +
  scale_color_discrete(name = 'Water Year') +
  theme_classic() +
  xlab('') +
  ylab('Cumulative N (kg)')
Nplot

resplot <- ggplot(data = flux_data, aes(x = dateTime_yearadj)) +
  geom_line(aes(y = res_cumul, group = wateryear, color = as_factor(wateryear))) +
  scale_x_datetime(date_breaks = '1 month', date_labels = '%b') +
  scale_color_discrete(name = 'Water Year') +
  theme_classic() +
  xlab('') +
  ylab('Cumulative Outflow (m3)')
resplot


##Plot all together
flux_facet <- flux_data %>%
  pivot_longer(c(Q_cumul, N_cumul, res_cumul), names_to = 'vars', values_to = 'loads')

ggplot(data = flux_facet, aes(x = dateTime_yearadj)) +
  geom_line(aes(y = loads, group = wateryear, color = as_factor(wateryear)), size = 1) +
  scale_x_datetime(date_breaks = '1 month', date_labels = '%b') +
  scale_color_discrete(name = 'Water Year') +
  theme_classic() +
  xlab('') +
  ylab('Cumulative Loads (m3 or kg)') +
  facet_wrap(~vars, ncol = 1, scales = 'free_y')

# FDC & Loads by flow percentile ---------------------------------------------------------------------

#Create Data
FDC_data <- all_hydro_15min %>%
  select(dateTime, wateryear, Q = riverQ_linterp, N) %>%
  mutate(Q_load = Q*60*15,
         N_load = N*Q*1000*60*15/1000000) %>%
  arrange(desc(Q)) %>%
  mutate(rank = seq(1:length(Q)),
         P = rank/(length(Q)+1)*100,
         percentile = 100 - P,
         percentile_bin = cut(percentile, breaks = c(0, 10, 40, 60, 90, 100), labels = c('0 - 10', '10 - 40', '40 - 60', '60 - 90', '90 - 100')))

FDC_data_byyear <- all_hydro_15min %>%
  select(dateTime, wateryear, Q = riverQ_linterp, N) %>%
  mutate(Q_load = Q*60*15,
         N_load = N*Q*1000*60*15/1000000) %>%
  group_by(wateryear) %>%
  arrange(desc(Q)) %>%
  mutate(rank = seq(1:length(Q)),
         P = rank/(length(Q)+1)*100,
         percentile = 100 - P,
         percentile_bin = cut(percentile, breaks = c(0, 10, 40, 60, 90, 100), labels = c('0 - 10', '10 - 40', '40 - 60', '60 - 90', '90 - 100')))

percentile_loads <- FDC_data %>%
  group_by(wateryear, percentile_bin) %>%
  summarise(load_within_percentile = sum(N_load, na.rm = TRUE)) %>%
  group_by(wateryear) %>%
  mutate(load_frac = load_within_percentile/sum(load_within_percentile))


##Full FDC
ggplot(data = FDC_data, aes(x = P, y = Q)) +
  geom_line()

##FDC per year
ggplot(data = FDC_data_byyear, aes(x = P, y = Q)) +
  geom_line(aes(color = as_factor(wateryear)))

##Loads by flow percentile and year
ggplot(data = percentile_loads) +
  geom_col(aes(x = factor(wateryear), y = load_frac, fill = percentile_bin), color = 'black', position = position_fill(reverse = TRUE))+
  scale_fill_manual(values = c('#282935','#455063','#607c94','#79acc4','#94def2'), guide = guide_legend(reverse = TRUE), name = "Flow Percentile") +
  xlab("") +
  ylab("Fraction of Annual Nitrate Load")


# Loads by event/nonevent, also number of event/nonevent observations -------------------------------------------------

#Calculate loads and summarize per year
N_yields <- all_hydro_15min %>%
  mutate(N_smooth = (lag(N, 4) + lag(N, 3) + lag(N, 2) + lag(N, 1) + N + lead(N, 1) + lead(N, 2) + lead(N, 3) + lead(N, 4))/9) %>%
  select(dateTime, wateryear, riverQ_linterp, N, N_smooth, event) %>%
  mutate(N_yield = (riverQ_linterp*1000)*(N_smooth/1000000)*60*15)

yield_summary <- N_yields %>%
  group_by(wateryear) %>%
  summarize(event_yield = sum(N_yield[!is.na(event)], na.rm = TRUE),
            nonevent_yield = sum(N_yield[is.na(event)], na.rm = TRUE)) %>%
  mutate(frac_event = event_yield/(event_yield+nonevent_yield),
         frac_nonevent = 1-frac_event)
yield_summary_long <- yield_summary %>%
  select(wateryear, event = frac_event, nonevent = frac_nonevent) %>%
  pivot_longer(-wateryear, names_to = 'type', values_to = 'frac')

##Loads by event/nonevent and year
ggplot(data = yield_summary_long) +
  geom_col(aes(x = factor(wateryear), y = frac, fill = type, alpha = type), color = 'black', position = position_fill(reverse = TRUE)) +
  scale_fill_manual(values = c('#607c94','white'), guide = guide_legend(reverse = TRUE), name = NULL) +
  scale_alpha_manual(values = c(100, 0), guide = 'none') +
  ggtitle('Fraction of Annual N Load during Events') +
  theme_classic()

#Calculate number of event/nonevent observations
event_obs <- all_hydro_15min %>%
  group_by(wateryear) %>%
  summarise(num_event = length(unique(dateTime[!is.na(event)])),
            num_nonevent = length(unique(dateTime[is.na(event)]))) %>%
  mutate(frac_event = num_event/(num_event + num_nonevent),
         frac_nonevent = 1-frac_event)
event_obs_long <- event_obs %>%
  select(wateryear, event = frac_event, nonevent = frac_nonevent) %>%
  pivot_longer(-wateryear, names_to = 'type', values_to = 'frac')

##Number of observations by event/nonevent and year
ggplot(data = event_obs_long) +
  geom_col(aes(x = factor(wateryear), y = frac, fill = type, alpha = type), color = 'black', position = position_fill(reverse = TRUE)) +
  scale_fill_manual(values = c('#607c94','white'), guide = guide_legend(reverse = TRUE), name = NULL) +
  scale_alpha_manual(values = c(100, 0), guide = 'none') +
  ggtitle('Fraction of observations during Events') +
  theme_classic()


# Monthly/Yearly Q by proportion of reservoir outflow -----------------------------

#calculate proportions of flow
Qs <- all_hydro_15min %>%
  select(dateTime, wateryear, ends_with('_linterp')) %>%
  select(!c(lawrenceQ_linterp, N_linterp)) %>%
  mutate(wateryear = factor(wateryear))

monthly_prop <- Qs %>%
  group_by(year = year(dateTime), month = month(dateTime)) %>%
  summarize(prop = sum(reservoir_sum_linterp)/sum(riverQ_linterp)) %>%
  mutate(wateryear = ifelse(month > 9, year+1, year))

mean(monthly_prop$prop)

monthly_prop_avg <- monthly_prop %>%
  group_by(month) %>%
  summarise(prop = mean(prop))
yearly_prop_avg <- monthly_prop %>%
  group_by(wateryear) %>%
  summarise(prop = mean(prop))

##Range of proportions of reservoir outflow
ggplot() +
  geom_histogram(data = monthly_prop, aes(x = prop), bins = 20) +
  xlim(0,1)

##Monthly proportions of reservoir outflow
ggplot() +
  geom_col(data = monthly_prop_avg, aes(x = factor(month), y = prop)) +
  ylim(c(0,1))

##Yearly proportions of reservoir outflow
ggplot() +
  geom_col(data = yearly_prop_avg, aes(x = factor(wateryear), y = prop)) +
  ylim(c(0,1))

# Reservoir Proportions ---------------------------------------------------

#gather Qs
Qs <- all_hydro_15min %>%
  select(dateTime, wateryear, ends_with('_linterp')) %>%
  select(!c(lawrenceQ_linterp, N_linterp)) %>%
  mutate(wateryear = factor(wateryear))

#sum total discharge from each source per year
Qs_year <- Qs %>%
  group_by(wateryear) %>%
  summarise(across(ends_with('_linterp'),sum)) %>%
  mutate(across(ends_with('_linterp'),~.*60*15)) %>%
  pivot_longer(!wateryear, names_to = 'res', values_to = 'Q')

#Total discharge from each source as proportion of River Q
Qs_year_props <- Qs %>%
  group_by(wateryear) %>%
  summarise(across(ends_with('_linterp'),sum)) %>%
  mutate(across(ends_with('_linterp'),~.*60*15)) %>%
  group_by(wateryear) %>%
  summarise(across(ends_with('_linterp'),~./riverQ_linterp)) %>%
  pivot_longer(!wateryear, names_to = 'res', values_to = 'prop')

group_by(Qs_year_props, res) %>%
  summarise(mean = mean(prop)) %>%
  arrange(desc(mean))

#Total discharge from each reservoir as proportion of all reservoir releases
Qs_year_props_resonly <- Qs %>%
  group_by(wateryear) %>%
  summarise(across(ends_with('_linterp'),sum)) %>%
  mutate(across(ends_with('_linterp'),~.*60*15)) %>%
  group_by(wateryear) %>%
  summarise(across(ends_with('_linterp'),~./reservoir_sum_linterp)) %>%
  pivot_longer(!wateryear, names_to = 'res', values_to = 'prop')

group_by(Qs_year_props_resonly, res) %>%
  summarise(mean = mean(prop)) %>%
  arrange(desc(mean))

##Yearly discharges from each reservoir
ggplot() +
  geom_col(data = subset(Qs_year, res == 'riverQ_linterp'), aes(x = wateryear, y = Q), color = 'black') +
  geom_col(data = subset(Qs_year, res != 'riverQ_linterp' & res != 'reservoir_sum_linterp'), aes(x = wateryear, y = Q, fill = res), color = 'black')

##Yearly proportions from each reservoir
ggplot() +
  geom_col(data = subset(Qs_year_props, res == 'riverQ_linterp'), aes(x = wateryear, y = prop), color = 'black') +
  geom_col(data = subset(Qs_year_props, res != 'riverQ_linterp' & res != 'reservoir_sum_linterp'), aes(x = wateryear, y = prop, fill = res), color = 'black')


# Monthly/Yearly precip across whole KSRB ---------------------------------
library(terra)

#Calculate precip totals
precip_monthly <- rast('C:/School/SAFE KAW/Data/DataFiles/climate_data/monthly_rasters/precip_ksrb.tiff') %>%
  terra::as.data.frame(., xy = TRUE) %>%
  pivot_longer(-c(x,y), names_to = 'date', values_to = 'precip') %>%
  group_by(date = ymd(date)) %>%
  summarize(total_precip = sum(precip)) %>%
  filter(date %within% interval('2013-10-01', '2022-09-01')) %>%
  mutate(month = factor(month(date)),
         wateryear = ifelse(month(date) > 9, year(date)+1, year(date)))

#monthly and yearly averages
monthly_precip_avg <- precip_monthly %>%
  group_by(month) %>%
  summarise(precip = mean(total_precip)) %>%
  mutate(prop = precip/sum(precip))
yearly_precip_avg <- precip_monthly %>%
  group_by(wateryear) %>%
  summarise(precip = mean(total_precip))

##Average monthly total precip
ggplot() +
  geom_col(data = monthly_precip_avg, aes(x = factor(month), y = precip))
##Average yearly total precip
ggplot() +
  geom_col(data = yearly_precip_avg, aes(x = factor(wateryear), y = precip))


# Avg reservoir release sched ---------------------------------------------
avg_releases_daily <- all_hydro_daily %>%
  filter(wateryear > 2013) %>%
  group_by(yday) %>%
  summarise(across(ends_with('Q'), mean, na.rm = TRUE)) %>%
  mutate(reservoir_sum = clintonQ + milfordQ + tuttleQ + perryQ + wilsonQ + kanopolisQ + wacondaQ,
         reservoir_ratio = reservoir_sum/riverQ)

ggplot(data = pivot_longer(avg_releases_daily, !c(yday, reservoir_ratio), names_to = 'location', values_to = 'Q')) +
  geom_line(aes(x = yday, y = Q, color = location))

ggplot(data = avg_releases_daily) +
  geom_point(aes(x = yday, y = reservoir_ratio)) +
  geom_smooth(aes(x = yday, y = reservoir_ratio)) +
  scale_y_continuous(limits = c(0,1.05))



# Q and N time series -----------------------------------------------------
all_hydro_weekly <- filter(all_hydro_daily, wateryear > 2013) %>%
  mutate(week = week(date),
         year = year(date)) %>%
  select(year, week, riverQ, N) %>%
  group_by(year, week) %>%
  summarize(weekQ = mean(riverQ, na.rm = TRUE),
            weekN = mean(N, na.rm = TRUE)) %>%
  mutate(date = ymd(paste0(year,'01-01')) + weeks(week - 1))

all_hydro_monthly <- filter(all_hydro_daily, wateryear > 2013) %>%
  mutate(datemonth = ymd(paste0(year(date),'-',month(date),'-01'))) %>%
  group_by(datemonth) %>%
  summarize(monthQ = mean(riverQ, na.rm = TRUE),
            monthN = mean(N, na.rm = TRUE))

all_hydro_yearly <- filter(all_hydro_daily, wateryear > 2013) %>%
  mutate(year = year(date)) %>%
  group_by(year) %>%
  summarize(yearQ = mean(riverQ, na.rm = TRUE),
            yearN = mean(N, na.rm = TRUE)) %>%
  filter(year > 2013)

all_hydro_yearly2 <- filter(all_hydro_15min, wateryear > 2013) %>%
  mutate(year = year(dateTime)) %>%
  group_by(year) %>%
  summarize(yearQ = mean(riverQ, na.rm = TRUE),
            yearN = mean(N, na.rm = TRUE)) %>%
  filter(year > 2013)

#15min
ggplot(data = all_hydro_15min) +
  geom_line(aes(x = dateTime, y = riverQ), color = 'black') +
  geom_line(aes(x = dateTime, y = N*100), color = 'red') +
  scale_y_continuous(sec.axis = sec_axis(trans = ~./100, name = 'Nitrate [mg/L]')) +
  ylab('Discharge [m^3/s]') +
  xlab('') +
  ggtitle('15 min Data') +
  theme_classic()

#daily
daily <- ggplot(data = subset(all_hydro_daily, wateryear > 2013)) +
  geom_line(aes(x = date, y = riverQ), color = 'black') +
  geom_line(aes(x = date, y = N*200), color = 'red', alpha = 0.7) +
  scale_y_continuous(sec.axis = sec_axis(trans = ~./200, name = 'Nitrate [mg/L]')) +
  ylab('Discharge [m^3/s]') +
  xlab('') +
  ggtitle('Daily Averages') +
  theme_classic() +
  theme(axis.title.y.right = element_text(color = 'red'),
        plot.title = element_text(size = 12))

#weekly
weekly <- ggplot(data = all_hydro_weekly) +
  geom_line(aes(x = date, y = weekQ), color = 'black') +
  geom_line(aes(x = date, y = weekN*200), color = 'red', alpha = 0.7) +
  scale_y_continuous(sec.axis = sec_axis(trans = ~./200, name = 'Nitrate [mg/L]')) +
  ylab('Discharge [m^3/s]') +
  xlab('') +
  ggtitle('Weekly Averages') +
  theme_classic() +
  theme(axis.title.y.right = element_text(color = 'red'),
        plot.title = element_text(size = 12))

#monthly
monthly <- ggplot(data = all_hydro_monthly) +
  geom_line(aes(x = datemonth, y = monthQ), color = 'black') +
  geom_line(aes(x = datemonth, y = monthN*200), color = 'red', alpha = 0.7)+
  scale_y_continuous(sec.axis = sec_axis(trans = ~./200, name = 'Nitrate [mg/L]')) +
  ylab('Discharge [m^3/s]') +
  xlab('') +
  ggtitle('Monthly Averages') +
  theme_classic() +
  theme(axis.title.y.right = element_text(color = 'red'),
        plot.title = element_text(size = 12))

#yearly
yearly <- ggplot(data = all_hydro_yearly) +
  geom_line(aes(x = year, y = yearQ), color = 'black') +
  geom_line(aes(x = year, y = yearN*200), color = 'red', alpha = 0.7) +
  scale_y_continuous(sec.axis = sec_axis(trans = ~./200, name = 'Nitrate [mg/L]')) +
  ylab('Discharge [m^3/s]') +
  xlab('') +
  ggtitle('Yearly Averages') +
  theme_classic() +
  theme(axis.title.y.right = element_text(color = 'red'),
        plot.title = element_text(size = 12))

ggarrange(daily, weekly, monthly, yearly, ncol = 1)
ggsave('./Figures/hydrographs/DeSoto_timescale_averages.tiff', device = 'tiff', width = 6.5, height = 8.5, units = 'in', compression = 'lzw')


# SPEI Time Series --------------------------------------------------------
#Load Data
spei1_ts <- read_csv('./DataFiles/hydro_data/drought/spei1_huc8.csv') %>%
  filter(date >= ymd('2013-10-01')) %>%
  mutate(scale = 1)
spei3_ts <- read_csv('./DataFiles/hydro_data/drought/spei3_huc8.csv') %>%
  filter(date >= ymd('2013-10-01')) %>%
  mutate(scale = 3)
spei6_ts <- read_csv('./DataFiles/hydro_data/drought/spei6_huc8.csv') %>%
  filter(date >= ymd('2013-10-01')) %>%
  mutate(scale = 6)
spei9_ts <- read_csv('./DataFiles/hydro_data/drought/spei9_huc8.csv') %>%
  filter(date >= ymd('2013-10-01')) %>%
  mutate(scale = 9)
spei12_ts <- read_csv('./DataFiles/hydro_data/drought/spei12_huc8.csv') %>%
  filter(date >= ymd('2013-10-01')) %>%
  mutate(scale = 12)
spei18_ts <- read_csv('./DataFiles/hydro_data/drought/spei18_huc8.csv') %>%
  filter(date >= ymd('2013-10-01')) %>%
  mutate(scale = 18)
spei24_ts <- read_csv('./DataFiles/hydro_data/drought/spei24_huc8.csv') %>%
  filter(date >= ymd('2013-10-01')) %>%
  mutate(scale = 24)


spei_ts <- rbind(spei1_ts, spei3_ts, spei6_ts, spei9_ts, spei12_ts, spei18_ts, spei24_ts) %>%
  mutate(scale = factor(scale)) %>%
  pivot_longer(!c(date,scale), names_to = 'location', values_to = 'spei')
spei_quants <- group_by(spei_ts, location, scale) %>%
  summarise(lower = quantile(spei, 0.25),
            upper = quantile(spei, 0.75),
            range = upper - lower)
      
loc_sel <- c('10270104', '10260008')
ggplot(data = subset(spei_ts, location %in% loc_sel)) +
  geom_line(aes(x = date, y = spei, color = location)) +
  geom_hline(data = subset(spei_quants, location %in% loc_sel), aes(yintercept = lower, color = location), linetype = 'dashed') +
  geom_hline(data = subset(spei_quants, location %in% loc_sel), aes(yintercept = upper, color = location), linetype = 'dashed') +
  facet_wrap(vars(scale))



spei_ts_mean <- spei_ts %>%
  group_by(date, scale) %>%
  summarise(spei = mean(spei))
spei_mean_quants <- group_by(spei_ts_mean,  scale) %>%
  summarise(lower = quantile(spei, 0.25),
            upper = quantile(spei, 0.75),
            range = upper - lower)

ggplot(data = spei_ts_mean) +
  geom_line(aes(x = date, y = spei)) +
  geom_hline(data = spei_mean_quants, aes(yintercept = lower), linetype = 'dashed') +
  geom_hline(data = spei_mean_quants, aes(yintercept = upper), linetype = 'dashed') +
  scale_color_manual(values = c('black')) +
  facet_wrap(vars(scale))


# Mean SPEI (basin-wide) per wateryear -------------------------------------------------
spei12 <- read_csv('./DataFiles/hydro_data/drought/spei9_huc8.csv') %>%
  filter(date >= ymd('2013-10-01') & date <= ymd('2022-09-30')) %>%
  mutate(wateryear = ifelse(month(date) > 9, year(date) + 1, year(date))) %>%
  pivot_longer(!c(date, wateryear), names_to = 'huc8', values_to = 'spei') %>%
  group_by(wateryear) %>%
  summarise(mean = mean(spei))

ggplot(spei12, aes(x = wateryear, y = mean)) +
  geom_line()


# Mean total precip (basin-wide) per wateryear ----------------------------
precip_mean <- read_csv('./DataFiles/hydro_data/precip/precip_huc8.csv') %>%
  filter(date >= ymd('2013-10-01') & date <= ymd('2022-09-30')) %>%
  mutate(wateryear = ifelse(month(date) > 9, year(date) + 1, year(date))) %>%
  pivot_longer(!c(date, wateryear), names_to = 'huc8', values_to = 'precip') %>%
  group_by(wateryear, huc8) %>%
  summarise(total_depth = sum(precip)) %>%
  filter(huc8 == '10270204') %>%
  group_by(wateryear) %>%
  summarise(mean_depth = mean(total_depth))

ggplot(precip_mean, aes(x = wateryear, y = mean_depth)) +
  geom_line()
# Double mass curve -------------------------------------------------------

dbl_mass <- tibble(
  date = all_hydro_daily$date,
  wateryear = factor(all_hydro_daily$wateryear),
  season = factor(month(floor_date(date, unit = 'season'))),
  q_cum = cumsum(all_hydro_daily$riverQ_linterp),
  p_cum = cumsum(all_hydro_daily$lawrence_precip)
)

ggplot(data = dbl_mass) +
  geom_line(aes(x = p_cum, y = q_cum, color = wateryear)) +
  geom_abline(slope = lm(q_cum ~ p_cum + 0, data = dbl_mass)$coefficients[1])


# CQ for all constituents -------------------------------------------------
dat <- read_csv(file.path('DataFiles','hydro_data','other_params','des_15min.csv'))

dat.daily <- group_by(dat, date(dateTime)) %>%
  summarise(q = mean(q, na.rm = T),
            n = mean(n, na.rm = T),
            c = mean(c, na.rm = T),
            t = mean(t, na.rm = T),)
dat.long <- pivot_longer(dat.daily, c(n,c,t), names_to = 'var', values_to = 'val')

ggplot(data = dat.long, aes(x = log10(q), y = log10(val))) +
  geom_point() +
  geom_smooth(method = 'lm', se = F, color = 'blue') +
  facet_wrap(vars(var), scales = 'free')


# FI and HI for all constituents ------------------------------------------
dat <- read_csv(file.path('DataFiles','hydro_data','other_params','CQ_summary.csv'))

fi.long <- pivot_longer(dat, starts_with('FI'), names_to = 'var', values_to = 'FI') %>%
  mutate(species = str_sub(var, -1)) %>%
  select(event, species, FI)
hi.long <- pivot_longer(dat, starts_with('HI'), names_to = 'var', values_to = 'HI') %>%
  mutate(species = str_sub(var, -1)) %>%
  select(event, species, HI)
dat.long <- left_join(fi.long, hi.long) %>%
  dplyr::rename(var = species)


ggplot(data = dat.long, aes(x = FI, y = HI)) +
  geom_point() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  facet_wrap(vars(var))

ggplot(data = dat.long, aes(y = HI, x = var)) +
  geom_boxplot()


# FI HI Scatter -----------------------------------------------------------
plot_dat <- filter(event_summary, !is.na(cluster)) %>%
  mutate(isflush = FI_n > 0,
         iscw = HI_n > 0,
         isres = reservoir_runoff_ratio > 0.5,
         res_season = paste0(season,'_',isres),
         res_quart = cut(reservoir_runoff_ratio, c(0,0.25,0.5,0.75,max(reservoir_runoff_ratio)),
                         labels = c(1,2,3,4)))
season_means <- group_by(plot_dat, season, isres) %>%
  summarise(FI_mean = mean(FI_n),
            HI_mean = mean(HI_n),
            FI_sd = var(FI_n)^0.5,
            HI_sd = var(HI_n)^0.5)

ggplot() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(data = plot_dat, aes(x = FI_n, y = HI_n, color = season)) +
  #geom_errorbar(data = season_means, aes(x = FI_mean, ymin = HI_mean - HI_sd, ymax = HI_mean + HI_sd, color = season)) +
  #geom_errorbar(data = season_means, aes(y = HI_mean, xmin = FI_mean - FI_sd, xmax = FI_mean + FI_sd, color = season)) +
  #geom_point(data = season_means, aes(x = FI_mean, y = HI_mean, color = season, shape = isres), size = 3) +
  scale_color_manual(values = c('#8DA9C4','#60992D','#C33C54','#F4B860')) +
  ylim(c(-1,1)) +
  xlim(c(-1,1))
ggplot(data = plot_dat, aes(y=season, x = FI_n, fill = season)) +
  geom_boxplot() +
  scale_fill_manual(values = c('#8DA9C4','#60992D','#C33C54','#F4B860')) +
  xlim(c(-1,1))

pairwise.wilcox.test(plot_dat$FI_n, plot_dat$season, p.adjust.method = 'fdr')
wilcox.test(plot_dat$FI_n[plot_dat$season == 'Winter'], plot_dat$FI_n[plot_dat$season == 'Spring'])
wilcox.test(event_summary$HI_n[event_summary$season == 'Winter'], event_summary$HI_n[event_summary$season == 'Summer'])

ggplot() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(data = plot_dat, aes(x = FI_n, y = HI_n, color = isres)) +
  ylim(c(-1,1)) +
  xlim(c(-1,1))
ggplot(data = plot_dat, aes(y=res_quart, x = FI_n, fill = res_quart, group = res_quart)) +
  geom_boxplot() +
  xlim(c(-1,1))

pairwise.wilcox.test(plot_dat$FI_n, plot_dat$res_quart, p.adjust.method = 'fdr')
wilcox.test(plot_dat$FI_n[plot_dat$res_quart == 1], plot_dat$FI_n[plot_dat$res_quart == 2])
wilcox.test(plot_dat$HI_n[plot_dat$res_quart == TRUE], plot_dat$HI_n[plot_dat$res_quart == FALSE])


ggplot(data = plot_dat, aes(y=season, x = FI_n, fill = season)) +
  geom_boxplot() +
  xlim(c(-1,1)) +
  scale_fill_manual(values = c('#8DA9C4','#60992D','#C33C54','#F4B860')) +
  facet_wrap(vars(isres), ncol = 1)
pairwise.wilcox.test(plot_dat$FI_n, plot_dat$res_season, p.adjust.method = 'bonferroni')

event_summary %>%
  filter(!is.na(cluster)) %>%
  mutate(isflush = FI_n > 0,
         iscw = HI_n > 0) %>%
  count(iscw)
85/(131+84)
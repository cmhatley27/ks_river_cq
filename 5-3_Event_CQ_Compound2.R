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
library(plotly)


# Create directories -----------------------------------------------------

if(!dir.exists('Figures/events')) {
  dir.create('Figures/events')
}

if(!dir.exists('Figures/events/raw')) {
  dir.create('Figures/events/raw')
}

# Load data --------------------------------------------

all_hydro_15min <- with_tz(read_csv('./DataFiles/hydro_data/all_hydro_15min.csv'), 'America/Chicago')
all_hydro_daily <- read_csv('./DataFiles/hydro_data/all_hydro_daily.csv')
events <- with_tz(read_csv('DataFiles/hydro_data/event_delineations/events_compound2.csv'), 'America/Chicago') %>%
  select(event_number, start_dateTime, end_dateTime)

# Calculate HI and FI -----------------------------------------------------

#Select only events, smooth N, and normalize
event_data <- all_hydro_15min %>%
  filter(!is.na(event_compound2)) %>% filter(is.na(riverQ_ingap)) %>%
  #smooth N
  mutate(N_smooth = (lag(N, 4) + lag(N, 3) + lag(N, 2) + lag(N, 1) + N + lead(N, 1) + lead(N, 2) + lead(N, 3) + lead(N, 4))/9,
         N_linterp_smooth = ifelse(is.na(N_ingap), (lag(N_linterp, 4) + lag(N_linterp, 3) + lag(N_linterp, 2) + lag(N_linterp, 1) + 
                                                      N_linterp + lead(N_linterp, 1) + lead(N_linterp, 2) + lead(N_linterp, 3) + lead(N_linterp, 4))/9, NA)) %>%
  #normalize Q and N for each event
  group_by(event_compound2) %>%
  mutate(Q_norm = (riverQ_linterp - min(riverQ_linterp))/(max(riverQ_linterp) - min(riverQ_linterp)),
         N_norm = ifelse(is.na(N_ingap), (N_linterp - min(N_linterp))/(max(N_linterp) - min(N_linterp)), NA),
         N_smooth_norm = ifelse(is.na(N_ingap), (N_linterp_smooth - min(N_linterp_smooth, na.rm = TRUE))/(max(N_linterp_smooth, na.rm = TRUE) - min(N_linterp_smooth, na.rm = TRUE)), NA))

#Initialize q_steps
q_step <- tibble(Q_norm = rep(seq(from = 0, to = 100)/100, nrow(events)),
                 event_compound2 = rep(events$event_number, each = 101))

#Split rising and falling limbs and merge q_steps onto each
rise_limbs <- event_data %>%
  group_by(event_compound2) %>%
  #rising limb = event start through max Q
  slice(1:which.max(Q_norm)) %>%
  #merge q_step and arrange by increasing Q_norm
  full_join(q_step)  %>%
  arrange(event_compound2, Q_norm) %>%
  select(dateTime, event_compound2, riverQ_linterp, Q_norm, N, N_smooth_norm) %>%
  #Interpolate between N_norms to fill in all q_steps with a concentration
  mutate(N_smooth_norm = na.approx(N_smooth_norm, x = Q_norm, xout = Q_norm, na.rm = FALSE)) %>%
  #select only q_steps
  filter(Q_norm %in% (seq(0:100)/100)) %>%
  select(event_compound2, dateTime, Q_norm, rise_N_norm = N_smooth_norm)


fall_limbs <- event_data %>%
  group_by(event_compound2) %>%
  slice(which.max(Q_norm):length(Q_norm)) %>%
  #merge q_step and arrange by increasing Q_norm
  full_join(q_step)  %>%
  arrange(event_compound2, Q_norm) %>%
  select(dateTime, event_compound2, riverQ_linterp, Q_norm, N, N_smooth_norm) %>%
  #Interpolate between N_norms to fill in all q_steps with a concentration
  mutate(N_smooth_norm = na.approx(N_smooth_norm, x = Q_norm, xout = Q_norm, na.rm = FALSE)) %>%
  #select only q_steps
  filter(Q_norm %in% (seq(0:100)/100)) %>%
  select(event_compound2, dateTime, Q_norm, fall_N_norm = N_smooth_norm)

hyst_merge <- rise_limbs %>%
  left_join(fall_limbs) %>%
  mutate(HI = rise_N_norm - fall_N_norm)

hyst_summary <- hyst_merge %>%
  group_by(event_compound2) %>%
  summarise(HI_mean = mean(HI, na.rm = TRUE),
            HI_sd = sqrt(var(HI, na.rm = TRUE)),
            HI_cv = HI_sd/HI_mean)

flush <- rise_limbs %>%
  group_by(event_compound2) %>%
  filter(!is.na(rise_N_norm)) %>%
  filter(Q_norm == max(Q_norm) | Q_norm == min(Q_norm)) %>%
  mutate(FI = rise_N_norm[2] - rise_N_norm[1])  %>%
  summarise(FI = mean(FI))

CQ_summary_compound2 <- hyst_summary %>%
  left_join(flush)
write_csv(CQ_summary_compound2, './DataFiles/hydro_data/CQ_summary_compound2.csv')




# Calculate Predictor Variables -------------------------------------------
CQ_summary_compound2 <- with_tz(read_csv('./DataFiles/hydro_data/CQ_summary_compound2.csv'), 'America/Chicago')

#Times/seasonality
event_temporal <- events %>%
  mutate(event_compound2 = event_number,
         wateryear = factor(ifelse(as.integer(month(start_dateTime)) >= 10, as.integer(year(start_dateTime)) + 1, year(start_dateTime))),
         month = factor(month(start_dateTime, label = TRUE)),
         season = factor(month(floor_date(start_dateTime, unit = 'season')), levels = c(12,3,6,9), labels = c('Winter', 'Spring', 'Summer', 'Fall')),
         duration = as.duration(interval(start = start_dateTime, end = end_dateTime))) %>%
  select(event_compound2, start_dateTime, end_dateTime, wateryear, month, season, duration)
#join to event summary
event_summary <- CQ_summary_compound2 %>%
  left_join(event_temporal)

#Event size/load
event_size <- event_data %>%
  group_by(event_compound2) %>%
  #runoff: integrate over time, then normalize by watershed area (in m2), then convert to mm
  #N: multiply concentration by runoff converted to L/s, integrate over time, normalize by watershed area (gives mg/m2 which is equivalent to kg/km2) 
  summarize(runoff = ifelse(any(!is.na(riverQ_ingap)), NA, sum(riverQ_linterp*60*15/(1.54767E11)*1000)),
            N_yield = ifelse(any(!is.na(N_ingap)), NA, sum(N_linterp_smooth*(riverQ_linterp*1000)*60*15/(1.54767E11))),
            reservoir_ratio = max(reservoir_sum_linterp, na.rm = TRUE)/max(riverQ_linterp, na.rm = TRUE),
            initial_Q = riverQ_linterp[1])
event_summary <- event_summary %>%
  left_join(event_size)
write_csv(event_summary, './DataFiles/hydro_data/event_summary_compound2.csv')


# Plotting ----------------------------------------------------------------

#Load Data
event_summary_compound2 <- with_tz(read_csv('./DataFiles/hydro_data/event_summary_compound2.csv'), 'America/Chicago') %>%
  mutate(wateryear = factor(ifelse(as.integer(month(start_dateTime)) >= 10, as.integer(year(start_dateTime)) + 1, year(start_dateTime))),
         month = factor(month(start_dateTime, label = TRUE)),
         season = factor(month(floor_date(start_dateTime, unit = 'season')), levels = c(12,3,6,9), labels = c('Winter', 'Spring', 'Summer', 'Fall')))

###Monthly averages
monthly_averages_compound2 <- event_summary_compound2 %>%
  group_by(month) %>%
  summarise(mean_FI = mean(FI, na.rm = TRUE),
            sd_FI = sd(FI, na.rm = TRUE),
            mean_HI = mean(HI_mean, na.rm = TRUE),
            sd_HI = sd(HI_mean, na.rm = TRUE))
###Seasonal averages
seasonal_averages_compound2 <- event_summary_compound2 %>%
  group_by(season) %>%
  summarise(mean_FI = mean(FI, na.rm = TRUE),
            sd_FI = sd(FI, na.rm = TRUE),
            mean_HI = mean(HI_mean, na.rm = TRUE),
            sd_HI = sd(HI_mean, na.rm = TRUE))
###Yearly averages
yearly_averages_compound2 <- event_summary_compound2 %>%
  group_by(wateryear) %>%
  summarise(mean_FI = mean(FI, na.rm = TRUE),
            sd_FI = sd(FI, na.rm = TRUE),
            mean_HI = mean(HI_mean, na.rm = TRUE),
            sd_HI = sd(HI_mean, na.rm = TRUE))
#HI histogram
ggplot(data = event_summary_compound2) +
  geom_histogram(aes(x = HI_mean), binwidth = 0.05) +
  scale_x_continuous(limits = c(-1,1))

#FI histogram
ggplot(data = event_summary_compound2) +
  geom_histogram(aes(x = FI)) +
  scale_x_continuous(limits = c(-1,1))

#HI/FI scatter
ggplot(data = event_summary_compound2) +
  geom_point(aes(x = FI, y = HI_mean)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  scale_x_continuous(limits = c(-1,1), name = 'Flushing Index') +
  scale_y_continuous(limits = c(-1,1), name = "Hysteresis Index")

#HI/FI scatter by month
ggplot(data = event_summary_compound2) +
  geom_point(aes(x = FI, y = HI_mean, color = season)) +
  geom_errorbar(data = monthly_averages_compound2, aes(y = mean_HI, xmin = mean_FI - sd_FI, xmax = mean_FI + sd_FI), alpha = 0.3) +
  geom_errorbar(data = monthly_averages_compound2, aes(x = mean_FI, ymin = mean_HI - sd_HI, ymax = mean_HI + sd_HI), alpha = 0.3) +
  geom_point(data = monthly_averages_compound2, aes(x = mean_FI, y = mean_HI), shape = 15, size = 2) +
  scale_color_manual(values = c('#8DA9C4','#60992D','#C33C54','#F4B860'), name = "Season") +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  scale_x_continuous(limits = c(-1,1)) +
  scale_y_continuous(limits = c(-1,1)) +
  theme_bw() +
  facet_wrap(vars(month), nrow = 4)

#HI/FI scatter by season
ggplot(data = event_summary_compound2) +
  geom_point(aes(x = FI, y = HI_mean, color = season)) +
  geom_errorbar(data = seasonal_averages_compound2, aes(y = mean_HI, xmin = mean_FI - sd_FI, xmax = mean_FI + sd_FI), alpha = 0.3) +
  geom_errorbar(data = seasonal_averages_compound2, aes(x = mean_FI, ymin = mean_HI - sd_HI, ymax = mean_HI + sd_HI), alpha = 0.3) +
  geom_point(data = seasonal_averages_compound2, aes(x = mean_FI, y = mean_HI, color = season), shape = 15, size = 2) +
  scale_color_manual(values = c('#8DA9C4','#60992D','#C33C54','#F4B860'), name = "Season") +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  scale_x_continuous(limits = c(-1,1)) +
  scale_y_continuous(limits = c(-1,1)) +
  facet_wrap(vars(season))

#HI/FI scatter by year
ggplot(data = event_summary_compound2) +
  geom_point(aes(x = FI, y = HI_mean, color = wateryear), alpha = 0.6) +
  geom_errorbar(data = yearly_averages_compound2, aes(y = mean_HI, xmin = mean_FI - sd_FI, xmax = mean_FI + sd_FI), alpha = 0.3) +
  geom_errorbar(data = yearly_averages_compound2, aes(x = mean_FI, ymin = mean_HI - sd_HI, ymax = mean_HI + sd_HI), alpha = 0.3) +
  geom_point(data = yearly_averages_compound2, aes(x = mean_FI, y = mean_HI, color = wateryear), shape = 15, size = 2) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  scale_x_continuous(limits = c(-1,1)) +
  scale_y_continuous(limits = c(-1,1)) +
  facet_wrap(vars(wateryear))

#HI/FI scatter by total runoff
ggplot(data = event_summary_compound2) +
  geom_point(aes(x = FI, y = HI_mean, color = log(runoff)), alpha = 1) +
  scale_color_binned(type = 'viridis') +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  scale_x_continuous(limits = c(-1,1)) +
  scale_y_continuous(limits = c(-1,1)) +
  facet_wrap(vars(month), nrow = 4)

#HI/FI scatter by initial discharge
ggplot(data = subset(event_summary_compound2, !is.nan(HI_mean))) +
  geom_point(aes(x = FI, y = HI_mean, color = log(initial_Q)), alpha = 1) +
  scale_color_gradient2(low = 'blue', mid = 'orange', high = 'red', midpoint = 5.5) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  scale_x_continuous(limits = c(-1,1)) +
  scale_y_continuous(limits = c(-1,1))

#HI/FI scatter by reservoir ratio
ggplot(data = event_summary_compound2) +
  geom_point(aes(x = FI, y = HI_mean, color = reservoir_ratio), alpha = 1) +
  scale_color_gradient2(low = 'blue', mid = 'orange', high = 'red', midpoint = 0.8) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  scale_x_continuous(limits = c(-1,1)) +
  scale_y_continuous(limits = c(-1,1))

?scale_color_gradient2

event_summary_compound2 %>%
  filter(!is.nan(HI_mean)) %>%
  mutate(log_runoff = log(runoff),
         log_initial_Q = log(initial_Q)) %>%
  select(HI_mean, FI, log_runoff, reservoir_ratio, log_initial_Q) %>%
  cor(., method = 'spearman')
# Graph single event ------------------------------------------------------

##127 - classic counterclockwise
##167 - counterclockwise with release
##139 - clockwise
##176 - clockwise
##149 - simple diluting, low hysteresis
##150 - simple flushing, low hysteresis
##140 - ???
event_select <- c(140)
#Graph event, time series
ts <- event_data %>%
  filter(event == event_select) %>%
  gather(key = "var", value = "val", c(reservoir_sum_linterp, riverQ, N_linterp_smooth)) %>%
  ggplot(aes(x = dateTime, y = val, color = dateTime)) +
  geom_point() +
  scale_color_datetime(low = 'blue', high = 'red') +
  facet_wrap(~var, ncol = 1, scales = "free_y")

#Graph event, CQ
cq <- event_data %>%
  filter(event == event_select) %>%
  ggplot(aes(x = riverQ_linterp, y = N_linterp_smooth, color = dateTime)) +
  geom_path() +
  scale_color_datetime(low = 'blue', high = 'red')

ggarrange(ts, cq, ncol = 2, legend = 'none')
ggsave(str_c('./Figures/events/raw/', as.character(event_select), '.tiff'), device = 'tiff', height = 5, width = 8, units =  'in', compression = 'lzw', dpi = 700)

cq_norm <- hyst_merge %>%
  filter(event == event_select) %>%
  ggplot(aes(x = Q_norm)) +
  geom_point(aes(y = rise_N_norm), color = 'blue') +
  geom_point(aes(y = fall_N_norm), color = 'red')
cq_norm

for(i in 1:nrow(events)) {
  event_select <- c(i)
  #Graph event, time series
  ts <- event_data %>%
    filter(event == event_select) %>%
    gather(key = "var", value = "val", c(reservoir_sum_linterp, riverQ, N_linterp_smooth)) %>%
    ggplot(aes(x = dateTime, y = val, color = dateTime)) +
    geom_point() +
    scale_color_datetime(low = 'blue', high = 'red') +
    facet_wrap(~var, ncol = 1, scales = "free_y")
  
  #Graph event, CQ
  cq <- event_data %>%
    filter(event == event_select) %>%
    ggplot(aes(x = riverQ_linterp, y = N_linterp_smooth, color = dateTime)) +
    geom_path() +
    scale_color_datetime(low = 'blue', high = 'red')
  
  ggarrange(ts, cq, ncol = 2, legend = 'none')
  ggsave(str_c('./Figures/events/raw/', as.character(i), '.tiff'), device = 'tiff', height = 5, width = 8, units =  'in', compression = 'lzw', dpi = 700)
  
}



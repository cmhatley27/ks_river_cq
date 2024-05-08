library(tidyverse)
library(ggpubr)
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

if(!dir.exists('DataFiles/hydro_data/event_data')) {
  dir.create('DataFiles/hydro_data/event_data')
}

# Load data --------------------------------------------

all_hydro_15min <- with_tz(read_csv('./DataFiles/hydro_data/all_hydro_15min.csv'), 'America/Chicago')
all_hydro_daily <- read_csv('./DataFiles/hydro_data/all_hydro_daily.csv')
events <- with_tz(read_csv('DataFiles/hydro_data/event_delineations/BFLOW_events_bfi_threshold_adjusted.csv'), 'America/Chicago') %>%
  select(event_number, start_dateTime, end_dateTime)
events_compound1 <- with_tz(read_csv('DataFiles/hydro_data/event_delineations/events_compound1.csv'), 'America/Chicago')
events_compound2 <- with_tz(read_csv('DataFiles/hydro_data/event_delineations/events_compound2.csv'), 'America/Chicago')

# Calculate HI and FI -----------------------------------------------------

#Select only events, smooth N, and normalize
event_data <- all_hydro_15min %>%
  filter(!is.na(event)) %>% filter(is.na(riverQ_ingap)) %>%
  #smooth N
  mutate(N_smooth = (lag(N, 4) + lag(N, 3) + lag(N, 2) + lag(N, 1) + N + lead(N, 1) + lead(N, 2) + lead(N, 3) + lead(N, 4))/9,
         N_linterp_smooth = ifelse(is.na(N_ingap), (lag(N_linterp, 4) + lag(N_linterp, 3) + lag(N_linterp, 2) + lag(N_linterp, 1) + 
                                                      N_linterp + lead(N_linterp, 1) + lead(N_linterp, 2) + lead(N_linterp, 3) + lead(N_linterp, 4))/9, NA)) %>%
  #normalize Q and N for each event
  group_by(event) %>%
  mutate(Q_norm = (riverQ_linterp - min(riverQ_linterp))/(max(riverQ_linterp) - min(riverQ_linterp)),
         N_norm = ifelse(is.na(N_ingap), (N_linterp - min(N_linterp))/(max(N_linterp) - min(N_linterp)), NA),
         N_smooth_norm = ifelse(is.na(N_ingap), (N_linterp_smooth - min(N_linterp_smooth, na.rm = TRUE))/(max(N_linterp_smooth, na.rm = TRUE) - min(N_linterp_smooth, na.rm = TRUE)), NA))
write_csv(event_data, './DataFiles/hydro_data/event_data/event_data.csv')

#Initialize q_steps
q_step <- tibble(Q_norm = rep(seq(from = 0, to = 100)/100, nrow(events)),
                 event = rep(1:nrow(events), each = 101))

#Split rising and falling limbs and merge q_steps onto each
rise_limbs <- event_data %>%
  group_by(event) %>%
  #rising limb = event start through max Q
  slice(1:which.max(Q_norm)) %>%
  #merge q_step and arrange by increasing Q_norm
  full_join(q_step)  %>%
  arrange(event, Q_norm) %>%
  select(dateTime, event, riverQ_linterp, Q_norm, N, N_smooth_norm) %>%
  #Interpolate between N_norms to fill in all q_steps with a concentration
  mutate(N_smooth_norm = na.approx(N_smooth_norm, x = Q_norm, xout = Q_norm, na.rm = FALSE)) %>%
  #select only q_steps
  filter(Q_norm %in% (seq(0:100)/100)) %>%
  select(event, dateTime, Q_norm, rise_N_norm = N_smooth_norm)
write_csv(rise_limbs, './DataFiles/hydro_data/event_data/rise_limbs.csv')

  
fall_limbs <- event_data %>%
  group_by(event) %>%
  slice(which.max(Q_norm):length(Q_norm)) %>%
  #merge q_step and arrange by increasing Q_norm
  full_join(q_step)  %>%
  arrange(event, Q_norm) %>%
  select(dateTime, event, riverQ_linterp, Q_norm, N, N_smooth_norm) %>%
  #Interpolate between N_norms to fill in all q_steps with a concentration
  mutate(N_smooth_norm = na.approx(N_smooth_norm, x = Q_norm, xout = Q_norm, na.rm = FALSE)) %>%
  #select only q_steps
  filter(Q_norm %in% (seq(0:100)/100)) %>%
  select(event, dateTime, Q_norm, fall_N_norm = N_smooth_norm)
write_csv(fall_limbs, './DataFiles/hydro_data/event_data/fall_limbs.csv')

hyst_merge <- rise_limbs %>%
  left_join(fall_limbs) %>%
  mutate(HI = rise_N_norm - fall_N_norm)

hyst_summary <- hyst_merge %>%
  group_by(event) %>%
  summarise(HI_mean = mean(HI, na.rm = TRUE),
            HI_sd = sqrt(var(HI, na.rm = TRUE)),
            HI_cv = HI_sd/HI_mean)

flush <- rise_limbs %>%
  group_by(event) %>%
  filter(!is.na(rise_N_norm)) %>%
  filter(Q_norm == max(Q_norm) | Q_norm == min(Q_norm)) %>%
  mutate(FI = rise_N_norm[2] - rise_N_norm[1])  %>%
  summarise(FI = mean(FI))

CQ_summary <- hyst_summary %>%
  left_join(flush)
write_csv(CQ_summary, './DataFiles/hydro_data/CQ_summary.csv')

# CQ slopes for rising, falling, and full events --------------------------
rising_slope <- event_data %>%
  group_by(event) %>%
  slice(1:which.max(Q_norm)) %>%
  select(event, dateTime, riverQ_linterp, N_smooth) %>%
  mutate(logQ = log10(riverQ_linterp), logN = log10(N_smooth)) %>%
  filter(!is.na(logN)) %>%
  summarise(rise_slope = lm(logN ~ logQ, na.action = 'na.omit')$coefficients[2])

falling_slope <- event_data %>%
  group_by(event) %>%
  slice(which.max(Q_norm):length(Q_norm)) %>%
  select(event, dateTime, riverQ_linterp, N_smooth) %>%
  mutate(logQ = log10(riverQ_linterp), logN = log10(N_smooth)) %>%
  filter(!is.na(logN)) %>%
  summarise(fall_slope = lm(logN ~ logQ, na.action = 'na.omit')$coefficients[2])

full_slope <- event_data %>%
  group_by(event) %>%
  select(event, dateTime, riverQ_linterp, N_smooth) %>%
  mutate(logQ = log10(riverQ_linterp), logN = log10(N_smooth)) %>%
  filter(!is.na(logN)) %>%
  summarise(full_slope = lm(logN ~ logQ, na.action = 'na.omit')$coefficients[2])

slope_summary <- full_slope %>%
  left_join(rising_slope, by = 'event') %>%
  left_join(falling_slope, by = 'event') %>%
  mutate(full_slope = if_else(is.na(rise_slope) | is.na(fall_slope), NA, full_slope))

#will have to grab event_summary_trim from next code file
#slope_and_predictors <- data.frame(slope_v_FI, event_summary_trim)
  
ggplot(data = slope_and_predictors) +
  geom_point(aes(x = rise_slope, y = fall_slope, color = season)) +
  xlim(c(-4,4)) +
  ylim(c(-4,4)) +
  geom_point(aes(x = 0.54, y = 0.54), shape = 9, size = 3, color = 'red') +
  facet_wrap(vars(season))
write_csv(slope_summary, './DataFiles/hydro_data/CQ_slope_summary.csv')

slope_v_FI <- CQ_summary %>%
  left_join(slope_summary, by = 'event')
ggplot(data = slope_v_FI, aes(x = rise_slope, y = FI)) +
  geom_point()
cor(x = slope_v_FI$rise_slope, y = slope_v_FI$FI, use = 'complete.obs')

# Plotting ----------------------------------------------------------------
CQ_summary <- read_csv('./DataFiles/hydro_data/CQ_summary.csv')
event_data <- with_tz(read_csv('./DataFiles/hydro_data/event_data/event_data.csv'), 'America/Chicago')
#* HI/FI scatter -----------------------------------------------------------
ggplot(data = subset(CQ_summary, event > 61)) +
  geom_point(aes(x = FI, y = HI_mean)) +
  geom_point(aes(x = mean(FI, na.rm = TRUE), y = mean(HI_mean, na.rm = TRUE)), color = 'red', shape = 15, size = 2) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  scale_x_continuous(limits = c(-1,1), name = 'Flushing Index') +
  scale_y_continuous(limits = c(-1,1), name = "Hysteresis Index")

plot_ly(data = CQ_summary, x = ~FI, y = ~HI_mean, text = ~paste(event))



#* Graph single event ------------------------------------------------------
##127 - classic counterclockwise
##167 - counterclockwise with release
##139 - clockwise
##176 - clockwise
##149 - simple diluting, low hysteresis
##150 - simple flushing, low hysteresis
##140 - ???
event_select <- c(115)

#Graph event, time series
ts <- event_data %>%
  filter(event == event_select) %>%
  gather(key = "var", value = "val", c(Flow = riverQ, Nitrate = N_linterp_smooth, Reservoir = reservoir_sum_linterp)) %>%
  ggplot(aes(x = dateTime, y = val, color = dateTime)) +
  geom_point() +
  scale_color_datetime(low = 'blue', high = 'red') +
  
  facet_wrap(~var, ncol = 1, scales = "free_y")
ts
#Graph event, CQ
cq <- event_data %>%
  filter(event == event_select) %>%
  ggplot(aes(x = riverQ_linterp, y = N_linterp_smooth, color = dateTime)) +
  geom_path() +
  xlab('River Discharge') +
  ylab('Nitrate Conc') +
  scale_color_datetime(low = 'blue', high = 'red', guide = 'none') +
  ggtitle(paste0('Event #',event_select))
cq
ggsave(str_c('./Figures/events/raw/cqloop/', as.character(event_select), '.tiff'), width = 6, height = 6, units = 'in', compression = 'lzw', dpi = 700)
#Show time series + CQ
ggarrange(ts, cq, ncol = 2, legend = 'none')


#* Loop and save all events ------------------------------------------------
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
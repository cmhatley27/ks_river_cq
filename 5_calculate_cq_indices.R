library(tidyverse)

library(zoo)
library(magrittr)
library(quantreg)

# Load data --------------------------------------------

all_hydro_15min <- with_tz(read_csv('./DataFiles/hydro_data/all_hydro_15min.csv'), 'America/Chicago')
events <- with_tz(read_csv('DataFiles/hydro_data/event_delineations/adjusted_events.csv'), 'America/Chicago') %>%
  select(event_number, start_dateTime, end_dateTime)

# Calculate HI and FI -----------------------------------------------------

#Select only events, smooth N, and normalize
event_data <- all_hydro_15min %>%
  filter(!is.na(event)) %>% filter(riverQ_ingap == 0) %>%
  #smooth N
  mutate(N_smooth = (lag(N, 4) + lag(N, 3) + lag(N, 2) + lag(N, 1) + N + lead(N, 1) + lead(N, 2) + lead(N, 3) + lead(N, 4))/9,
         N_linterp_smooth = ifelse(N_ingap == 0, (lag(N_linterp, 4) + lag(N_linterp, 3) + lag(N_linterp, 2) + lag(N_linterp, 1) + 
                                                      N_linterp + lead(N_linterp, 1) + lead(N_linterp, 2) + lead(N_linterp, 3) + lead(N_linterp, 4))/9, NA)) %>%
  #normalize Q and N for each event
  group_by(event) %>%
  mutate(Q_norm = (riverQ_linterp - min(riverQ_linterp))/(max(riverQ_linterp) - min(riverQ_linterp)),
         N_norm = ifelse(N_ingap == 0, (N_linterp - min(N_linterp))/(max(N_linterp) - min(N_linterp)), NA),
         N_smooth_norm = ifelse(N_ingap == 0, (N_linterp_smooth - min(N_linterp_smooth, na.rm = TRUE))/(max(N_linterp_smooth, na.rm = TRUE) - min(N_linterp_smooth, na.rm = TRUE)), NA))
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
  left_join(flush) %>%
  select(event, FI_n = FI, HI_n = HI_mean)
write_csv(CQ_summary, './DataFiles/CQ_summary.csv')
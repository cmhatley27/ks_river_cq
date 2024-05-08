library(tidyverse)
library(dataRetrieval)
library(zoo)


# Download data from NWIS -------------------------------------------------

#we want params:
#00060 - Q
#99133 - N
#00095 - Cond
#63680 - Turb

whatNWISsites(stateCd = 'KS', parameterCd = '63680')
dat <- readNWISuv("06892350", c("00060", "99133", "00095", "63680"), "2013-10-01", "2022-09-30", tz = 'America/Chicago')

dateTimes <- data.frame(dateTime = seq.POSIXt(as.POSIXlt('2013-10-01 00:00:00'), as.POSIXlt('2022-09-30 23:45:00'), by = '15 min'))

gage <- left_join(dateTimes, dat) %>%
  select(dateTime,
         q = X_00060_00000,
         n = X_99133_00000,
         c = X_.YSI.EXO._00095_00000,
         t = X_.YSI.EXO._63680_00000) %>%
  mutate(q = q/35.314666212661)

write_csv(gage, file.path('DataFiles','hydro_data','other_params','des_15min.csv'))

# gage_daily <- gage %>%
#   group_by(date = date(gage$dateTime)) %>%
#   summarise(across(c(q,n,c,t),mean))
# 
# ggplot(data = gage_daily) +
#   geom_point(aes(x = q, y = c))

# Load data --------------------------------------------

all_hydro_15min <- with_tz(read_csv('./DataFiles/hydro_data/all_hydro_15min.csv'), 'America/Chicago')
des_15min <- with_tz(read_csv('./DataFiles/hydro_data/other_params/des_15min.csv'), 'America/Chicago')
events <- with_tz(read_csv('DataFiles/hydro_data/event_delineations/BFLOW_events_bfi_threshold_adjusted.csv'), 'America/Chicago') %>%
  select(event_number, start_dateTime, end_dateTime)


# Assign events -----------------------------------------
event_list <- list()
for(i in seq_along(events$start_dateTime)){
event_list[[i]] <- interval(events$start_dateTime[i], events$end_dateTime[i])
}

event_assignment <- function(x, list){
  in_event <- x %within% list
  for(i in seq_along(list)){
    in_event <- ifelse(x %within% list[[i]], i, in_event)
  }
  return(in_event)
}

des_15min <- mutate(des_15min, 
                    in_event = dateTime %within% event_list,
                    event = event_assignment(dateTime, event_list))


# Flag NAs ----------------------------------------------------------------
library(FluMoDL)

na_length <- function(x, threshold){
  rle_obj <- rle(is.na(x))
  vec <- rep(rle_obj$lengths*rle_obj$values, times = rle_obj$lengths)
  return(ifelse(vec > threshold, TRUE, FALSE))
}

fillterp <- function(x){
  if(is.na(x[1])){
    x[1] <- unique(x)[2]
  }
  if(is.na(x[length(x)])){
    x[length(x)] <- unique(x)[length(unique(x))]
  }
  return(linterp(x, max_allow = NULL))
}

gap_threshold <- 12 #3 hours of NA

des_15min <- mutate(des_15min,
                    across(c(q,n,c,t), 
                           list(gap = ~na_length(.x,gap_threshold), linterp = ~fillterp(.x)),
                           .names = '{.col}_{.fn}'))

# Prep for index calculation -----------------------------------------------------
#Smooth out time series of constituents
smooth_9 <- function(x){
  smoothed <- (lag(x, 4) + lag(x, 3) + lag(x, 2) + lag(x, 1) + x + lead(x, 1) + lead(x, 2) + lead(x, 3) + lead(x, 4))/9
  return(smoothed)
}

#Normalize event values to 0:1
normalizer <- function(x){
  return((x - min(x))/(max(x) - min(x)))
}

#Select only events, smooth conc, and normalize
event_data <- des_15min %>%
  filter(event != 0) %>% filter(!q_gap) %>%
  #smooth
  mutate(across(c(n_linterp, c_linterp, t_linterp),
                ~smooth_9(.x),
                .names = '{.col}_smooth')) %>%
  #normalize by event
  group_by(event) %>%
  mutate(across(c(q_linterp, ends_with('_smooth')),
                ~normalizer(.x),
                .names = '{.col}_norm')) %>%
  #Replace gaps with NA
  mutate(n_linterp_smooth_norm = ifelse(n_gap, NA, n_linterp_smooth_norm),
         c_linterp_smooth_norm = ifelse(c_gap, NA, c_linterp_smooth_norm),
         t_linterp_smooth_norm = ifelse(t_gap, NA, t_linterp_smooth_norm))
write_csv(event_data, './DataFiles/hydro_data/other_params/event_data.csv')

#event viewer
plot_dat <- pivot_longer(event_data, c(n_linterp_smooth_norm, c_linterp_smooth_norm, t_linterp_smooth_norm),
                         names_to = 'var', values_to = 'val')
ggplot(data = subset(plot_dat, event == 182), aes(x = q_linterp_norm, y = val, color = dateTime)) +
  geom_path() +
  scale_color_gradient(low = 'blue', high = 'red') +
  facet_wrap(vars(var))
# Calculate FI and HI -----------------------------------------------------

#Initialize q_steps
q_step <- tibble(q_linterp_norm = rep(seq(from = 0, to = 100)/100, nrow(events)),
                 event = rep(1:nrow(events), each = 101))



toss_gaps <- function(c, q){
  if(sum(is.na(c)) > 0.95*length(c)){
    out <- NA
  }
  else{
    out <- na.approx(c, x = q, xout = q, na.rm = FALSE)
  }
  return(out)
}

#Split rising and falling limbs and merge q_steps onto each
rise_limbs <- event_data %>%
  group_by(event) %>%
  #rising limb = event start through max Q
  slice(1:which.max(q_linterp_norm)) %>%
  #merge q_step and arrange by increasing q_linterp_norm
  full_join(q_step)  %>%
  arrange(event, q_linterp_norm) %>%
  select(dateTime, event, q_linterp, q_linterp_norm, n, n_linterp_smooth_norm, c, c_linterp_smooth_norm, t, t_linterp_smooth_norm) %>%
  #Interpolate between N_norms to fill in all q_steps with a concentration
  mutate(n_linterp_smooth_norm = na.approx(n_linterp_smooth_norm, x = q_linterp_norm, xout = q_linterp_norm, na.rm = FALSE),
         c_linterp_smooth_norm = toss_gaps(c_linterp_smooth_norm, q_linterp_norm),
         t_linterp_smooth_norm = toss_gaps(t_linterp_smooth_norm, q_linterp_norm)) %>%
  #select only q_steps
  filter(q_linterp_norm %in% (seq(0:100)/100)) %>%
  select(event, dateTime, q_norm = q_linterp_norm, rise_n_norm = n_linterp_smooth_norm, rise_c_norm = c_linterp_smooth_norm, rise_t_norm = t_linterp_smooth_norm)


fall_limbs <- event_data %>%
  group_by(event) %>%
  slice(which.max(q_linterp_norm):length(q_linterp_norm)) %>%
  #merge q_step and arrange by increasing q_linterp_norm
  full_join(q_step)  %>%
  arrange(event, q_linterp_norm) %>%
  select(dateTime, event, q_linterp, q_linterp_norm, n, n_linterp_smooth_norm, c, c_linterp_smooth_norm, t, t_linterp_smooth_norm) %>%
  #Interpolate between N_norms to fill in all q_steps with a concentration
  mutate(n_linterp_smooth_norm = na.approx(n_linterp_smooth_norm, x = q_linterp_norm, xout = q_linterp_norm, na.rm = FALSE),
         c_linterp_smooth_norm = toss_gaps(c_linterp_smooth_norm, q_linterp_norm),
         t_linterp_smooth_norm = toss_gaps(t_linterp_smooth_norm, q_linterp_norm)) %>%
  #select only q_steps
  filter(q_linterp_norm %in% (seq(0:100)/100)) %>%
  select(event, dateTime, q_norm = q_linterp_norm, fall_n_norm = n_linterp_smooth_norm, fall_c_norm = c_linterp_smooth_norm, fall_t_norm = t_linterp_smooth_norm)


hyst_merge <- rise_limbs %>%
  left_join(fall_limbs) %>%
  mutate(HI_n = rise_n_norm - fall_n_norm,
         HI_c = rise_c_norm - fall_c_norm,
         HI_t = rise_t_norm - fall_t_norm)

hyst_summary <- hyst_merge %>%
  group_by(event) %>%
  summarise(HI_n = mean(HI_n, na.rm = TRUE),
            HI_c = mean(HI_c, na.rm = TRUE),
            HI_t = mean(HI_t, na.rm = TRUE))


flush_n <- rise_limbs %>%
  group_by(event) %>%
  filter(!is.na(rise_n_norm)) %>%
  filter(q_norm == max(q_norm) | q_norm == min(q_norm)) %>%
  mutate(FI = rise_n_norm[2] - rise_n_norm[1])  %>%
  summarise(FI_n = mean(FI))

flush_c <- rise_limbs %>%
  group_by(event) %>%
  filter(!is.na(rise_c_norm)) %>%
  filter(q_norm == max(q_norm) | q_norm == min(q_norm)) %>%
  mutate(FI = rise_c_norm[2] - rise_c_norm[1])  %>%
  summarise(FI_c = mean(FI))

flush_t <- rise_limbs %>%
  group_by(event) %>%
  filter(!is.na(rise_t_norm)) %>%
  filter(q_norm == max(q_norm) | q_norm == min(q_norm)) %>%
  mutate(FI = rise_t_norm[2] - rise_t_norm[1])  %>%
  summarise(FI_t = mean(FI))


CQ_summary <- hyst_summary %>%
  left_join(flush_n) %>%
  left_join(flush_c) %>%
  left_join(flush_t)
write_csv(CQ_summary, './DataFiles/hydro_data/other_params/CQ_summary.csv')


plot_dat <- pivot_longer(CQ_summary,
                         starts_with('HI'), names_to = 'HI_var', values_to = 'HI') %>%
  pivot_longer(., starts_with('FI'), names_to = 'FI_var', values_to = 'FI') %>%
  filter(str_sub(FI_var, -1) == str_sub(HI_var, -1)) %>%
  mutate(var = str_sub(FI_var, -1)) %>%
  select(event, var, FI, HI)

ggplot(data = plot_dat, aes(x = FI, y = HI)) +
  geom_point() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  facet_wrap(vars(var))

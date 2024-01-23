library(tidyverse)
library(reshape)
library(corrplot)




# Load all relevant data --------------------------------------------------

all_hydro_15min <- with_tz(read_csv('./DataFiles/hydro_data/all_hydro_15min.csv'), 'America/Chicago')
all_hydro_daily <- read_csv('./DataFiles/hydro_data/all_hydro_daily.csv')
events <- with_tz(read_csv('DataFiles/hydro_data/event_delineations/BFLOW_events_bfi_threshold_adjusted.csv'), 'America/Chicago') %>%
  select(event_number, start_dateTime, end_dateTime)
CQ_summary <- with_tz(read_csv('./DataFiles/hydro_data/CQ_summary.csv'), 'America/Chicago')
event_data <- with_tz(read_csv('./DataFiles/hydro_data/event_data/event_data.csv'), 'America/Chicago')
rise_limbs <- with_tz(read_csv('./DataFiles/hydro_data/event_data/rise_limbs.csv'), 'America/Chicago')
fall_limbs <- with_tz(read_csv('./DataFiles/hydro_data/event_data/fall_limbs.csv'), 'America/Chicago')


# Calculate Predictor Variables -------------------------------------------

#* Times/seasonality -------------------------------------------------------
event_temporal <- events %>%
  mutate(
    #Event number
    event = row_number(start_dateTime),
    #Water Year
    wateryear = factor(ifelse(as.integer(month(start_dateTime)) >= 10, as.integer(year(start_dateTime)) + 1, year(start_dateTime))),
    #Month
    month = factor(month(start_dateTime, label = TRUE)),
    #Season (Winter = Dec - Feb, Spring = Mar - May, etc.)
    season = factor(month(floor_date(start_dateTime, unit = 'season')), levels = c(12,3,6,9), labels = c('Winter', 'Spring', 'Summer', 'Fall')),
    #Julian day
    yday = yday(start_dateTime),
    #Event duration (in seconds)
    duration = time_length(as.duration(interval(start = start_dateTime, end = end_dateTime))),
    #Duration since previous event (in seconds)
    duration_since_last = time_length(as.duration(interval(start = lag(end_dateTime), end = start_dateTime)))
    ) %>%
  select(-event, -event_number, start_dateTime, end_dateTime, wateryear, month, season, yday, duration, duration_since_last)



#* Event size/load ---------------------------------------------------------
event_size <- event_data %>%
  group_by(event) %>%
  summarize(
    #N Yield: multiply concentration by runoff converted to L/s, integrate over time, normalize by watershed area (gives mg/m2 which is equivalent to kg/km2)
    N_yield = ifelse(any(!is.na(N_ingap)), NA, sum(N_linterp_smooth*(riverQ_linterp*1000)*60*15/(1.54767E11))),
    #Initial discharge
    initial_Q = riverQ_linterp[1],
    #Maximum discharge
    max_Q = max(riverQ_linterp, na.rm = TRUE),
    #Initial N
    initial_N = N_linterp_smooth[1],
    #Maximum N
    max_N = max(N_linterp_smooth, na.rm = TRUE),
    #Change in N from start to max
    delta_N = max_N - initial_N,
    #Change in N as ratio from initial N
    delta_ratio_N = delta_N/initial_N) %>%
  select(-event)


  

#* Precip and reservoir windows --------------------------------------------


#** Precip windows ----------------------------------------------------------

#Function for adding up precip in a window from event start
cumul_precip <- function(date_start, window_length, dates, precip_series){
  df <- data.frame(date = dates, precip = precip_series)
  df <- df %>%
    mutate(cumul_precip = cumsum(precip))
  precip_total <- df$cumul_precip[df$date == date_start] - max(df$cumul_precip[df$date == date_start - window_length - 1], 0)
  return(ifelse(date_start %within% interval(dates[1], dates[1] + window_length), NA, precip_total))
}

#**Loops ante-event precip function over every gage+window length combo
before_event_precip <- function(precip_frame, window_lengths, event_frame = event_data){
  #extract dates
  dates <- precip_frame$date
  #extract just precip series
  precips <- precip_frame %>%
    select(ends_with('precip'))
  #initialize list of output vectors
  window_totals <- list()
  #Loop through each combo of gage+window length and calculate cumulative precip
  for(gage in 1:ncol(precips)){
    for(window in 1:length(window_lengths)){
      summary_frame <- event_frame %>%
        group_by(event) %>%
        summarize(date_start = date(dateTime[1]),
                  date_max = date(dateTime[which.max(riverQ_linterp)]),
                  precip_in_window = cumul_precip(date_start = date_start, window_length = window_lengths[window], dates = precip_frame$date, precip_series = precips[[gage]])
        )
      window_totals[[paste('total_prior_', colnames(precips[gage]), '_', window_lengths[window], 'day_window', sep = "")]] <- summary_frame$precip_in_window
    }
  }
  return(window_totals)
}

precip_windows_ante <- before_event_precip(precip_frame = all_hydro_daily, window_lengths = c(30,60,365))

#Function for adding up precip in a window from event start minus a buffer to event max
event_precip <- function(date_start, date_max, buffer_length, dates, precip_series){
  df <- data.frame(date = dates, precip = precip_series)
  df <- df %>%
    mutate(cumul_precip = cumsum(precip))
  precip_total <- df$cumul_precip[df$date == date_max] - max(df$cumul_precip[df$date == date_start - buffer_length - 1], 0)
  return(ifelse(date_start %within% interval(dates[1], dates[1] + buffer_length), NA, precip_total))
}

#**Loops within-event precip function over every gage+buffer length combo
within_event_precip <- function(precip_frame, buffer_lengths, event_frame = event_data){
  #extract dates
  dates <- precip_frame$date
  #extract just precip series
  precips <- precip_frame %>%
    select(ends_with('precip'))
  #initialize list of output vectors
  window_totals <- list()
  #Loop through each combo of gage+buffer length and calculate cumulative precip
  for(gage in 1:ncol(precips)){
    for(buffer in 1:length(buffer_lengths)){
      summary_frame <- event_frame %>%
        group_by(event) %>%
        summarize(date_start = date(dateTime[1]),
                  date_max = date(dateTime[which.max(riverQ_linterp)]),
                  precip_in_window = event_precip(date_start = date_start, date_max = date_max, buffer_length = buffer_lengths[buffer], dates = precip_frame$date, precip_series = precips[[gage]])
        )
      window_totals[[paste('total_event_', colnames(precips[gage]), '_', buffer_lengths[buffer], 'day_buffer', sep = "")]] <- summary_frame$precip_in_window
    }
  }
  return(window_totals)
}

precip_windows_within <- within_event_precip(precip_frame = all_hydro_daily, buffer_lengths = c(0,1,2))

#Function for counting number of wet days in a window from event start
count_wet <- function(date_start, window_length, dates, precip_series){
  df <- data.frame(date = dates, precip = precip_series)
  df <- df %>%
    filter(date %within% interval(date_start - window_length, date_start)) %>%
    filter(precip != 0)
  return(ifelse(date_start %within% interval(dates[1], dates[1] + window_length), NA, nrow(df)))
}

#**Loops wet days function over every gage+window length combo
count_wet_days <- function(precip_frame, window_lengths, event_frame = event_data){
  #extract dates
  dates <- precip_frame$date
  #extract just precip series
  precips <- precip_frame %>%
    select(ends_with('precip'))
  #initialize list of output vectors
  window_totals <- list()
  #Loop through each combo of gage+window length and calculate cumulative precip
  for(gage in 1:ncol(precips)){
    for(window in 1:length(window_lengths)){
      summary_frame <- event_frame %>%
        group_by(event) %>%
        summarize(date_start = date(dateTime[1]),
                  date_max = date(dateTime[which.max(riverQ_linterp)]),
                  precip_count = count_wet(date_start = date_start, window_length = window_lengths[window], dates = precip_frame$date, precip_series = precips[[gage]])
        )
      window_totals[[paste('number_wet_', colnames(precips[gage]), '_', window_lengths[window], 'day_window', sep = "")]] <- summary_frame$precip_count
    }
  }
  return(window_totals)
}

precip_wet_days <- count_wet_days(precip_frame = all_hydro_daily, window_lengths = c(30,90,365))

#Combine all into dataframe
precip_summary <- data.frame(precip_windows_ante, precip_windows_within, precip_wet_days)
write_csv(precip_summary, './DataFiles/hydro_data/precip_summary.csv')


#** Reservoir windows ----------------------------------------------------------

#Function for adding up 15min flows between event start minus a buffer to event end
event_flow <- function(dateTime_start, dateTime_end, buffer_length, dateTimes, flow_series){
  df <- data.frame(dateTime = dateTimes, cum_flows = cumsum(flow_series))
  flow_total <- df$cum_flows[df$dateTime == dateTime_end] - df$cum_flows[df$dateTime == dateTime_start - days(buffer_length) - minutes(15)]
  return(ifelse(date(dateTime_start) < ymd('2013-10-01'), NA, flow_total))
}

#**Loops within-event flow function over every gage and buffer length combo
within_event_flows <- function(flows_frame, buffer_lengths, event_frame = event_data){
  #extract dates
  dateTimes <- flows_frame$dateTime
  #extract just precip series
  flows <- flows_frame %>%
    select(ends_with('_linterp'))
  #initialize list of output vectors
  window_totals <- list()
  #Loop through each combo of gage+buffer length and calculate cumulative precip
  for(gage in 1:ncol(flows)){
    for(buffer in 1:length(buffer_lengths)){
      summary_frame <- event_frame %>%
        group_by(event) %>%
        summarize(dateTime_start = dateTime[1],
                  dateTime_end = max(dateTime),
                  flow_in_window = event_flow(dateTime_start = dateTime_start, dateTime_end = dateTime_end, buffer_length = buffer_lengths[buffer], dateTimes = flows_frame$dateTime, flow_series = flows[[gage]])
        )
      window_totals[[paste('total_event_flow_', word(colnames(flows[gage]),1, sep = "_"), '_', buffer_lengths[buffer], 'day_buffer', sep = "")]] <- summary_frame$flow_in_window
    }
  }
  return(window_totals)
}

flow_windows_within <- within_event_flows(flows_frame = all_hydro_15min, buffer_lengths = c(0,1,2))

#Function for calculting max - initial flow for each gage
delta_flows <- function(flows_frame = all_hydro_15min, buffer_lengths){
  #extract just flow series
  flows <- flows_frame %>%
    select(dateTime, ends_with('_linterp'))
  max_flows <- flows_frame %>%
    group_by(event) %>%
    summarize(across(ends_with('_linterp'), max)) %>%
    filter(!is.na(event))
  deltas <- list()
  ratios <- list()
  for(buffer in 1:length(buffer_lengths)){
    start_dateTimes <- flows_frame %>%
      group_by(event) %>%
      summarize(dateTime = dateTime[1] - days(buffer_lengths[buffer]) + minutes(90)) %>%
      filter(!is.na(event))
    start_flows <- start_dateTimes %>%
      left_join(flows) %>%
      select(-dateTime) %>%
      filter(!is.na(event))
    delta_frame <- max_flows - start_flows
    colnames(delta_frame) <- paste('delta_', word(colnames(delta_frame),1,sep = '_'), '_', buffer_lengths[buffer], 'day_buffer', sep = '')
    ratio_frame <- delta_frame / start_flows
    colnames(ratio_frame) <- paste('delta_ratio_', word(colnames(start_flows),1,sep = '_'), '_', buffer_lengths[buffer], 'day_buffer', sep = '')
    deltas[[buffer]] <- delta_frame
    ratios[[buffer]] <- ratio_frame
  }
  results <- list(deltas, ratios)
  return(results)
}

flow_deltas <- data.frame(delta_flows(buffer_lengths = c(0,1,2))) %>%
  select(!contains('event'))

#Combine all into data frame
flow_summary <- data.frame(flow_windows_within, flow_deltas) %>%
  #Convert to mm
  mutate(across(!starts_with('delta_'), ~.x*60*15/(1.54767E11)*1000))

#Split into runoff and reservoir flows for easier calculations
runoff_summary <- flow_summary %>%
  select(contains('riverQ'))
write_csv(runoff_summary, './DataFiles/hydro_data/runoff_summary.csv')

reservoir_summary <- flow_summary %>%
  select(!contains('river') & !contains('lawrence') & !contains('upstream') &!contains('_N', ignore.case = FALSE)) 
write_csv(reservoir_summary, './DataFiles/hydro_data/reservoir_summary.csv')

#Synthesize predictor sets --------------------------------------------

#Join all predictors and save
precip_summary <- read_csv('./DataFiles/hydro_data/precip_summary.csv')
runoff_summary <- read_csv('./DataFiles/hydro_data/runoff_summary.csv')
reservoir_summary <- read_csv('./DataFiles/hydro_data/reservoir_summary.csv')
event_summary <- data.frame(CQ_summary,
                            event_temporal,
                            event_size,
                            precip_summary,
                            runoff_summary,
                            reservoir_summary) %>%
  mutate(runoff_precip_ratio = total_event_flow_riverQ_0day_buffer/total_event_lawrence_precip_2day_buffer,
         runoff_precip_ratio_norm = ifelse(is.infinite(runoff_precip_ratio), 1, runoff_precip_ratio/max(subset(runoff_precip_ratio, !is.infinite(runoff_precip_ratio)))),
         reservoir_precip_ratio = total_event_flow_reservoir_2day_buffer/total_event_lawrence_precip_2day_buffer,
         reservoir_precip_ratio_norm = ifelse(is.infinite(reservoir_precip_ratio), 1, reservoir_precip_ratio/max(subset(reservoir_precip_ratio, !is.infinite(reservoir_precip_ratio)))),
         reservoir_runoff_ratio = total_event_flow_reservoir_0day_buffer/total_event_flow_riverQ_0day_buffer) %>%
  mutate(runoff_precip_ratio_norm = log10(runoff_precip_ratio_norm),
         reservoir_precip_ratio_norm = log10(reservoir_precip_ratio_norm))

write_csv(event_summary, './DataFiles/hydro_data/event_summary.csv')
event_summary <- with_tz(read_csv('./DataFiles/hydro_data/event_summary.csv'), 'America/Chicago')


# Extraneous/irrelevant predictors to drop: -------------------------------
extraneous_predictors <- event_summary %>%
  select(c(
event,
start_dateTime,
end_dateTime,
HI_sd,
HI_cv,
N_yield,
max_N,
delta_N,
delta_ratio_N,
total_event_lawrence_precip_1day_buffer,
total_event_lawrence_precip_2day_buffer,
total_event_manhattan_precip_0day_buffer,
total_event_manhattan_precip_1day_buffer,
total_event_salina_precip_0day_buffer,
total_event_salina_precip_1day_buffer,
total_event_flow_riverQ_1day_buffer,
total_event_flow_riverQ_2day_buffer,
delta_riverQ_1day_buffer,
delta_riverQ_2day_buffer,
delta_ratio_riverQ_1day_buffer,
delta_ratio_riverQ_2day_buffer,
total_event_flow_clintonQ_0day_buffer,
total_event_flow_clintonQ_2day_buffer,
total_event_flow_milfordQ_0day_buffer,
total_event_flow_milfordQ_1day_buffer,
total_event_flow_perryQ_0day_buffer,
total_event_flow_perryQ_2day_buffer,
total_event_flow_tuttleQ_0day_buffer,
total_event_flow_tuttleQ_1day_buffer,
total_event_flow_reservoir_0day_buffer,
total_event_flow_reservoir_1day_buffer,
delta_clintonQ_0day_buffer,
delta_clintonQ_2day_buffer,
delta_milfordQ_0day_buffer,
delta_milfordQ_1day_buffer,
delta_perryQ_0day_buffer,
delta_perryQ_2day_buffer,
delta_tuttleQ_0day_buffer,
delta_tuttleQ_1day_buffer,
delta_reservoir_0day_buffer,
delta_reservoir_1day_buffer,
delta_ratio_clintonQ_0day_buffer,
delta_ratio_clintonQ_2day_buffer,
delta_ratio_milfordQ_0day_buffer,
delta_ratio_milfordQ_1day_buffer,
delta_ratio_perryQ_0day_buffer,
delta_ratio_perryQ_2day_buffer,
delta_ratio_tuttleQ_0day_buffer,
delta_ratio_tuttleQ_1day_buffer,
delta_ratio_reservoir_0day_buffer,
delta_ratio_reservoir_1day_buffer,
runoff_precip_ratio,
reservoir_precip_ratio,
contains('kanopolis'), 
contains('waconda'), 
contains('wilson'), 
contains('belvue'), 
contains('topeka')))

event_summary_trim <- event_summary %>%
  select(!c(names(extraneous_predictors)))


# Update names and reorder ------------------------------------------------
names_update <- read_csv('./DataFiles/hydro_data/var_names.csv')

colnames(event_summary_trim) <- names_update$new_name

event_summary_trim <- event_summary_trim %>%
  select(c(
    FI, HI, wateryear, month, season, yday,
    initial_Q, initial_N, duration_since_last,
    max_Q, duration, total_Q, delta_Q, delta_Q_ratio,
    starts_with('precip_L'), starts_with('precip_M'), starts_with('precip_S'),
    starts_with('res_C'), starts_with('res_P'), starts_with('res_T'), starts_with('res_M'), starts_with('res_sum'),
    runoff_precip_ratio, reservoir_precip_ratio, reservoir_runoff_ratio
  ))

# Check correlations ------------------------------------------------------

#* High correlations between predictors ---------------------------------------------
r_cutoff <- 0.85
high_corr_trim <- event_summary_trim %>%
  select(where(is.numeric)) %>%
  cor(., method = 'pearson', use = 'pairwise.complete.obs') %>%
  melt(.) %>%
  filter(abs(value) >= r_cutoff) %>%
  filter(X1 != X2)

corr_trim <- event_summary_trim %>%
  select(where(is.numeric)) %>%
  cor(., method = 'pearson', use = 'pairwise.complete.obs')

tiff(filename = './Figures/corrplot_trimmed.tiff', width = 10, height = 10, units = 'in', res = 800, compression = 'lzw')
corrplot(corr_trim, method = 'color', tl.cex = 0.5, tl.col = 'black', type = 'lower', diag = FALSE,
         col = COL2('RdBu', n = 13), col.lim = c(min(corr_trim), 1))
dev.off()



#Correlated columns to drop:
#Salina precip, highly correlated with manhattan
#Reservoir sums, highly correlated with tuttle or milford
event_summary_trim <- event_summary_trim %>%
  select(!c(
    max_Q,
    runoff_precip_ratio,
    contains('_S_'),
    res_sum_total,
    res_sum_delta,
    res_sum_delta_ratio))

write_csv(event_summary_trim, './DataFiles/hydro_data/event_summary_trim.csv')
event_summary_trim <- read_csv('./DataFiles/hydro_data/event_summary_trim.csv')

corr_final <- event_summary_trim %>%
  select(where(is.numeric)) %>%
  cor(., method = 'pearson', use = 'pairwise.complete.obs')

tiff(filename = './Figures/corrplot_final.tiff', width = 10, height = 10, units = 'in', res = 800, compression = 'lzw')
corrplot(corr_final, method = 'color', tl.cex = 0.6, tl.col = 'black', type = 'lower', diag = FALSE,
         col = COL2('RdBu', n = 20), col.lim = c(-1, 1))
dev.off()


# Correlations with CQ metrics --------------------------------------------

CQ_metric_corr <- melt(corr_final) %>%
  filter(X1 == 'FI' | X1 == 'HI') %>%
  arrange(X2)
write_csv(CQ_metric_corr, './DataFiles/hydro_data/CQ_metric_corr.csv')

cor.sig <- 0.1386
ggplot(data = subset(CQ_metric_corr, X1 == 'FI' & abs(value) >= 0), aes(x = X2, y = value)) +
  geom_col() +
  geom_hline(yintercept = cor.sig, color = 'blue', linetype = 'dashed') +
  geom_hline(yintercept = -cor.sig, color = 'blue', linetype = 'dashed') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab('') +
  ylab('r')
ggplot(data = subset(CQ_metric_corr, X1 == 'HI' & abs(value) >= 0), aes(x = X2, y = value)) +
  geom_col() +
  geom_hline(yintercept = cor.sig, color = 'blue', linetype = 'dashed') +
  geom_hline(yintercept = -cor.sig, color = 'blue', linetype = 'dashed') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab('') +
  ylab('r')


# Correlations between water yield and N yield -------------------------------------------------------------------------
yield_slopes <- event_summary %>%
  select(event, wateryear, month, season, Q_yield = total_event_flow_riverQ_0day_buffer, N_yield) %>%
  filter(!is.na(N_yield)) %>%
  group_by(season) %>%
  do(slope = as.numeric(lm(N_yield ~ Q_yield, data = .)$coefficients[2]))

yield_correlations <- event_summary %>%
  select(event, wateryear, month, season, Q_yield = total_event_flow_riverQ_0day_buffer, N_yield) %>%
  filter(!is.na(N_yield)) %>%
  group_by(season) %>%
  summarise(cor = cor(Q_yield, N_yield))

yield_stats <- yield_slopes %>% left_join(yield_correlations)

ggplot(data = event_summary, aes(x = total_event_flow_riverQ_0day_buffer, y = N_yield, color = season)) +
  geom_point() +
  scale_color_manual(values = c('#8DA9C4','#60992D','#C33C54','#F4B860'), name = "Season") +
  geom_smooth(se = FALSE, method = 'lm', color = 'black', linetype = 'dashed') +
  geom_text(data = yield_stats, aes(x = 13, y = 15, label = paste0('slope: ', round(as.numeric(slope),3))), color = 'black') +
  geom_text(data = yield_stats, aes(x = 13, y = 12, label = paste0('r2: ', round(cor, 3))), color = 'black') +
  xlab('Event runoff (m3/km2)') +
  ylab('Event N Yield (kg/km2)') +
  facet_wrap(vars(season))
round((yield_stats$slope), 2)

# Plotting ----------------------------------------------------------------

#Load Data
event_summary <- with_tz(read_csv('./DataFiles/hydro_data/event_summary.csv'), 'America/Chicago') %>%
  mutate(wateryear = factor(ifelse(as.integer(month(start_dateTime)) >= 10, as.integer(year(start_dateTime)) + 1, year(start_dateTime))),
         month = factor(month(start_dateTime, label = TRUE)),
         season = factor(month(floor_date(start_dateTime, unit = 'season')), levels = c(12,3,6,9), labels = c('Winter', 'Spring', 'Summer', 'Fall')))


# * Scatters --------------------------------------------------------------

#Single predictor vs single CQ metric, by month
ggplot(data = event_summary, aes(x = total_event_flow_riverQ_0day_buffer, y = N_yield, color = season)) +
  geom_point() +
  scale_color_manual(values = c('#8DA9C4','#60992D','#C33C54','#F4B860'), name = "Season") +
  facet_wrap(vars(month), nrow = 4)
cor(x = event_summary$total_event_flow_riverQ_0day_buffer, y = event_summary$N_yield, use = 'complete.obs')

###Monthly averages
monthly_averages <- event_summary %>%
  group_by(month) %>%
  summarise(mean_FI = mean(FI, na.rm = TRUE),
            sd_FI = sd(FI, na.rm = TRUE),
            mean_HI = mean(HI_mean, na.rm = TRUE),
            sd_HI = sd(HI_mean, na.rm = TRUE),
            season = season)
###Seasonal averages
seasonal_averages <- event_summary %>%
  group_by(season) %>%
  summarise(mean_FI = mean(FI, na.rm = TRUE),
            sd_FI = sd(FI, na.rm = TRUE),
            mean_HI = mean(HI_mean, na.rm = TRUE),
            sd_HI = sd(HI_mean, na.rm = TRUE))
###Yearly averages
yearly_averages <- event_summary %>%
  group_by(wateryear) %>%
  summarise(mean_FI = mean(FI, na.rm = TRUE),
            sd_FI = sd(FI, na.rm = TRUE),
            mean_HI = mean(HI_mean, na.rm = TRUE),
            sd_HI = sd(HI_mean, na.rm = TRUE))

#HI/FI scatter by month
ggplot(data = event_summary) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(aes(x = FI, y = HI_mean, color = season)) +
  geom_point(data = monthly_averages, aes(x = mean_FI, y = mean_HI), color = 'red', shape = 15, size = 2) +
  scale_color_manual(values = c('#8DA9C4','#60992D','#C33C54','#F4B860'), name = "Season") +
  scale_x_continuous(limits = c(-1,1)) +
  scale_y_continuous(limits = c(-1,1)) +
  facet_wrap(vars(month), nrow = 4)

#HI/FI scatter by season
ggplot(data = event_summary) +
  geom_point(aes(x = FI, y = HI_mean, color = season)) +
  geom_errorbar(data = seasonal_averages, aes(y = mean_HI, xmin = mean_FI - sd_FI, xmax = mean_FI + sd_FI), alpha = 0.3) +
  geom_errorbar(data = seasonal_averages, aes(x = mean_FI, ymin = mean_HI - sd_HI, ymax = mean_HI + sd_HI), alpha = 0.3) +
  geom_point(data = seasonal_averages, aes(x = mean_FI, y = mean_HI), color = 'red', shape = 15, size = 2) +
  scale_color_manual(values = c('#8DA9C4','#60992D','#C33C54','#F4B860'), name = "Season") +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  scale_x_continuous(limits = c(-1,1)) +
  scale_y_continuous(limits = c(-1,1)) +
  facet_wrap(vars(season))

#HI/FI scatter by year
ggplot(data = event_summary) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(aes(x = FI, y = HI_mean, color = season)) +
  scale_color_manual(values = c('#8DA9C4','#60992D','#C33C54','#F4B860'), name = "Season") +
  scale_x_continuous(limits = c(-1,1)) +
  scale_y_continuous(limits = c(-1,1)) +
  facet_wrap(vars(wateryear))

#HI/FI scatter by total runoff
ggplot(data = event_summary) +
  geom_point(aes(x = FI, y = HI_mean, color = log(runoff)), alpha = 1) +
  scale_color_binned(type = 'viridis') +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  scale_x_continuous(limits = c(-1,1)) +
  scale_y_continuous(limits = c(-1,1)) +
  facet_wrap(vars(month), nrow = 4)

#HI/FI scatter by initial discharge
ggplot(data = subset(event_summary, !is.nan(HI_mean))) +
  geom_point(aes(x = FI, y = HI_mean, color = log(initial_Q)), alpha = 1) +
  scale_color_gradient2(low = 'blue', mid = 'orange', high = 'red', midpoint = 5.5) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  scale_x_continuous(limits = c(-1,1)) +
  scale_y_continuous(limits = c(-1,1))

#HI/FI scatter by reservoir ratio
ggplot(data = event_summary) +
  geom_point(aes(x = FI, y = HI_mean, color = reservoir_ratio), alpha = 1) +
  scale_color_gradient2(low = 'blue', mid = 'orange', high = 'red', midpoint = 0.8) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  scale_x_continuous(limits = c(-1,1)) +
  scale_y_continuous(limits = c(-1,1))


# * Histograms ------------------------------------------------------------
#HI histogram
ggplot(data = event_summary) +
  geom_histogram(aes(x = HI_mean), binwidth = 0.05) +
  scale_x_continuous(limits = c(-1,1))

#FI histogram
ggplot(data = event_summary) +
  geom_histogram(aes(x = FI), bins = 20) +
  scale_x_continuous(limits = c(-1,1))


# * Box Plots -------------------------------------------------------------
ggplot(data = event_summary, aes(x = FI, y = 1)) +
  geom_boxplot(width = 0.1) +
  geom_jitter(height = 0.01) +
  theme_classic() +
  xlim(c(-1,1)) +
  ylab('')


ggplot(data = event_summary, aes(x = FI, y = factor(season, levels = c('Fall','Summer','Spring','Winter')), fill = season)) +
  geom_vline(xintercept = c(-1,-.5,0,.5,1), color = 'grey50', linetype = 'dashed') +
  geom_boxplot(width = 0.5) +
  scale_fill_manual(values = c('#8DA9C4','#60992D','#C33C54','#F4B860'), name = "Season", guide = 'none') +
  geom_jitter(height = 0.2) +
  theme_classic() +
  xlim(c(-1,1)) +
  ylab('')

# **CQ stats by month and year --------------------------------------------

#FI by month
ggplot(data = event_summary, aes(x = month, y = FI, fill = season)) +
  geom_boxplot(width = 0.5) +
  scale_fill_manual(values = c('#8DA9C4','#60992D','#C33C54','#F4B860'), name = "Season") +
  geom_jitter(width = 0.2) +
  ylim(c(-1,1))
#HI by month
ggplot(data = event_summary, aes(x = month, y = HI_mean, fill = season)) +
  geom_boxplot(width = 0.5) +
  scale_fill_manual(values = c('#8DA9C4','#60992D','#C33C54','#F4B860'), name = "Season") +
  geom_jitter(width = 0.2) +
  ylim(c(-1,1))

#FI by year
ggplot(data = event_summary, aes(x = wateryear, y = FI, fill = wateryear)) +
  geom_boxplot(width = 0.5) +
  geom_jitter(width = 0.2) +
  ylim(c(-1,1))
#HI by year
ggplot(data = event_summary, aes(x = wateryear, y = HI_mean, fill = wateryear)) +
  geom_boxplot(width = 0.5) +
  geom_jitter(width = 0.2) +
  ylim(c(-1,1))

#N yield by month
ggplot(data = event_summary, aes(x = month, y = N_yield, fill = season)) +
  geom_boxplot(width = 0.5) +
  scale_fill_manual(values = c('#8DA9C4','#60992D','#C33C54','#F4B860'), name = "Season") +
  geom_jitter(width = 0.2)
#Runoff by month
ggplot(data = event_summary, aes(x = month, y = runoff, fill = season)) +
  geom_boxplot(width = 0.5) +
  scale_fill_manual(values = c('#8DA9C4','#60992D','#C33C54','#F4B860'), name = "Season") +
  geom_jitter(width = 0.2)
#N/Q by month
ggplot(data = event_summary, aes(x = month, y = N_yield/runoff, fill = season)) +
  geom_boxplot(width = 0.5) +
  scale_fill_manual(values = c('#8DA9C4','#60992D','#C33C54','#F4B860'), name = "Season") +
  geom_jitter(width = 0.2)

# ** Q and N by year and month --------------------------------------------

#Q by month
ggplot(data = all_hydro_daily, aes(x = factor(month(date, label = TRUE)), y = log10(riverQ_linterp), fill = factor(month(floor_date(date, unit = 'season')), levels = c(12,3,6,9), labels = c('Winter', 'Spring', 'Summer', 'Fall')))) +
  geom_boxplot(width = 0.5) +
  scale_fill_manual(values = c('#8DA9C4','#60992D','#C33C54','#F4B860'), name = "Season") +
  geom_jitter(width = 0.2, size = 0.2) +
  xlab('Month')
#N by month
ggplot(data = all_hydro_daily, aes(x = factor(month(date, label = TRUE)), y = N, fill = factor(month(floor_date(date, unit = 'season')), levels = c(12,3,6,9), labels = c('Winter', 'Spring', 'Summer', 'Fall')))) +
  geom_boxplot(width = 0.5) +
  scale_fill_manual(values = c('#8DA9C4','#60992D','#C33C54','#F4B860'), name = "Season") +
  geom_jitter(width = 0.2, size = 0.2) +
  xlab('Month')
#P by month
all_hydro_daily %>%
  group_by(month = factor(month(date, label = TRUE))) %>%
  summarize(P_mean = mean(lawrence_precip)) %>%
ggplot(aes(x = month, y = P_mean)) +
  geom_col(width = 0.5, fill = 'steelblue1', color = 'black') +
  xlab('Month')


#Q by year
ggplot(data = all_hydro_daily, aes(x = factor(wateryear), y = log10(riverQ_linterp), fill = factor(wateryear))) +
  geom_boxplot(width = 0.5) +
  geom_jitter(width = 0.2, size = 0.2) +
  xlab('year') +
  theme(legend.position = 'none')
#N by year
ggplot(data = all_hydro_daily, aes(x = factor(wateryear), y = N, fill = factor(wateryear))) +
  geom_boxplot(width = 0.5) +
  geom_jitter(width = 0.2, size = 0.2) +
  xlab('year') +
  theme(legend.position = 'none')
#P by year
all_hydro_daily %>%
  group_by(wateryear = factor(wateryear)) %>%
  summarize(P_mean = mean(lawrence_precip)) %>%
  ggplot(aes(x = wateryear, y = P_mean)) +
  geom_col(width = 0.5, fill = 'steelblue1', color = 'black') +
  xlab('Year')

# * Autocorrelation -----------------------------------------------------------
Q_acf <- acf(all_hydro_daily$riverQ_linterp, lag.max = 365, na.action = na.pass)
N_acf <- acf(all_hydro_daily$N, lag.max = 365, na.action = na.pass)

QN_ccf <- ccf(all_hydro_15min$riverQ_linterp, all_hydro_15min$perryQ_linterp, lag.max = 672, na.action = na.pass)
(which.max(QN_ccf$acf) - 672) / 96



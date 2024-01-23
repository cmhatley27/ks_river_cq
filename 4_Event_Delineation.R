library(tidyverse)

library(dataRetrieval)
library(hydroEvents)
library(rnoaa)
library(ggpubr)
library(FluMoDL)
library(zoo)
library(magrittr)
library(quantreg)
library(plotly)


# Create folders and load data --------------------------------------------

if(!dir.exists('DataFiles/hydro_data/event_delineations')) {
  dir.create('DataFiles/hydro_data/event_delineations')
}
if(!dir.exists('Figures/hydrographs')) {
  dir.create('Figures/hydrographs')
}

all_hydro_15min <- read_csv('./DataFiles/hydro_data/all_hydro_15min.csv')
all_hydro_daily <- read_csv('./DataFiles/hydro_data/all_hydro_daily.csv')


# Separate Baseflow using Lyne-Hollick Filter ------------------------------------------------------

#Lyne-Holick filter function
baseflow_BFLOW <- function(Q, beta=0.925, passes=3){
  # R implementation of BFLOW baseflow separation algorithm as
  # described in Arnold & Allen (1999). This is the same as the 
  # original digital filter proposed by Lyne & Holick (1979) and
  # tested in Nathan & McMahon (1990). 
  #
  # It is called BFLOW because of this website: 
  #   http://www.envsys.co.kr/~swatbflow/USGS_GOOGLE/display_GoogleMap_for_SWAT_BFlow.cgi?state_name=indiana
  #
  # This is effectively the same as the 'BaseflowSeparation' function 
  # in the EcoHydRology package but with slightly different handling of 
  # start/end dates.
  #
  # Inputs:
  #   Q = discharge timeseries (no missing data) (any units are OK)
  #   beta = filter parameter; recommended value 0.925 (Nathan & McMahon, 1990); 0.9-0.95 reasonable range
  #   passes = how many times to go through the data (3=default=forward/backward/forward)
  #       
  # Output:
  #   bf = baseflow timeseries, same length and units as Q
  
  ## package dependencies
  require(zoo)
  require(dplyr)
  require(magrittr)
  
  # Q for use in calculations
  bfP <- Q
  
  for (p in 1:passes){
    # figure out start and end
    if ((p %% 2)==0){
      # backward pass
      i.start <- length(Q)-1
      i.end   <- 1
      i.fill  <- length(Q)
      ts      <- -1
    } else {
      # forward pass
      i.start <- 2
      i.end   <- length(Q)
      i.fill  <- 1
      ts      <- 1
    }
    
    # make empty vector
    qf <- rep(NaN, length=length(Q))
    
    # fill in value for timestep that will be ignored by filter
    if (p==1){
      qf[i.fill] <- bfP[1]*0.5
    } else {
      qf[i.fill] <- max(c(0, (Q[i.fill]-bfP[i.fill])))
    }
    
    # go through rest of timeseries
    for (i in i.start:i.end){
      qf[i] <- 
        (beta*qf[i-ts] + ((1+beta)/2)*(bfP[i]-bfP[i-ts]))
      
      # check to make sure not too high/low
      if (qf[i] > bfP[i]) qf[i] <- bfP[i]
      if (qf[i] < 0) qf[i] <- 0
    }
    
    # calculate bf for this pass
    bfP <- bfP-qf
    
    # when p==passes, return bfP
    if (p==passes){
      bf <- bfP
    }
    
  } # end of passes loop
  
  return(bf)
}

#Filter parameters
alpha <- 0.992
passes <- 12

#Run filter
baseflow_15min <- data.frame(
  dateTime = all_hydro_15min$dateTime,
  riverQ_linterp = all_hydro_15min$riverQ_linterp,
  BFLOW = baseflow_BFLOW(all_hydro_15min$riverQ_linterp, alpha, passes)) %>%
  mutate(BFLOW_bfi = BFLOW/riverQ_linterp)


# Delineate events based on BFI threshold ----------------------------------------------------

#Select bfi threshold and compare bfi distributions
bfi_threshold <- 0.8

ggplot(data = baseflow_15min) +
  geom_histogram(aes(x = BFLOW_bfi), binwidth = 0.01, fill = NA, color = 'blue') +
  scale_x_continuous(limits = c(0,1))

#Create events data frame using threshold
BFLOW_events_bfi_threshold <- data.frame(
  start_dateTime = baseflow_15min$dateTime[baseflow_15min$BFLOW_bfi< bfi_threshold & lag(baseflow_15min$BFLOW_bfi) >= bfi_threshold],
  end_dateTime = baseflow_15min$dateTime[baseflow_15min$BFLOW_bfi< bfi_threshold & lead(baseflow_15min$BFLOW_bfi) >= bfi_threshold],
  event_length = rle(baseflow_15min$BFLOW_bfi< bfi_threshold)$lengths[rle(baseflow_15min$BFLOW_bfi< bfi_threshold)$values]
) %>%
  #remove miniscule events
  filter(event_length > 3)

#Fix first and last rows of events data frame, which get set to NA from the above 
BFLOW_events_bfi_threshold$start_dateTime[1] <- baseflow_15min$dateTime[1]
BFLOW_events_bfi_threshold$end_dateTime[nrow(BFLOW_events_bfi_threshold)] <- baseflow_15min$dateTime[nrow(baseflow_15min)]

#Calculate max flows in each event
event.max <- function(event_frame, start_column = 'start_dateTime', end_column = 'end_dateTime', length_column = 'event_length',
                      Q_frame , Q_column, Q_dateTime_column) {
  
  event_frame$event_max <- double(length = nrow(event_frame))
  
  for(i in 1:nrow(event_frame)) {
    event_interval <- interval(event_frame[[start_column]][i], event_frame[[end_column]][i])
    Q_slice <- Q_frame[[Q_column]][Q_frame[[Q_dateTime_column]] %within% event_interval]
    event_frame$event_max[i] <- max(Q_slice)
  }
  return(event_frame)
}

BFLOW_events_bfi_threshold <- event.max(BFLOW_events_bfi_threshold, Q_frame = all_hydro_15min, Q_column = 'riverQ_linterp', Q_dateTime_column = 'dateTime')

#Write results
write_csv(BFLOW_events_bfi_threshold, './DataFiles/hydro_data/event_delineations/BFLOW_events_bfi_threshold.csv')


##Event statistics
ggplot(data = BFLOW_events_bfi_threshold) +
  geom_histogram(aes(x = event_max), binwidth = 5) +
  geom_vline(xintercept = 45)
ggplot(data = BFLOW_events_bfi_threshold) +
  geom_histogram(aes(x = event_length), binwidth = 12) +
  geom_vline(xintercept = 12)


# Plot events for manual adjustment ------------------------------------------

BFLOW_events_bfi_threshold <- read_csv('./DataFiles/hydro_data/event_delineations/BFLOW_events_bfi_threshold.csv')
events <- read_csv('DataFiles/hydro_data/event_delineations/BFLOW_events_bfi_threshold_adjusted.csv') %>%
  select(start_dateTime, end_dateTime)

ggplot(data = all_hydro_15min) +
  ##Event limits and labels
  #raw
  #geom_rect(data = BFLOW_events_bfi_threshold, aes(xmin = start_dateTime, xmax = end_dateTime), 
            #ymin = -200, ymax = 1.05*max(all_hydro_15min$riverQ_linterp), fill = 'grey90', color = 'grey50', linetype = 'dashed') +
  #geom_text(data = BFLOW_events_bfi_threshold, aes(x = start_dateTime, label = row_number(start_dateTime)), y = 0.3*max(all_hydro_15min$riverQ_linterp)) +
  #adjusted
  geom_rect(data = events, aes(xmin = start_dateTime, xmax = end_dateTime), 
            ymin = -200, ymax = 1.05*max(all_hydro_15min$riverQ_linterp), fill = 'grey90', color = 'grey50', linetype = 'dashed') +
  geom_text(data = events, aes(x = start_dateTime, label = row_number(start_dateTime)), y = 0.3*max(all_hydro_15min$riverQ_linterp)) +
  
  ##Precip
  #geom_col(data = all_hydro_daily, aes(x = as_datetime(date), y = lawrence_precip*10), fill = 'steelblue3', alpha = 0.5) +
  #geom_col(data = all_hydro_daily, aes(x = as_datetime(date), y = manhattan_precip*10), color = 'purple', fill = NA) +
  #geom_col(data = all_hydro_daily, aes(x = as_datetime(date), y = topeka_precip*10), color = 'darkblue', fill = NA) +
  #geom_col(data = all_hydro_daily, aes(x = as_datetime(date), y = salina_precip*10), color = 'tomato1', fill = NA) +
  #geom_col(aes(x = dateTime, y = belvue_precip*40), color = 'blue', alpha = 0.5) +
  
  ##River Q
  #geom_line(aes(x = dateTime, y = riverQ)) +
  geom_line(aes(x = dateTime, y = riverQ_linterp), linetype = 'longdash') +
  #geom_line(aes(x=dateTime, y = riverQ_ingap), linetype = 'longdash', color = 'red') +
  #geom_line(aes(x = dateTime, y = lawrenceQ_linterp), linetype = 'longdash', color = 'grey80') +
  
  ## Nitrate
  #geom_line(aes(x = dateTime, y = N*300), color = 'red', alpha = 0.7) +
  
  ##Reservoir Outflows
  #geom_line(aes(x = dateTime, y = clintonQ), color = 'red') +
  geom_line(aes(x = dateTime, y = milfordQ_linterp), color = 'gold') +
  #geom_line(aes(x = dateTime, y = perryQ), color = 'lightblue') +
  geom_line(aes(x = dateTime, y = tuttleQ_linterp), color = 'purple') +
  geom_line(aes(x = dateTime, y = reservoir_sum_linterp), color = 'darkgreen') +
  #geom_line(aes(x = dateTime, y = reservoir_sum_upstream_linterp), color = 'pink') +
  
  ##Baseflows
  #geom_line(aes(x = dateTime, y = HYSEP_fixed), linetype = 'dashed', color = 'red') +
  #geom_line(aes(x = dateTime, y = HYSEP_sliding), linetype = 'dashed', color = 'blue') +
  #geom_line(aes(x = dateTime, y = HYSEP_local), linetype = 'dashed', color = 'purple') +
  #geom_line(aes(x = dateTime, y = UKIH), linetype = 'dashed', color = 'orange') +
  #geom_line(data = baseflow_15min, aes(x = dateTime, y = BFLOW), linetype = 'dashed', color = 'black') +
  #geom_line(data = baseflow_15min, aes(x = dateTime, y = Eckhardt), linetype = 'dashed', color = 'black') +
  
  ##BFI target line
  #geom_line(data = baseflow_15min,aes(x = dateTime, y = riverQ_linterp*bfi_threshold), linetype = 'dashed', color = 'red') +

  ##Theme and scale
  theme_classic() +
  scale_y_continuous(limits = c(0, 3000), name = 'Q [m3/s]') +
                     #sec.axis = sec_axis(trans = ~. / 75, name = 'N [mg/L]')) +
  scale_x_datetime(limits = as_datetime(c('2017-10-01', '2018-09-30')),
                   date_breaks = '1 month', date_labels = '%b', name = NULL)
ggsave("./Figures/hydrographs/DeSoto_2019_half2_delineation.tiff", device = "tiff", height = 3.5, width = 10.5, units = "in", compression = "lzw", dpi = 700)
  


# Zoom in on individual event ---------------------------------------------

event_select <- 455
time_padding <- days(7)

ggplot(data = all_hydro_15min) +
  #Event limits and labels
  geom_rect(data = BFLOW_events_bfi_threshold, aes(xmin = start_dateTime, xmax = end_dateTime), 
            ymin = -200, ymax = 1.05*max(all_hydro_15min$riverQ_linterp), fill = 'grey90', color = 'grey50', linetype = 'dashed', alpha = 0.5) +
  geom_text(data = BFLOW_events_bfi_threshold, aes(x = start_dateTime, label = row_number(start_dateTime)), y = 1.1*BFLOW_events_bfi_threshold$event_max[event_select]) +
  #Precip
  geom_col(data = all_hydro_daily, aes(x = as_datetime(date), y = lawrence_precip*10), fill = 'steelblue3', alpha = 0.3) +
  geom_col(data = all_hydro_daily, aes(x = as_datetime(date), y = manhattan_precip*10), color = 'purple', fill = NA) +
  geom_col(data = all_hydro_daily, aes(x = as_datetime(date), y = topeka_precip*10), color = 'darkblue', fill = NA) +
  geom_col(data = all_hydro_daily, aes(x = as_datetime(date), y = salina_precip*10), color = 'tomato1', fill = NA) +
  geom_col(aes(x = dateTime, y = belvue_precip*40), color = 'blue', alpha = 0.5) +
  #River Q
  geom_line(aes(x = dateTime, y = riverQ)) +
  geom_line(aes(x = dateTime, y = riverQ_linterp), linetype = 'longdash') +
  geom_line(aes(x=dateTime, y = riverQ_ingap), linetype = 'longdash', color = 'red') +
  geom_line(aes(x = dateTime, y = lawrenceQ_linterp), linetype = 'longdash', color = 'grey70') +
  ## Nitrate
  geom_line(aes(x = dateTime, y = N*100), color = 'red', alpha = 0.7) +
  ##Reservoir Outflows
  geom_line(aes(x = dateTime, y = clintonQ), color = 'red') +
  geom_line(aes(x = dateTime, y = milfordQ), color = 'gold') +
  geom_line(aes(x = dateTime, y = perryQ), color = 'lightblue') +
  geom_line(aes(x = dateTime, y = tuttleQ), color = 'purple') +
  #geom_line(aes(x = dateTime, y = wilsonQ), color = 'pink') +
  geom_line(aes(x = dateTime, y = reservoir_sum_linterp), color = 'darkgreen') +
  #geom_line(aes(x = dateTime, y = reservoir_sum_upstream_linterp), color = 'pink') +
  ##Baseflows
  #geom_line(aes(x = dateTime, y = HYSEP_fixed), linetype = 'dashed', color = 'red') +
  #geom_line(aes(x = dateTime, y = HYSEP_sliding), linetype = 'dashed', color = 'blue') +
  #geom_line(aes(x = dateTime, y = HYSEP_local), linetype = 'dashed', color = 'purple') +
  #geom_line(aes(x = dateTime, y = UKIH), linetype = 'dashed', color = 'orange') +
  #geom_line(data = baseflow_15min, aes(x = dateTime, y = BFLOW), linetype = 'dashed', color = 'black') +
  #geom_line(data = baseflow_15min, aes(x = dateTime, y = Eckhardt), linetype = 'dashed', color = 'black') +
  ##BFI target line
  #geom_line(data = baseflow_15min,aes(x = dateTime, y = riverQ_linterp*bfi_threshold), linetype = 'dashed', color = 'red') +
  #Theme and scale
  theme_classic() +
  theme(panel.grid.major.x = element_line(color = 'grey60'), panel.grid.minor.x = element_line(color = 'grey80')) +
  scale_y_continuous(limits = c(0, max(1.25*BFLOW_events_bfi_threshold$event_max[event_select], 400))) +
  scale_x_datetime(limits = c(BFLOW_events_bfi_threshold$start_dateTime[event_select]-5*time_padding, BFLOW_events_bfi_threshold$end_dateTime[event_select]+time_padding),
                   date_breaks = '24 hours', date_labels = '%m/%d %R', date_minor_breaks = '6 hours') +
  coord_cartesian(xlim = c(BFLOW_events_bfi_threshold$start_dateTime[event_select]-time_padding, BFLOW_events_bfi_threshold$end_dateTime[event_select]+time_padding))

#Interactive version
ggplotly(
  ggplot(data = all_hydro_15min) +
  geom_line(aes(x = dateTime, y = riverQ)) +
  theme_classic() +
  scale_y_continuous(limits = c(0, max(1.25*BFLOW_events_bfi_threshold$event_max[event_select], 500))) +
  scale_x_datetime(limits = c(BFLOW_events_bfi_threshold$start_dateTime[event_select]-5*time_padding, BFLOW_events_bfi_threshold$end_dateTime[event_select]+5*time_padding),
                   date_breaks = '24 hours', date_labels = '%m/%d %R', date_minor_breaks = '4 hours') +
  coord_cartesian(xlim = c(BFLOW_events_bfi_threshold$start_dateTime[event_select]-time_padding, BFLOW_events_bfi_threshold$end_dateTime[event_select]+time_padding))
)


# AFTER MANUAL ADJUSTMENTS ------------------------------------------------
##Adjusted events created as a separate .csv @ DataFiles/hydro_data/event_delineations/BFLOW_events_bfi_threshold_adjusted.csv

# Load data --------------------------------------------

all_hydro_15min <- with_tz(read_csv('./DataFiles/hydro_data/all_hydro_15min.csv'), 'America/Chicago')
all_hydro_daily <- read_csv('./DataFiles/hydro_data/all_hydro_daily.csv')
events <- with_tz(read_csv('DataFiles/hydro_data/event_delineations/BFLOW_events_bfi_threshold_adjusted.csv'), 'America/Chicago') %>%
  select(event_number, start_dateTime, end_dateTime)


# Plot long time series ---------------------------------------------------

ggplot(data = all_hydro_15min) +
  ##Event limits and labels
  #raw
  #geom_rect(data = BFLOW_events_bfi_threshold, aes(xmin = start_dateTime, xmax = end_dateTime), 
  #ymin = -200, ymax = 1.05*max(all_hydro_15min$riverQ_linterp), fill = 'grey90', color = 'grey50', linetype = 'dashed') +
  #geom_text(data = BFLOW_events_bfi_threshold, aes(x = start_dateTime, label = row_number(start_dateTime)), y = 0.3*max(all_hydro_15min$riverQ_linterp)) +
  #adjusted
  geom_rect(data = events, aes(xmin = start_dateTime, xmax = end_dateTime), 
            ymin = -200, ymax = 1.05*max(all_hydro_15min$riverQ_linterp), fill = 'grey90', color = 'grey50', linetype = 'dashed') +
  geom_text(data = events, aes(x = start_dateTime, label = row_number(start_dateTime)), y = 600) +
  
  ##Precip
  #geom_col(data = all_hydro_daily, aes(x = as_datetime(date), y = lawrence_precip*10), fill = 'steelblue3', alpha = 0.5) +
  #geom_col(data = all_hydro_daily, aes(x = as_datetime(date), y = manhattan_precip*10), color = 'purple', fill = NA) +
  #geom_col(data = all_hydro_daily, aes(x = as_datetime(date), y = topeka_precip*10), color = 'darkblue', fill = NA) +
  #geom_col(data = all_hydro_daily, aes(x = as_datetime(date), y = salina_precip*10), color = 'tomato1', fill = NA) +
  #geom_col(aes(x = dateTime, y = belvue_precip*40), color = 'blue', alpha = 0.5) +
  
  ##River Q
  geom_line(aes(x = dateTime, y = riverQ)) +
  geom_line(aes(x = dateTime, y = riverQ_linterp), linetype = 'longdash') +
  geom_line(aes(x=dateTime, y = riverQ_ingap), linetype = 'longdash', color = 'red') +
  #geom_line(aes(x = dateTime, y = lawrenceQ_linterp), linetype = 'longdash', color = 'grey80') +
  
  ## Nitrate
  geom_line(aes(x = dateTime, y = N*100), color = 'red', alpha = 0.7) +
  
  ##Reservoir Outflows
  #geom_line(aes(x = dateTime, y = clintonQ), color = 'red') +
  #geom_line(aes(x = dateTime, y = milfordQ), color = 'gold') +
  #geom_line(aes(x = dateTime, y = perryQ), color = 'lightblue') +
  #geom_line(aes(x = dateTime, y = tuttleQ), color = 'purple') +
  #geom_line(aes(x = dateTime, y = reservoir_sum_linterp), color = 'darkgreen') +
  #geom_line(aes(x = dateTime, y = reservoir_sum_upstream_linterp), color = 'pink') +
  
  ##Baseflows
  #geom_line(aes(x = dateTime, y = HYSEP_fixed), linetype = 'dashed', color = 'red') +
  #geom_line(aes(x = dateTime, y = HYSEP_sliding), linetype = 'dashed', color = 'blue') +
  #geom_line(aes(x = dateTime, y = HYSEP_local), linetype = 'dashed', color = 'purple') +
  #geom_line(aes(x = dateTime, y = UKIH), linetype = 'dashed', color = 'orange') +
  #geom_line(data = baseflow_15min, aes(x = dateTime, y = BFLOW), linetype = 'dashed', color = 'black') +
  #geom_line(data = baseflow_15min, aes(x = dateTime, y = Eckhardt), linetype = 'dashed', color = 'black') +
  
  ##BFI target line
#geom_line(data = baseflow_15min,aes(x = dateTime, y = riverQ_linterp*bfi_threshold), linetype = 'dashed', color = 'red') +

##Theme and scale
theme_classic() +
  scale_y_continuous(limits = c(0, 3000), name = 'Q [m3/s]',
  sec.axis = sec_axis(trans = ~. / 300, name = 'N [mg/L]')) +
  scale_x_datetime(limits = as_datetime(c('2014-10-01', '2015-09-30')),
                   date_breaks = '1 month', date_labels = '%b', name = NULL)

ggsave("./Figures/hydrographs/DeSoto_2018_nitrate_rescale.tiff", device = "tiff", height = 3.5, width = 10.5, units = "in", compression = "lzw", dpi = 700)


# Zoom in on individual event ---------------------------------------------

event_select <- 5
time_padding <- days(14)

ggplot(data = all_hydro_15min) +
  #Event limits and labels
  geom_rect(data = events, aes(xmin = start_dateTime, xmax = end_dateTime), 
            ymin = -200, ymax = 1.05*max(all_hydro_15min$riverQ_linterp), fill = 'grey90', color = 'grey50', linetype = 'dashed') +
  geom_text(data = events, aes(x = start_dateTime, label = row_number(start_dateTime)), y = 0.13*max(all_hydro_15min$riverQ_linterp)) +
  #Precip
  geom_col(data = all_hydro_daily, aes(x = as_datetime(date), y = lawrence_precip*10), fill = 'steelblue3', alpha = 0.3) +
  geom_col(data = all_hydro_daily, aes(x = as_datetime(date), y = manhattan_precip*10), color = 'purple', fill = NA) +
  geom_col(data = all_hydro_daily, aes(x = as_datetime(date), y = topeka_precip*10), color = 'darkblue', fill = NA) +
  geom_col(data = all_hydro_daily, aes(x = as_datetime(date), y = salina_precip*10), color = 'tomato1', fill = NA) +
  geom_col(aes(x = dateTime, y = belvue_precip*40), color = 'blue', alpha = 0.5) +
  #River Q
  geom_line(aes(x = dateTime, y = riverQ)) +
  geom_line(aes(x = dateTime, y = riverQ_linterp), linetype = 'longdash') +
  geom_line(aes(x=dateTime, y = riverQ_ingap), linetype = 'longdash', color = 'red') +
  geom_line(aes(x = dateTime, y = lawrenceQ_linterp), linetype = 'longdash', color = 'grey70') +
  ## Nitrate
  geom_line(aes(x = dateTime, y = N*100), color = 'red', alpha = 0.7) +
  ##Reservoir Outflows
  geom_line(aes(x = dateTime, y = clintonQ), color = 'red') +
  geom_line(aes(x = dateTime, y = milfordQ), color = 'gold') +
  geom_line(aes(x = dateTime, y = perryQ), color = 'lightblue') +
  geom_line(aes(x = dateTime, y = tuttleQ), color = 'purple') +
  #geom_line(aes(x = dateTime, y = wilsonQ), color = 'pink') +
  geom_line(aes(x = dateTime, y = reservoir_sum_linterp), color = 'darkgreen') +
  #geom_line(aes(x = dateTime, y = reservoir_sum_upstream_linterp), color = 'pink') +
  ##Baseflows
  #geom_line(aes(x = dateTime, y = HYSEP_fixed), linetype = 'dashed', color = 'red') +
  #geom_line(aes(x = dateTime, y = HYSEP_sliding), linetype = 'dashed', color = 'blue') +
  #geom_line(aes(x = dateTime, y = HYSEP_local), linetype = 'dashed', color = 'purple') +
  #geom_line(aes(x = dateTime, y = UKIH), linetype = 'dashed', color = 'orange') +
  #geom_line(data = baseflow_15min, aes(x = dateTime, y = BFLOW), linetype = 'dashed', color = 'black') +
  #geom_line(data = baseflow_15min, aes(x = dateTime, y = Eckhardt), linetype = 'dashed', color = 'black') +
  ##BFI target line
  #geom_line(data = baseflow_15min,aes(x = dateTime, y = riverQ_linterp*bfi_threshold), linetype = 'dashed', color = 'red') +
  #Theme and scale
theme_classic() +
  theme(panel.grid.major.x = element_line(color = 'grey60'), panel.grid.minor.x = element_line(color = 'grey80')) +
  scale_y_continuous(limits = c(0, 400)) +
  scale_x_datetime(limits = c(events$start_dateTime[event_select]-5*time_padding, events$end_dateTime[event_select]+time_padding),
                   date_breaks = '24 hours', date_labels = '%m/%d %R', date_minor_breaks = '6 hours') +
  coord_cartesian(xlim = c(events$start_dateTime[event_select]-time_padding, events$end_dateTime[event_select]+time_padding))

#Interactive version
ggplotly(
  ggplot(data = all_hydro_15min) +
    geom_line(aes(x = dateTime, y = riverQ)) +
    theme_classic() +
    scale_y_continuous(limits = c(0, max(1.25*events$event_max[event_select], 500))) +
    scale_x_datetime(limits = c(events$start_dateTime[event_select]-5*time_padding, events$end_dateTime[event_select]+5*time_padding),
                     date_breaks = '24 hours', date_labels = '%m/%d %R', date_minor_breaks = '4 hours') +
    coord_cartesian(xlim = c(events$start_dateTime[event_select]-time_padding, events$end_dateTime[event_select]+time_padding))
)

# Assign each observation as event/non-event -------------------------------------------------

all_hydro_15min <- with_tz(read_csv('./DataFiles/hydro_data/all_hydro_15min.csv'), 'America/Chicago')
all_hydro_daily <- read_csv('./DataFiles/hydro_data/all_hydro_daily.csv')
events <- with_tz(read_csv('DataFiles/hydro_data/event_delineations/BFLOW_events_bfi_threshold_adjusted.csv'), 'America/Chicago') %>%
  select(event_number, start_dateTime, end_dateTime)

#Function that creates a list of interval for each event
#Allows usage of lubridate::%within% to check if a data point is inside an event
interval_list <- function(data){
  
  list <- list()
  
  for(i in 1:nrow(data)) {
    list[[i]] <- interval(data[['start_dateTime']][i], data[['end_dateTime']][i])
  }
  names(list) <- data[['event_number']]
  return(list)
}

event_list <- interval_list(events)

#Assigns each data point to its resepective event
event_assignment <- function(data, events, event_output){
  
  data$event_output <- double(length = nrow(data))
  
  for(i in 1:length(events)){
    data$event_output <- ifelse(data[['dateTime']] %within% events[[i]], names(events)[[i]], data$event_output)
    
  }
  return(data)
}

all_hydro_15min <- event_assignment(all_hydro_15min, event_list) %>%
  mutate(event = ifelse(event_output == 0, NA, event_output)) %>%
  select(!event_output)

write_csv(all_hydro_15min, './DataFiles/hydro_data/all_hydro_15min.csv')


#Assign Compound Events

events <- with_tz(read_csv('DataFiles/hydro_data/event_delineations/BFLOW_events_bfi_threshold_adjusted.csv'), 'America/Chicago') %>%
  select(event_number, start_dateTime, end_dateTime, compound_start_rank1, compound_start_rank2) %>%
  mutate(compound_start_rank1 = ifelse(is.na(compound_start_rank1), event_number, compound_start_rank1),
         compound_start_rank2 = ifelse(is.na(compound_start_rank2), event_number, compound_start_rank2))

events_compound1 <- events %>%
  group_by(event_number = compound_start_rank1) %>%
  summarise(start_dateTime = min(start_dateTime),
            end_dateTime = max(end_dateTime))
write_csv(events_compound1, './DataFiles/hydro_data/event_delineations/events_compound1.csv')

events_compound2 <- events %>%
  group_by(event_number = compound_start_rank2) %>%
  summarise(start_dateTime = min(start_dateTime),
            end_dateTime = max(end_dateTime))
write_csv(events_compound2, './DataFiles/hydro_data/event_delineations/events_compound2.csv')

events_list_compound1 <- interval_list(events_compound1)
all_hydro_15min <- event_assignment(all_hydro_15min, events_list_compound1) %>%
  mutate(event_compound1 = ifelse(event_output == 0, NA, event_output)) %>%
  select(!event_output)

events_list_compound2 <- interval_list(events_compound2)
all_hydro_15min <- event_assignment(all_hydro_15min, events_list_compound2) %>%
  mutate(event_compound2 = ifelse(event_output == 0, NA, event_output)) %>%
  select(!event_output)

write_csv(all_hydro_15min, './DataFiles/hydro_data/all_hydro_15min.csv')        




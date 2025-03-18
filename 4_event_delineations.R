# Load libraries and data -------------------------------------------------

library(tidyverse)
library(zoo)

# load timeseries with at minimum a dateTime column named 'dateTime' and a 
# discharge column named 'q'. I don't think NAs are allowed
dat <- read_csv('./DataFiles/hydro_data/all_hydro_15min.csv') %>%
  select(dateTime, q = riverQ_linterp)

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


# Run filter --------------------------------------------------------------

# Filter parameters. More passes for more frequent data. Default is 3 for daily,
# ~11 recommended for 15-min. Beta can be anywhere > 0.9, just tune to something
# that looks good. More info in Ladson et al. (2013), DOI 10.7158/W12-028.2013.17.1
beta <- 0.992
passes <- 11

#Run filter
dat_bf <- data.frame(
  dateTime = dat$dateTime,
  q = dat$q,
  bf = baseflow_BFLOW(dat$q, beta, passes)) %>%
  mutate(bfi = bf/q)
median(dat_bf$bf)

# dateTime range for checking parameters
check_year <- 2017
dat_check <- filter(dat_bf, year(dateTime) == check_year)
# plot and check parameters
ggplot(data = dat_check) +
  geom_line(aes(x = dateTime, y = q)) +
  geom_line(aes(x = dateTime, y = bf), color = 'red') +
  theme_classic()


# Delineate events based on BFI threshold ----------------------------------------------------

# Select bfi threshold. I just picked a value that seemed good based on the histogram
bfi_threshold <- 0.8

ggplot(data = dat_bf) +
  geom_histogram(aes(x = bfi), binwidth = 0.01) +
  geom_vline(xintercept = bfi_threshold, color = 'red') +
  scale_x_continuous(limits = c(0,1))

# Create events data frame using threshold. Min_length argument filters out events
# less than the specified length in days. By default this function will 
delineator <- function(x, bfi_threshold, min_length = NA){
  dateTime <- x[['dateTime']]
  bfi <- x[['bfi']]

  events <- data.frame(
    start = dateTime[bfi < bfi_threshold & lag(bfi >= bfi_threshold)],
    end = dateTime[bfi < bfi_threshold & lead(bfi >= bfi_threshold)]
    )
  
  events$start[1] <- dateTime[1]
  events$end[length(events$end)] <- dateTime[length(dateTime)]
  
  events$length <- as.numeric(as.duration(interval(events$start, events$end)),'days')
  
  if(!is.na(min_length)){events <- filter(events, length >= min_length)}
  
  events[['event_id']] <- seq(1:nrow(events))
  
  return(events)
}

events <- delineator(dat_bf, bfi_threshold, min_length = 3/96)
write_csv(events, './DataFiles/hydro_data/event_delineations/automated_delineations.csv')

# plot it up --------------------------------------------------------------
# Use these plots to make manual adjustments. Zoom into individual events, check
# for missing values, etc. All adjustments and justifications logged in 
# /DataFiles/hydro_data/event_delineations/adjustment_log.xlsx

check_year <- 2017
dat_check <- filter(dat_bf, year(dateTime) == check_year)
events_check <- filter(events, year(start) == check_year)

ggplot(data = dat_check) +
  geom_line(aes(x = dateTime, y = q)) +
  geom_line(aes(x = dateTime, y = bf), color = 'red') +
  geom_rect(data = events_check,
            aes(xmin = start, xmax = end, ymin = 0, ymax = 1.05*max(dat_check$q)),
            color = 'grey50', alpha = 0.3, linetype = 'dashed') +
  geom_text(data = events_check,
            aes(x = start, label = event_id),
            y = 0.95*max(dat_check$q)) +
  theme_classic()


# After all manual adjustments:
# assign each observation to event ----------------------------------------
all_hydro_15min <- with_tz(read_csv('./DataFiles/hydro_data/all_hydro_15min.csv'), 'America/Chicago')
events <- with_tz(read_csv('./DataFiles/hydro_data/event_delineations/adjusted_events.csv'), 'America/Chicago') %>%
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

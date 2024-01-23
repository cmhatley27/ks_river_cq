library(tidyverse)
library(lubridate)
library(dataRetrieval)
library(rnoaa)
library(hydroEvents)
library(FluMoDL)
library(jsonlite)
library(curl)


# Create folders and load data --------------------------------------------

all_hydro_15min <- with_tz(read_csv('./DataFiles/hydro_data/all_hydro_15min.csv'), tz = 'America/Chicago')
all_hydro_daily <- read_csv('./DataFiles/hydro_data/all_hydro_daily.csv')

if(!dir.exists('DataFiles/hydro_data/reservoir_outflows/na_events')) {
  dir.create('DataFiles/hydro_data/reservoir_outflows/na_events')
}


# NA Handling -----------------------------------------------

##Discharge: 
##Assigning >12 consecutive NAs (3 hours) as a gap

##Gaps are not removed from the dataset before calculating daily averages or delineating events
##Event delineation essentially interpolates between gaps anyway
##This is just for visualization when adjusting events so I can throw out events that occur during gaps


#*** DeSoto gage -------------------------------------------------------------


#*** *** Discharge ---------------------------------------------------------------


##Detect NA events and number of consecutive NAs in each
DeSoto_15min_NAevents <- data.frame(
  start_dateTime = all_hydro_15min$dateTime[is.na(all_hydro_15min$riverQ) & !is.na(lag(all_hydro_15min$riverQ))],
  end_dateTime = all_hydro_15min$dateTime[is.na(all_hydro_15min$riverQ) & !is.na(lead(all_hydro_15min$riverQ))],
  na_length = rle(is.na(all_hydro_15min$riverQ))$lengths[rle(is.na(all_hydro_15min$riverQ))$values]
) %>%
  # Select gap length here (15 minute increments, e.g. gap length of 12 = 3 hours)
  mutate(gap = na_length > 12)
write_csv(DeSoto_15min_NAevents, './DataFiles/hydro_data/DeSoto_gage/DeSoto_15min_NAevents_N.csv')

##Summary stats of NA events
summary(DeSoto_15min_NAevents$na_length)
sum(DeSoto_15min_NAevents$na_length[DeSoto_15min_NAevents$gap])
sum(DeSoto_15min_NAevents$na_length[DeSoto_15min_NAevents$gap])/nrow(all_hydro_15min)



#*** *** Concentration -----------------------------------------------------------

##Detect NA events and number of consecutive NAs in each
DeSoto_15min_NAevents_N <- data.frame(
  start_dateTime = all_hydro_15min$dateTime[is.na(all_hydro_15min$N) & !is.na(lag(all_hydro_15min$N))],
  end_dateTime = all_hydro_15min$dateTime[is.na(all_hydro_15min$N) & !is.na(lead(all_hydro_15min$N))],
  na_length = rle(is.na(all_hydro_15min$N))$lengths[rle(is.na(all_hydro_15min$N))$values]
) %>%
  # Select gap length here (15 minute increments, e.g. gap length of 12 = 3 hours)
  mutate(gap = na_length > 12)
write_csv(DeSoto_15min_NAevents_N, './DataFiles/hydro_data/DeSoto_gage/DeSoto_15min_NAevents_N.csv')

##Summary stats of NA events
summary(DeSoto_15min_NAevents_N$na_length)
sum(DeSoto_15min_NAevents_N$na_length[DeSoto_15min_NAevents_N$gap])
sum(DeSoto_15min_NAevents_N$na_length[DeSoto_15min_NAevents_N$gap])/nrow(all_hydro_15min)



#*** Reservoir outflow gages -------------------------------------------------------------

#*** *** Clinton -----------------------------------------------------------------

##Detect NA events and number of consecutive NAs in each
clinton_15min_NAevents <- data.frame(
  start_dateTime = all_hydro_15min$dateTime[is.na(all_hydro_15min$clintonQ) & !is.na(lag(all_hydro_15min$clintonQ))],
  end_dateTime = all_hydro_15min$dateTime[is.na(all_hydro_15min$clintonQ) & !is.na(lead(all_hydro_15min$clintonQ))],
  na_length = rle(is.na(all_hydro_15min$clintonQ))$lengths[rle(is.na(all_hydro_15min$clintonQ))$values]
) %>%
  # Select gap length here (15 minute increments, e.g. gap length of 12 = 3 hours)
  mutate(gap = na_length > 12)
write_csv(clinton_15min_NAevents, './DataFiles/hydro_data/reservoir_outflows/NA_events/clinton_15min_NAevents.csv')

##Summary stats of NA events
summary(clinton_15min_NAevents$na_length)
##Total length of all NAs, and as fraction of full data set
sum(clinton_15min_NAevents$na_length)
sum(clinton_15min_NAevents$na_length)/nrow(all_hydro_15min)
##Total length of all gaps, and as fraction of full data set
sum(clinton_15min_NAevents$na_length[clinton_15min_NAevents$gap])
sum(clinton_15min_NAevents$na_length[clinton_15min_NAevents$gap])/nrow(all_hydro_15min)


#*** *** Milford -----------------------------------------------------------------

##Detect NA events and number of consecutive NAs in each
milford_15min_NAevents <- data.frame(
  start_dateTime = all_hydro_15min$dateTime[is.na(all_hydro_15min$milfordQ) & !is.na(lag(all_hydro_15min$milfordQ))],
  end_dateTime = all_hydro_15min$dateTime[is.na(all_hydro_15min$milfordQ) & !is.na(lead(all_hydro_15min$milfordQ))],
  na_length = rle(is.na(all_hydro_15min$milfordQ))$lengths[rle(is.na(all_hydro_15min$milfordQ))$values]
) %>%
  # Select gap length here (15 minute increments, e.g. gap length of 12 = 3 hours)
  mutate(gap = na_length > 12)
write_csv(milford_15min_NAevents, './DataFiles/hydro_data/reservoir_outflows/NA_events/milford_15min_NAevents.csv')

##Summary stats of NA events
summary(milford_15min_NAevents$na_length)
##Total length of all NAs, and as fraction of full data set
sum(milford_15min_NAevents$na_length)
sum(milford_15min_NAevents$na_length)/nrow(all_hydro_15min)
##Total length of all gaps, and as fraction of full data set
sum(milford_15min_NAevents$na_length[milford_15min_NAevents$gap])
sum(milford_15min_NAevents$na_length[milford_15min_NAevents$gap])/nrow(all_hydro_15min)


#*** *** Perry -----------------------------------------------------------------

##Detect NA events and number of consecutive NAs in each
perry_15min_NAevents <- data.frame(
  start_dateTime = all_hydro_15min$dateTime[is.na(all_hydro_15min$perryQ) & !is.na(lag(all_hydro_15min$perryQ))],
  end_dateTime = all_hydro_15min$dateTime[is.na(all_hydro_15min$perryQ) & !is.na(lead(all_hydro_15min$perryQ))],
  na_length = rle(is.na(all_hydro_15min$perryQ))$lengths[rle(is.na(all_hydro_15min$perryQ))$values]
) %>%
  # Select gap length here (15 minute increments, e.g. gap length of 12 = 3 hours)
  mutate(gap = na_length > 12)
write_csv(perry_15min_NAevents, './DataFiles/hydro_data/reservoir_outflows/NA_events/perry_15min_NAevents.csv')

##Summary stats of NA events
summary(perry_15min_NAevents$na_length)
##Total length of all NAs, and as fraction of full data set
sum(perry_15min_NAevents$na_length)
sum(perry_15min_NAevents$na_length)/nrow(all_hydro_15min)
##Total length of all gaps, and as fraction of full data set
sum(perry_15min_NAevents$na_length[perry_15min_NAevents$gap])
sum(perry_15min_NAevents$na_length[perry_15min_NAevents$gap])/nrow(all_hydro_15min)


#*** *** Tuttle -----------------------------------------------------------------

##Detect NA events and number of consecutive NAs in each
tuttle_15min_NAevents <- data.frame(
  start_dateTime = all_hydro_15min$dateTime[is.na(all_hydro_15min$tuttleQ) & !is.na(lag(all_hydro_15min$tuttleQ))],
  end_dateTime = all_hydro_15min$dateTime[is.na(all_hydro_15min$tuttleQ) & !is.na(lead(all_hydro_15min$tuttleQ))],
  na_length = rle(is.na(all_hydro_15min$tuttleQ))$lengths[rle(is.na(all_hydro_15min$tuttleQ))$values]
) %>%
  # Select gap length here (15 minute increments, e.g. gap length of 12 = 3 hours)
  mutate(gap = na_length > 12)
write_csv(tuttle_15min_NAevents, './DataFiles/hydro_data/reservoir_outflows/NA_events/tuttle_15min_NAevents.csv')

##Summary stats of NA events
summary(tuttle_15min_NAevents$na_length)
##Total length of all NAs, and as fraction of full data set
sum(tuttle_15min_NAevents$na_length)
sum(tuttle_15min_NAevents$na_length)/nrow(all_hydro_15min)
##Total length of all gaps, and as fraction of full data set
sum(tuttle_15min_NAevents$na_length[tuttle_15min_NAevents$gap])
sum(tuttle_15min_NAevents$na_length[tuttle_15min_NAevents$gap])/nrow(all_hydro_15min)



# Create Q columns of only gap measurements -------------------------------

##For visualization purposes

#*** Generate *lists* of gap events as interval objects -------------------------------------------

DeSoto_15min_NAevents <- read_csv('./DataFiles/hydro_data/DeSoto_gage/DeSoto_15min_NAevents.csv')
DeSoto_15min_NAevents_N <- read_csv('./DataFiles/hydro_data/DeSoto_gage/DeSoto_15min_NAevents_N.csv')
clinton_15min_NAevents <- read_csv('./DataFiles/hydro_data/reservoir_outflows/NA_events/clinton_15min_NAevents.csv')
milford_15min_NAevents <- read_csv('./DataFiles/hydro_data/reservoir_outflows/NA_events/milford_15min_NAevents.csv')
perry_15min_NAevents <- read_csv('./DataFiles/hydro_data/reservoir_outflows/NA_events/perry_15min_NAevents.csv')
tuttle_15min_NAevents <- read_csv('./DataFiles/hydro_data/reservoir_outflows/NA_events/tuttle_15min_NAevents.csv')

gap.list <- function(na_frame, start_column = 'start_dateTime', end_column = 'end_dateTime', length_column = 'na_length', gap_length = 12) {
  
  gap_frame <- na_frame %>%
    filter(na_frame[[length_column]] > gap_length)
  interval_list <- list()
  
  for(i in 1:nrow(gap_frame)) {
    interval_list[[i]] <- interval(gap_frame[[start_column]][i], gap_frame[[end_column]][i])
  }
  return(interval_list)
}

DeSoto_gap_list <- gap.list(DeSoto_15min_NAevents)
DeSoto_N_gap_list <- gap.list(DeSoto_15min_NAevents_N)
clinton_gap_list <- gap.list(clinton_15min_NAevents)
milford_gap_list <- gap.list(milford_15min_NAevents)
perry_gap_list <- gap.list(perry_15min_NAevents)
tuttle_gap_list <- gap.list(tuttle_15min_NAevents)


#*** Create columns for Q within gaps, update master hydro data frame ---------------------------------------------

all_hydro_15min <- all_hydro_15min %>%
  mutate(riverQ_ingap = ifelse(dateTime %within% DeSoto_gap_list, riverQ_linterp, NA),
         N_ingap = ifelse(dateTime %within% DeSoto_N_gap_list, N_linterp, NA),
         clintonQ_ingap = ifelse(dateTime %within% clinton_gap_list, clintonQ_linterp, NA),
         milfordQ_ingap = ifelse(dateTime %within% milford_gap_list, milfordQ_linterp, NA),
         perryQ_ingap = ifelse(dateTime %within% perry_gap_list, perryQ_linterp, NA),
         tuttleQ_ingap = ifelse(dateTime %within% tuttle_gap_list, tuttleQ_linterp, NA))
write_csv(all_hydro_15min, './DataFiles/hydro_data/all_hydro_15min.csv')


# Showing that trimming data set of gaps and then delineating has no effect  --------

##Trim data set of detected gaps
DeSoto_15min_gaptrimmed <- all_hydro_15min %>%
  mutate(in_gap = ifelse(dateTime %within% DeSoto_gap_list, 1, 0)) %>%
  filter(in_gap == 0) %>% select(!in_gap)

DeSoto_daily_gaptrimmed <- DeSoto_15min_gaptrimmed %>%
  group_by(dateTime = date(dateTime)) %>%
  summarise(riverQ_linterp = mean(riverQ_linterp))


##Run event delineation on de-gapped daily data
bfi_threshold <- 0.85
mindiff <- 1

DeSoto_daily_gaptrimmed.bfi <- eventBaseflow(DeSoto_daily_gaptrimmed$riverQ_linterp, 
                                      BFI_Th = bfi_threshold, min.diff = mindiff)
DeSoto_daily_gaptrimmed.bfi <- DeSoto_daily_gaptrimmed.bfi %>%
  mutate(event_start = date(DeSoto_daily_gaptrimmed$dateTime[DeSoto_daily_gaptrimmed.bfi$srt]),
         event_end = date(DeSoto_daily_gaptrimmed$dateTime[DeSoto_daily_gaptrimmed.bfi$end]),
         event_max_date = date(DeSoto_daily_gaptrimmed$dateTime[DeSoto_daily_gaptrimmed.bfi$which.max])) %>%
  select(event_start, event_end, event_max_date, event_max_Q = max)

startyear <- 2015
endyear <- 2015

ggplot() +
  #Event Highlights & Labels
  geom_rect(data = DeSoto_daily_gaptrimmed.bfi, aes(xmin = event_start, xmax = event_end), 
            ymin = -200, ymax = 1.05*max(DeSoto_daily_gaptrimmed$riverQ_linterp), fill = 'grey95', color = 'grey50', linetype = 'dashed') +
  geom_text(data = DeSoto_daily_gaptrimmed.bfi, aes(x = event_start + 1, label = row_number(event_start)), y = 0.9*max(DeSoto_daily_gaptrimmed$riverQ_linterp)) +
  #Gap Highlights
  geom_rect(data = subset(DeSoto_NAs, gap == TRUE), aes(xmin = date(start_dateTime), xmax = date(end_dateTime)), 
            ymin = -200, ymax = 1.05*max(DeSoto_daily_gaptrimmed$riverQ_linterp), fill = 'red', color = 'darkred', linetype = 'dashed', alpha = 0.1) +
  #River Q and baseflow
  geom_line(data = DeSoto_daily_gaptrimmed, aes(x = dateTime, y = riverQ_linterp)) +
  theme_classic() +
  scale_x_date(limits = as_date(c(str_c(startyear-1, '-10-01'), str_c(endyear, '-09-30'))),
               date_breaks = str_c(endyear-startyear+1,' months'),
               date_labels = '%b %y')

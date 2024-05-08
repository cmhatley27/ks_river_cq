library(tidyverse)
log <- read_csv(file.path('DataFiles','hydro_data','event_delineations','BFLOW_events_bfi_threshold_adjustment_log.csv'))

starts <- filter(log, action %in% c('adjust_start', 'adjust_both')) %>%
  mutate(start_diff = abs(as.numeric((adjusted_start - start_dateTime)/60/60))) %>%
  select(start_diff)

ends <- filter(log, action %in% c('adjust_end', 'adjust_both')) %>%
          mutate(end_diff = abs(as.numeric((adjusted_end - end_dateTime)/60/60))) %>%
          select(end_diff)

diffs <- c(starts$start_diff, ends$end_diff)

summary(diffs)

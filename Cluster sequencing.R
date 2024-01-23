library(tidyverse)
library(cluster)  
library(factoextra)
library(reshape)


# Load data ---------------------------------------------------------------


# * Basic characteristics --------------------------------------------------

CQ_summary <- read_csv('./DataFiles/hydro_data/event_char/CQ_summary.csv') %>%
  select(FI, HI = HI_mean)

simple_temporal <- read_csv('./DataFiles/hydro_data/event_char/simple_temporal.csv')
simple_size <- read_csv('./DataFiles/hydro_data/event_char/simple_size.csv')

N_yield <- read_csv('./DataFiles/hydro_data/event_char/N_yield.csv') %>%
  cbind(., CQ_summary)

# * Gridded Precip Predictors -----------------------------------------------
precip_inevent_gridded <- read_csv('./DataFiles/hydro_data/event_char/precip_inevent_gridded.csv') %>%
  select(contains('1d')) %>%
  select(contains(c('10260008', '10260010', '10260015', '102701')))

precip_anteevent_gridded <- read_csv('./DataFiles/hydro_data/event_char/precip_anteevent_gridded.csv') %>%
  select(contains(c('10260008', '10260010', '10260015', '102701')))

precip_unintercepted_volumes <- read_csv('./DataFiles/hydro_data/event_char/precip_unintercepted_volumes.csv')

max_precip_coords <- read_csv('./DataFiles/hydro_data/event_char/max_precip_coords.csv')

avg_precip_intensity <- read_csv('./DataFiles/hydro_data/event_char/avg_precip_intensity.csv') %>%
  select(contains('1d')) %>%
  select(contains(c('10260008', '10260010', '10260015', '102701')))

max_precip_intensity <- read_csv('./DataFiles/hydro_data/event_char/max_precip_intensity.csv') %>%
  select(contains(c('10260008', '10260010', '10260015', '102701')))


# * Outflows Predictors ---------------------------------------------------
outflows_inevent_totals <- read_csv('./DataFiles/hydro_data/event_char/outflows_inevent_totals.csv') %>%
  select(contains('1d'))

outflows_deltas <- read_csv('./DataFiles/hydro_data/event_char/outflows_deltas.csv') %>%
  select(contains('1d'))

flow_ratios <- read_csv('./DataFiles/hydro_data/event_char/flow_ratios.csv')


# * Drought Predictors ------------------------------------------------------
spei1 <- read_csv('./DataFiles/hydro_data/event_char/spei1.csv') %>%
  select(contains(c('10260015', '10270104')))

spei3 <- read_csv('./DataFiles/hydro_data/event_char/spei3.csv') %>%
  select(contains(c('10260015', '10270104')))

spei6 <- read_csv('./DataFiles/hydro_data/event_char/spei6.csv') %>%
  select(contains(c('10260015', '10270104')))

spei12 <- read_csv('./DataFiles/hydro_data/event_char/spei12.csv') %>%
  select(contains(c('10260015', '10270104')))

# * Correlated Predictors -------------------------------------------------

correlated_predictors <- unlist(as.vector(read_csv('./DataFiles/hydro_data/event_char/correlated_predictors.csv'))) 

# * Combine ----------------------------------------------

event_summary <- tibble(
  CQ_summary, 
  simple_temporal, 
  simple_size,
  precip_inevent_gridded,
  precip_anteevent_gridded,
  precip_unintercepted_volumes,
  max_precip_coords,
  avg_precip_intensity,
  max_precip_intensity,
  outflows_inevent_totals,
  outflows_deltas,
  flow_ratios,
  spei1,
  spei3,
  spei6,
  spei12
) %>%
  mutate(wateryear = factor(wateryear),
         month = factor(month),
         season = factor(season)) %>%
  select(!contains(correlated_predictors))

rm(list = ls()[ls() != 'event_summary' & ls() != 'N_yield'])


# Cluster by CQ Indices --------------------------------------------------

# * create clusters -------------------------------------------------------

hclust_data <- event_summary %>%
  filter(complete.cases(.)) %>%
  select(FI, HI)

dist <- daisy(hclust_data, metric = 'euclidean')

cls <- agnes(dist, method = "ward")

k <- 8

cls_id <- cutree(cls, k = k)

#Add cluster ID back in
event_summary_cluster <- event_summary %>%
  filter(complete.cases(.)) %>%
  mutate(cluster = cls_id,
         season = factor(season, levels = c('Winter', 'Spring', 'Summer', 'Fall')))
ggplot(data = event_summary_cluster) +
  geom_point(aes(x = FI, y = HI, color = factor(cluster))) +
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)

# Markov chain ------------------------------------------------------------

clus_t <- cls_id[1:(length(cls_id)-1)]
clus_t1 <- cls_id[-1]

freq <- matrix(data = NA,k,k)
ex <- matrix(data = NA,k,k)
tp <- matrix(data = NA,k,k)

for(i in seq(k)){
  for(j in seq(k)){
    is_t <- clus_t == i
    is_t1 <- clus_t1 == j
    
    freq[i,j] <- sum(is_t & is_t1)
    ex[i,j] <- sum(is_t)*sum(is_t1)/length(clus_t)
    tp[i,j] <- sum(is_t & is_t1)/sum(is_t)
  }
}

chisq <- sum(((freq-ex)^2)/ex)
df <- (k-1)^2
p <- pchisq(chisq, df)



# Visualizing -------------------------------------------------------------
for(k in seq(3,8)){
  k <- k
  cls_id <- cutree(cls, k = k)
  event_summary_cluster <- event_summary %>%
    filter(complete.cases(.)) %>%
    mutate(cluster = cls_id,
           season = factor(season, levels = c('Winter', 'Spring', 'Summer', 'Fall')))
  ggplot(data = event_summary_cluster) +
    geom_point(aes(x = FI, y = HI, color = factor(cluster))) +
    geom_hline(yintercept = 0)+
    geom_vline(xintercept = 0)

  add_cluster <- left_join(event_summary, event_summary_cluster)

  all_hydro_15min <- with_tz(read_csv('./DataFiles/hydro_data/all_hydro_15min.csv'), tz = 'America/Chicago')
  event_data <- with_tz(read_csv('./DataFiles/hydro_data/event_delineations/BFLOW_events_bfi_threshold_adjusted.csv'), tz = 'America/Chicago')%>%
    select(start_dateTime, end_dateTime) %>%
    cbind(., cluster = add_cluster$cluster) %>%
    mutate(wateryear = ifelse(month(start_dateTime) > 9, year(start_dateTime) + 1, year(start_dateTime)))

  for(years in seq(2014,2022)){
    year_select <- years
    ggplot() +
      geom_rect(data = subset(event_data, wateryear == year_select), aes(xmin = start_dateTime, xmax = end_dateTime, ymin = 0, ymax = 3000,
                                         fill = factor(cluster)), alpha = 0.5) +
      geom_text(data = subset(event_data, wateryear == year_select), aes(x = start_dateTime + as.duration(interval(start_dateTime, end_dateTime))/2, 
                                                                         y = max(all_hydro_15min$riverQ_linterp[all_hydro_15min$wateryear == year_select]), label = cluster)) +
      geom_line(data = subset(all_hydro_15min, wateryear == year_select), aes(x = dateTime, y = riverQ_linterp)) +
      scale_fill_discrete(guide = 'none') +
      scale_x_datetime(date_breaks = '1 month', date_labels = '%b %y') +
      xlab(year_select) +
      ylab('Q [m3/s]')
    ggsave(paste0('./Figures/events/cluster_sequences/k',k,'_',year_select,'.tiff'), device = 'tiff', compression = 'lzw',
          height = 4, width = 8, units = 'in')
  }
}

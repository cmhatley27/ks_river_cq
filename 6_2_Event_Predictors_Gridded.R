library(tidyverse)
library(sf)
library(terra)


# Handy functions ---------------------------------------------------------

#raster to data frame
r2df <- function(raster, value_name = 'value'){
  df <- terra::as.data.frame(raster, xy = TRUE) %>%
    pivot_longer(-c(x,y), names_to = 'date', values_to = value_name) %>%
    mutate(date = ymd(date))}
  
# Load in data ------------------------------------------------------------

#Data frames
all_hydro_15min <- with_tz(read_csv('./DataFiles/hydro_data/all_hydro_15min.csv'), 'America/Chicago')
events <- with_tz(read_csv('DataFiles/hydro_data/event_delineations/BFLOW_events_bfi_threshold_adjusted.csv'), 'America/Chicago') %>%
  select(event_number, start_dateTime, end_dateTime)
CQ_summary <- with_tz(read_csv('./DataFiles/hydro_data/CQ_summary.csv'), 'America/Chicago')

#Shapefiles
eksrb_boundary <- st_read('C:/School/SAFE KAW/Data/DataFiles/shapefiles/watershed_bndry.shp')
ksrb_boundary <- st_read('C:/School/SAFE KAW/Data/DataFiles/shapefiles/created/ksrb_boundary.shp')
huc8_boundaries <- st_read('C:/School/SAFE KAW/Data/DataFiles/shapefiles/created/huc8_boundaries.shp')
huc12_boundaries_eksrb <- st_read('C:/School/SAFE KAW/Data/DataFiles/shapefiles/created/huc12_boundaries_eksrb.shp')
reservoirs <- st_read('C:/School/SAFE KAW/Data/DataFiles/shapefiles/created/reservoirs.shp')
rivers_ksrb <- st_read('C:/School/SAFE KAW/Data/DataFiles/shapefiles/created/rivers_ksrb.shp')

#Rasters
precip_ksrb <- rast('C:/School/SAFE KAW/Data/DataFiles/climate_data/daily_rasters/precip_ksrb.tiff')
precip_ksrb_df <- r2df(precip_ksrb, 'precip')

spei1_ksrb <- rast('C:/School/SAFE KAW/Data/DataFiles/climate_data/monthly_rasters/spei1_ksrb.tiff')
spei1_ksrb_df <- r2df(spei1_ksrb, 'spei1')

spei3_ksrb <- rast('C:/School/SAFE KAW/Data/DataFiles/climate_data/monthly_rasters/spei3_ksrb.tiff')
spei3_ksrb_df <- r2df(spei3_ksrb, 'spei3')

spei6_ksrb <- rast('C:/School/SAFE KAW/Data/DataFiles/climate_data/monthly_rasters/spei6_ksrb.tiff')
spei6_ksrb_df <- r2df(spei6_ksrb, 'spei6')

spei12_ksrb <- rast('C:/School/SAFE KAW/Data/DataFiles/climate_data/monthly_rasters/spei12_ksrb.tiff')
spei12_ksrb_df <- r2df(spei12_ksrb, 'spei12')

speis <- sds(list(spei1_ksrb, spei3_ksrb, spei6_ksrb, spei12_ksrb))


# Mean daily precip by HUC8 -----------------------------------------------

precip_huc8 <- tibble(
  date = as_date(names(precip_ksrb))
)

for(huc8 in 1:nrow(huc8_boundaries)){

  huc8_ids <- huc8_boundaries$HUC8
  
  local_boundary <- huc8_boundaries %>%
    filter(HUC8 == huc8_ids[huc8])
  
  local_precip <- mask(crop(precip_ksrb, local_boundary), vect(local_boundary)) %>%
    r2df(., 'precip')
  
  mean_precip <- local_precip %>%
    group_by(date) %>%
    summarise(mean_precip = mean(precip))
  
  precip_huc8[,huc8+1] <- mean_precip[,2]
}

names(precip_huc8)[-1] <- huc8_ids


# Adjustment for distance -------------------------------------------------

desoto_centroid <- huc8_boundaries %>%
  filter(HUC8 == '10270104') %>%
  st_centroid(.)

centroid_distances <- tibble(
  huc8 = huc8_boundaries$HUC8,
  distance = as.numeric(st_distance(desoto_centroid, st_centroid(huc8_boundaries))[,1:nrow(huc8_boundaries)])
) %>%
  arrange(huc8)

# Total daily precip by HUC8 per event ------------------------------------

# * 0 day buffer ----------------------------------------------------------

event_precip_huc8_0day <- matrix(nrow = nrow(events), ncol = nrow(huc8_boundaries))

for(event in 1:nrow(events)){
event_precip <- precip_huc8 %>%
  filter(date %within% interval(date(events$start_dateTime[event]), date(events$end_dateTime[event]))) %>%
  pivot_longer(-date, names_to = 'huc8', values_to = 'precip') %>%
  group_by(huc8) %>%
  summarise(event_precip = round(sum(precip), 3))
event_precip_huc8_0day[event,] <- event_precip$event_precip
}

colnames(event_precip_huc8_0day) <- paste0('precip_',event_precip$huc8)

event_precip_huc8_0day <- cbind(events[,1], event_precip_huc8_0day)

# * 1 day buffer ----------------------------------------------------------

event_precip_huc8_1day <- matrix(nrow = nrow(events), ncol = nrow(huc8_boundaries))

for(event in 1:nrow(events)){
  event_precip <- precip_huc8 %>%
    filter(date %within% interval(date(events$start_dateTime[event]) - 1, date(events$end_dateTime[event]) - 1)) %>%
    pivot_longer(-date, names_to = 'huc8', values_to = 'precip') %>%
    group_by(huc8) %>%
    summarise(event_precip = round(sum(precip), 3))
  event_precip_huc8_1day[event,] <- event_precip$event_precip
}

colnames(event_precip_huc8_1day) <- paste0('precip_',event_precip$huc8)

event_precip_huc8_1day <- cbind(events[,1], event_precip_huc8_1day)


# Correlations for HUC8 precips -------------------------------------------

library(reshape)
library(corrplot)
precip_corr <- event_precip_huc8 %>%
  select(-event_number) %>%
  cor(., method = 'pearson', use = 'pairwise.complete.obs')
corrplot(precip_corr, method = 'pie', tl.cex = 0.7, tl.col = 'black', type = 'lower', diag = FALSE, col.lim = c(0, 1))
precip_corr_melt <- melt(precip_corr) %>%
  filter(X1 != X2)
?corrplot
desoto_corrs <- melt(precip_corr) %>%
  filter(X1 == 'precip_10270104')

mean_corrs <- precip_corr_melt %>%
  group_by(X1) %>%
  summarise(mean_corr = mean(value),
            median_corr = median(value))



# Precip corr and distancce from desoto -----------------------------------

desoto_combine <- centroid_distances %>%
  mutate(corr = desoto_corrs$value)
ggplot(data = desoto_combine) +
  geom_point(aes(x = distance, y = corr))

ggplot(data = desoto_combine) +
  geom_point(aes(x = distance, y = distance/75000 + 1))

# Map ---------------------------------------------------------------------

#precip, single date
ggplot() +
  geom_raster(data = subset(precip_ksrb_df, date == '2016-05-27'), aes(x=x, y=y, fill = precip)) +
  geom_sf(data = huc8_boundaries, fill = NA) +
  geom_sf(data = rivers_ksrb, color = 'lightblue') +
  geom_sf(data = reservoirs, fill = 'lightblue') +
  scale_fill_gradientn(colours = terrain.colors(n=12)) +
  coord_sf()

#spei12, single month
ggplot() +
  geom_raster(data = subset(spei12_ksrb_df, date == '2019-05-01'), aes(x=x, y=y, fill = spei12)) +
  geom_sf(data = huc8_boundaries, fill = NA) +
  geom_sf(data = rivers_ksrb, color = 'lightblue') +
  geom_sf(data = reservoirs, fill = 'lightblue') +
  scale_fill_gradient2() +
  coord_sf()


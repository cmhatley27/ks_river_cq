
# Load --------------------------------------------------------------------
library(tidyverse)
library(sf)
library(tmap)
library(terra)
library(FedData)
source('Theme+Settings.R')


# Load common shapefiles --------------------------------------------------
ks <- st_read('./DataFiles/shapefiles/created/state_borders.shp')
ksrb <- st_read('./DataFiles/shapefiles/created/ksrb_boundary.shp')
huc8s <- st_read('./DataFiles/shapefiles/created/huc8_boundaries.shp')
streams <- st_read('./DataFiles/shapefiles/created/rivers_ksrb.shp')
reservoirs <- st_read('./DataFiles/shapefiles/created/reservoirs.shp')
rivers <- filter(streams, str_detect(streams$GNIS_NAME, 'River'))


# Calculate MAP -----------------------------------------------------------
precip_ksrb <- rast('C:/School/SAFE KAW/Data/DataFiles/climate_data/daily_rasters/precip_ksrb.tiff')
precip_ksrb_df <- terra::as.data.frame(precip_ksrb, xy = TRUE) %>%
  pivot_longer(-c(x,y), names_to = 'date', values_to = 'precip') %>%
  mutate(date = ymd(date))

ksrb_map <- precip_ksrb_df %>%
  group_by(year(date), x, y) %>%
  summarise(tap = sum(precip)) %>%
  group_by(x,y) %>%
  summarise(map = mean(tap))

ksrb_map_rast <- rast(ksrb_map, crs = 'EPSG:4269')


tm_shape(ksrb_map_rast) +
  tm_raster(palette = '-viridis', style = 'cont', legend.reverse = TRUE) +
  tm_shape(ksrb) + tm_borders(col = 'black') +
  tm_shape(rivers) + tm_lines(col = 'lightblue1') +
  tm_shape(huc8s) + tm_borders() +
  tm_shape(reservoirs) + tm_dots(col = 'red', size = 0.1)


# Find length of KS river past DeSoto -------------------------------------------------

ks_river <- filter(streams, str_detect(GNIS_NAME, 'Kansas'))
desoto_long <- -94.9646893
desoto_bbox <- st_bbox(c(xmin = desoto_long,
                         st_bbox(ks_river)['xmax'],
                         st_bbox(ks_river)['ymax'],
                         st_bbox(ks_river)['ymin']),
                       crs = st_crs(4269))
ks_river_crop <- st_crop(ks_river, desoto_bbox)

tm_shape(ksrb) + tm_borders() +
  tm_shape(ks_river) + tm_lines() +
  tm_shape(ks_river_crop) +tm_lines(col = 'red')

sum(st_length(ks_river_crop))


# Find length of big blue river downstream of gage ------------------------
bb_river <- filter(streams, str_detect(GNIS_NAME, 'Big Blue'))
tuttle_lat <- 39.2372186
tuttle_bbox <- st_bbox(c(st_bbox(bb_river)['xmin'],
                         st_bbox(bb_river)['xmax'],
                         ymax = tuttle_lat,
                         st_bbox(bb_river)['ymin']),
                       crs = st_crs(4269))
bb_river_crop <- st_crop(bb_river, tuttle_bbox)

tm_shape(ksrb) + tm_borders() +
  tm_shape(bb_river) + tm_lines() +
  tm_shape(bb_river_crop) +tm_lines(col = 'red')

sum(st_length(bb_river_crop))



# Monthly precip by huc8 --------------------------------------------------

precip_huc8 <- read_csv(file.path('DataFiles','hydro_data','precip','precip_huc8.csv'))

monthly_precip <- pivot_longer(precip_huc8, !date, names_to = 'huc8', values_to = 'precip') %>%
  group_by(huc8, year = year(date), month = month(date)) %>%
  summarise(month_precip = sum(precip)) %>%
  group_by(huc8, month) %>%
  summarise(month_mean = mean(month_precip))

ggplot(data = monthly_precip, aes(x = month, y = month_mean)) +
  geom_col() +
  facet_wrap(vars(huc8))

monthly_precip %>%
  group_by(month) %>%
  summarise(p = mean(month_mean)) %>%
  mutate(prop = p/sum(p))


# LULC --------------------------------------------------------------------
lulc <- get_nlcd(ksrb, 'ksrb_lulc', extraction.dir = file.path('DataFiles','land_use')) %>%
  terra::project(., 'epsg:4269') %>%
  terra::mask(., ksrb)
lulc <- rast(file.path('DataFiles', 'land_use', 'ksrb_lulc_NLCD_Land_Cover_2019.tif'))

lulc_summary <- freq(lulc) %>%
  mutate(pct = count/sum(count))


lulc_agg <- terra::aggregate(lulc, fact=10, fun = 'modal')

freq(lulc_agg) %>%
  mutate(pct = count/sum(count))

lulc_map <- tm_shape(lulc_agg) + tm_raster(col = 'Class', drop.levels = TRUE, palette = pal_nlcd()$Color) +
  tm_shape(huc8s) + tm_borders()
lulc_map
tmap_save(lulc_map, filename = file.path('Figures','maps','lulc.tiff'),
          height = 4, width = 9.7, units = 'in', compression = 'lzw')


lulc_agg_df <- terra::as.data.frame(lulc_agg, xy = TRUE)
lulc_xsum <- lulc_agg_df %>%
  mutate(band = cut(x, 12, labels = seq(1:12)))
table(select(lulc_xsum, band, Class))
ggplot(lulc_agg_df) +
  geom_bar(aes(x = x, y = Class))


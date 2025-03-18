
# Load --------------------------------------------------------------------
library(tidyverse)
library(sf)
library(terra)
library(tidyterra)
library(FedData)
library(ggspatial)
source('Theme+Settings.R')


# Load common shapefiles --------------------------------------------------
ksrb <- st_read('./DataFiles/shapefiles/created/ksrb_boundary.shp')
huc8s <- st_read('./DataFiles/shapefiles/created/huc8_boundaries.shp')
streams <- st_read('./DataFiles/shapefiles/created/rivers_ksrb.shp')
reservoirs <- st_read('./DataFiles/shapefiles/created/reservoirs.shp')
rivers <- filter(streams, str_detect(streams$GNIS_NAME, 'River'))
states <- st_read('./DataFiles/shapefiles/States_shapefile.shp') %>%
  filter(!(State_Code %in% c('AK', 'HI')))


# Get extent --------------------------------------------------------------
ggplot() +
  geom_sf(data = states) +
  # geom_spatraster(data = ksrb_map_rast) + 
  # scale_fill_viridis_c(name = 'Mean Annual Precipitation [mm]', na.value = NA, direction = -1, guide = 'none') +
  geom_sf(data = ksrb, color = 'black', fill = NA, linewidth = 1)  +
  # geom_sf(data = extent, color = 'red', fill = NA, linewidth = 1) +
  theme_void()

ggsave(file.path('Figures','final','study area','US_extent.png'), height = 3, width = 6.5, units = 'in')


# Hydrology ---------------------------------------------------------------
river_sel <- c('Kansas|Big Blue|Delaware|Wakarusa|Smoky Hill|Saline|Solomon|Republican')
rivers_fil <- filter(rivers, str_detect(GNIS_NAME, river_sel))
ggplot() +
  geom_sf(data = ksrb, color = 'black', fill = 'grey90')  +
  geom_sf(data = rivers_fil, color = '#466b9f', linewidth = 0.2)  +
  geom_sf(data = reservoirs, fill = '#466b9f', color = 'black', linewidth = 0.1) 
  # annotation_scale() +
  theme_void()
ggsave(file.path('Figures','final','study area','hydro.png'), height = 3, width = 6.5, units = 'in')


# LULC --------------------------------------------------------------------
#download LULC data if not already downloaded
# lulc <- get_nlcd(ksrb, 'ksrb_lulc', extraction.dir = file.path('DataFiles','land_use')) %>%
#   terra::project(., 'epsg:4269') %>%
#   terra::mask(., ksrb)

lulc <- rast(file.path('DataFiles', 'land_use', 'ksrb_lulc_NLCD_Land_Cover_2019.tif')) %>%
  project(., ksrb) %>%
  mask(., ksrb)

lulc_agg <- terra::aggregate(lulc, fact=10, fun = 'modal')
rm(lulc)

huc8s_terra <- vect(huc8s) %>%
  project(., lulc_agg)

freq_table <- freq(lulc_agg, zones = huc8s_terra)
frac <- freq_table %>%
  group_by(zone) %>%
  summarise(crop = sum(count[value == 'Cultivated Crops'])/sum(count),
            grass = sum(count[value %in% c('Pasture/Hay', 'Grassland/Herbaceous')])/sum(count),
            developed = sum(count[str_detect(value, 'Developed')])/sum(count),
            forest = sum(count[str_detect(value, 'Forest')])/sum(count))

huc8s_frac <- cbind(huc8s, frac)

ggplot(data = lulc_agg, aes(fill = value)) +
  geom_spatraster(data = lulc_agg) +
  geom_sf(data = ksrb, fill = NA, color = 'black') +
  theme_void()

levs1 <- levels(lulc_agg)[[1]]
tab1<-coltab(lulc_agg)[[1]]

tab2 <- tab1
#forest
tab2[42:44,2:4] <- rep(c(0,156,0),each = 3)
#grassland
tab2[82,2:4] <- tab1[72,2:4]
#urban
tab2[22:25,2:4] <- tab1[24,2:4]
#water
tab2[12,2:4] <- c(70,107,159)

l2 <- lulc_agg
coltab(l2) <- tab2

ggplot(data = l2, aes(fill = value)) +
  geom_spatraster(data = l2) +
  geom_sf(data = ksrb, fill = NA, color = 'black') +
  geom_sf(data = huc8s, fill = NA, color = 'black', linewidth = 0.3) +
  # geom_sf(data = rivers, fill = NA, color = 'lightblue3', linewidth = 0.3)  +
  # geom_sf(data = reservoirs, fill = 'lightblue3', color = 'black', linewidth = 0.1) +
  theme_void() +
  theme(legend.position = 'none')
ggsave(file.path('Figures','final','study area','lulc_simple.png'), height = 3, width = 6.5, units = 'in')
# ggsave(file.path('Figures','maps','lulc_noleg.png'), height = 2.238, width = 5.035, units = 'in')
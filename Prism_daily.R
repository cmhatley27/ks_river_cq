library(tidyverse)
library(sf)
library(terra)
library(prism)
library(SPEI)


# Download and prepare daily precip data -----------------------------------------------

#get dates vector
all_hydro_daily <- read_csv('./DataFiles/hydro_data/all_hydro_daily.csv')
dates <- all_hydro_daily$date

#start download
prism_set_dl_dir('C:/School/SAFE KAW/Data/DataFiles/climate_data/prism/daily_raw')
get_prism_dailys(type = 'ppt', dates = dates, keepZip = FALSE)


#get file names
dirs <- list.dirs('C:/School/SAFE KAW/Data/DataFiles/climate_data/prism/daily_raw', recursive = FALSE)
files <- list.files(dirs, pattern = '.bil$', full.names = TRUE)

#assemble raster stack
precip_all <- rast(files)

#Crop to full KSRB
ksrb_boundary <- st_read('C:/School/SAFE KAW/Data/DataFiles/shapefiles/created/ksrb_boundary.shp')
precip_crop <- crop(precip_all, ksrb_boundary)
precip_ksrb <- mask(precip_crop, ksrb_boundary)
#shorten field names to just date
names(precip_ksrb) <- ymd(str_sub(names(precip_ksrb), 24, 31))
#Save new raster
writeRaster(precip_ksrb, 'C:/School/SAFE KAW/Data/DataFiles/climate_data/daily_rasters/precip_ksrb.tiff', overwrite = TRUE)

#crop to EKSRB
eksrb_boundary <- st_read('C:/School/SAFE KAW/Data/DataFiles/shapefiles/watershed_bndry.shp')
precip_crop <- crop(precip_all, eksrb_boundary)
precip_eksrb <- mask(precip_crop, eksrb_boundary)
#shorten field names to just date
names(precip_eksrb) <- str_sub(names(precip_eksrb), 24, 31)
#Save new raster
writeRaster(precip_eksrb, 'C:/School/SAFE KAW/Data/DataFiles/climate_data/daily_rasters/precip_eksrb.tiff', overwrite = TRUE)

#crop to EKSRB extension
eksrb_extended_boundary <- st_read('C:/School/SAFE KAW/Data/DataFiles/shapefiles/created/eksrb_extended_boundary.shp')
precip_crop <- crop(precip_all, eksrb_extended_boundary)
precip_eksrb_extended <- mask(precip_crop, eksrb_extended_boundary)
#shorten field names to just date
names(precip_eksrb_extended) <- str_sub(names(precip_eksrb_extended), 24, 31)
#Save new raster
writeRaster(precip_eksrb_extended, 'C:/School/SAFE KAW/Data/DataFiles/climate_data/daily_rasters/precip_eksrb_extended.tiff', overwrite = TRUE)

# Look at data (Full KSRB) ------------------------------------------------

#Load and convert to data frame
precip_ksrb <- rast('C:/School/SAFE KAW/Data/DataFiles/climate_data/daily_rasters/precip_ksrb.tiff')
precip_ksrb_df <- terra::as.data.frame(precip_ksrb, xy = TRUE) %>%
  pivot_longer(-c(x,y), names_to = 'date', values_to = 'precip') %>%
  mutate(date = ymd(date))

#Load map shapefiles
state_borders <- st_read('C:/School/SAFE KAW/Data/DataFiles/shapefiles/created/state_borders.shp')
eksrb_boundary <- st_read('C:/School/SAFE KAW/Data/DataFiles/shapefiles/watershed_bndry.shp')
eksrb_extended_boundary <- st_read('C:/School/SAFE KAW/Data/DataFiles/shapefiles/created/eksrb_extended_boundary.shp')
ksrb_boundary <- st_read('C:/School/SAFE KAW/Data/DataFiles/shapefiles/created/ksrb_boundary.shp')
reservoirs <- st_read('C:/School/SAFE KAW/Data/DataFiles/shapefiles/created/reservoirs.shp')
rivers_ksrb <- st_read('C:/School/SAFE KAW/Data/DataFiles/shapefiles/created/rivers_ksrb.shp')
huc8_boundaries <- st_read('C:/School/SAFE KAW/Data/DataFiles/shapefiles/created/huc8_boundaries.shp')

#Add up all precip for each grid cell
total_precip_bycell <- precip_ksrb_df %>%
  group_by(x,y) %>%
  summarise(avg_annual = sum(precip)/10)
ggplot() +
  geom_raster(data = total_precip_bycell, aes(x=x,y=y, fill = avg_annual)) +
  geom_sf(data = ksrb_boundary, fill = NA) +
  #geom_sf(data = huc8_boundaries, fill = NA) +
  geom_sf(data = rivers_ksrb, color = 'grey70', linewidth = 0.2) +
  geom_sf(data = reservoirs, fill = 'lightblue', color = 'lightblue') +
  scale_fill_viridis_c(direction = -1, name = 'Avg Annual Precip [mm]') +
  theme_void() +
  coord_sf()
ggsave('C:/School/SAFE KAW/Data/Figures/maps/ksrb_annual_precip.tiff', device = 'tiff', width = 10, height = 5.8, units = 'in', compression = 'lzw', dpi = 700)

US_borders <- st_transform(st_read('C:/School/SAFE KAW/Data/DataFiles/shapefiles/States_shapefile.shp'), crs(ksrb_boundary)) %>%
  filter(State_Code != 'AK' &
           State_Code != 'HI')
ksrb_extent <- st_as_sfc(st_bbox(ksrb_boundary))

ggplot() +
  geom_sf(data = US_borders) +
  geom_sf(data = ksrb_extent, color = 'red', fill = NA, linewidth = 1) +
  geom_raster(data = total_precip_bycell, aes(x=x,y=y, fill = avg_annual)) +
  geom_sf(data = ksrb_boundary, fill = NA) +
  scale_fill_viridis_c(direction = -1, guide = 'none') +
  theme_void()
ggsave('C:/School/SAFE KAW/Data/Figures/maps/us_location.tiff', device = 'tiff', width = 10, height = 5.8, units = 'in', compression = 'lzw', dpi = 700)

state_borders <- st_read('C:/School/SAFE KAW/Data/DataFiles/shapefiles/States_shapefile.shp') %>%
  st_transform(., crs(ksrb_boundary)) %>%
  st_cast(., 'MULTILINESTRING') %>%
  st_crop(., xmin = -103.79108, ymin = 38.38057, xmax = -94.54938, ymax = 41.28698)

# Calculate mean daily precip per huc8 ---------------------------------------------------------

#raster to data frame function
r2df <- function(raster, value_name = 'value'){
  df <- terra::as.data.frame(raster, xy = TRUE) %>%
    pivot_longer(-c(x,y), names_to = 'date', values_to = value_name) %>%
    mutate(date = ymd(date))}

# Load in data

#Data frames
all_hydro_15min <- with_tz(read_csv('./DataFiles/hydro_data/all_hydro_15min.csv'), 'America/Chicago')
events <- with_tz(read_csv('DataFiles/hydro_data/event_delineations/BFLOW_events_bfi_threshold_adjusted.csv'), 'America/Chicago') %>%
  select(event_number, start_dateTime, end_dateTime)
event_data <- with_tz(read_csv('./DataFiles/hydro_data/event_data/event_data.csv'), 'America/Chicago')
event_summary_trim <- with_tz(read_csv('./DataFiles/hydro_data/event_summary_trim.csv'), 'America/Chicago') %>%
  mutate(wateryear = factor(wateryear),
         month = factor(month),
         season = factor(season))
event_summary <- with_tz(read_csv('./DataFiles/hydro_data/event_summary.csv'), 'America/Chicago')

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

ggplot() +
  geom_raster(data = subset(precip_ksrb_df, date == 'YYYY-MM-DD'),
              aes(x = x, y = y, fill = precip))

# Mean daily precip by HUC8

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
write_csv(precip_huc8, './DataFiles/hydro_data/precip_huc8.csv')



# Look at data (EKSRB) ------------------------------------------------------------

#Load and convert to data frame
precip_eksrb <- rast('C:/School/SAFE KAW/Data/DataFiles/climate_data/daily_rasters/precip_eksrb.tiff')
precip_eksrb_df <- terra::as.data.frame(precip_eksrb, xy = TRUE) %>%
  pivot_longer(-c(x,y), names_to = 'date', values_to = 'precip') %>%
  mutate(date = ymd(date))

#Load map shapefiles
eksrb_boundary <- st_read('C:/School/SAFE KAW/Data/DataFiles/shapefiles/watershed_bndry.shp')
reservoirs <- st_read('C:/School/SAFE KAW/Data/DataFiles/shapefiles/created/reservoirs.shp')
rivers_eksrb <- st_read('C:/School/SAFE KAW/Data/DataFiles/shapefiles/created/rivers_eksrb.shp')

#Add up all precip for each grid cell
total_precip_bycell <- precip_eksrb_df %>%
  group_by(x,y) %>%
  summarise(total = sum(precip))
ggplot() +
  geom_raster(data = total_precip_bycell, aes(x=x,y=y, fill = total)) +
  geom_sf(data = eksrb_boundary, fill = NA) +
  geom_sf(data = reservoirs, fill = 'lightblue') +
  geom_sf(data = rivers_eksrb, color = 'lightblue') +
  coord_sf()
  

#Daily sequence 
start_date <- '2021-05-15'
end_date <- '2021-05-22'
ggplot() +
  geom_raster(data = subset(precip_eksrb_df, date %within% interval(start_date, end_date)), aes(x=x, y=y, fill = precip)) +
  geom_sf(data = eksrb_boundary, fill = NA) +
  geom_sf(data = reservoirs, fill = 'lightblue') +
  geom_sf(data = rivers_eksrb, color = 'lightblue') +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) +
  coord_sf() +
  facet_wrap(vars(date), nrow = 2)


#my bday
ggplot() +
  geom_raster(data = subset(precip_eksrb_df, month(date) == 5 & day(date) == 27), aes(x=x, y=y, fill = precip)) +
  geom_sf(data = eksrb_boundary, fill = NA) +
  geom_sf(data = reservoirs, fill = 'lightblue') +
  geom_sf(data = rivers_eksrb, color = 'lightblue') +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) +
  coord_sf() +
  facet_wrap(vars(date), nrow = 2)


#rainest day?
total_precip_byday <- precip_eksrb_df %>%
  group_by(date) %>%
  summarise(total = sum(precip))

#precip totals
precip_eksrb_extended <- crop(precip_ksrb, eksrb_extended_boundary) %>%
  mask(., eksrb_extended_boundary) %>%
  terra::as.data.frame(., xy = TRUE) %>%
  pivot_longer(-c(x,y), names_to = 'date', values_to = 'precip') %>%
  mutate(date = ymd(date))

eksrb_extended_total <- sum(precip_eksrb_extended$precip)  
eksrb_total <- sum(precip_eksrb_df$precip)
ksrb_total <- sum(precip_ksrb_df$precip)

eksrb_total/ksrb_total
eksrb_extended_total/ksrb_total



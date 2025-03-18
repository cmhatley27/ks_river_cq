library(tidyverse)
library(sf)
library(terra)
library(prism)


# Download and prepare daily precip data -----------------------------------------------

#get dates vector
all_hydro_daily <- read_csv('./DataFiles/hydro_data/all_hydro_daily.csv')
dates <- all_hydro_daily$date

#download precip
prism_set_dl_dir('./DataFiles/prism_data/prism/daily_raw/precip')
get_prism_dailys(type = 'ppt', dates = dates, keepZip = FALSE)

#get file names
dirs <- list.dirs('./DataFiles/prism_data/prism/daily_raw/precip', recursive = FALSE)
files <- list.files(dirs, pattern = '.bil$', full.names = TRUE)

#assemble raster stack
precip_all <- rast(files)

#Crop to full KSRB
ksrb_boundary <- st_read('./DataFiles/shapefiles/created/ksrb_boundary.shp')
precip_crop <- crop(precip_all, ksrb_boundary)
precip_ksrb <- mask(precip_crop, ksrb_boundary)
#shorten field names to just date
names(precip_ksrb) <- ymd(str_sub(names(precip_ksrb), 24, 31))
#Save new raster
writeRaster(precip_ksrb, './DataFiles/prism_data/daily_rasters/precip_ksrb.tiff', overwrite = TRUE)


# Download and prepare daily mean temperature data ------------------------

#get dates vector
all_hydro_daily <- read_csv('./DataFiles/hydro_data/all_hydro_daily.csv')
dates <- all_hydro_daily$date

#download temp
prism_set_dl_dir('./DataFiles/prism_data/prism/daily_raw/temp')
get_prism_dailys(type = 'tmean', dates = dates, keepZip = FALSE)

#get file names
dirs <- list.dirs('./DataFiles/prism_data/prism/daily_raw/temp', recursive = FALSE)
files <- list.files(dirs, pattern = '.bil$', full.names = TRUE)

#assemble raster stack
temp_all <- rast(files)

#Crop to full KSRB
ksrb_boundary <- st_read('./DataFiles/shapefiles/created/ksrb_boundary.shp')
temp_crop <- crop(temp_all, ksrb_boundary)
temp_ksrb <- mask(temp_crop, ksrb_boundary)
#shorten field names to just date
names(temp_ksrb) <- ymd(str_sub(names(temp_ksrb), 26, 33))
#Save new raster
writeRaster(temp_ksrb, './DataFiles/prism_data/daily_rasters/temp_ksrb.tiff', overwrite = TRUE)


# global mean daily precip ------------------------------------------------
precip <- rast('./DataFiles/prism_data/daily_rasters/precip_ksrb.tiff')
global_mean_precip <- global(precip, 'mean', na.rm = T) %>%
  mutate(date = rownames(.)) %>%
  select(date, precip = mean)
write_csv(global_mean_precip, './DataFiles/hydro_data/precip/precip_global_mean.csv')


# Calculate mean daily precip per huc8 ---------------------------------------------------------

#raster to data frame function
r2df <- function(raster, value_name = 'value'){
  df <- terra::as.data.frame(raster, xy = TRUE) %>%
    pivot_longer(-c(x,y), names_to = 'date', values_to = value_name) %>%
    mutate(date = ymd(date))}

#load huc8 boundaries and precip data
huc8_boundaries <- st_read('./DataFiles/shapefiles/created/huc8_boundaries.shp')

precip_ksrb <- rast('./DataFiles/prism_data/daily_rasters/precip_ksrb.tiff')
precip_ksrb_df <- r2df(precip_ksrb, 'precip')

# calculate means
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
write_csv(precip_huc8, './DataFiles/hydro_data/precip/precip_huc8.csv')
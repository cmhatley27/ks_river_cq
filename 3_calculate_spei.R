# Load libraries ----------------------------------------------------------
library(tidyverse)
library(sf)
library(terra)
library(prism)
library(SPEI) # SPEI package uses monthly climate data, so need to download that


# Set file directories -----------------------------------------------------

#point this to the directory where raw PRISM data should be downloaded into, or
#where the raw files already are
prism_dir <- file.path('DataFiles','prism_data','prism','monthly_raw')

#point this to the shapefile that outlines the study area
boundary_dir <- file.path('DataFiles','shapefiles','created','ksrb_boundary.shp')

#point this to where you want all resulting rasters to save
save_dir <- file.path('DataFiles','prism_data','monthly_rasters')

# download PRISM data -----------------------------------------------------
###This is monthly but can easily be adapted for yearly/daily values by changing
###the prism::get... functions to their yearly/daily counterparts

year_range <- c(1981:2022)

#Set download directory
prism_set_dl_dir(file.path(prism_dir))

#Download needed data - precip, min temp, and max temp
get_prism_monthlys(type = 'ppt', years = year_range, mon = c(1:12), keepZip = FALSE)

get_prism_monthlys(type = 'tmin', years = year_range, mon = c(1:12), keepZip = FALSE)

get_prism_monthlys(type = 'tmax', years = year_range, mon = c(1:12), keepZip = FALSE)


# prep raw data for SPEI calc ---------------------------------------------

#get file names
dirs <- list.dirs(prism_dir, recursive = FALSE)
files <- list.files(dirs, pattern = '.bil$', full.names = TRUE)

#Assemble raster stacks of all dates
precip_all <- rast(subset(files, str_detect(files, 'ppt')))
tmin_all <- rast(subset(files, str_detect(files, 'tmin')))
tmax_all <- rast(subset(files, str_detect(files, 'tmax')))

#Crop to selected boundary shapfile. Boundary needs to be in lat/long crs (epsg:4269)
#to match with PRISM output
boundary <- st_read(boundary_dir) %>% st_transform(., 'epsg:4269')
precip_crop <- mask(crop(precip_all, boundary), boundary)
tmin_crop <- mask(crop(tmin_all, boundary), boundary)
tmax_crop <- mask(crop(tmax_all, boundary), boundary)

#shorten field names to just date to make future steps easier. The date is contained
#in the second to last 'word' of each file's name
names(precip_crop) <- str_sub(word(names(precip_crop), -2, sep = '_'))
names(tmin_crop) <- str_sub(word(names(tmin_crop), -2, sep = '_'))
names(tmax_crop) <- str_sub(word(names(tmax_crop), -2, sep = '_'))

#convert to data frames
precip_crop_df <- terra::as.data.frame(precip_crop, xy = TRUE)  %>%
  pivot_longer(-c(x,y), names_to = 'date', values_to = 'precip') %>%
  filter(!is.na(date))
tmin_crop_df <- terra::as.data.frame(tmin_crop, xy = TRUE)  %>%
  pivot_longer(-c(x,y), names_to = 'date', values_to = 'tmin') %>%
  filter(!is.na(date))
tmax_crop_df <- terra::as.data.frame(tmax_crop, xy = TRUE)  %>%
  pivot_longer(-c(x,y), names_to = 'date', values_to = 'tmax') %>%
  filter(!is.na(date))


# Calculate Hargreaves PET for SPEI ---------------------------------------
#Built-in Hargreaves PET function in the SPEI package outputs data in a weird
#format, so using the 'hydrosystems' package instead. Allows for more flexible
#data inputting which simplifies the whole process

#Comment out the following 2 install lines if already installed
# install.packages("remotes")
# remotes::install_github("tanerumit/hydrosystems")
library(hydrosystems)


#Combine climate data & calculate PET and water balance
climate <- precip_crop_df %>%
  left_join(tmin_crop_df) %>%
  left_join(tmax_crop_df) %>%
  mutate(date = ym(date),
         pet = hargreavesPET(date, (tmin + tmax)/2, tmax - tmin, y),
         clim_balance = precip - pet) %>%
  group_by(x,y) %>% arrange(date, .by_group = TRUE)

#Calculate SPEI. Select desired timescales on the following vector (in months)
scales <- c(12)
speis_out <- tibble(date = climate$date)
for(scale_sel in 1:length(scales)){
  spei <- climate %>%
    group_by(x,y) %>%
    summarise(val = spei(data = clim_balance, scale = scales[scale_sel], verbose = FALSE)$fitted)
  colnames(spei)[3] <- paste0('spei',scales[scale_sel])
  speis_out <- cbind(speis_out, spei)
  print(paste0('scale ', scales[scale_sel], ' finished'))
}

climate_w_spei <- left_join(climate, select(speis_out, 'x','y',date,starts_with('spei')))


# Quick check of results - annual mean SPEI -------------------------------

spei_year_means <- group_by(climate_w_spei, date = year(date)) %>%
  summarise(mean_spei = mean(spei12, na.rm = TRUE))
ggplot(spei_year_means, aes(x = date, y = mean_spei)) +
  geom_hline(yintercept = 0) +
  geom_point() +
  geom_line() +
  ylim(c(-2,2))

# Save SPEI rasters -----------------------------------------------------------

#convert each variable in the climate data frame to a raster and save
for(vars in seq(4,ncol(climate_w_spei))){
  var_df <- climate_w_spei %>%
    select(x,y,date,vars) %>%
    pivot_wider(id_cols = c(x,y), names_from = date, values_from = 4)
  
  var_rast <- rast(var_df, crs = crs(boundary), digits = 4)
  writeRaster(var_rast, paste0(file.path(save_dir),'/',names(climate[vars]),'.tiff'), overwrite = TRUE)
}

# Calculate mean SPEI for each HUC8 ---------------------------------------
#raster to data frame function
r2df <- function(raster, value_name = 'value'){
  df <- terra::as.data.frame(raster, xy = TRUE) %>%
    pivot_longer(-c(x,y), names_to = 'date', values_to = value_name) %>%
    mutate(date = ymd(date))}

#load huc8 boundaries and spei12 data
huc8_boundaries <- st_read('./DataFiles/shapefiles/created/huc8_boundaries.shp')

spei12_ksrb <- rast('./DataFiles/prism_data/monthly_rasters/spei12_ksrb.tiff')
spei12_ksrb_df <- r2df(spei12_ksrb, 'spei12')

# calculate means
spei12_huc8 <- tibble(
  date = as_date(names(spei12_ksrb))
)

for(huc8 in 1:nrow(huc8_boundaries)){
  
  huc8_ids <- huc8_boundaries$HUC8
  
  local_boundary <- huc8_boundaries %>%
    filter(HUC8 == huc8_ids[huc8])
  
  local_spei12 <- mask(crop(spei12_ksrb, local_boundary), vect(local_boundary)) %>%
    r2df(., 'spei12')
  
  mean_spei12 <- local_spei12 %>%
    group_by(date) %>%
    summarise(mean_spei12 = mean(spei12))
  
  spei12_huc8[,huc8+1] <- mean_spei12[,2]
}

names(spei12_huc8)[-1] <- huc8_ids
write_csv(spei12_huc8, './DataFiles/hydro_data/drought/spei12_huc8.csv')
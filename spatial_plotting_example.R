library(tidyverse)
library(sf)
library(terra)

#function to turn rasters with a time dimension into plottable dataframes
#if no time dimension, you don't need this function and can just use terra::as.data.frame()
r2df <- function(raster, value_name = 'value'){
  df <- terra::as.data.frame(raster, xy = TRUE) %>%
    pivot_longer(-c(x,y), names_to = 'date', values_to = value_name) %>%
    mutate(date = ymd(date))}

#load in a raster and convert to data frame
spei12_ksrb <- rast('C:/School/SAFE KAW/Data/DataFiles/climate_data/monthly_rasters/spei12_ksrb.tiff')
spei12_ksrb_df <- r2df(spei12_ksrb, 'spei12')

#load in a shapefile
ksrb_boundary <- st_read('C:/School/SAFE KAW/Data/DataFiles/shapefiles/created/ksrb_boundary.shp')

#plot raster + shapefile
#raster is subsetted because there is a time dimension - can only plot one date at a time
#shapefile is plotted after raster so it appears on top, but needs its fill set to NA so it doesnt conver raster
ggplot() +
  geom_raster(data = subset(spei12_ksrb_df, date = '2016-05-01'), aes(x = x, y = y, fill = spei12)) +
  geom_sf(data = ksrb_boundary, fill = NA, color = 'black')

library(tidyverse)
library(sf)
library(terra)
library(prism)
library(SPEI)
library(trend)


# Download monthly climate data -------------------------------------------


#start download
prism_set_dl_dir('C:/School/SAFE KAW/Data/DataFiles/climate_data/prism/monthly_raw')

get_prism_monthlys(type = 'ppt', years = c(1981:2022), mon = c(1:12), keepZip = FALSE)

get_prism_monthlys(type = 'tmin', years = c(1981:2022), mon = c(1:12), keepZip = FALSE)

get_prism_monthlys(type = 'tmax', years = c(1981:2022), mon = c(1:12), keepZip = FALSE)


#get file names
dirs <- list.dirs('C:/School/SAFE KAW/Data/DataFiles/climate_data/prism/monthly_raw', recursive = FALSE)
files <- list.files(dirs, pattern = '.bil$', full.names = TRUE)

#Assemble raster stacks
precip_all <- rast(subset(files, str_detect(files, 'ppt')))
tmin_all <- rast(subset(files, str_detect(files, 'tmin')))
tmax_all <- rast(subset(files, str_detect(files, 'tmax')))

#Crop to full KSRB
ksrb_boundary <- st_read('C:/School/SAFE KAW/Data/DataFiles/shapefiles/created/ksrb_boundary.shp')
precip_ksrb <- mask(crop(precip_all, ksrb_boundary), ksrb_boundary)
tmin_ksrb <- mask(crop(tmin_all, ksrb_boundary), ksrb_boundary)
tmax_ksrb <- mask(crop(tmax_all, ksrb_boundary), ksrb_boundary)

#shorten field names to just date
names(precip_ksrb) <- str_sub(names(precip_ksrb), 24, 29)
names(tmin_ksrb) <- str_sub(names(tmin_ksrb), 25, 30)
names(tmax_ksrb) <- str_sub(names(tmax_ksrb), 25, 30)

#convert to data frames
precip_ksrb_df <- terra::as.data.frame(precip_ksrb, xy = TRUE)  %>%
  pivot_longer(-c(x,y), names_to = 'date', values_to = 'precip') %>%
  mutate(date = ym(date)) %>%
  filter(!is.na(date))

tmin_ksrb_df <- terra::as.data.frame(tmin_ksrb, xy = TRUE)  %>%
  pivot_longer(-c(x,y), names_to = 'date', values_to = 'tmin') %>%
  mutate(date = ym(date)) %>%
  filter(!is.na(date))
  
tmax_ksrb_df <- terra::as.data.frame(tmax_ksrb, xy = TRUE)  %>%
  pivot_longer(-c(x,y), names_to = 'date', values_to = 'tmax') %>%
  mutate(date = ym(date)) %>%
  filter(!is.na(date))

#Calculate Hargreaves PET
tmin_wide <- tmin_ksrb_df %>%
  pivot_wider(id_cols = date, names_from = c(x,y), values_from = tmin) %>%
  select(!date)

tmax_wide <- tmax_ksrb_df %>%
  pivot_wider(id_cols = date, names_from = c(x,y), values_from = tmax) %>%
  select(!date)
unique(tmin_ksrb_df$date)
lats <- as.numeric(word(names(tmin_wide), 2, sep = '_'))

pet <- hargreaves(tmin_wide, tmax_wide, lat = lats)
pet_df <- as_data_frame(pet) %>%
  mutate(date = unique(tmin_ksrb_df$date)) %>%
  pivot_longer(!date, names_to = 'coord', values_to = 'pet') %>%
  select(coord,date,pet)


#Combine climate data & calculate spei
climate_monthly_ksrb <- precip_ksrb_df %>%
  left_join(tmin_ksrb_df) %>%
  left_join(tmax_ksrb_df) %>%
  mutate(pet = pet_df$pet) %>%
  mutate(clim_balance = precip - pet) %>%
  group_by(x,y) %>%
  mutate(spei1 = spei(data = clim_balance, scale = 1, verbose = FALSE)$fitted,
         spei3 = spei(data = clim_balance, scale = 3, verbose = FALSE)$fitted,
         spei6 = spei(data = clim_balance, scale = 6, verbose = FALSE)$fitted,
         spei9 = spei(data = clim_balance, scale = 9, verbose = FALSE)$fitted,
         spei12 = spei(data = clim_balance, scale = 12, verbose = FALSE)$fitted,
         spei18 = spei(data = clim_balance, scale = 18, verbose = FALSE)$fitted,
         spei24 = spei(data = clim_balance, scale = 24, verbose = FALSE)$fitted) %>%
  ungroup()

#convert to rasters
eksrb_boundary <- st_read('C:/School/SAFE KAW/Data/DataFiles/shapefiles/watershed_bndry.shp')

for(vars in seq(4,ncol(climate_monthly_ksrb))){
var_df <- climate_monthly_ksrb %>%
  select(x,y,date,vars) %>%
  pivot_wider(id_cols = c(x,y), names_from = date, values_from = 4)

var_rast <- rast(var_df, crs = crs(eksrb_boundary), digits = 4)
writeRaster(var_rast, paste0('C:/School/SAFE KAW/Data/DataFiles/climate_data/monthly_rasters/',names(climate_monthly_ksrb[vars]),'_ksrb.tiff'), overwrite = TRUE)
}



# Average each SPEI per HUC8 ----------------------------------------------

#raster to data frame function
r2df <- function(raster, value_name = 'value'){
  df <- terra::as.data.frame(raster, xy = TRUE) %>%
    pivot_longer(-c(x,y), names_to = 'date', values_to = value_name) %>%
    mutate(date = ymd(date))}

#Load data
huc8_boundaries <- st_read('C:/School/SAFE KAW/Data/DataFiles/shapefiles/created/huc8_boundaries.shp')

spei1_ksrb <- rast('C:/School/SAFE KAW/Data/DataFiles/climate_data/monthly_rasters/spei1_ksrb.tiff')
spei1_ksrb_df <- r2df(spei1_ksrb, 'spei1')

spei3_ksrb <- rast('C:/School/SAFE KAW/Data/DataFiles/climate_data/monthly_rasters/spei3_ksrb.tiff')
spei3_ksrb_df <- r2df(spei3_ksrb, 'spei3')

spei6_ksrb <- rast('C:/School/SAFE KAW/Data/DataFiles/climate_data/monthly_rasters/spei6_ksrb.tiff')
spei6_ksrb_df <- r2df(spei6_ksrb, 'spei6')

spei9_ksrb <- rast('C:/School/SAFE KAW/Data/DataFiles/climate_data/monthly_rasters/spei9_ksrb.tiff')
spei9_ksrb_df <- r2df(spei9_ksrb, 'spei9')

spei12_ksrb <- rast('C:/School/SAFE KAW/Data/DataFiles/climate_data/monthly_rasters/spei12_ksrb.tiff')
spei12_ksrb_df <- r2df(spei12_ksrb, 'spei12')

spei18_ksrb <- rast('C:/School/SAFE KAW/Data/DataFiles/climate_data/monthly_rasters/spei18_ksrb.tiff')
spei18_ksrb_df <- r2df(spei18_ksrb, 'spei18')

spei24_ksrb <- rast('C:/School/SAFE KAW/Data/DataFiles/climate_data/monthly_rasters/spei24_ksrb.tiff')
spei24_ksrb_df <- r2df(spei24_ksrb, 'spei24')

# Mean SPEI1

spei1_huc8 <- tibble(
  date = as_date(names(spei1_ksrb))
)

for(huc8 in 1:nrow(huc8_boundaries)){
  
  huc8_ids <- huc8_boundaries$HUC8
  
  local_boundary <- huc8_boundaries %>%
    filter(HUC8 == huc8_ids[huc8])
  
  local_spei1 <- mask(crop(spei1_ksrb, local_boundary), vect(local_boundary)) %>%
    r2df(., 'spei1')
  
  mean_spei1 <- local_spei1 %>%
    group_by(date) %>%
    summarise(mean_spei1 = mean(spei1))
  
  spei1_huc8[,huc8+1] <- mean_spei1[,2]
}

names(spei1_huc8)[-1] <- huc8_ids
write_csv(spei1_huc8, './DataFiles/hydro_data/drought/spei1_huc8.csv')

# Mean SPEI3

spei3_huc8 <- tibble(
  date = as_date(names(spei3_ksrb))
)

for(huc8 in 1:nrow(huc8_boundaries)){
  
  huc8_ids <- huc8_boundaries$HUC8
  
  local_boundary <- huc8_boundaries %>%
    filter(HUC8 == huc8_ids[huc8])
  
  local_spei3 <- mask(crop(spei3_ksrb, local_boundary), vect(local_boundary)) %>%
    r2df(., 'spei3')
  
  mean_spei3 <- local_spei3 %>%
    group_by(date) %>%
    summarise(mean_spei3 = mean(spei3))
  
  spei3_huc8[,huc8+1] <- mean_spei3[,2]
}

names(spei3_huc8)[-1] <- huc8_ids
write_csv(spei3_huc8, './DataFiles/hydro_data/drought/spei3_huc8.csv')

# Mean SPEI6

spei6_huc8 <- tibble(
  date = as_date(names(spei6_ksrb))
)

for(huc8 in 1:nrow(huc8_boundaries)){
  
  huc8_ids <- huc8_boundaries$HUC8
  
  local_boundary <- huc8_boundaries %>%
    filter(HUC8 == huc8_ids[huc8])
  
  local_spei6 <- mask(crop(spei6_ksrb, local_boundary), vect(local_boundary)) %>%
    r2df(., 'spei6')
  
  mean_spei6 <- local_spei6 %>%
    group_by(date) %>%
    summarise(mean_spei6 = mean(spei6))
  
  spei6_huc8[,huc8+1] <- mean_spei6[,2]
}

names(spei6_huc8)[-1] <- huc8_ids
write_csv(spei6_huc8, './DataFiles/hydro_data/drought/spei6_huc8.csv')

# Mean SPEI9

spei9_huc8 <- tibble(
  date = as_date(names(spei9_ksrb))
)

for(huc8 in 1:nrow(huc8_boundaries)){
  
  huc8_ids <- huc8_boundaries$HUC8
  
  local_boundary <- huc8_boundaries %>%
    filter(HUC8 == huc8_ids[huc8])
  
  local_spei9 <- mask(crop(spei9_ksrb, local_boundary), vect(local_boundary)) %>%
    r2df(., 'spei9')
  
  mean_spei9 <- local_spei9 %>%
    group_by(date) %>%
    summarise(mean_spei9 = mean(spei9))
  
  spei9_huc8[,huc8+1] <- mean_spei9[,2]
}

names(spei9_huc8)[-1] <- huc8_ids
write_csv(spei9_huc8, './DataFiles/hydro_data/drought/spei9_huc8.csv')

# Mean SPEI12

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

# Mean SPEI18

spei18_huc8 <- tibble(
  date = as_date(names(spei18_ksrb))
)

for(huc8 in 1:nrow(huc8_boundaries)){
  
  huc8_ids <- huc8_boundaries$HUC8
  
  local_boundary <- huc8_boundaries %>%
    filter(HUC8 == huc8_ids[huc8])
  
  local_spei18 <- mask(crop(spei18_ksrb, local_boundary), vect(local_boundary)) %>%
    r2df(., 'spei18')
  
  mean_spei18 <- local_spei18 %>%
    group_by(date) %>%
    summarise(mean_spei18 = mean(spei18))
  
  spei18_huc8[,huc8+1] <- mean_spei18[,2]
}

names(spei18_huc8)[-1] <- huc8_ids
write_csv(spei18_huc8, './DataFiles/hydro_data/drought/spei18_huc8.csv')

# Mean SPEI24

spei24_huc8 <- tibble(
  date = as_date(names(spei24_ksrb))
)

for(huc8 in 1:nrow(huc8_boundaries)){
  
  huc8_ids <- huc8_boundaries$HUC8
  
  local_boundary <- huc8_boundaries %>%
    filter(HUC8 == huc8_ids[huc8])
  
  local_spei24 <- mask(crop(spei24_ksrb, local_boundary), vect(local_boundary)) %>%
    r2df(., 'spei24')
  
  mean_spei24 <- local_spei24 %>%
    group_by(date) %>%
    summarise(mean_spei24 = mean(spei24))
  
  spei24_huc8[,huc8+1] <- mean_spei24[,2]
}

names(spei24_huc8)[-1] <- huc8_ids
write_csv(spei24_huc8, './DataFiles/hydro_data/drought/spei24_huc8.csv')

# Plotting ----------------------------------------------------------------

#May PET across KSRB for 2018 and 2019
ggplot() +
  geom_raster(data = subset(climate_monthly_ksrb, date == '2018-05-01'), aes(x=x, y=y, fill = pet)) +
  geom_sf(data = ksrb_boundary, fill = NA, color = 'red') +
  scale_fill_gradientn(colours = terrain.colors(n = 12)) +
  coord_sf()
ggplot() +
  geom_raster(data = subset(climate_monthly_ksrb, date == '2019-05-01'), aes(x=x, y=y, fill = pet)) +
  geom_sf(data = ksrb_boundary, fill = NA, color = 'red') +
  scale_fill_gradientn(colours = terrain.colors(n = 12)) +
  coord_sf()

#Avg monthly PET across KSRB 2013-2022
avg_monthly_pet <- climate_monthly_ksrb %>%
  group_by(x,y,month = month(date)) %>%
  summarise(mean_pet = mean(pet))
ggplot() +
  geom_raster(data = avg_monthly_pet, aes(x=x,y=y,fill = mean_pet)) +
  geom_sf(data = ksrb_boundary, fill = NA, color = 'black') +
  scale_fill_gradientn(colours = terrain.colors(n = 12)) +
  coord_sf() +
  facet_wrap(vars(month))


#SPEI for a single grid cell (lawrence)
lawrence_climate_monthly <- subset(climate_monthly_ksrb, near(x,-95.2353, tol = 0.02) & near(y, 38.9717, tol = 0.02)) %>%
  pivot_longer(!c(x,y,date), names_to = 'var', values_to = 'value')

ggplot(data = subset(lawrence_climate_monthly, str_detect(var, 'spei'))) +
  geom_hline(yintercept = 0) +
  geom_line(aes(x = date, y = value, color = factor(var, levels = c('spei1', 'spei3', 'spei6', 'spei12'))), linewidth = 1) +
  scale_color_discrete(guide = 'none') +
  scale_x_date(date_breaks = '4 years', date_labels = '%Y') +
  coord_cartesian(ylim = c(-3,3)) +
  facet_wrap(vars(factor(var, levels = c('spei1', 'spei3', 'spei6', 'spei12'))))

colby_climate_monthly <- subset(climate_monthly_ksrb, near(x,-101.0525, tol = 0.02085) & near(y, 39.3958, tol = 0.02085)) %>%
  pivot_longer(!c(x,y,date), names_to = 'var', values_to = 'value')

ggplot(data = subset(colby_climate_monthly, str_detect(var, 'spei'))) +
  geom_hline(yintercept = 0) +
  geom_line(aes(x = date, y = value, color = factor(var, levels = c('spei1', 'spei3', 'spei6', 'spei12'))), linewidth = 1) +
  scale_color_discrete(guide = 'none') +
  scale_x_date(date_breaks = '4 years', date_labels = '%Y') +
  coord_cartesian(ylim = c(-3,3)) +
  facet_wrap(vars(factor(var, levels = c('spei1', 'spei3', 'spei6', 'spei12'))))


#SPEI trend
spei12_trend <- climate_monthly_ksrb %>%
  filter(!is.na(spei12)) %>%
  group_by(x,y) %>%
  summarize(trend = mk.test(spei12)$statistic)

ggplot(data = spei12_trend) +
  geom_raster(aes(x=x,y=y,fill=trend)) +
  scale_fill_gradient2() +
  coord_quickmap()

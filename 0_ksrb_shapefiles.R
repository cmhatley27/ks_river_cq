library(tidyverse)
library(sf)

#HUC 8 Boundaries
huc8_all <- st_read('C:/School/SAFE KAW/Data/DataFiles/shapefiles/HUC8_US.shp') %>%
  st_transform(., st_crs('EPSG:4269'))
huc8_boundaries <- huc8_all %>%
  filter(str_detect(huc8_all$HUC8, '^1025') | str_detect(huc8_all$HUC8, '^1026') | str_detect(huc8_all$HUC8, '^1027')) %>%
  mutate(HUC6 = str_sub(HUC8,1,6),
         HUC4 = str_sub(HUC8,1,4))
st_write(huc8_boundaries, 'C:/School/SAFE KAW/Data/DataFiles/shapefiles/created/huc8_boundaries.shp', delete_layer = TRUE)

#KS River Basin Boundary
ksrb_boundary <- huc8_boundaries %>%
  st_union(.)
st_write(ksrb_boundary, 'C:/School/SAFE KAW/Data/DataFiles/shapefiles/created/ksrb_boundary.shp', delete_layer = TRUE)

#State Borders
state_borders <- st_read('C:/School/SAFE KAW/Data/DataFiles/shapefiles/States_shapefile.shp') %>%
  st_transform(., crs(ksrb_boundary)) %>%
  st_cast(., 'MULTILINESTRING') %>%
  st_crop(., xmin = -103.79108, ymin = 38.38057, xmax = -94.54938, ymax = 41.28698) %>%
  filter(State_Code == 'KS')
st_write(state_borders, 'C:/School/SAFE KAW/Data/DataFiles/shapefiles/created/state_borders.shp', delete_layer = TRUE)

#Reservoirs
waterbodies <- st_read('C:/School/SAFE KAW/Data/DataFiles/shapefiles/NHDWaterbody.shp')

reservoirs <- waterbodies %>%
  #filter(AREASQKM > 29) %>%
  #st_intersection(., ksrb_boundary)
  filter(GNIS_NAME == 'Clinton Lake' |
         GNIS_NAME == 'Perry Lake' & AREASQKM > 1 |
         GNIS_NAME == 'Tuttle Creek Lake' |
         GNIS_NAME == 'Milford Lake' |
           GNIS_NAME == 'Waconda Lake' |
           GNIS_NAME == 'Wilson Lake' & AREASQKM > 1|
           GNIS_NAME == 'Kanopolis Lake')

st_write(reservoirs, 'C:/School/SAFE KAW/Data/DataFiles/shapefiles/created/reservoirs.shp', delete_layer = TRUE)

#River
ks_streams <- st_read('C:/School/SAFE KAW/Data/DataFiles/shapefiles/kansas_streams_order3_9_webmerc.shp') %>%
  st_transform(., crs(ksrb_boundary))
ks_rivers <- ks_streams %>%
  filter(STRAHLER > 3) %>%
  st_intersection(., ksrb_boundary)

ne_rivers <- st_read('C:/School/SAFE KAW/Data/DataFiles/shapefiles/HYDRO_MajorStream_DNR.shp') %>%
  st_transform(., crs(ksrb_boundary)) %>%
  st_intersection(., ksrb_boundary) %>%
  filter(HUC_LVL < 17)

ks_join <- ks_rivers %>%
  select(GNIS_NAME, REACHCODE, geometry)
ne_join <- ne_rivers %>%
  select(GNIS_NAME = GNIS_Name, REACHCODE = ReachCode, geometry)
rivers_ksrb <- rbind(ks_join, ne_join) %>%
  st_transform(., crs(ksrb_boundary))

st_write(rivers_ksrb, 'C:/School/SAFE KAW/Data/DataFiles/shapefiles/created/rivers_ksrb.shp', delete_layer = TRUE)

#KS River outlet:
outlet <- st_sfc(st_point(c(-94.6091232, 39.0900045)), crs = crs(ksrb_boundary))
st_write(outlet, 'C:/School/SAFE KAW/Data/DataFiles/shapefiles/created/ksrb_outlet.shp', delete_layer = TRUE)

outlet <- st_read('C:/School/SAFE KAW/Data/DataFiles/shapefiles/created/ksrb_outlet.shp')
library(tidyverse)
library(sf)
library(terra)


#EKSRB Boundary 
eksrb_boundary <- st_read('C:/School/SAFE KAW/Data/DataFiles/shapefiles/watershed_bndry.shp')

#HUC 8 Boundaries
huc8_all <- st_read('C:/School/SAFE KAW/Data/DataFiles/shapefiles/HUC8_CONUS/HUC8_US.shp')  %>%
  st_transform(., crs(eksrb_boundary))
huc8_boundaries <- huc8_all %>%
  filter(str_detect(huc8_all$HUC8, '^1025') | str_detect(huc8_all$HUC8, '^1026') | str_detect(huc8_all$HUC8, '^1027')) %>%
  mutate(HUC6 = str_sub(HUC8,1,6),
         HUC4 = str_sub(HUC8,1,4))
st_write(huc8_boundaries, 'C:/School/SAFE KAW/Data/DataFiles/shapefiles/created/huc8_boundaries.shp', delete_layer = TRUE)

huc8_names <- tibble(huc8 = huc8_boundaries$HUC8,
                     huc8name = huc8_boundaries$NAME)

huc8_centroids <- tibble(
  huc8 = huc8_boundaries$HUC8,
  x = st_coordinates(st_centroid(huc8_boundaries))[,1],
  y = st_coordinates(st_centroid(huc8_boundaries))[,2]) %>%
  arrange(x)
write_csv(huc8_centroids, 'C:/School/SAFE KAW/Data/DataFiles/shapefiles/created/huc8_centroids.csv')


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

#EKSRB extension
huc8_selection <- as.character(c(10270104, 10270103, 10270102, 10270101, 10260008, 10250017, 10270205, 10260006, 10260010, 10260015, 10250016, 10270207, 10270206, 10270202, 10270204, 10270203, 	10270201))
eksrb_extended_boundary <- huc8_boundaries %>%
  filter(HUC8 %in% huc8_selection) %>%
  st_union(.)
st_write(eksrb_extended_boundary, 'C:/School/SAFE KAW/Data/DataFiles/shapefiles/created/eksrb_extended_boundary.shp', delete_layer = TRUE)


#HUC 12 Boundaries
huc12_ks <- st_read('C:/School/SAFE KAW/Data/DataFiles/shapefiles/WBDHU12mod.shp') %>%
  select(!huc12num)

huc12_ne <- st_read('C:/School/SAFE KAW/Data/DataFiles/shapefiles/HYDRO_WBDHU12_USGS.shp') %>%
  st_transform(., crs(huc12_ks)) %>%
  rename(referenceg = GNIS_ID) %>%
  select(!OBJECTID)

names(huc12_ks) <- tolower(names(huc12_ks))
names(huc12_ne) <- tolower(names(huc12_ne))

huc12_all <- rbind(huc12_ks, huc12_ne)

huc12_boundaries <- huc12_all %>%
  filter(str_detect(huc12_all$huc12, '^1025') | str_detect(huc12_all$huc12, '^1026') | str_detect(huc12_all$huc12, '^1027')) %>%
  mutate(huc8 = str_sub(huc12, 1, 8)) %>%
  left_join(huc8_names, by = 'huc8')


ggplot() +
  geom_sf(data = ksrb_boundary, fill = NA, color = 'yellow') +
  geom_sf(data = huc12_boundaries, fill = NA) +
  geom_sf(data = huc8_boundaries, fill = NA, linewidth = 1)


st_write(huc12_boundaries, 'C:/School/SAFE KAW/Data/DataFiles/shapefiles/created/huc12_boundaries.shp', delete_layer = TRUE)

huc12_boundaries_eksrb <- huc12_boundaries[st_intersects(eksrb_boundary, st_centroid(huc12_boundaries))[[1]],]
st_write(huc12_boundaries_eksrb, 'C:/School/SAFE KAW/Data/DataFiles/shapefiles/created/huc12_boundaries_eksrb.shp', delete_layer = TRUE)

huc12_boundaries_eksrb_extended <- huc12_boundaries[st_intersects(eksrb_extended_boundary, st_centroid(huc12_boundaries))[[1]],]
st_write(huc12_boundaries_eksrb_extended, 'C:/School/SAFE KAW/Data/DataFiles/shapefiles/created/huc12_boundaries_eksrb_extended.shp', delete_layer = TRUE)

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

rivers_eksrb <- rivers_ksrb %>%
  st_intersection(., eksrb_boundary)
st_write(rivers_eksrb, 'C:/School/SAFE KAW/Data/DataFiles/shapefiles/created/rivers_eksrb.shp', delete_layer = TRUE)

#Point
#KS River outlet:
outlet <- st_sfc(st_point(c(-94.6091232, 39.0900045)), crs = crs(ksrb_boundary))
st_write(outlet, 'C:/School/SAFE KAW/Data/DataFiles/shapefiles/created/ksrb_outlet.shp', delete_layer = TRUE)

outlet <- st_read('C:/School/SAFE KAW/Data/DataFiles/shapefiles/created/ksrb_outlet.shp')

# Load all in and map!! -----------------------------------------------------
state_borders <- st_read('C:/School/SAFE KAW/Data/DataFiles/shapefiles/created/state_borders.shp')
eksrb_boundary <- st_read('C:/School/SAFE KAW/Data/DataFiles/shapefiles/watershed_bndry.shp')
ksrb_boundary <- st_read('C:/School/SAFE KAW/Data/DataFiles/shapefiles/created/ksrb_boundary.shp')
eksrb_extended_boundary <- st_read('C:/School/SAFE KAW/Data/DataFiles/shapefiles/created/eksrb_extended_boundary.shp')
huc8_boundaries <- st_read('C:/School/SAFE KAW/Data/DataFiles/shapefiles/created/huc8_boundaries.shp')
huc12_boundaries <- st_read('C:/School/SAFE KAW/Data/DataFiles/shapefiles/created/huc12_boundaries.shp')
huc12_boundaries_eksrb <- st_read('C:/School/SAFE KAW/Data/DataFiles/shapefiles/created/huc12_boundaries_eksrb.shp')
huc12_boundaries_eksrb_extended <- st_read('C:/School/SAFE KAW/Data/DataFiles/shapefiles/created/huc12_boundaries_eksrb_extended.shp')
reservoirs <- st_read('C:/School/SAFE KAW/Data/DataFiles/shapefiles/created/reservoirs.shp')
rivers_ksrb <- st_read('C:/School/SAFE KAW/Data/DataFiles/shapefiles/created/rivers_ksrb.shp')
rivers_eksrb <- st_read('C:/School/SAFE KAW/Data/DataFiles/shapefiles/created/rivers_eksrb.shp')

st_area(ksrb_boundary)/1000000
sum(huc8_boundaries$AREASQKM)
#Full KSRB
ggplot() +
  geom_sf(data = huc8_boundaries, aes(fill = HUC6)) +
  #geom_sf_label(data = huc8_boundaries, aes(label = HUC8)) +
  geom_sf(data = ksrb_boundary, fill = NA, color = 'black', linewidth = 1.2) +
  #geom_sf(data = rivers_ksrb, color = 'lightblue') +
  geom_sf(data = reservoirs, fill = 'lightblue') +
  #geom_sf(data = st_centroid(huc8_boundaries)) +
  geom_sf(data = outlet, color = 'orange') +
  coord_sf()

#EKSRB
ggplot() +
  geom_sf(data = huc12_boundaries_eksrb, aes(fill = huc8name)) +
  geom_sf(data = rivers_eksrb, color = 'lightblue') +
  geom_sf(data = reservoirs, fill = 'lightblue') +
  scale_fill_discrete(name = 'Catchment Name') +
  coord_sf()

#EKSRB_Extended
ggplot() +
  geom_sf(data = huc12_boundaries_eksrb_extended, aes(fill = huc8name), color = NA) +
  geom_sf(data = st_intersection(rivers_ksrb, eksrb_extended_boundary), color = 'lightblue') +
  geom_sf(data = reservoirs, fill = 'lightblue') +
  scale_fill_discrete(name = 'Catchment Name') +
  coord_sf()

#Area of EKSRB
(sum(huc12_boundaries_eksrb$areasqkm))
(sum(huc12_boundaries_eksrb$areasqkm))/(sum(huc12_boundaries$areasqkm))

sum(huc12_boundaries_eksrb_extended$areasqkm)
sum(huc12_boundaries_eksrb_extended$areasqkm)/sum(huc12_boundaries$areasqkm)

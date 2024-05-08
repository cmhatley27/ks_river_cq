
# Load libraries and data ----------------------------------------------------------
library(tidyverse)
library(cluster)  
library(factoextra)
library(reshape)
library(brunnermunzel)
library(ggpubr)
source('load_and_combine.R')


quantize <- function(x, q){
  if(q[1] == q[2]){
    return(cut(unlist(x), 
               quantile(unlist(x), c(0, q[1], 1)), 
               labels = c('low', 'high'),
               include.lowest = TRUE))
  } else{
    return(cut(unlist(x), 
               quantile(unlist(x), c(0, q, 1)), 
               labels = c('low', 'mid', 'high'),
               include.lowest = TRUE))
  }
}

# Calculate water balance data -----------------------------------------------------

#EKSRB area in m^2
library(sf)
huc8_boundaries <- st_read(file.path('DataFiles','shapefiles','created','huc8_boundaries.shp')) 
huc8_areas <- tibble(HUC8 = huc8_boundaries$HUC8,
                     area_sqkm = huc8_boundaries$AREASQKM) %>%
  filter(HUC8 %in% c('10260008', '10260010', '10260015', '10270101', '10270102', '10270104'))
area <- sum(huc8_areas$area_sqkm)*1000000

#water balance in mm
water_bal <- event_summary_cluster %>%
  #select(FI, HI, cluster, total_Q, p_unint_vol, r_sum_1d_total, reservoir_runoff_ratio, reservoir_precip_ratio, delta_N, loadrate) %>%
  mutate(p_depth = p_unint_vol/area*1000,
         bal = (r_sum_1d_total + p_unint_vol - total_Q)/area*1000,
         bal_onlyprecip = (p_unint_vol - total_Q)/area*1000,
         bal_overprecip = bal/(p_unint_vol/area*1000),
         bal_overinputs = bal/((p_unint_vol+r_sum_1d_total)/area*1000),
         runoff_ratio = event_summary_cluster$runoff_precip_ratio)
           


# Load libraries ----------------------------------------------------------
library(tidyverse)
library(sf)
library(tmap)
source('8_combine_all_event_data.R')
source('Theme+Settings.R')

event_summary <- select(event_summary, !contains('global'))
# calculate correlations --------------------------------------------------------
huc8s <- st_read('./DataFiles/shapefiles/created/huc8_boundaries.shp')

climate_vars <- colnames(select(event_summary, 
                                starts_with('d_') & ends_with(c('spei12')),
                                # ends_with('_ante') & contains(c('30d')), #remove short term climate, insignificant
                                starts_with('p_') & ends_with('1d_total')))
climate_var_order <- c('1d_total', 'spei12')

res_var <- 'reservoir_runoff_ratio'
res_name <- 'Res/Runoff'

response_vars <- c('FI_n')

sig_level <- 0.05

dat <- select(event_summary, c(all_of(c(response_vars, res_var, climate_vars)))) %>%
  mutate(res_end = cut(.[[res_var]], breaks = c(0,0.5,1.5), labels = c('low','high'))) %>%
  pivot_longer(all_of(climate_vars), names_to = 'climate_var', values_to = 'climate_val')

cor_calcs <- function(x, clim, sig, res = NULL, res_side = c('low','mid','high')){
  if(is.null(res)){
    cor_stat <- cor.test(x,clim, use = 'pairwise.complete.obs', method = 'spearman')$estimate
    cor_p <- cor.test(x,clim, use = 'pairwise.complete.obs', method = 'spearman')$p.value
  }
  else{
    cor_stat <- cor.test(x[res == res_side],clim[res == res_side], use = 'pairwise.complete.obs', method = 'spearman')$estimate
    cor_p <- cor.test(x[res == res_side],clim[res == res_side], use = 'pairwise.complete.obs', method = 'spearman')$p.value
  }
  return(data.frame(stat = cor_stat,
                    sig = ifelse(cor_p <= sig, TRUE, FALSE)))
}

cors <- dat %>%
  group_by(climate_var) %>%
  summarise(across(all_of(response_vars),
                   ~cor_calcs(.x, climate_val, sig_level),
                   .names = '{.col}', .unpack = TRUE),
            across(all_of(response_vars),
                   ~cor_calcs(.x, climate_val, sig_level, res_end, res_side = 'low'),
                   .names = '{.col}_precip', .unpack = TRUE),
            across(all_of(response_vars),
                   ~cor_calcs(.x, climate_val, sig_level, res_end, res_side = 'high'),
                   .names = '{.col}_res', .unpack = TRUE)) %>%
  mutate(HUC8 = str_sub(climate_var, 3, 10),
         var = str_sub(climate_var, 12, -1)) %>%
  select(climate_var, HUC8, var, everything())

cor_stats <- select(cors, c(climate_var, HUC8, var, ends_with('stat')))
cor_sigs <- select(cors, c(climate_var, HUC8, var, ends_with('sig')))

huc8_dat <- list()
for(element in 1:(length(colnames(cor_stats))-3)){
  huc8_dat[[element]] <- pivot_wider(cor_stats, id_cols = HUC8, names_from = var, values_from = colnames(cor_stats)[element+3]) %>%
    left_join(huc8s, .) %>%
    pivot_longer(!colnames(huc8s), names_to = 'var', values_to = 'val') %>%
    mutate(var = factor(var, levels = climate_var_order))
}
names(huc8_dat) <- colnames(cor_stats)[-c(1:3)]

huc8_sigs <- list()
for(element in 1:(length(colnames(cor_sigs))-3)){
  huc8_sigs[[element]] <- pivot_wider(cor_sigs, id_cols = HUC8, names_from = var, values_from = colnames(cor_sigs)[element+3]) %>%
    left_join(huc8s, .) %>%
    pivot_longer(!colnames(huc8s), names_to = 'var', values_to = 'val') %>%
    mutate(var = factor(var, levels = climate_var_order)) %>%
    st_centroid(.)
}
names(huc8_sigs) <- colnames(cor_sigs)[-c(1:3)]

# Map ---------------------------------------------------------------------

for(i in 1:length(huc8_dat)){
  map <- tm_shape(huc8_dat[[i]]) +
    tm_polygons('val', palette = 'div', breaks = seq(-0.32,0.32, by = 0.04),
                legend.show = TRUE, legend.reverse = TRUE, title = 'rho') +
    tm_facets('var', nrow = 1, ncol = 2) +
    tm_layout(panel.show = FALSE) +
    tm_shape(huc8_sigs[[i]]) + tm_dots(col = 'val', palette = c(NA, 'blue'), size = 0.2,
                                       legend.show = FALSE) + tm_facets('var', nrow = 1, ncol = 2)
  map
  tmap_save(map, filename = file.path('Figures','final','correlation map',paste0(names(huc8_dat)[i],'.png')),
            height = 3, width = 12, units = 'in')
  print(paste0('Figure ',i,'/',length(huc8_dat),' Completed!'))
}


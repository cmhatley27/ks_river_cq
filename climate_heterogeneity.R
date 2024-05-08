
# Load libraries ----------------------------------------------------------
library(tidyverse)
library(cluster)  
library(factoextra)
library(reshape)
library(brunnermunzel)
library(ggpubr)
library(sf)
library(tmap)
source('load_and_combine.R')
source('Theme+Settings.R')


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


# BM differences --------------------------------------------
huc8s <- st_read('./DataFiles/shapefiles/created/huc8_boundaries.shp')
streams <- st_read('./DataFiles/shapefiles/created/rivers_ksrb.shp')
reservoirs <- st_read('./DataFiles/shapefiles/created/reservoirs.shp')

climate_vars <- colnames(select(event_summary, starts_with('d_'), ends_with('_ante'), ends_with('1d_total')))
climate_var_order <- c('15d_ante', '30d_ante', '60d_ante', '90d_ante', '180d_ante', '365d_ante',
                       'spei1', 'spei3', 'spei6', 'spei9', 'spei12', 'spei18', 'spei24', '1d_total')

climate_q <- c(0.25,0.75)

res_var <- 'reservoir_runoff_ratio'
res_name <- 'Res/Runoff'
res_q <- c(0.5,0.5)

response_vars <- c('FI_n', 'HI_n', 'log_loadrate', 'initial_N', 'initial_Q', 'delta_N')

sig_level <- 0.05

dat <- select(event_summary, all_of(c(response_vars, res_var, climate_vars))) %>%
  mutate(res_end = quantize(select(.,res_var), res_q)) %>%
  mutate(across(all_of(climate_vars), ~quantize(.x, climate_q))) %>%
  pivot_longer(all_of(climate_vars), names_to = 'climate_var', values_to = 'climate_end')

bm_calcs <- function(x, clim, sig, res = NULL, res_side = c('low', 'high')){
  if(is.null(res)){
    bm_stat <- brunnermunzel.test(x[clim == 'low'], x[clim == 'high'])$statistic
    bm_p <- brunnermunzel.test(x[clim == 'low'], x[clim == 'high'])$statistic
  }
  else{
    bm_stat <- brunnermunzel.test(x[clim == 'low' & res == res_side], x[clim == 'high' & res == res_side])$statistic
    bm_p <- brunnermunzel.test(x[clim == 'low' & res == res_side], x[clim == 'high' & res == res_side])$statistic
  }
  return(ifelse(bm_p <= sig, bm_stat/abs(bm_stat), NA))
}

bm <- dat %>%
  group_by(climate_var) %>%
  summarise(across(all_of(response_vars),
                   ~bm_calcs(.x, climate_end, sig_level),
                   .names = '{.col}_stat'),
            across(all_of(response_vars),
                   ~bm_calcs(.x, climate_end, sig_level, res_end, res_side = 'low'),
                   .names = '{.col}_lowres_stat'),
            across(all_of(response_vars),
                   ~bm_calcs(.x, climate_end, sig_level, res_end, res_side = 'high'),
                   .names = '{.col}_highres_stat')) %>%
  mutate(HUC8 = str_sub(climate_var, 3, 10),
         var = str_sub(climate_var, 12, -1)) %>%
  select(climate_var, HUC8, var, everything())

huc8_dat <- list()
for(element in 1:(length(colnames(bm))-3)){
  huc8_dat[[element]] <- pivot_wider(bm, id_cols = HUC8, names_from = var, values_from = colnames(bm)[element+3]) %>%
                                   left_join(huc8s, .) %>%
                                   pivot_longer(!colnames(huc8s), names_to = 'var', values_to = 'val') %>%
                                   mutate(var = factor(var, levels = climate_var_order))
  if(str_detect(colnames(bm)[element+3], '_p')){
    huc8_dat[[element]] <- mutate(huc8_dat[[element]], sig = ifelse(val <= sig_level, TRUE, FALSE))
  }
}
names(huc8_dat) <- colnames(bm)[-c(1:3)]


# Map ---------------------------------------------------------------------
#basic
tm_shape(huc8_dat$fi_climate_p) +
  tm_polygons('sig', palette = c('grey85', 'coral2')) +
  tm_facets('var')

#fancy
tm_shape(huc8_dat$fi_climate_lowres_p) +
  tm_polygons('sig', palette = c('grey85', 'coral2')) +
  tm_shape(streams) + tm_lines(col = 'grey75', lwd = 0.5) +
  tm_shape(reservoirs) + tm_fill(col = 'skyblue')


# correlations --------------------------------------------------------
huc8s <- st_read('./DataFiles/shapefiles/created/huc8_boundaries.shp')
streams <- st_read('./DataFiles/shapefiles/created/rivers_ksrb.shp')
reservoirs <- st_read('./DataFiles/shapefiles/created/reservoirs.shp')

climate_vars <- colnames(select(event_summary, 
                                #starts_with('d_') & ends_with(c('spei12')), 
                                ends_with('_ante') & contains(c('15d', '30d','90d','365d')),
                                starts_with('p_') & ends_with('1d_total')))
climate_var_order <- c('1d_total', '15d_ante', '30d_ante', '90d_ante', '365d_ante')

res_var <- 'reservoir_runoff_ratio'
res_name <- 'Res/Runoff'
res_q <- c(0.5,0.5)

response_vars <- c('FI_n', 'HI_n')

sig_level <- 0.05

dat <- select(event_summary, c(all_of(c(response_vars, res_var, climate_vars)))) %>%
  mutate(res_end = quantize(select(.,res_var), res_q)) %>%
  pivot_longer(all_of(climate_vars), names_to = 'climate_var', values_to = 'climate_val')

cor_calcs <- function(x, clim, sig, res = NULL, res_side = c('low','high')){
  if(is.null(res)){
  cor_stat <- cor.test(x,clim, use = 'pairwise.complete.obs')$estimate
  cor_p <- cor.test(x,clim, use = 'pairwise.complete.obs')$p.value
  }
  else{
    cor_stat <- cor.test(x[res == res_side],clim[res == res_side], use = 'pairwise.complete.obs')$estimate
    cor_p <- cor.test(x[res == res_side],clim[res == res_side], use = 'pairwise.complete.obs')$p.value
  }
  return(ifelse(cor_p <= sig, cor_stat, NA))
}

cors <- dat %>%
  group_by(climate_var) %>%
  summarise(across(all_of(response_vars),
                   ~cor_calcs(.x, climate_val, sig_level),
                   .names = '{.col}_stat'),
            across(all_of(response_vars),
                   ~cor_calcs(.x, climate_val, sig_level, res_end, res_side = 'low'),
                   .names = '{.col}_lowres_stat'),
            across(all_of(response_vars),
                   ~cor_calcs(.x, climate_val, sig_level, res_end, res_side = 'high'),
                   .names = '{.col}_highres_stat')) %>%
  mutate(HUC8 = str_sub(climate_var, 3, 10),
         var = str_sub(climate_var, 12, -1)) %>%
  select(climate_var, HUC8, var, everything())
    
huc8_dat <- list()
for(element in 1:(length(colnames(cors))-3)){
  huc8_dat[[element]] <- pivot_wider(cors, id_cols = HUC8, names_from = var, values_from = colnames(cors)[element+3]) %>%
    left_join(huc8s, .) %>%
    pivot_longer(!colnames(huc8s), names_to = 'var', values_to = 'val') %>%
    mutate(var = factor(var, levels = climate_var_order))
  if(str_detect(colnames(cors)[element+3], '_p')){
    huc8_dat[[element]] <- mutate(huc8_dat[[element]], sig = ifelse(val <= sig_level, TRUE, FALSE))
  }
}
names(huc8_dat) <- colnames(cors)[-c(1:3)]

# Map ---------------------------------------------------------------------
for(i in 1:length(huc8_dat)){
  map <- tm_shape(huc8_dat[[i]]) +
    tm_polygons('val', palette = 'div', breaks = seq(-0.4,0.4, by = 0.05), legend.show = TRUE) +
    tm_facets('var', nrow = 1, ncol = 5)
  tmap_save(map, filename = file.path('Figures','maps','correlations', 'reduced2',paste0(names(huc8_dat)[i],'.tiff')),
             height = 1.5, width = 12, units = 'in', compression = 'lzw')
  print(paste0('Figure ',i,'/',length(huc8_dat),' Completed!'))
}

test_map <- tm_shape(huc8_dat$FI_n_stat) +
  tm_polygons('val', palette = 'div', breaks = seq(-0.4,0.4, by = 0.05), legend.show = FALSE) +
  tm_facets('var')
tmap_save(test_map, filename = file.path('Figures','maps','cluster_summaries',paste0('test.tiff')),
          height = 6.67, width = 8, units = 'in', compression = 'lzw')
test_map <- tm_shape(huc8_dat$FI_n_lowres_stat) +
  tm_polygons('val', palette = 'div', breaks = seq(-0.4,0.4, by = 0.05), legend.show = FALSE) +
  tm_facets('var')
tmap_save(test_map, filename = file.path('Figures','maps','cluster_summaries',paste0('test_lowres.tiff')),
          height = 6.67, width = 8, units = 'in', compression = 'lzw')
test_map <- tm_shape(huc8_dat$FI_n_highres_stat) +
  tm_polygons('val', palette = 'div', breaks = seq(-0.4,0.4, by = 0.05), legend.show = FALSE) +
  tm_facets('var')
tmap_save(test_map, filename = file.path('Figures','maps','cluster_summaries',paste0('test_highres.tiff')),
          height = 6.67, width = 8, units = 'in', compression = 'lzw')



# regressions --------------------------------------------------------------
huc8s <- st_read('./DataFiles/shapefiles/created/huc8_boundaries.shp')
streams <- st_read('./DataFiles/shapefiles/created/rivers_ksrb.shp')
reservoirs <- st_read('./DataFiles/shapefiles/created/reservoirs.shp')

climate_vars <- colnames(select(event_summary, 
                                starts_with('d_') & ends_with(c('1','3','6','12')), 
                                ends_with('_ante') & contains(c('30d','90d','180d','365d')),
                                ends_with('1d_total')))
climate_var_order <- c('30d_ante', '90d_ante', '180d_ante', '365d_ante',
                       'pet1', 'pet3', 'pet6', 'pet12',
                       'bal1', 'bal3', 'bal6', 'bal12',
                       'spei1', 'spei3', 'spei6', 'spei12',
                       '1d_total')

res_var <- 'reservoir_runoff_ratio'
res_name <- 'Res/Runoff'
res_q <- c(0.5,0.5)

response_vars <- c('FI_n', 'HI_n')

sig_level <- 0.05

dat <- select(event_summary, all_of(c(response_vars, res_var, climate_vars))) %>%
  mutate(res_end = quantize(select(.,res_var), res_q)) %>%
  pivot_longer(all_of(climate_vars), names_to = 'climate_var', values_to = 'climate_val')

regresso <- function(x, clim, res = NULL, sig){
  frame <- tibble(x,clim,res) %>%
    filter(complete.cases(.)) %>%
    mutate(across(where(is.numeric),  ~(.x - mean(.x))/sqrt(var(.x))))
  if(is.null(res)){
    mod <- lm(x ~ clim, data = frame)
  } else{
  mod <- lm(x ~ clim + res, data = frame)
  }
  return(ifelse(summary(mod)$coefficients[2,4] <= sig, mod$coefficients[2], NA))
}

regs <- dat %>%
  group_by(climate_var) %>%
  summarise(across(all_of(response_vars),
                   ~regresso(.x, climate_val, res = reservoir_runoff_ratio, sig = sig_level),
                   .names = '{.col}_coef'),
            across(all_of(response_vars),
                   ~regresso(.x, climate_val, sig = sig_level),
                   .names = '{.col}_coef_nores')) %>%
  mutate(HUC8 = str_sub(climate_var, 3, 10),
         var = str_sub(climate_var, 12, -1)) %>%
  select(climate_var, HUC8, var, everything())

huc8_dat <- list()
for(element in 1:(length(colnames(regs))-3)){
  huc8_dat[[element]] <- pivot_wider(regs, id_cols = HUC8, names_from = var, values_from = colnames(regs)[element+3]) %>%
    left_join(huc8s, .) %>%
    pivot_longer(!colnames(huc8s), names_to = 'var', values_to = 'val') %>%
    mutate(var = factor(var, levels = climate_var_order))
  if(str_detect(colnames(regs)[element+3], '_p')){
    huc8_dat[[element]] <- mutate(huc8_dat[[element]], sig = ifelse(val <= sig_level, TRUE, FALSE))
  }
}
names(huc8_dat) <- colnames(regs)[-c(1:3)]

# Map ---------------------------------------------------------------------

test_map <- tm_shape(huc8_dat$FI_n_coef) +
  tm_polygons('val', palette = 'div', breaks = seq(-0.3,0.3, by = 0.05)) +
  tm_facets('var')
tmap_save(test_map, filename = file.path('Figures','maps','cluster_summaries',paste0('test.tiff')),
          height = 6.67, width = 8, units = 'in', compression = 'lzw')

tm_shape(huc8_dat$FI_n_coef_nores) +
  tm_polygons('val', palette = 'div', breaks = seq(-0.3,0.3, by = 0.05)) +
  tm_facets('var')


# clusters ----------------------------------------------------------------
huc8s <- st_read('./DataFiles/shapefiles/created/huc8_boundaries.shp')
streams <- st_read('./DataFiles/shapefiles/created/rivers_ksrb.shp')
reservoirs <- st_read('./DataFiles/shapefiles/created/reservoirs.shp')

climate_vars <- colnames(select(event_summary, 
                                starts_with('d_') & ends_with(c('1','3','6','12')), 
                                ends_with('_ante') & contains(c('7d', '15d', '30d','90d','180d','270d', '365d')),
                                ends_with('1d_total')))
climate_var_order <- c('7d_ante', '15d_ante', '30d_ante', '90d_ante', '180d_ante', '270d_ante', '365d_ante',
                       'pet1', 'pet3', 'pet6', 'pet12',
                       'bal1', 'bal3', 'bal6', 'bal12',
                       'spei1', 'spei3', 'spei6', 'spei12',
                       '1d_total')

res_var <- 'reservoir_runoff_ratio'
res_name <- 'Res/Runoff'
res_q <- c(0.5,0.5)

response_vars <- c('cluster')

sig_level <- 0.05

dat <- select(event_summary, all_of(c(response_vars, res_var, climate_vars))) %>%
  filter(!is.na(cluster)) %>%
  mutate(res_end = quantize(select(.,res_var), res_q),
         cluster = paste0('clus',cluster)) 

clus_summary <- dat %>%
  group_by(cluster) %>%
  summarise(across(all_of(climate_vars),
                   ~median(.x))) %>%
  pivot_longer(!cluster, names_to = 'climate_var', values_to = 'climate_val') %>%
  pivot_wider(names_from = cluster, values_from = 'climate_val') %>%
  mutate(HUC8 = str_sub(climate_var, 3, 10),
         var = str_sub(climate_var, 12, -1)) %>%
  select(climate_var, HUC8, var, everything())

huc8_dat <- tibble()
for(element in 1:(length(colnames(clus_summary))-3)){
  huc8_dat_element <- pivot_wider(clus_summary, id_cols = HUC8, names_from = var, values_from = colnames(clus_summary)[element+3]) %>%
    left_join(huc8s, .) %>%
    pivot_longer(!colnames(huc8s), names_to = 'var', values_to = 'val') %>%
    mutate(var = factor(var, levels = climate_var_order),
           cluster = paste0('clus',element))
  huc8_dat <- rbind(huc8_dat, huc8_dat_element)
}


# Map ---------------------------------------------------------------------
#Summary maps of each cluster
for(clus in 1:(length(colnames(clus_summary))-3)){
  huc8_dat_fil <- filter(huc8_dat, cluster == paste0('clus',clus))
  d <- tm_shape(subset(huc8_dat_fil, str_detect(var, 'spei'))) +
    tm_polygons('val', palette = 'div', breaks = seq(-1,1, by = 0.05), legend.show = FALSE) +
    tm_facets(by = 'var', ncol = 1)
  p <- tm_shape(subset(huc8_dat_fil, str_detect(var, '1d_total'))) +
    tm_polygons('val', breaks = seq(0,20, by = 1), legend.show = FALSE, palette = 'Blues')
  comb <- tmap_arrange(d,p, nrow = 1)
  tmap_save(comb, filename = file.path('Figures','maps','cluster_summaries',paste0('clus',clus,'.tiff')),
            height = 6, width = 9, units = 'in', compression = 'lzw')
}

summary(huc8_dat$val[huc8_dat$var == 'spei1'])

#Maps of all vars together
vars <- c('7d_ante', '15d_ante', '30d_ante', '90d_ante', '180d_ante', '270d_ante', '365d_ante')
vars <- c('pet1$', 'pet3', 'pet6', 'pet12')

pal <- 'Oranges'
var1 <- tm_shape(subset(huc8_dat, str_detect(var, vars[1]))) +
  tm_polygons('val', n=20, legend.show = F, palette = pal) +
  tm_facets('cluster', free.scales = F, nrow = 1)
var2 <- tm_shape(subset(huc8_dat, str_detect(var, vars[2]))) +
  tm_polygons('val', n=20, legend.show = F, palette = pal) +
  tm_facets('cluster', free.scales = F, nrow = 1)
var3 <- tm_shape(subset(huc8_dat, str_detect(var, vars[3]))) +
  tm_polygons('val', n=20, legend.show = F, palette = pal) +
  tm_facets('cluster', free.scales = F, nrow = 1)
var4 <- tm_shape(subset(huc8_dat, str_detect(var, vars[4]))) +
  tm_polygons('val', n=20, legend.show = F, palette = pal) +
  tm_facets('cluster', free.scales = F, nrow = 1)
var5 <- tm_shape(subset(huc8_dat, str_detect(var, vars[5]))) +
  tm_polygons('val', n=20, legend.show = F, palette = pal) +
  tm_facets('cluster', free.scales = F, nrow = 1)
var6 <- tm_shape(subset(huc8_dat, str_detect(var, vars[6]))) +
  tm_polygons('val', n=20, legend.show = F, palette = pal) +
  tm_facets('cluster', free.scales = F, nrow = 1)
var7 <- tm_shape(subset(huc8_dat, str_detect(var, vars[7]))) +
  tm_polygons('val', n=20, legend.show = F, palette = pal) +
  tm_facets('cluster', free.scales = F, nrow = 1)
all_vars <- tmap_arrange(var1,var2,var3,var4, nrow = 4)
tmap_save(all_vars, filename = file.path('Figures','maps','cluster_summaries',paste0('pet_all.tiff')),
          height = 6, width = 27.4, units = 'in', compression = 'lzw')

d <- tm_shape(subset(huc8_dat, str_detect(var, 'spei'))) +
  tm_polygons('val', palette = 'div', breaks = seq(-1,1, by = 0.05), legend.show = FALSE) +
  tm_facets(by = c('var','cluster'))
tmap_save(d, filename = file.path('Figures','maps','cluster_summaries','k8',paste0('d_all.tiff')),
          height = 6, width = 27.4, units = 'in', compression = 'lzw')

p <- tm_shape(subset(huc8_dat, str_detect(var, '1d_total'))) +
  tm_polygons('val', breaks = seq(0,40, by = 2), legend.show = FALSE, palette = 'Blues') +
  tm_facets('cluster', nrow = 2)
tmap_save(p, filename = file.path('Figures','maps','cluster_summaries','k8', paste0('p_all.tiff')),
          height = 6, width = 24, units = 'in', compression = 'lzw')


# season and res summaries -----------------------------------------------
seasons <- select(event_summary, season, cluster) %>%
  filter(!is.na(cluster))

ggplot(data = seasons, aes(x = season, fill = season)) +
  geom_bar(color = 'black') +
  scale_fill_manual(values = c('#8DA9C4','#60992D','#C33C54','#F4B860')) +
  facet_wrap(vars(cluster), nrow = 2)

months <- select(event_summary, month, season, cluster) %>%
  filter(!is.na(cluster)) %>%
  mutate(month = factor(month, levels = c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun' ,'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec')))

ggplot(data = months, aes(x = month, fill = season)) +
  geom_bar(color = 'black') +
  scale_fill_manual(values = c('#8DA9C4','#60992D','#C33C54','#F4B860')) +
  facet_wrap(vars(cluster), nrow = 2)
  
res <- select(event_summary, cluster, reservoir_runoff_ratio, r_t_share, r_m_share, r_p_share) %>%
  filter(!is.na(cluster)) %>%
  pivot_longer(!cluster, names_to = 'var', values_to = 'val') %>%
  mutate(var = factor(var, levels = c('reservoir_runoff_ratio', 'r_t_share', 'r_m_share', 'r_p_share')))
ggplot(data = res, aes(x=factor(cluster), y = val, fill = factor(cluster))) +
  geom_hline(yintercept = c(0.25,0.5,0.75), linetype = 'dashed') +
  geom_boxplot() +
  ylim(c(0,1)) +
  facet_wrap(vars(var), nrow = 2)

res_totals <- select(event_summary, cluster, starts_with('r_') & ends_with('1d_total')) %>%
  filter(!is.na(cluster)) %>%
  pivot_longer(!cluster, names_to = 'var', values_to = 'val')
ggplot(data = subset(res_totals, val <= 10^9), aes(x=factor(cluster), y = val, fill = factor(cluster))) +
  geom_boxplot() +
  #ylim(c(0,10^9)) +
  facet_wrap(vars(var), nrow = 2, scales = 'free_y')

res_deltas <- select(event_summary, cluster, starts_with('r_') & ends_with('1d_delta_rat')) %>%
  filter(!is.na(cluster)) %>%
  pivot_longer(!cluster, names_to = 'var', values_to = 'val')
ggplot(data = subset(res_deltas, val <= 10), aes(x=factor(cluster), y = val, fill = factor(cluster))) +
  geom_boxplot() +
  #ylim(c(0,10^9)) +
  facet_wrap(vars(var), nrow = 2, scales = 'free_y')

ggplot(data = res, aes(x=factor(cluster), y = (val), fill = factor(cluster))) +
  geom_hline(yintercept = c(0.25,0.5,0.75),linetype = 'dashed') +
  ylim(c(0,1)) +
  geom_boxplot() +
  facet_wrap(vars(var))

sqrt(var(event_summary$r_c_share))
size <- select(event_summary, cluster, initial_Q, max_Q) %>%
  filter(!is.na(cluster)) %>%
  pivot_longer(!cluster, names_to = 'var', values_to = 'val') %>%
  mutate(var = factor(var, levels = c('initial_Q', 'max_Q')))

ggplot(data = size, aes(x=var, y = (val), fill = var)) +
  geom_hline(yintercept = seq(250,1750,by=250),linetype = 'dashed') +
  ylim(c(0,2000)) +
  geom_boxplot() +
  facet_wrap(vars(cluster), nrow = 2)

ggplot(data = size, aes(x=factor(cluster), y = (val), fill = factor(cluster))) +
  geom_hline(yintercept = seq(250,1750,by=250),linetype = 'dashed') +
  ylim(c(0,2000)) +
  geom_boxplot() +
  facet_wrap(vars(var))

ggplot(data = event_summary, aes(x=factor(cluster), y = initial_N, fill = factor(cluster))) +
  geom_boxplot()


ggplot(data = event_summary, aes(x = factor(cluster), y = total_Q, fill = factor(cluster))) +
  geom_boxplot()


ggplot(data = subset(event_summary, !is.na(cluster)), aes(x = FI_n, y = HI_n, color = factor(cluster))) +
  geom_point() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  xlim(c(-1,1)) +
  ylim(c(-1,1))

ggplot(data = subset(event_summary, !is.na(cluster)), aes(x = FI_t, y = HI_t, color = factor(cluster))) +
  geom_point() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  xlim(c(-1,1)) +
  ylim(c(-1,1)) +
  facet_wrap(vars(cluster))


ggsave(file.path('Figures','maps','cluster_summaries','scatter.tiff'),
       height = 6, width =8, units = 'in', compression = 'lzw')


# misc --------------------------------------------------------------------
summary <- select(event_summary, starts_with('p_') & ends_with('365d_ante')) %>%
  summarise(across(everything(), ~summary(.x)))


asdf <- cor.test(event_summary$HI_t, event_summary$HI_n, use = 'pairwise.complete.obs')
?cor.test



ggplot(data = event_summary, aes(x = yday, y = N_load/total_Q)) +
  geom_point() +
  geom_smooth()
ggplot(data = event_summary, aes(x = yday, y = FI_n)) +
  geom_point() +
  geom_smooth()

ggplot(data = event_summary, aes(x = total_Q, y = N_load)) +
  geom_point() +
  geom_smooth()

ggplot(data = event_summary, aes(x = season, y = FI_n)) +
  geom_boxplot()

ggplot(data = event_summary) +
  geom_boxplot(aes(x = factor(cluster), y = FI_c, fill = factor(cluster)))

ggplot(data = event_summary) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  ylim(c(-1,1)) +
  xlim(c(-1,1)) +
  geom_point(aes(x = FI_n, y = HI_n)) 

asdf <- event_summary %>%
  filter(!is.na(cluster)) %>%
  select(cluster, FI_n) %>%
  mutate(FI_n_before = lag(FI_n),
         cluster_before = lag(cluster)) %>%
  filter(cluster == 7)

table(asdf$cluster_before)

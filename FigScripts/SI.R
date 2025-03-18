library(tidyverse)
source('Theme+Settings.R')

# Figure S1: precip distributions and snow percent -------------------------------------------------------
library(terra)

#load daily precip and mean temp
precip <- rast('./DataFiles/climate_data/daily_rasters/precip_ksrb.tiff')
precip_mean <- global(precip, 'mean', na.rm = T) %>%
  mutate(date = rownames(.)) %>%
  select(date, precip = mean)

temp <- rast('./DataFiles/climate_data/daily_rasters/temp_ksrb.tiff')
temp_mean <- global(temp, 'mean', na.rm = T) %>%
  mutate(date = rownames(.)) %>%
  select(date, temp = mean)

#combine and roughly estimate snow as precip on days when mean temp < 0
dat <- left_join(precip_mean, temp_mean, by = 'date') %>%
  mutate(date = ymd(date),
         wateryear = ifelse(month(date) > 9, year(date) + 1, year(date)),
         freeze = temp < 0,
         snow = precip*freeze)

#aggregate to month, calculate mean monthly total precip
monthly_precip <- dat %>%
  group_by(wateryear, month = month(date)) %>%
  summarise(total_precip = sum(precip)) %>%
  group_by(month) %>%
  summarise(avg_precip = mean(total_precip))

#plot monthly bars
ggplot(monthly_precip, aes(x = factor(month), y = avg_precip)) +
  geom_col(color = 'black') +
  ylab('Mean Total Precip [mm]') +
  xlab('Month')
ggsave('./Figures/final/SI/1_precip/month_bars.png', height = 3, width = 5, units = 'in')

#aggregate snow and precip to water year
yearly_snow <- dat %>%
  group_by(wateryear) %>%
  summarise(precip = sum(precip, na.rm = T),
            snow = sum(snow, na.rm = T),
            snow_frac = snow/precip) %>%
  pivot_longer(!wateryear, names_to = 'var', values_to = 'val') %>%
  filter(var != 'snow_frac')

#plot yearly snow bars
ggplot(yearly_snow, aes(x = factor(wateryear), y = val, fill = var)) +
  geom_col(color = 'black', position = position_jitter(width = 0, height = 0)) +
  ylab('Mean Total Precip [mm]') +
  xlab('Water Year') +
  scale_fill_manual(values = c('grey80', 'black'), labels = c('Rain','Snow'), name = 'Precip')
ggsave('./Figures/final/SI/1_precip/year_snow_bars.png', height = 3, width = 5, units = 'in')

# Table S1: Reservoir flow stats ---------------------------------------------
dat <- read_csv('./DataFiles/hydro_data/all_hydro_daily.csv') %>%
  group_by(wateryear) %>%
  summarise(q = sum(riverQ_linterp),
            res = sum(reservoir_sum_linterp),
            t = sum(tuttleQ_linterp),
            m = sum(milfordQ_linterp),
            p = sum(perryQ_linterp),
            c = sum(clintonQ_linterp),
            k = sum(kanopolisQ_linterp),
            wi = sum(wilsonQ_linterp),
            wa = sum(wacondaQ_linterp))  %>%
  mutate(res_q = res/q,
         across(!c(wateryear, q, res, res_q), ~./res, .names = '{.col}_pct')) %>%
  filter(wateryear != 2013)

means <- dat %>%
  summarise(across(everything(), mean))


# Table S3: HUC 8 precip stats ----------------------------------------------------------------
source('8_combine_all_event_data.R')

dat <- select(event_summary,
              c(initial_Q, max_Q, delta_Q, delta_rat_Q, duration, duration_since_last,
                starts_with('p_'), ends_with('spei12'), 
                reservoir_runoff_ratio, ends_with('share')))

dat_summary <- dat %>%
  summarise(across(everything(), list(min = min, median = median, mean = mean, max = max))) %>%
  pivot_longer(everything(), values_to = 'val', names_to = 'var') %>%
  mutate(stat = word(var, -1, -1, sep = '_'),
         var = word(var, 1, -2, sep = '_'),
         huc_res = word(var, 2, 2, sep = '_')) %>%
  pivot_wider(id_cols = c(var, huc_res), names_from = stat, values_from = val) %>%
  mutate(across(c(min,median,mean,max), ~round(.x, digits = 2)))
write_csv(dat_summary, './DataFiles/hydro_data/event_char/summary_stats.csv')

# Figure S2: Climate conditions in each year ------------------------------------
library(sf)
library(ggspatial)

hucs_sel <- c('10270104','10260008','10270201','10250002','10260003','10260012')
spei12_huc8 <- read_csv('./DataFiles/hydro_data/drought/spei12_huc8.csv') %>%
  pivot_longer(!date, names_to = 'huc', values_to = 'spei') %>%
  filter(huc %in% hucs_sel)

ggplot(spei12_huc8, aes(x = date, y = spei, color = huc)) +
  geom_line() +
  scale_x_date(limits = ymd(c('2013-10-01','2022-09-30')), name = NULL) +
  scale_color_discrete(name = 'HUC 8 Subbasin') +
  ylab('1-Year SPEI [-]') +
  theme(legend.position = 'bottom')
ggsave('./Figures/final/SI/2_spei_timeseries/spei12.png', width = 5, height = 3, units = 'in')

huc8s <- st_read('./DataFiles/shapefiles/created/huc8_boundaries.shp') %>%
  mutate(select = ifelse(HUC8 %in% c('10270104','10260008','10270201','10250002','10260003','10260012'),T,F))
ggplot() +
  geom_sf(data = huc8s, color = 'grey45') +
  geom_sf(data = subset(huc8s, select == T), aes(fill = HUC8)) +
  geom_sf_label(data = subset(huc8s, select == T), aes(label = HUC8)) +
  scale_fill_discrete(guide = NULL) +
  theme_void()
ggsave('./Figures/final/SI/2_spei_timeseries/basemap.png', width = 5, height = 2, units = 'in')

# Figure S3: Non-spatial event characteristics by season ---------------------
source('8_combine_all_event_data.R')
source('Theme+Settings.R')

szn_dat_static <- select(event_summary, c(season, initial_Q, delta_Q, max_Q, delta_rat_Q, duration, duration_since_last,
                  reservoir_runoff_ratio, all_of(ends_with('_share')))) %>%
  pivot_longer(!season, names_to = 'var', values_to = 'val') %>%
  mutate(var = factor(var, 
                      levels = c('initial_Q', 'delta_Q', 'max_Q', 'delta_rat_Q',
                                                     'duration', 'duration_since_last',
                                                     'reservoir_runoff_ratio', 'r_t_share', 'r_m_share', 'r_p_share', 'r_c_share',
                                                     'r_k_share','r_wi_share','r_wa_share'),
                      labels = c('Initial Q [m^3/s]', 'Delta Q [m^3/s]', 'Max Q [m^3/s]', 'Delta/Initial Q [-]', 
                                 'Duration [days]', 'Duration Since Last [days]',
                                 'Summed Outflow Ratio [-]', 'Tuttle Ratio [-]', 'Milford Ratio [-]', 'Perry Ratio [-]', 'Clinton Ratio [-]',
                                 'Kanopolis Ratio [-]','Wilson Ratio [-]','Waconda Ratio [-]')))

ggplot(szn_dat_static, aes(x = season, y = val, fill = season)) +
  geom_boxplot() +
  facet_wrap(vars(var), scales = 'free_y', ncol = 3) +
  scale_fill_manual(values = season_palette, name = 'Season') +
  ylab('') + xlab('') +
  theme(text = element_text(size = 9), legend.position = 'bottom')
ggsave('./Figures/final/SI/3_seasonal_static/3_seasonal_static.png', width = 6.5, height = 9, units = 'in')

vars <- unique(szn_dat_static$var)
for(test_var in 1:length(vars)){
  dat_sub <- subset(szn_dat_static, var == vars[test_var])
  print(vars[test_var])
  conover.test(dat_sub$val, dat_sub$season)
}

# Figure S4: Spatial event characteristics by season ------------------------
source('8_combine_all_event_data.R')
source('Theme+Settings.R')
library(sf)
library(conover.test)

huc_sel <- c('10270104','10260008','10270201','10250002','10260003','10260012', 'Basin_Avg')

szn_dat_precip <- select(event_summary, c(season, starts_with('p') & ends_with('1d_total')))
colnames(szn_dat_precip) <- c('season', word(colnames(szn_dat_precip)[-1], 2, sep = '_'))
szn_dat_precip$Basin_Avg <- rowMeans(select(szn_dat_precip, !season))
for(huc in 1:length(huc_sel)){
  plot_dat <- select(szn_dat_precip, c(season, huc_sel[huc]))
  plot <- ggplot(plot_dat, aes(x = season, fill = season, y = .data[[huc_sel[huc]]])) +
    geom_boxplot() +
    scale_fill_manual(values = season_palette, name = 'Season', guide = 'none') +
    scale_x_discrete(labels = NULL, name = NULL) +
    scale_y_continuous(name = NULL) +
    theme(text = element_text(size = 9)) +
    ggtitle('Event Precip [mm]')
  ggsave(paste0('./Figures/final/SI/4_seasonal_spatial/',huc_sel[huc],'_event_precip.png'), plot = plot,
         width = 1.67, height = 1.67, units = 'in')
}
for(test_var in 1:length(huc_sel)){
  print(huc_sel[test_var])
  conover.test(szn_dat_precip[[huc_sel[test_var]]], szn_dat_precip$season)
}


szn_dat_ap30 <- select(event_summary, c(season, starts_with('p') & ends_with('30d_ante')))
colnames(szn_dat_ap30) <- c('season', word(colnames(szn_dat_ap30)[-1], 2, sep = '_'))
szn_dat_ap30$Basin_Avg <- rowMeans(select(szn_dat_ap30, !season))
for(huc in 1:length(huc_sel)){
  plot_dat <- select(szn_dat_ap30, c(season, huc_sel[huc]))
  plot <- ggplot(plot_dat, aes(x = season, fill = season, y = .data[[huc_sel[huc]]])) +
    geom_boxplot() +
    scale_fill_manual(values = season_palette, name = 'Season', guide = 'none') +
    scale_x_discrete(labels = NULL, name = NULL) +
    scale_y_continuous(name = NULL) +
    theme(text = element_text(size = 9)) +
    ggtitle('30-day AP [mm]')
  ggsave(paste0('./Figures/final/SI/4_seasonal_spatial/',huc_sel[huc],'_ante_precip.png'), plot = plot,
         width = 1.67, height = 1.67, units = 'in')
}
vars <- c('10270104','10260008','10270201','10250002','10260003','10260012', 'Basin_Avg')
for(test_var in 1:length(huc_sel)){
  print(huc_sel[test_var])
  conover.test(szn_dat_ap30[[huc_sel[test_var]]], szn_dat_ap30$season)
}


szn_dat_spei <- select(event_summary, c(season, starts_with('d') & ends_with('spei12')))
colnames(szn_dat_spei) <- c('season', word(colnames(szn_dat_spei)[-1], 2, sep = '_'))
szn_dat_spei$Basin_Avg <- rowMeans(select(szn_dat_spei, !season))
for(huc in 1:length(huc_sel)){
  plot_dat <- select(szn_dat_spei, c(season, huc_sel[huc]))
  plot <- ggplot(plot_dat, aes(x = season, fill = season, y = .data[[huc_sel[huc]]])) +
    geom_boxplot() +
    scale_fill_manual(values = season_palette, name = 'Season', guide = 'none') +
    scale_x_discrete(labels = NULL, name = NULL) +
    scale_y_continuous(name = NULL) +
    theme(text = element_text(size = 9)) +
    ggtitle('1-Year SPEI [-]')
  ggsave(paste0('./Figures/final/SI/4_seasonal_spatial/',huc_sel[huc],'_spei12.png'), plot = plot,
         width = 1.67, height = 1.67, units = 'in')
}
for(test_var in 1:length(huc_sel)){
  print(huc_sel[test_var])
  conover.test(szn_dat_spei[[huc_sel[test_var]]], szn_dat_spei$season)
}


huc8s <- st_read('./DataFiles/shapefiles/created/huc8_boundaries.shp') %>%
  mutate(select = ifelse(HUC8 %in% c('10270104','10260008','10270201','10250002','10260003','10260012'),T,F))
ggplot() +
  geom_sf(data = huc8s, aes(fill = select), color = 'grey45') +
  scale_fill_manual(guide = 'none', values = c('grey90','steelblue')) +
  geom_sf_label(data = subset(huc8s, select == T), aes(label = HUC8)) +
  theme_void()
ggsave('./Figures/final/SI/4_seasonal_spatial/basemap.png', width = 9, height = 4.5, units = 'in')



# Figure S5: Non-spatial event characteristics by precip/res -----------------
source('8_combine_all_event_data.R')
source('Theme+Settings.R')

event_summary <- event_summary %>%
  mutate(res_end = cut(reservoir_runoff_ratio, breaks = c(0,0.5,1.5), labels = c('Precip','Reservoir')))

res_dat_static <- select(event_summary, c(res_end, initial_Q, delta_Q, max_Q, delta_rat_Q, duration, duration_since_last,
                                          all_of(ends_with('_share')))) %>%
  pivot_longer(!res_end, names_to = 'var', values_to = 'val') %>%
  mutate(var = factor(var, 
                      levels = c('initial_Q', 'delta_Q', 'max_Q', 'delta_rat_Q',  
                                 'duration', 'duration_since_last',
                                 'r_t_share', 'r_m_share', 'r_p_share', 'r_c_share',
                                 'r_k_share','r_wi_share','r_wa_share'),
                      labels = c('Initial Q [m^3/s]', 'Delta Q [m^3/s]', 'Max Q [m^3/s]', 'Delta/Initial Q [-]',  
                                 'Duration [days]', 'Duration Since Last [days]',
                                 'Tuttle Ratio [-]', 'Milford Ratio [-]', 'Perry Ratio [-]', 'Clinton Ratio [-]',
                                 'Kanopolis Ratio [-]','Wilson Ratio [-]','Waconda Ratio [-]')))

ggplot(res_dat_static, aes(x = res_end, y = val, fill = res_end)) +
  geom_boxplot() +
  facet_wrap(vars(var), scales = 'free_y', ncol = 3) +
  scale_fill_manual(values = res_palette, name = 'Grouping') +
  ylab('') + xlab('') +
  theme(text = element_text(size = 9), legend.position = 'bottom')
ggsave('./Figures/final/SI/5_res_static/res_static.png', width = 6.5, height = 9, units = 'in')

vars <- unique(res_dat_static$var)
for(test_var in 1:length(vars)){
  dat_sub <- subset(res_dat_static, var == vars[test_var])
  print(vars[test_var])
  print(wilcox.test(dat_sub$val[dat_sub$res_end == 'Precip'], dat_sub$val[dat_sub$res_end == 'Reservoir'])$p.value)
}

res_seasons <- as.data.frame(table(select(event_summary, season, res_end)))
ggplot(res_seasons, aes(x = season, y = Freq, fill = res_end)) +
  geom_col(position = position_dodge2(), color = 'black') +
  scale_fill_manual(values = res_palette, guide = 'none') +
  ylab('n') +
  scale_x_discrete(name = NULL) +
  ggtitle('Occurrence by season')
ggsave('./Figures/final/SI/5_res_static/res_seasons.png', width = 4, height = 9/5, units = 'in')

# Figure S6: Spatial event characteristics by precip/res ---------------------
source('8_combine_all_event_data.R')
source('Theme+Settings.R')
library(sf)
library(conover.test)

event_summary <- event_summary %>%
  mutate(res_end = cut(reservoir_runoff_ratio, breaks = c(0,0.5,1.5), labels = c('Precip','Reservoir')))

huc_sel <- c('10270104','10260008','10270201','10250002','10260003','10260012', 'Basin_Avg')

res_dat_precip <- select(event_summary, c(res_end, starts_with('p') & ends_with('1d_total')))
colnames(res_dat_precip) <- c('res_end', word(colnames(res_dat_precip)[-1], 2, sep = '_'))
res_dat_precip$Basin_Avg <- rowMeans(select(res_dat_precip, !res_end))
for(huc in 1:length(huc_sel)){
  plot_dat <- select(res_dat_precip, c(res_end, huc_sel[huc]))
  plot <- ggplot(plot_dat, aes(x = res_end, fill = res_end, y = .data[[huc_sel[huc]]])) +
    geom_boxplot() +
    scale_fill_manual(values = res_palette, name = 'Season', guide = 'none') +
    scale_x_discrete(labels = NULL, name = NULL) +
    scale_y_continuous(name = NULL) +
    theme(text = element_text(size = 9)) +
    ggtitle('Event Precip [mm]')
  ggsave(paste0('./Figures/final/SI/6_res_spatial/',huc_sel[huc],'_event_precip.png'), plot = plot,
         width = 1.67, height = 1.67, units = 'in')
}

for(test_var in 1:length(huc_sel)){
  print(huc_sel[test_var])
  print(wilcox.test(res_dat_precip[[huc_sel[test_var]]][res_dat_precip$res_end == 'Precip'],
              res_dat_precip[[huc_sel[test_var]]][res_dat_precip$res_end == 'Reservoir'])$p.value)
}


res_dat_ap30 <- select(event_summary, c(res_end, starts_with('p') & ends_with('30d_ante')))
colnames(res_dat_ap30) <- c('res_end', word(colnames(res_dat_ap30)[-1], 2, sep = '_'))
res_dat_ap30$Basin_Avg <- rowMeans(select(res_dat_ap30, !res_end))
for(huc in 1:length(huc_sel)){
  plot_dat <- select(res_dat_ap30, c(res_end, huc_sel[huc]))
  plot <- ggplot(plot_dat, aes(x = res_end, fill = res_end, y = .data[[huc_sel[huc]]])) +
    geom_boxplot() +
    scale_fill_manual(values = res_palette, name = 'Season', guide = 'none') +
    scale_x_discrete(labels = NULL, name = NULL) +
    scale_y_continuous(name = NULL) +
    theme(text = element_text(size = 9)) +
    ggtitle('30-Day AP [mm]')
  ggsave(paste0('./Figures/final/SI/6_res_spatial/',huc_sel[huc],'_ante_precip.png'), plot = plot,
         width = 1.67, height = 1.67, units = 'in')
}
for(test_var in 1:length(huc_sel)){
  print(huc_sel[test_var])
  print(wilcox.test(res_dat_ap30[[huc_sel[test_var]]][res_dat_ap30$res_end == 'Precip'],
                    res_dat_ap30[[huc_sel[test_var]]][res_dat_ap30$res_end == 'Reservoir'])$p.value)
}


res_dat_spei <- select(event_summary, c(res_end, starts_with('d') & ends_with('spei12')))
colnames(res_dat_spei) <- c('res_end', word(colnames(res_dat_spei)[-1], 2, sep = '_'))
res_dat_spei$Basin_Avg <- rowMeans(select(res_dat_spei, !res_end))
for(huc in 1:length(huc_sel)){
  plot_dat <- select(res_dat_spei, c(res_end, huc_sel[huc]))
  plot <- ggplot(plot_dat, aes(x = res_end, fill = res_end, y = .data[[huc_sel[huc]]])) +
    geom_boxplot() +
    scale_fill_manual(values = res_palette, name = 'Season', guide = 'none') +
    scale_x_discrete(labels = NULL, name = NULL) +
    scale_y_continuous(name = NULL) +
    theme(text = element_text(size = 9)) +
    ggtitle('1-Year SPEI [-]')
  ggsave(paste0('./Figures/final/SI/6_res_spatial/',huc_sel[huc],'_spei.png'), plot = plot,
         width = 1.67, height = 1.67, units = 'in')
}
for(test_var in 1:length(huc_sel)){
  print(huc_sel[test_var])
  print(wilcox.test(res_dat_spei[[huc_sel[test_var]]][res_dat_spei$res_end == 'Precip'],
                    res_dat_spei[[huc_sel[test_var]]][res_dat_spei$res_end == 'Reservoir'])$p.value)
}

huc8s <- st_read('./DataFiles/shapefiles/created/huc8_boundaries.shp') %>%
  mutate(select = ifelse(HUC8 %in% c('10270104','10260008','10270201','10250002','10260003','10260012'),T,F))
ggplot() +
  geom_sf(data = huc8s, aes(fill = select), color = 'grey45') +
  scale_fill_manual(guide = 'none', values = c('grey90','steelblue')) +
  geom_sf_label(data = subset(huc8s, select == T), aes(label = HUC8)) +
  theme_void()
ggsave('./Figures/final/SI/6_res_spatial/basemap.png', width = 9, height = 4.5, units = 'in')





# Figure S7: Spatial correlations for 30d AP and HI -----------------------
library(tidyverse)
library(sf)
library(tmap)
source('8_combine_all_event_data.R')
source('Theme+Settings.R')

event_summary <- select(event_summary, !contains('global'))
#calculate correlations
huc8s <- st_read('./DataFiles/shapefiles/created/huc8_boundaries.shp')

climate_vars <- colnames(select(event_summary, 
                                starts_with('d_') & ends_with(c('spei12')),
                                ends_with('_ante') & contains(c('30d')),
                                starts_with('p_') & ends_with('1d_total')))
climate_var_order <- c('1d_total', '30d_ante', 'spei12')

res_var <- 'reservoir_runoff_ratio'
res_name <- 'Res/Runoff'

response_vars <- c('FI_n', 'HI_n')

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

#make figs
for(i in 1:length(huc8_dat)){
  map <- tm_shape(huc8_dat[[i]]) +
    tm_polygons('val', palette = 'div', breaks = seq(-0.32,0.32, by = 0.04),
                legend.show = TRUE, legend.reverse = TRUE, title = 'rho') +
    tm_facets('var', nrow = 1, ncol = 3) +
    tm_layout(panel.show = F) +
    tm_shape(huc8_sigs[[i]]) + tm_dots(col = 'val', palette = c(NA, 'blue'), size = 0.2,
                                       legend.show = FALSE) + tm_facets('var', nrow = 1, ncol = 3)
  # map
  tmap_save(map, filename = file.path('Figures','final','SI','7_correlation_maps',paste0(names(huc8_dat)[i],'.png')),
            height = 3, width = 12, units = 'in')
  print(paste0('Figure ',i,'/',length(huc8_dat),' Completed!'))
}



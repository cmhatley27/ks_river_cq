library(tidyverse)
library(cluster)  
library(factoextra)


# Load data ---------------------------------------------------------------


# Basic characteristics --------------------------------------------------

CQ_summary <- read_csv('./DataFiles/hydro_data/event_char/CQ_summary.csv') %>%
  select(FI, HI = HI_mean)

simple_temporal <- read_csv('./DataFiles/hydro_data/event_char/simple_temporal.csv')
simple_size <- read_csv('./DataFiles/hydro_data/event_char/simple_size.csv')

N_yield <- read_csv('./DataFiles/hydro_data/event_char/N_yield.csv') %>%
  cbind(., CQ_summary)

#* Gridded Precip Predictors -----------------------------------------------
precip_inevent_gridded <- read_csv('./DataFiles/hydro_data/event_char/precip_inevent_gridded.csv') %>%
  select(contains('1d')) %>%
  select(contains(c('10260008', '10260010', '10260015', '102701')))

precip_anteevent_gridded <- read_csv('./DataFiles/hydro_data/event_char/precip_anteevent_gridded.csv') %>%
  select(contains(c('10260008', '10260010', '10260015', '102701')))

precip_unintercepted_volumes <- read_csv('./DataFiles/hydro_data/event_char/precip_unintercepted_volumes.csv')

max_precip_coords <- read_csv('./DataFiles/hydro_data/event_char/max_precip_coords.csv')

avg_precip_intensity <- read_csv('./DataFiles/hydro_data/event_char/avg_precip_intensity.csv') %>%
  select(contains('1d')) %>%
  select(contains(c('10260008', '10260010', '10260015', '102701')))

max_precip_intensity <- read_csv('./DataFiles/hydro_data/event_char/max_precip_intensity.csv') %>%
  select(contains(c('10260008', '10260010', '10260015', '102701')))


# * Outflows Predictors ---------------------------------------------------
outflows_inevent_totals <- read_csv('./DataFiles/hydro_data/event_char/outflows_inevent_totals.csv') %>%
  select(contains('1d'))

outflows_deltas <- read_csv('./DataFiles/hydro_data/event_char/outflows_deltas.csv') %>%
  select(contains('1d'))

flow_ratios <- read_csv('./DataFiles/hydro_data/event_char/flow_ratios.csv')


# Drought Predictors ------------------------------------------------------
spei1 <- read_csv('./DataFiles/hydro_data/event_char/spei1.csv') %>%
  select(contains(c('10260015', '10270104')))

spei3 <- read_csv('./DataFiles/hydro_data/event_char/spei3.csv') %>%
  select(contains(c('10260015', '10270104')))

spei6 <- read_csv('./DataFiles/hydro_data/event_char/spei6.csv') %>%
  select(contains(c('10260015', '10270104')))

spei9 <- read_csv('./DataFiles/hydro_data/event_char/spei9.csv') %>%
  select(contains(c('10260015', '10270104')))

spei12 <- read_csv('./DataFiles/hydro_data/event_char/spei12.csv') %>%
  select(contains(c('10260015', '10270104')))

spei18 <- read_csv('./DataFiles/hydro_data/event_char/spei18.csv') %>%
  select(contains(c('10260015', '10270104')))

spei24 <- read_csv('./DataFiles/hydro_data/event_char/spei24.csv') %>%
  select(contains(c('10260015', '10270104')))

# * Correlated Predictors -------------------------------------------------

correlated_predictors <- unlist(as.vector(read_csv('./DataFiles/hydro_data/event_char/correlated_predictors.csv'))) 

#* Combine ----------------------------------------------

event_summary <- tibble(
  CQ_summary, 
  simple_temporal, 
  simple_size,
  precip_inevent_gridded,
  precip_anteevent_gridded,
  precip_unintercepted_volumes,
  max_precip_coords,
  avg_precip_intensity,
  max_precip_intensity,
  outflows_inevent_totals,
  outflows_deltas,
  flow_ratios,
  spei1,
  spei3,
  spei6,
  spei9,
  spei12,
  spei18,
  spei24
) %>%
  mutate(wateryear = factor(wateryear),
         month = factor(month),
         season = factor(season)) %>%
  select(!contains(correlated_predictors[correlated_predictors != 'reservoir_runoff_ratio']))

rm(list = ls()[ls() != 'event_summary' & ls() != 'N_yield'])


# Cluster by CQ Indices --------------------------------------------------

# * create clusters -------------------------------------------------------

hclust_data <- event_summary %>%
  filter(complete.cases(.)) %>%
  select(FI, HI)

dist <- daisy(hclust_data, metric = 'euclidean')

cls <- agnes(dist, method = "ward")

k <- 8

cls_id <- cutree(cls, k = k)

#Add cluster ID back in and add N yield data
event_summary_cluster <- event_summary %>%
  filter(complete.cases(.)) %>%
  mutate(cluster = cls_id,
         season = factor(season, levels = c('Winter', 'Spring', 'Summer', 'Fall'))) %>%
  left_join(., filter(N_yield, !is.na(initial_N) & !is.na(HI))) %>%
  mutate(loadrate = N_load/duration,
         log_loadrate = log10(loadrate),
         initial_loadrate = initial_Q*initial_N*60*60*24/1000,
         log_initial_loadrate = log10(initial_loadrate),
         delta_N = max_N - initial_N)

find_event <- event_summary %>%
  mutate(event_no = seq(1,nrow(event_summary))) %>%
  select(event_no, FI, HI) %>%
  left_join(., event_summary_cluster, by = c('FI', 'HI')) %>%
  select(event_no, cluster, FI, HI, reservoir_precip_ratio)


# FI HI scatter with large precip events/loads flagged ------------------------

#FI HI scatter with N loads and large precip events flagged
##Disturbing - highest load events and highest precip events do not show a pattern
ggplot(data = filter(event_summary_cluster, initial_N > 0)) +
  geom_point(aes(x = FI, y = HI, color = reservoir_precip_ratio)) +
  #geom_point(data = filter(event_summary_cluster, p_10270104_1d_total > quantile(event_summary_cluster$p_10270104_1d_total, 0.8)),
             #aes(x = FI, y = HI), shape = 1, color = 'red', size = 2) +
  # geom_point(data = filter(event_summary_cluster, d_10260015_spei12 > quantile(event_summary_cluster$d_10260015_spei12, 0.9)),
  #            aes(x = FI, y = HI), shape = 1, color = 'red', size = 2) +
  scale_color_viridis_b(end = 0.95, breaks = seq(0,0.5, by = 0.1)) +
  ylab('') +
  xlim(-1,1) +
  ylim(-1,1) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0)

ggplot(data = event_summary_cluster) +
  geom_point(aes(x = FI, y = HI, color = factor(cluster))) +
  xlim(-1,1) +
  ylim(-1,1) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0)

  


# box plots by cluster ----------------------------------------------------

##distributions of loads per cluster
ggplot(data = event_summary_cluster) +
  geom_boxplot(aes(x = cluster, y = log_loadrate, group = cluster, fill = factor(cluster))) +
  scale_fill_discrete(guide = 'none') +
  scale_x_continuous(breaks = seq(8)) +
  xlab('Cluster') +
  ylab('Load / Day [log(kg)/d]')

##distributions of N dynamics (initial, max, and diff)
event_summary_cluster %>%
  select(cluster, initial_N, max_N, delta_N) %>%
  filter(!is.na(.[,4])) %>%
  pivot_longer(!cluster, names_to = 'var', values_to = 'value') %>%
  group_by(var) %>%
  filter(value <= IQR(value)+quantile(value, 0.75) &
           value >= quantile(value, 0.25)-IQR(value)) %>%
  ggplot() +
  geom_boxplot(aes(y = value, x = cluster, fill = factor(cluster))) +
  facet_wrap(vars(var)) +
  scale_fill_discrete(guide = 'none') +
  scale_x_continuous(breaks = seq(8)) +
  scale_y_continuous() +
  xlab('Cluster') +
  ylab('Conc')
ggsave('./Figures/boxplots/k8/spei12_west.tiff', device = 'tiff', compression = 'lzw', height = 3, width = 3, units = 'in')

##Single Var box plot
event_summary_cluster %>%
  select(cluster, value = d_10260015_spei12) %>%
  filter(!is.na(.[,2])) %>%
  filter(value <= IQR(value)+quantile(value, 0.75) &
           value >= quantile(value, 0.25)-IQR(value)) %>%
  ggplot() +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_boxplot(aes(y = value, x = cluster, fill = factor(cluster))) +
  scale_fill_discrete(guide = 'none') +
  scale_x_continuous(breaks = seq(8)) +
  scale_y_continuous() +
  xlab('Cluster') +
  ylab('SPEI 12')
#ggsave('./Figures/boxplots/k8/spei12_west.tiff', device = 'tiff', compression = 'lzw', height = 3, width = 3, units = 'in')


variable <- 'p_10260015_365d_ante'
clus1 <- 1
clus2 <- 3
wilcox.test(as_vector(subset(event_summary_cluster[,variable], event_summary_cluster$cluster == clus1)), 
            as_vector(subset(event_summary_cluster[,variable], event_summary_cluster$cluster == clus2)))

##distributions of spei per cluster
event_summary_cluster %>%
  select(cluster, starts_with('d_')) %>%
  pivot_longer(!cluster, names_to = 'var', values_to = 'spei') %>%
  ggplot() +
  geom_boxplot(aes(x = cluster, y = spei, group = cluster, fill = factor(cluster))) +
  facet_wrap(vars(var))



## Number high load/precip events per cluster -------------------------------------
#interesting characteristics:
#log_loadrate, initial_N, max_N, d_10260015_spei12
characteristic <- event_summary_cluster$d_10260015_spei12
p <- 0.25
high_events <- mutate(event_summary_cluster,
                      high_value = characteristic <= quantile(characteristic, p, na.rm = TRUE)) %>%
               filter(!is.na(characteristic)) %>%
               group_by(cluster) %>%
               summarise(n = n(),
                         n_high = sum(high_value)) %>%
               mutate(pct = n_high/n)

ggplot(data = high_events) +
  geom_col(aes(x = cluster, y = n), fill = 'skyblue', color = 'black') +
  geom_col(aes(x = cluster, y = n_high), fill = 'coral2', color = 'black')

ggplot(data = high_events) +
  geom_col(aes(x = cluster, y = 1), fill = 'skyblue', color = 'black') +
  geom_col(aes(x = cluster, y = pct), fill = 'coral2', color = 'black') +
  ylab('% of events in the upper percentile')


# Distributions of characteristics per load/conc regimes ------------------
group_characteristic <- 'd_10260015_spei12'
group_characteristic_name <- 'SPEI Antecedent Moisture Conditions'
p_ends <- c(0.1,0.9)
end_dists <- mutate(event_summary_cluster,
                    end = ifelse(event_summary_cluster[,group_characteristic] <= quantile(event_summary_cluster[,group_characteristic], p_ends[1], na.rm = TRUE), 'low', 'mid'),
                    end = ifelse(event_summary_cluster[,group_characteristic] >= quantile(event_summary_cluster[,group_characteristic], p_ends[2], na.rm = TRUE), 'high', end))

select(end_dists, end,
       'Flushing Index' = FI,
       'Reservoir Precip Ratio' = reservoir_precip_ratio,
       'Reservoir Runoff Ratio' = reservoir_runoff_ratio,
       'N Load / Day' = log_loadrate) %>%
  pivot_longer(!end, names_to = 'var', values_to = 'value') %>%
  group_by(var) %>%
  filter(!is.na(value)) %>%
  filter(value <= IQR(value)+quantile(value, 0.75) &
           value >= quantile(value, 0.25)-IQR(value)) %>%
  ggplot() +
    geom_boxplot(aes(x = end, y = value, fill = end)) +
    scale_x_discrete(limits = c('Drought' = 'low', 'Flood' = 'high')) +
    scale_fill_manual(limits = c('low', 'high'), values = c('coral2', 'skyblue')) +
    xlab(group_characteristic_name) +
    #xlab(paste0(group_characteristic_name, ' Regime (q = ', p_ends[1],', ',p_ends[2],')')) +
    facet_wrap(vars(var), scales = 'free_y')
                                 
ggplot(data = event_summary_cluster) +
  geom_point(aes(x = initial_N, y = max_N, color = FI)) +
  scale_color_gradient2()


# Other plots w relationships between FI, antecedent, loads ---------------
##What are relationships between FI, antecedent conditions, and loads?

#load vs initial conditions
ggplot(data = event_summary_cluster, aes(x = log10(initial_Q), y = log_loadrate)) +
  geom_point()

ggplot(data = event_summary_cluster, aes(x = log10(initial_Q), y = log10(total_Q))) +
  geom_point()
ggplot(data = event_summary_cluster, aes(x = log10(initial_N), y = log10(N_load))) +
  geom_point()


ggplot(data = event_summary_cluster, aes(x = log_initial_loadrate, y = log_loadrate)) +
  geom_point()
ggplot(data = event_summary_cluster, aes(x = d_10260015_spei12, y = log_loadrate)) +
  geom_point()
select(event_summary_cluster, log_loadrate, d_10260015_spei12) %>%
  cor(., use = 'pairwise.complete.obs')

#load vs FI
##for dry initial conditions, FI has weak (r = 0.33) correlation with load
## for wet initial conditions, loads are high no matter FI
ggplot(data = event_summary_cluster, aes(x = FI, y = log_loadrate, color = log10(initial_Q))) +
  geom_point() +
  scale_color_viridis_b(end = 0.7, breaks = median(log10(event_summary_cluster$initial_Q)))
ggplot(data = event_summary_cluster, aes(x = FI, y = log_loadrate, color = d_10260015_spei12)) +
  geom_point() +
  scale_color_steps2(breaks = c(-0.5,0,0.5))

#FI vs initial conditions
ggplot(data = event_summary_cluster, aes(x = log10(initial_Q), y = FI)) +
  geom_point()
ggplot(data = event_summary_cluster, aes(x = d_10260015_spei12, y = FI)) +
  geom_point()


#N load most correlated with initial conditions (initial daily load broken by season)
ggplot(data = event_summary_cluster) +
  geom_point(aes(x = initial_loadrate, y = loadrate, color = season)) +
  geom_smooth(aes(x = initial_loadrate, y = loadrate, color = season), method = 'lm', se = FALSE) +
  scale_color_manual(values = c('#8DA9C4','#60992D','#C33C54','#F4B860'), guide = 'none')
ggplot(data = filter(event_summary_cluster, initial_N > 0)) +
  geom_point(aes(x = initial_loadrate, y = loadrate, color = wateryear)) +
  geom_smooth(aes(x = initial_loadrate, y = loadrate, color = wateryear), method = 'lm', se = FALSE)
ggplot(data = filter(event_summary_cluster, initial_N > 0)) +
  geom_point(aes(x = log10(initial_Q*initial_N*60*60*24/1000), y = log10(N_load/duration), color = wateryear)) +
  geom_smooth(aes(x = log10(initial_Q*initial_N*60*60*24/1000), y = log10(N_load/duration), color = wateryear), method = 'lm', se = FALSE)

#Difference between concentrating and diluting events on Load vs. Initial Load relationship
ggplot(data = filter(event_summary_cluster, initial_N > 0)) +
  geom_point(aes(x = initial_Q*initial_N*60*60*24/1000, y = (N_load/duration), color = FI)) +
  scale_color_steps2(breaks = c(-0.5,0,0.5), low = 'blue', high = 'red') +
  geom_abline(slope = 1) +
  geom_smooth(data = filter(event_summary_cluster, FI > 0.5), aes(x = initial_Q*initial_N*60*60*24/1000, y = (N_load/duration)), method = 'lm', se = FALSE, color = 'red') +
  geom_smooth(data = filter(event_summary_cluster, FI < -0.5), aes(x = initial_Q*initial_N*60*60*24/1000, y = (N_load/duration)), method = 'lm', se = FALSE, color = 'blue')
#Same but logged
ggplot(data = filter(event_summary_cluster, initial_N > 0)) +
  geom_point(aes(x = log10(initial_Q*initial_N*60*60*24/1000), y = log10(N_load/duration), color = FI)) +
  scale_color_steps2(breaks = c(-0.5,0,0.5), low = 'blue', high = 'red') +
  geom_abline(slope = 1) +
  geom_smooth(data = filter(event_summary_cluster, FI > 0.5), aes(x = log10(initial_Q*initial_N*60*60*24/1000), y = log10(N_load/duration)), method = 'lm', se = FALSE, color = 'red') +
  geom_smooth(data = filter(event_summary_cluster, FI < -0.5), aes(x = log10(initial_Q*initial_N*60*60*24/1000), y = log10(N_load/duration)), method = 'lm', se = FALSE, color = 'blue')

#Load vs FI for each year
## Events during wet periods have high loads and variable FI
## Events during dry periods have lower loads and generally? higher FI
ggplot(data = filter(event_summary_cluster, initial_N > 0)) +
  geom_point(aes(x = FI, y = log10(N_load/duration), color = d_10260015_spei12)) +
  scale_color_steps2(breaks = c(-2,-1,0,1,2), low = 'red', high = 'blue') +
  facet_wrap(vars(wateryear)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 4)
#Swap spei and FI
ggplot(data = event_summary_cluster) +
  geom_point(aes(x = d_10270104_spei6, y = log10(N_load/duration), color = FI)) +
  scale_color_steps2(breaks = c(-0.5,0,0.5), low = 'blue', high = 'red') +
  #scale_color_viridis_b() +
  facet_wrap(vars(wateryear)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 4)

#Comparing different SPEI measurements
select(event_summary_cluster, initial_N, initial_Q, FI, starts_with('d_')) %>%
  pivot_longer(!c(initial_N,initial_Q,FI), names_to = 'spei', values_to = 'val') %>%
  ggplot(aes(x = val, y = log10(initial_N*initial_Q*60*60*24/1000), color = FI)) +
  scale_color_steps2(breaks = c(-0.5,0,0.5), low = 'blue', high = 'red') +
  geom_point() +
  facet_wrap(vars(spei))


#Average concentration vs FI for each wateryear
## similar pattern as loads bc loads correlated with concentration
ggplot(data = event_summary_cluster) +
  geom_point(aes(x = FI, y = (N_load/total_Q*1000), color = d_10260015_spei6)) +
  #scale_color_steps2(breaks = c(-0.5,0,0.5), low = 'blue', high = 'red') +
  scale_color_steps2(midpoint = c(0), breaks = c(-1,0,1)) +
  facet_wrap(vars(wateryear)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 1.13)


#Load vs event size, broken by diluting vs flushing events
ggplot(data = event_summary_cluster, aes(x = log10(total_Q/duration), y = log10(N_load/duration))) +
  geom_point(aes(color = FI)) +
  scale_color_steps2(breaks = 0, high = 'red', low = 'blue') +
  geom_smooth(data = subset(event_summary_cluster, FI < -0.17), method = 'lm', se = FALSE) +
  geom_smooth(data = subset(event_summary_cluster, FI > 0.5), method = 'lm', se = FALSE, color = 'red') +
  facet_wrap(vars(season))
#Broken by flood vs drought events
ggplot(data = event_summary_cluster, aes(x = log10(total_Q/duration), y = log10(N_load/duration))) +
  geom_point(aes(color = d_10270104_spei12)) +
  scale_color_steps2(breaks = c(-1,0,1,2), high = 'blue', low = 'red') +
  geom_smooth(data = subset(event_summary_cluster, d_10270104_spei12 < -0.812), method = 'lm', se = FALSE, color = 'red') +
  geom_smooth(data = subset(event_summary_cluster, d_10270104_spei12 > 1.821), method = 'lm', se = FALSE, color = 'blue')
lm(log10(N_load/duration)~log10(total_Q/duration), data = subset(event_summary_cluster, d_10270104_spei12 < -0.812))
lm(log10(N_load/duration)~log10(total_Q/duration), data = subset(event_summary_cluster, d_10270104_spei12 > 1.821))

#Concentration vs Event size (normal CQ stuff)
##Wet events much more chemostatic than drought events
ggplot(data = event_summary_cluster, aes(x = log10(total_Q/duration), y = log10(N_load/total_Q))) +
  geom_point(aes(color = d_10270104_spei12)) +
  scale_color_steps2(breaks = c(-1,0,1,2), high = 'blue', low = 'red') +
  geom_smooth(data = subset(event_summary_cluster, d_10270104_spei12 < -0.812), method = 'lm', se = FALSE, color = 'red') +
  geom_smooth(data = subset(event_summary_cluster, d_10270104_spei12 > 1.821), method = 'lm', se = FALSE, color = 'blue')
lm(log10(N_load/total_Q)~log10(total_Q/duration), data = subset(event_summary_cluster, d_10270104_spei12 < -0.812))
lm(log10(N_load/total_Q)~log10(total_Q/duration), data = subset(event_summary_cluster, d_10270104_spei12 > 1.821))


ggplot(data = event_summary_cluster, aes(x = (total_Q/duration), y = (N_load/total_Q))) +
  geom_point(aes(color = d_10270104_spei12)) +
  scale_color_steps2(breaks = c(-1,0,1,2), high = 'blue', low = 'red') +
  geom_smooth(data = subset(event_summary_cluster, d_10270104_spei12 < -0.812), method = 'lm', se = FALSE, color = 'red') +
  geom_smooth(data = subset(event_summary_cluster, d_10270104_spei12 > 1.821), method = 'lm', se = FALSE, color = 'blue')


#Average concentration vs FI
##Diluting events do indeed result in lower average concentrations than the initial
ggplot(data = event_summary_cluster) +
  geom_point(aes(x = FI, y = ((N_load/total_Q*1000)/initial_N)-1, color = initial_N)) +
  ylim(-1,1) +
  geom_hline(yintercept = 0)


ggplot(data = event_summary_cluster) +
  geom_boxplot(aes(x = factor(cluster), fill = factor(cluster), y = loadrate))
ggplot(data = event_summary_cluster) +
  geom_boxplot(aes(x = factor(cluster), fill = factor(cluster), y = initial_N))







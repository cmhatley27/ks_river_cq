
# Load libraries and data ----------------------------------------------------------
library(tidyverse)
library(cluster)  
library(factoextra)
library(reshape)
library(brunnermunzel)
library(ggpubr)
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

# Climate Reservoir Box Plot ---------------
climate_var <- 'p_10260010_30d_ante'
climate_name <- climate_var
climate_q <- c(0.25,0.75)
climate_quant <- quantile(unlist(event_summary[,climate_var]), probs = c(climate_q))

res_var <- 'reservoir_runoff_ratio'
res_name <- 'Res/Runoff'
res_q <- c(0.5,0.5)
res_quant <- quantile(unlist(event_summary[,res_var]), probs = c(res_q))

response_vars <- c('FI_n', 'max_N')

mid_group <- FALSE

#Event Characteristics
# response_vars <- colnames(event_summary_cluster %>%
#                             select(p_10270104_1d_total, duration_since_last, duration, initial_Q, delta_rat_Q, HI, 
#                                    r_t_1d_delta, r_m_1d_delta, r_t_share, r_m_share, initial_N, p_10260015_1d_total, r_c_1d_delta, r_c_share,
#                                    r_p_share, r_p_1d_delta, reservoir_runoff_ratio,
#                                    FI))


plot_dat <- event_summary %>%
  mutate(climate_end = quantize(select(.,all_of(climate_var)), climate_q),
         res_end = quantize(select(.,all_of(res_var)), res_q)) %>%
  pivot_longer(all_of(response_vars), names_to = 'var', values_to = 'val') %>%
  select(var, val, climate_end, res_end)

if(mid_group == FALSE) {
  plot_dat <- filter(plot_dat, climate_end != 'mid')
  #Climate boxes
  clim <- ggplot(data = plot_dat, aes(x = climate_end, y = val)) +
    geom_boxplot(aes(fill = climate_end)) +
    scale_fill_manual(limits = c('low', 'high'), values = c('coral2', 'skyblue'), na.value = NA,
                      labels = c(paste0('< ',climate_q[1]*100, 'th'), paste0('> ',climate_q[2]*100, 'th')),
                      name = paste0(climate_name, ' Percentile')) +
    scale_x_discrete(limits = c('low', 'high'), labels = c(paste0('< ',climate_q[1]*100, 'th'), paste0('> ',climate_q[2]*100, 'th')),
                     name = paste(climate_name, ' Percentile')) +
    facet_wrap(vars(var), scales = 'free_y') +
    ylab(NULL)
  #Reservoir boxes
  res <- ggplot(data = subset(plot_dat, res_end != 'mid'), aes(x = res_end, y = val)) +
    geom_boxplot(aes(fill = climate_end)) +
    geom_vline(xintercept = 1.5, linetype = 'dashed') +
    scale_fill_manual(limits = c('low', 'high'), values = c('coral2', 'skyblue'), na.value = NA,
                      labels = c(paste0('< ',climate_q[1]*100, 'th'), paste0('> ',climate_q[2]*100, 'th')),
                      name = paste0(climate_name, ' Percentile')) +
    scale_x_discrete(limits = c('low', 'high'), labels = c('Precip-Driven', 'Reservoir-Driven'),
                     name = paste(climate_name, ' Percentile')) +
    facet_wrap(vars(var), scales = 'free_y') +
    ylab(NULL)
}

if(mid_group == TRUE){
  #Plot only climate boxes
  if(climate_q[1] != 0.5) {
    clim <- ggplot(data = plot_dat, aes(x = climate_end, y = val)) +
      geom_boxplot(aes(fill = climate_end)) +
      scale_fill_manual(limits = c('low', 'mid', 'high'), values = c('coral2', '#E9ECEE', 'skyblue'), na.value = NA,
                        labels = c(paste0('< ',climate_q[1]*100, 'th'), paste0(climate_q[1]*100,'th - ',climate_q[2]*100,'th'), paste0('> ',climate_q[2]*100, 'th')),
                        name = paste0(climate_name, ' Percentile')) +
      scale_x_discrete(limits = c('low', 'mid', 'high'), labels = c(paste0('< ',climate_q[1]*100, 'th'), paste0(climate_q[1]*100,'th - ',climate_q[2]*100,'th'), paste0('> ',climate_q[2]*100, 'th')),
                       name = paste(climate_name, ' Percentile')) +
      facet_wrap(vars(var), scales = 'free_y') +
      ylab(NULL)
  } else {
    clim <- ggplot(data = plot_dat, aes(x = climate_end, y = val)) +
      geom_boxplot(aes(fill = climate_end)) +
      scale_fill_manual(limits = c('low', 'mid', 'high'), values = c('coral2', '#E9ECEE', 'skyblue'), na.value = NA,
                        labels = c(paste0('< ',climate_q[1]*100, 'th'), paste0(climate_q[1]*100,'th - ',climate_q[2]*100,'th'), paste0('> ',climate_q[2]*100, 'th')),
                        name = paste0(climate_name, ' Percentile')) +
      scale_x_discrete(limits = c('low', 'high'), labels = c(paste0('< ',climate_q[1]*100, 'th'), paste0('> ',climate_q[2]*100, 'th')),
                       name = paste0(climate_name, ' Percentile')) +
      facet_wrap(vars(var), scales = 'free_y')
  }
  
  #Plot reservoir boxes
  if(climate_q[1] != 0.5) {
    res <- ggplot(data = plot_dat, aes(x = res_end, y = val)) +
      geom_boxplot(aes(fill = climate_end)) +
      scale_fill_manual(limits = c('low', 'mid', 'high'), values = c('coral2', '#E9ECEE', 'skyblue'), na.value = NA,
                        labels = c(paste0('< ',climate_q[1]*100, 'th'), paste0(climate_q[1]*100,'th - ',climate_q[2]*100,'th'), paste0('> ',climate_q[2]*100, 'th')),
                        name = paste0(climate_name, ' Percentile')) +
      scale_x_discrete(limits = c('low', 'mid', 'high'), labels = c('Precip', 'Mixed', 'Res'),
                       name = paste(climate_name, ' Percentile')) +
      geom_vline(xintercept = c(1.5,2.5), linetype = 'dashed') +
      facet_wrap(vars(var), scales = 'free_y') +
      ylab(NULL)
  } else {
    res <- ggplot(data = plot_dat, aes(x = res_end, y = val)) +
      geom_boxplot(aes(fill = climate_end)) +
      scale_fill_manual(limits = c('low', 'mid', 'high'), values = c('coral2', '#E9ECEE', 'skyblue'), na.value = NA,
                        labels = c(paste0('< ',climate_q[1]*100, 'th'), paste0(climate_q[1]*100,'th - ',climate_q[2]*100,'th'), paste0('> ',climate_q[2]*100, 'th')),
                        name = paste0(climate_name, ' Percentile')) +
      scale_x_discrete(limits = c('low', 'high'), labels = c('Precip', 'Res'),
                       name = paste0(climate_name, ' Percentile')) +
      geom_vline(xintercept = 1.5, linetype = 'dashed') +
      facet_wrap(vars(var), scales = 'free_y')
  }
}
clim
res

#Calculate significant differences
climate_combs <- combn(unique(plot_dat$climate_end), 2)
climate_sig <- tibble(
  response_var = rep(response_vars, each = ncol(climate_combs)),
  climate1 = rep(climate_combs[1,], length.out = length(response_var)),
  climate2 = rep(climate_combs[2,], length.out = length(response_var)),
  p = NA,
  climate_var = climate_var
)

for(combo in 1:nrow(climate_sig)){
  dat_fil <- filter(plot_dat, var == unlist(climate_sig[combo, 'response_var']))
  
  climate_sig[combo, 'p'] <- round(brunnermunzel.test(filter(dat_fil, climate_end == unlist(climate_sig[combo,'climate1']))$val,
                                                      filter(dat_fil, climate_end == unlist(climate_sig[combo,'climate2']))$val)$p.value,
                                   digits = 3)
}


res_combs <- combn(unique(plot_dat$res_end), 2)
res_sig <- tibble(
  response_var = rep(response_vars, each = ncol(res_combs)*length(unique(plot_dat$climate_end))),
  climate = rep(rep(unique(plot_dat$climate_end), each = ncol(res_combs)), length.out = length(response_var)),
  res1 = rep(res_combs[1,], length.out = length(response_var)),
  res2 = rep(res_combs[2,], length.out = length(response_var)),
  p = NA,
  climate_var = climate_var
)

for(combo in 1:nrow(res_sig)){
  dat_fil <- filter(plot_dat, var == unlist(res_sig[combo, 'response_var']) & climate_end == unlist(res_sig[combo, 'climate']))
  
  res_sig[combo, 'p'] <- round(brunnermunzel.test(filter(dat_fil, res_end == unlist(res_sig[combo,'res1']))$val,
                                                  filter(dat_fil, res_end == unlist(res_sig[combo,'res2']))$val)$p.value,
                               digits = 3)
}

climate_combs <- combn(unique(plot_dat$climate_end), 2)
climate_res_sig <- tibble(
  response_var = rep(response_vars, each = ncol(climate_combs)*length(unique(plot_dat$climate_end))),
  res = rep(rep(unique(plot_dat$res_end), each = ncol(climate_combs)), length.out = length(response_var)),
  clim1 = rep(climate_combs[1,], length.out = length(response_var)),
  clim2 = rep(climate_combs[2,], length.out = length(response_var)),
  p = NA,
  climate_var = climate_var
)

for(combo in 1:nrow(climate_res_sig)){
  dat_fil <- filter(plot_dat, var == unlist(res_sig[combo, 'response_var']) & res_end == unlist(climate_res_sig[combo, 'res']))
  
  climate_res_sig[combo, 'p'] <- round(brunnermunzel.test(filter(dat_fil, climate_end == unlist(climate_res_sig[combo,'clim1']))$val,
                                                  filter(dat_fil, climate_end == unlist(climate_res_sig[combo,'clim2']))$val)$p.value,
                               digits = 3)
}

#Print distribution summary
print(climate_quant)
print(res_quant)
print(table(select(plot_dat,climate_end))/length(response_vars))
print(table(select(plot_dat, climate_end, res_end))/length(response_vars))

# FI HI scatter -----------------------------------------------------------
climate_var <- 'p_10260015_30d_ante'
climate_name <- '12 Month SPEI (West & East Mean)'
climate_q <- c(0.5,0.5)

res_var <- 'reservoir_runoff_ratio'
res_name <- 'Res/Runoff'
res_q <- c(0.5,0.5)

plot_dat <- event_summary_cluster %>%
  mutate(climate_end = quantize(select(.,all_of(climate_var)), climate_q),
         res_end = quantize(select(.,all_of(res_var)), res_q),
         grouping = factor(paste0(climate_end,'-',res_end),
                           levels = c('low-high', 'mid-high', 'high-high',
                                      'low-mid', 'mid-mid', 'high-mid',
                                      'low-low', 'mid-low', 'high-low'),
                           labels = c('dry-reservoir', 'avg-reservoir', 'wet-reservoir',
                                      'dry-mix', 'avg-mix', 'wet-mix',
                                      'dry-precip', 'avg-precip', 'wet-precip'))) %>%
  select(FI,HI,climate_end,res_end,grouping)


group_means <- group_by(plot_dat, grouping) %>%
  summarise(FI_mean = mean(FI),
            HI_mean = mean(HI))

ggplot(data = plot_dat, aes(x = FI, y = HI, color = grouping)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point() +
  geom_point(data = group_means, aes(x = FI_mean, y = HI_mean, color = grouping), shape = 15, size = 3) +
  facet_wrap(vars(grouping))

ggplot(data = event_summary_cluster, aes(x = FI, y = HI, color = quantize(select(event_summary_cluster,all_of(climate_var)), climate_q))) +
  geom_point() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  facet_wrap(vars(quantize(reservoir_runoff_ratio, res_q), season), nrow = 2)



# Significant Differences for each climate category -----------------------
climate_indicators <- colnames(select(event_summary_cluster, starts_with('d_'), ends_with('_ante'), initial_Q))

climate_summary <- tibble()
res_summary <- tibble()

for(indicator in 1:length(climate_indicators)){
  climate_var <- climate_indicators[indicator]
  climate_q <- c(0.25,0.75)
  
  res_var <- 'reservoir_runoff_ratio'
  res_q <- c(0.25,0.75)
  
  response_vars <- c('FI', 'HI', 'log_loadrate', 'delta_N')
  
  plot_dat <- event_summary_cluster %>%
    mutate(climate_end = quantize(select(.,all_of(climate_var)), climate_q),
           res_end = quantize(select(.,all_of(res_var)), res_q)) %>%
    pivot_longer(all_of(response_vars), names_to = 'var', values_to = 'val') %>%
    select(var, val, climate_end, res_end)

  #Calculate significant differences
  climate_combs <- combn(unique(plot_dat$climate_end), 2)
  climate_sig <- tibble(
    response_var = rep(response_vars, each = ncol(climate_combs)),
    climate1 = rep(climate_combs[1,], length.out = length(response_var)),
    climate2 = rep(climate_combs[2,], length.out = length(response_var)),
    p = NA,
    stat = NA,
    climate_var = climate_var
  )
  
  for(combo in 1:nrow(climate_sig)){
    dat_fil <- filter(plot_dat, var == unlist(climate_sig[combo, 'response_var']))
    test <- brunnermunzel.test(filter(dat_fil, climate_end == unlist(climate_sig[combo,'climate1']))$val,
                               filter(dat_fil, climate_end == unlist(climate_sig[combo,'climate2']))$val)
    
    climate_sig[combo, 'p'] <- round(test$p.value, digits = 3)
    climate_sig[combo, 'stat'] <- round(test$statistic, digits = 3)
  }
  
  climate_summary <- rbind(climate_summary, climate_sig)
  
  
  res_combs <- combn(unique(plot_dat$res_end), 2)
  res_sig <- tibble(
    response_var = rep(response_vars, each = ncol(res_combs)*length(unique(plot_dat$climate_end))),
    climate = rep(rep(unique(plot_dat$climate_end), each = ncol(res_combs)), length.out = length(response_var)),
    res1 = rep(res_combs[1,], length.out = length(response_var)),
    res2 = rep(res_combs[2,], length.out = length(response_var)),
    p = NA,
    stat = NA,
    climate_var = climate_var
  )
  
  for(combo in 1:nrow(res_sig)){
    dat_fil <- filter(plot_dat, var == unlist(res_sig[combo, 'response_var']) & climate_end == unlist(res_sig[combo, 'climate']))
    test <- brunnermunzel.test(filter(dat_fil, res_end == unlist(res_sig[combo,'res1']))$val,
                               filter(dat_fil, res_end == unlist(res_sig[combo,'res2']))$val)
    res_sig[combo, 'p'] <- round(test$p.value, digits = 3)
    res_sig[combo, 'stat'] <- round(test$statistic, digits = 3)
  }
  
  res_summary <- rbind(res_summary, res_sig)
}

filter(res_summary, p <= 0.05) %>%
  group_by(climate_var) %>%
  summarise(n = n()) %>%
  arrange(desc(n)) %>%
  mutate(combos = ((length(unique(climate_q))+1)*(length(unique(res_q))+1)*length(response_vars)),
         pct = n/combos)



# Box Plots for ALL climate variables -------------------------------------
climate_vars <- colnames(select(event_summary_cluster, starts_with('d_'), ends_with('_ante'), initial_Q))
climate_q <- c(0.25,0.75)

res_var <- 'reservoir_runoff_ratio'
res_name <- 'Res/Runoff'
res_q <- c(0.25,0.75)

response_var <- 'FI'

plot_dat <- select(water_bal, all_of(c(response_var, res_var, climate_vars))) %>%
  dplyr::rename(response_var = response_var) %>%
  mutate(res_end = quantize(select(.,res_var), res_q)) %>%
  mutate(across(all_of(climate_vars), ~quantize(.x, climate_q))) %>%
  pivot_longer(c(all_of(climate_vars)), names_to = 'climate_var', values_to = 'climate_end')

ggplot(data = subset(plot_dat, climate_end != 'mid'), aes(x = climate_end, y = response_var)) +
  geom_boxplot(aes(fill = climate_end)) +
  facet_wrap(vars(climate_var)) +
  ylab(response_var)

ggplot(data = subset(plot_dat, climate_end != 'mid' & res_end != 'mid'), aes(x = climate_end, y = response_var)) +
  geom_boxplot(aes(fill = climate_end, alpha = res_end)) +
  facet_wrap(vars(climate_var)) +
  ylab(response_var)


# Median response variable value across ALL climate variables --------------
climate_vars <- colnames(select(event_summary_cluster, starts_with('d_'), 
                                ends_with('_ante'), initial_Q))


climate_q <- c(0.25,0.75)

res_var <- 'reservoir_runoff_ratio'
res_name <- 'Res/Runoff'
res_q <- c(0.25,0.75)

response_var <- 'delta_N'

plot_dat <- select(water_bal, all_of(c(response_var, res_var, climate_vars))) %>%
  dplyr::rename(response_var = response_var) %>%
  mutate(res_end = quantize(select(.,res_var), res_q)) %>%
  mutate(across(all_of(climate_vars), ~quantize(.x, climate_q))) %>%
  pivot_longer(c(all_of(climate_vars)), names_to = 'climate_var', values_to = 'climate_end')


ggplot(data = subset(plot_dat, climate_end != 'mid'), aes(x = climate_end, y = response_var)) +
  geom_boxplot(aes(fill = climate_end))

ggplot(data = subset(plot_dat, climate_end != 'mid' & res_end != 'mid'), aes(x = climate_end, y = response_var)) +
  geom_boxplot(aes(fill = climate_end, alpha = res_end)) +
  ylab(response_var)


climate_medians <- group_by(plot_dat, climate_end) %>%
  summarise(median = median(response_var, na.rm = TRUE))
res_medians <- group_by(plot_dat, climate_end, res_end) %>%
  summarise(median = median(response_var, na.rm = TRUE))


#B-M tests probably not legit bc pooling events from multiple climate indicators
#means that the groupings are not independent - e.g. an event in the 'dry' group
#according to SPEI 12 will most likely also be put in the 'dry' group by 12-month AP.
#Essentially inflating the number of samples in each group to reduce the p-value.
brunnermunzel.test(plot_dat$response_var[plot_dat$climate_end == 'low'], plot_dat$response_var[plot_dat$climate_end == 'high'])$p.value
brunnermunzel.test(plot_dat$response_var[plot_dat$climate_end == 'low' & plot_dat$res_end == 'low'], 
                   plot_dat$response_var[plot_dat$climate_end == 'low' & plot_dat$res_end == 'high'])$p.value
brunnermunzel.test(plot_dat$response_var[plot_dat$climate_end == 'high' & plot_dat$res_end == 'low'], 
                   plot_dat$response_var[plot_dat$climate_end == 'high' & plot_dat$res_end == 'high'])$p.value

# Compare co-occurence of climate categories ------------------------------
ggplot(data = water_bal) +
  geom_jitter(aes(x = spei3_w, y = spei3_e, color = cut(reservoir_runoff_ratio, breaks = quantile(reservoir_runoff_ratio, c(0,0.25,0.75,1)))), width = 0.1) +
  geom_abline(slope = 1) +
  geom_hline(yintercept = quantile(water_bal$spei3_w, c(0.25,0.75))) +
  geom_vline(xintercept = quantile(water_bal$spei3_e, c(0.25,0.75))) +
  scale_color_discrete(name = 'res_cat')

table(
  water_bal %>%
    mutate(res_cat = cut(water_bal$reservoir_runoff_ratio, breaks = quantile(water_bal$reservoir_runoff_ratio, c(0,0.25,0.75,1))),
           clim_cat = cut(water_bal$spei18_w, breaks = quantile(water_bal$spei18_w, c(0,0.25,0.75,1)))) %>%
    # filter(spei18_w < quantile(spei18_w, 0.75),
    #        spei1_w > quantile(spei1_w, 0.75)) %>%
    select(res_cat, clim_cat, wateryear)
)

view(
  water_bal %>%
    mutate(res_cat = cut(water_bal$reservoir_runoff_ratio, breaks = quantile(water_bal$reservoir_runoff_ratio, c(0,0.25,0.75,1)),
                         labels = c('precip','mix','res')),
           clim_cat = cut(water_bal$spei18_w, breaks = quantile(water_bal$spei18_w, c(0,0.25,0.75,1)),
                          labels = c('dry', 'avg', 'wet'))) %>%
    # filter(spei18_w < quantile(spei18_w, 0.75),
    #        spei1_w > quantile(spei1_w, 0.75)) %>%
    select(res_cat, clim_cat, wateryear, month, season, initial_Q, total_Q)
)

water_bal %>%
  mutate(res_cat = cut(water_bal$reservoir_runoff_ratio, breaks = quantile(water_bal$reservoir_runoff_ratio, c(0,0.25,0.75,1)),
                       labels = c('precip','mix','res')),
         clim_cat = cut(water_bal$spei1_w, breaks = quantile(water_bal$spei1_w, c(0,0.25,0.75,1)),
                        labels = c('dry', 'avg', 'wet'))) %>%
  # filter(spei18_w < quantile(spei18_w, 0.75),
  #        spei1_w > quantile(spei1_w, 0.75)) %>%
  select(res_cat, clim_cat, wateryear, month, season, initial_Q, total_Q, delta_Q, delta_rat_Q, duration) %>%
  group_by(clim_cat, res_cat) %>%
  summarise(across(c(ends_with('_Q'), duration), median))

select(water_bal, spei3_w, spei3_e) %>% cor(.)
summary(lm(spei3_e ~ spei3_w, data = water_bal))

select(water_bal, ends_with('_ante')) %>%
  select(contains('10270104'), contains('10260015')) %>%
  mutate(event = seq(1:nrow(water_bal))) %>%
  pivot_longer(!event, names_to = 'var', values_to = 'ap') %>%
  ggplot() +
  geom_histogram(aes(x = ap)) +
  facet_wrap(vars(var))



# rational box plots ------------------------------------------------------
quant_vars <- c('reservoir_precip_ratio', 'd_10260015_spei12', 'p_10260015_365d_ante')
q <- c(0.25,0.75)

dat <- event_summary_cluster %>%
  mutate(across(quant_vars, ~quantize(.,q), .names = '{.col}'))

ggplot(subset(dat, d_10260015_spei12_end != 'mid')) +
  geom_boxplot(aes(x = d_10260015_spei12_end, y = FI))

ggplot(subset(dat, p_10260015_365d_ante_end != 'mid')) +
  geom_boxplot(aes(x = p_10260015_365d_ante_end, y = FI))

ggplot(subset(dat, p_10260015_365d_ante_end != 'mid')) +
  geom_boxplot(aes(x = p_10260015_365d_ante_end, y = delta_N)) +
  scale_y_continuous(limits = c(-2,2))


plot_dat <- event_summary %>%
  select(FI_n, season, reservoir_runoff_ratio) %>%
  mutate(reservoir_end = quantize(reservoir_runoff_ratio, c(0.5,0.5))) %>%
  filter(!is.na(FI_n)) %>%
  group_by(season) %>%
  summarise(p = brunnermunzel.test(FI_n[reservoir_end == 'low'], FI_n[reservoir_end == 'high'])$p.value,
            stat = brunnermunzel.test(FI_n[reservoir_end == 'low'], FI_n[reservoir_end == 'high'])$statistic)

brunnermunzel.test(plot_dat$FI_n[plot_dat$reservoir_end == 'low'], plot_dat$FI_n[plot_dat$reservoir_end == 'high'])


ggplot(data = event_summary, aes(x = FI_n, y = max_N)) +
  geom_point()
cor(event_summary$initial_N, event_summary$max_N, use = 'pairwise.complete.obs', method = 'spearman')

event_summary_fil <- filter(event_summary, !is.na(cluster))
table(event_summary_fil$wateryear)
summary(event_summary_fil$duration_since_last)
summary(count(event_summary_fil, wateryear)$n)
ggplot(event_summary_fil, aes(x = reservoir_runoff_ratio)) + geom_histogram()

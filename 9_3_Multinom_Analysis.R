library(tidyverse)
library(nnet)


# Load Data ---------------------------------------------------------------


# Basic Characteristics

CQ_summary <- read_csv('./DataFiles/hydro_data/event_char/CQ_summary.csv') %>%
  select(FI, HI = HI_mean)

simple_temporal <- read_csv('./DataFiles/hydro_data/event_char/simple_temporal.csv')
simple_size <- read_csv('./DataFiles/hydro_data/event_char/simple_size.csv')

# Gridded Precip Predictors
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


# Outflow Predictors
outflows_inevent_totals <- read_csv('./DataFiles/hydro_data/event_char/outflows_inevent_totals.csv') %>%
  select(contains('1d'))

outflows_deltas <- read_csv('./DataFiles/hydro_data/event_char/outflows_deltas.csv') %>%
  select(contains('1d'))

flow_ratios <- read_csv('./DataFiles/hydro_data/event_char/flow_ratios.csv')


# Cluster Assignments

cls_id <- read_csv('C:/School/SAFE KAW/Data/DataFiles/cluster_data/cluster_assignment.csv')

# Correlated Predicors

correlated_predictors <- unlist(as.vector(read_csv('./DataFiles/hydro_data/event_char/correlated_predictors.csv'))) 

# Combine

event_summary_cluster <- tibble(
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
  flow_ratios
) %>%
  mutate(wateryear = factor(wateryear),
         month = factor(month),
         season = factor(season),
         spring = factor(ifelse(season == 'Spring', 1, 0))) %>%
  filter(complete.cases(.)) %>%
  mutate(cluster = cls_id$value) %>%
  #remove variables that are correlated > 0.8
  select(!contains(correlated_predictors)) %>%
  #remove vars of non-clinical importance
  select(!c(wateryear, month, yday))


#Set Reference categories
event_summary_cluster$cluster <- relevel(factor(event_summary_cluster$cluster, levels = seq(1:8)), ref = 8)
# event_summary_cluster$wateryear <- relevel(factor(event_summary_cluster$wateryear), ref = '2017')
event_summary_cluster$season <- relevel(event_summary_cluster$season, ref = 'Spring')


rm(list = ls()[ls() != 'event_summary_cluster'])



# Gather final model ------------------------------------------------------

sig_preds <- c('season', 'delta_Q', 'delta_rat_Q', 'duration', 'duration_since_last', 
               'initial_N', 'initial_Q', 'p_10260015_1d_total', 'p_10260015_maxint',
               'Eastern Precip Total' = 'p_10270104_1d_total', 'p_10270104_90d_ante', 'r_c_1d_delta', 'r_m_share',
               'r_t_1d_delta', 'reservoir_precip_ratio')

mod_frame <- select(event_summary_cluster, cluster, sig_preds) %>%
  mutate(initial_Q_log = log(ifelse(initial_Q <= 0, 0.01, initial_Q)),
         r_t_1d_delta_log = log(ifelse(r_t_1d_delta <= 0, 0.01, r_t_1d_delta))) %>%
  select(!c(initial_Q, r_t_1d_delta))

mod_frame$cluster <- relevel(factor(mod_frame$cluster, levels = seq(1:8)), ref = 8)
mod <- multinom(cluster ~ 1 + . + 
                  delta_rat_Q*season + delta_rat_Q*duration + p_10260015_1d_total*reservoir_precip_ratio,
                data = mod_frame)

mod_coef <- (exp(summary(mod)$coefficients))

mod_std <- summary(mod)$standard.errors
mod_wald <- abs(summary(mod)$coefficients/mod_std)
mod_wald_p <- 2*(1 - pnorm(mod_wald))
mod_sig <- ifelse(mod_coef > 0, 1, -1) * ifelse(abs(mod_wald_p) > 0.05, 0, 1)



# Pairwise differences in full data v model vars --------------------------
event_summary_cluster$cluster <- relevel(factor(event_summary_cluster$cluster, levels = seq(1:8)), ref = 1)
# Get sig pairwise differences
#Calculate KS statistics
event_cdfs <- event_summary_cluster %>%
  select(!c(FI, HI, season, spring)) %>%
  pivot_longer(!cluster, names_to = 'var', values_to = 'value') %>%
  group_by(cluster, var) %>%
  arrange((value)) %>%
  mutate(rank = seq(1:length(value)),
         P = rank/(length(value) + 1)*100) %>%
  ungroup()

k <- 8
combs <- combn(k, 2)

cluster_ks <- tibble(
  var = sort(unique(event_cdfs$var))
)

for(comb in seq(ncol(combs))){
  comb_ks <- event_cdfs %>%
    group_by(var) %>%
    summarize(ks = ks.test(value[cluster == combs[1,comb]], value[cluster == combs[2,comb]])$statistic,
              sig = ks.test(value[cluster == combs[1,comb]], value[cluster == combs[2,comb]])$p.value)
  
  cluster_ks[,comb+1] <- comb_ks[,3]
  colnames(cluster_ks)[comb+1] <- paste0('sig_',combs[1,comb],'_',combs[2,comb])
}

#Filter for significance
diff_vars <- cluster_ks %>%
  pivot_longer(!var, names_to = 'combo', values_to = 'sig') %>%
  filter(sig <= 0.05)

#Number of significant differences per var
diff_vars_obs_all <- diff_vars %>%
  count(var) %>%
  arrange(desc(n))

diff_vars_obs <- tibble(var = unique(cluster_ks$var)) %>% arrange(var)
for(cluster in seq(k)) {
  cluster_char <- as.character(cluster)
  diff_vars_clust <- cluster_ks %>%
    select(var, c(contains(paste0('_',cluster_char,'_')), ends_with(paste0('_',cluster_char)))) %>%
    pivot_longer(-var, names_to = 'combo', values_to = 'p') %>%
    mutate(sig = ifelse(p <= 0.05, 1, 0)) %>%
    group_by(var) %>%
    summarise(n = sum(sig))
  diff_vars_obs[,cluster+1] <- diff_vars_clust$n
}
colnames(diff_vars_obs)[-1] <- as.character(seq(1:8))

pivot_longer(diff_vars_obs, !var, names_to = 'cluster', values_to = 'sig_n') %>%
  ggplot() +
  geom_col(aes(x = var, y = sig_n, group = cluster)) +
  facet_wrap(vars(cluster)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



#Repeat for vars that are in final model
mod_diff_vars <- filter(diff_vars, var %in% sig_preds)

mod_diff_vars_obs_all <- mod_diff_vars %>%
  count(var) %>%
  arrange(desc(n))

mod_diff_vars_obs <- filter(diff_vars_obs, var %in% sig_preds)

pivot_longer(mod_diff_vars_obs, !var, names_to = 'cluster', values_to = 'sig_n') %>%
  ggplot() +
  geom_col(aes(x = var, y = sig_n, group = cluster, fill = var)) +
  facet_wrap(vars(cluster)) +
  scale_fill_manual(values = c(rep('#DB995A', 3), rep('#654236', 3), rep('#8EA4D2', 4), rep('#DA7635', 4)), guide = 'none') +
  scale_x_discrete(labels = c('Delta Q', 'Delta Q/Initial Q', 'Duration', 'Duration Since Last', 'Initial N', 'Initial Q',
                              'Western Precip Total', 'Western Max Intensity', 'Eastern Precip Total', 'Eastern 90d Antecedent',
                              'Clinton Delta', 'Milford % of Outflows', 'Tuttle Delta', 'Reservoir Precip Ratio')) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#box plots of significantly different vars, remove outliers
diff_cutoff <- 4
event_summary_cluster %>%
  select(all_of(diff_vars_obs_all$var[diff_vars_obs_all$n >= diff_cutoff]), cluster) %>%
  pivot_longer(!cluster, names_to = 'var', values_to = 'value') %>%
  group_by(var) %>%
  filter(value <= IQR(value)+quantile(value, 0.75) &
           value >= quantile(value, 0.25)-IQR(value)) %>%
  ggplot() +
  geom_boxplot(aes(y = value, x = cluster, fill = factor(cluster))) +
  scale_fill_discrete(guide = 'none') +
  facet_wrap(vars(var), scales = 'free_y')

#box plots of vars only in final model
event_summary_cluster %>%
  select(all_of(sig_preds[sig_preds != 'season']), cluster) %>%
  pivot_longer(!cluster, names_to = 'var', values_to = 'value') %>%
  group_by(var) %>%
  filter(value <= IQR(value)+quantile(value, 0.75) &
           value >= quantile(value, 0.25)-IQR(value)) %>%
  ggplot() +
  geom_boxplot(aes(y = value, x = cluster, fill = factor(cluster))) +
  scale_fill_discrete(guide = 'none') +
  facet_wrap(vars(var), scales = 'free_y')

#Var difference by season
event_summary_cluster %>%
  select(season, value = delta_rat_Q, cluster) %>%
  filter(value <= IQR(value)+quantile(value, 0.75) &
           value >= quantile(value, 0.25)-IQR(value)) %>%
  ggplot() +
  geom_boxplot(aes(y = value, x = season, fill = factor(cluster))) +
  scale_fill_discrete(guide = 'none') +
  facet_wrap(vars(cluster))


# Collect significant coefficients for each reference level ---------------


# Effect of rescaling variables on model coefficients ---------------------




# Effect of reference cluster ---------------------------------------------

predictions_byref <- list()
for(cluster in unique(event_summary_cluster$cluster)){
  predictions_byref[[cluster]] <- tibble(event = seq(1:length(event_summary_cluster$cluster)))
  for(reference in seq(1:8)){
    mod_frame$cluster <- relevel(factor(mod_frame$cluster, levels = seq(1:8)), ref = reference)
    mod <- multinom(cluster ~ 1 + . + 
                      delta_rat_Q*season + delta_rat_Q*duration + p_10260015_1d_total*reservoir_precip_ratio,
                    data = mod_frame)
    predictions_byref[[cluster]][,as.character(reference)] <- predict(mod, type = 'probs')[,as.character(cluster)]
  }
}

predictions_byref_summary <- tibble(
  cluster = unique(event_summary_cluster$cluster),
  range = NA,
  std = NA
)
for(cluster in unique(event_summary_cluster$cluster)){
  summary <- predictions_byref[[as.character(cluster)]] %>% pivot_longer(!event, names_to = 'reference', values_to = 'prob') %>%
    group_by(event) %>%
    summarise(prob_max = max(prob),
              prob_min = min(prob),
              prob_range = prob_max - prob_min,
              prob_mean = mean(prob),
              prob_std = var(prob)^0.5)
  predictions_byref_summary[cluster,'range'] <- mean(summary$prob_range)
  predictions_byref_summary[cluster, 'std'] <- mean(summary$prob_std)
}



# Visualizing probs across 1 variable -------------------------------------

mod_frame$cluster <- relevel(factor(mod_frame$cluster, levels = seq(1:8)), ref = 8)
mod <- multinom(cluster ~ 1 + . + 
                  delta_rat_Q*season + delta_rat_Q*duration + p_10260015_1d_total*reservoir_precip_ratio,
                data = mod_frame)

#old method - synthesize data using median values of all other variables
multi_logit_viz_single <- function(predictor_select, low_quant = 0.1, high_quant = 0.9, split_season = FALSE){
  
  #Get predictor values as vector
  predictor_vals <- unlist(as.vector(mod_frame[,predictor_select]))
  
  #Set up new data frame and model
  mod_frame_replacement <- mod_frame
  colnames(mod_frame_replacement)[colnames(mod_frame_replacement) == predictor_select] <- 'predictor'
  
  mod_synth <- multinom(cluster ~ 1 + . + 
                          delta_rat_Q*season + delta_rat_Q*duration + p_10260015_1d_total*reservoir_precip_ratio,
                        data = mod_frame_replacement)
  
  #Synthesize data from medians
  var_quant_low <- quantile(predictor_vals, probs = low_quant)
  var_quant_high <- quantile(predictor_vals, probs = high_quant)
  var_quant_steps <- seq(var_quant_low, var_quant_high, by = (var_quant_high - var_quant_low)/100)
  
  predictor_medians <- apply(select(mod_frame_replacement, !c(cluster, predictor, season)), 2, median)
  predictor_medians <- t(as.data.frame(predictor_medians))
  rownames(predictor_medians) <- NULL
  predictor_medians <- predictor_medians[rep(1, each = 101*4),]
  
  synth_data <- tibble(season = rep(unique(mod_frame_replacement$season), each = 101),
                       predictor = rep(var_quant_steps, 4)) %>%
    cbind(., predictor_medians)
  
  synth_predict <- cbind(synth_data, predict(mod_synth, newdata = synth_data, 'probs')) %>%
    select(season, predictor, as.character(seq(8))) %>%
    pivot_longer(!c(season, predictor), names_to = 'cluster', values_to = 'probability') %>%
    group_by(season, predictor, cluster) %>%
    summarize(p = mean(probability)) 
  
  # if(str_detect(predictor_select, '_log') == TRUE){
  #   synth_predict$predictor <- exp(synth_predict$predictor)
  # }
  
  plot <- ggplot(data = synth_predict) +
    geom_line(aes(x = predictor, y = p, color = cluster), linewidth = 1) +
    xlab(predictor_select) +
    facet_wrap(vars(season))
  
  if(split_season == FALSE){
    synth_predict <- synth_predict %>%
      group_by(predictor, cluster) %>%
      summarise(p = mean(p))
      
    plot <- ggplot(data = synth_predict) +
      geom_line(aes(x = predictor, y = p, color = cluster), linewidth = 1) +
      xlab(predictor_select)
  }
  
  return(plot)
}

multi_logit_viz_single('p_10260015_maxint', low_quant = 0.1, high_quant = 0.9, split_season = FALSE)







# Tidy up box plots -------------------------------------------------------
event_summary_cluster$cluster <- relevel(factor(event_summary_cluster$cluster, levels = seq(1:8)), ref = 1)
#Name all vars
all_vars <- c(
            #Event Size
              'Duration' = 'duration', 'Delta Q' = 'delta_Q', 'Delta Q/Initial Q' = 'delta_rat_Q', 'Total Runoff' = 'total_Q', 'Max Q' = 'max_Q',
            #Initial Conditions
              'Duration Since Last' = 'duration_since_last', 'Initial Q' = 'initial_Q', 'Initial N' = 'initial_N', 
            #Total in-event precip
              'Western Total' = 'p_10260015_1d_total', 'Eastern Total' = 'p_10270104_1d_total', 
            #Long term antecedent precip
              'Western 365d Antecedent' = 'p_10260015_365d_ante', 'Eastern 365d Antecedent' = 'p_10270104_365d_ante',
            #Intermediate term antecedent precip, only in east
              'Eastern 180d Antecedent' = 'p_10270104_180d_ante','Eastern 90d Antecedent' = 'p_10270104_90d_ante', 'Eastern 60d Antecedent' = 'p_10270104_60d_ante','Eastern 30d Antecedent' = 'p_10270104_30d_ante',
            #Short term antecedent precip
              'Western 15d Antecedent' = 'p_10260015_15d_ante','Eastern 15d Antecedent' = 'p_10270104_15d_ante', 
            #Average precip Intensities
              'Western Average Intensity' = 'p_10260015_1d_avgint', 'Eastern Average Intensity' = 'p_10270104_1d_avgint',
            #Max precip intensities
              'Western Max Intensity' = 'p_10260015_maxint', 'Eastern Max Intensity' = 'p_10270104_maxint', 
            #Coordinates of max precip
              'Latitude of Max Total' = 'p_max_y', 'Longitude of Max Total' = 'p_max_x',
            #Reservoir outflows
              'Tuttle Total' = 'r_t_1d_total', 'Tuttle Delta' = 'r_t_1d_delta', 'Tuttle Delta/Initial' = 'r_t_1d_delta_rat', 'Tuttle % of All Outflows' = 'r_t_share',
              'Clinton Total' = 'r_c_1d_total', 'Clinton Delta' = 'r_c_1d_delta', 'Clinton Delta/Initial' = 'r_c_1d_delta_rat', 'Clinton % of All Outflows' = 'r_c_share',
              'Milford Total' = 'r_m_1d_total', 'Milford Delta' = 'r_m_1d_delta', 'Milford Delta/Initial' = 'r_m_1d_delta_rat', 'Milford % of All Outflows' = 'r_m_share',
              'Perry Total' = 'r_p_1d_total', 'Perry Delta' = 'r_p_1d_delta', 'Perry Delta/Initial' = 'r_p_1d_delta_rat', 'Perry % of All Outflows' = 'r_p_share',
            #Outflow ratios
              'Reservoir Outflows/Total EKSRB Precip' = 'reservoir_precip_ratio', 'Reservoir Outflows/Total Runoff' = 'reservoir_runoff_ratio'
            )

#Get order from sig differences
order_vars <- factor(all_vars, levels = c(diff_vars_obs_all$var, subset(all_vars, !(all_vars %in% diff_vars_obs_all$var))))
order_names <- vector(length = length(order_vars))
for(varnum in seq(length(order_vars))){
  order_names[varnum] <- names(order_vars[order_vars == levels(order_vars)[varnum]])
}


#box plots of all vars

event_summary_cluster %>%
  select(all_of(all_vars), cluster) %>%
  pivot_longer(!cluster, names_to = 'var', values_to = 'value') %>%
  group_by(var) %>%
  filter(value <= IQR(value)+quantile(value, 0.75) &
           value >= quantile(value, 0.25)-IQR(value)) %>%
  ggplot() +
  geom_boxplot(aes(y = value, x = cluster, fill = factor(cluster))) +
  scale_fill_discrete(guide = 'none') +
  facet_wrap(vars(factor(var, levels = order_names)), scales = 'free_y')

#Only Precip Vars
event_summary_cluster %>%
  select(all_of(all_vars[str_detect(all_vars, '^p_')]), cluster) %>%
  pivot_longer(!cluster, names_to = 'var', values_to = 'value') %>%
  group_by(var) %>%
  filter(value <= IQR(value)+quantile(value, 0.75) &
           value >= quantile(value, 0.25)-IQR(value)) %>%
  ggplot() +
  geom_boxplot(aes(y = value, x = cluster, fill = factor(cluster))) +
  scale_fill_discrete(guide = 'none') +
  xlab('') +
  ylab('Precip [mm]') +
  facet_wrap(vars(factor(var, levels = order_names[order_names %in% names(all_vars[str_detect(all_vars, '^p_')])])), scales = 'free_y')

#Only Reservoir Vars
event_summary_cluster %>%
  select(all_of(all_vars[str_detect(all_vars, '^r_')]), cluster) %>%
  pivot_longer(!cluster, names_to = 'var', values_to = 'value') %>%
  group_by(var) %>%
  filter(value <= IQR(value)+quantile(value, 0.75) &
           value >= quantile(value, 0.25)-IQR(value)) %>%
  ggplot() +
  geom_boxplot(aes(y = value, x = cluster, fill = factor(cluster))) +
  scale_fill_discrete(guide = 'none') +
  xlab('') +
  ylab('Outflow [m3, m3/s]') +
  facet_wrap(vars(factor(var, levels = c(
    'Tuttle Total', 'Tuttle Delta', 'Tuttle Delta/Initial', 'Tuttle % of All Outflows',
    'Clinton Total', 'Clinton Delta', 'Clinton Delta/Initial', 'Clinton % of All Outflows',
    'Milford Total', 'Milford Delta', 'Milford Delta/Initial', 'Milford % of All Outflows',
    'Perry Total', 'Perry Delta', 'Perry Delta/Initial', 'Perry % of All Outflows'
  ))), scales = 'free_y')
  #facet_wrap(vars(factor(var, levels = order_names[order_names %in% names(all_vars[str_detect(all_vars, '^r_')])])), scales = 'free_y')


#Only Event Size Vars
event_summary_cluster %>%
  select(all_of(all_vars[all_vars %in% c('duration', 'delta_Q', 'delta_rat_Q', 'max_Q', 'total_Q')]), cluster) %>%
  pivot_longer(!cluster, names_to = 'var', values_to = 'value') %>%
  group_by(var) %>%
  filter(value <= IQR(value)+quantile(value, 0.75) &
           value >= quantile(value, 0.25)-IQR(value)) %>%
  ggplot() +
  geom_boxplot(aes(y = value, x = cluster, fill = factor(cluster))) +
  scale_fill_discrete(guide = 'none') +
  xlab('') +
  ylab('') +
  facet_wrap(vars(factor(var, levels = order_names[order_names %in% names(all_vars[all_vars %in% c('duration', 'delta_Q', 'delta_rat_Q', 'max_Q', 'total_Q')])])), scales = 'free_y')


#Only Initial Condition Vars
event_summary_cluster %>%
  select(all_of(all_vars[all_vars %in% c('duration_since_last', 'initial_Q', 'initial_N')]), cluster) %>%
  pivot_longer(!cluster, names_to = 'var', values_to = 'value') %>%
  group_by(var) %>%
  filter(value <= IQR(value)+quantile(value, 0.75) &
           value >= quantile(value, 0.25)-IQR(value)) %>%
  ggplot() +
  geom_boxplot(aes(y = value, x = cluster, fill = factor(cluster))) +
  scale_fill_discrete(guide = 'none') +
  xlab('') +
  ylab('') +
  facet_wrap(vars(factor(var, levels = order_names[order_names %in% names(all_vars[all_vars %in% c('duration_since_last', 'initial_Q', 'initial_N')])])), scales = 'free_y')


#Only Ratio vars
event_summary_cluster %>%
  select(all_of(all_vars[all_vars %in% c('reservoir_precip_ratio', 'reservoir_runoff_ratio')]), cluster) %>%
  pivot_longer(!cluster, names_to = 'var', values_to = 'value') %>%
  group_by(var) %>%
  filter(value <= IQR(value)+quantile(value, 0.75) &
           value >= quantile(value, 0.25)-IQR(value)) %>%
  ggplot() +
  geom_boxplot(aes(y = value, x = cluster, fill = factor(cluster))) +
  scale_fill_discrete(guide = 'none') +
  xlab('') +
  ylab('') +
  facet_wrap(vars(factor(var, levels = order_names[order_names %in% names(all_vars[all_vars %in% c('reservoir_precip_ratio', 'reservoir_runoff_ratio')])])), scales = 'free_y')




#box plots of vars only in final model
event_summary_cluster %>%
  select(all_of(all_vars[all_vars %in% sig_preds]), cluster) %>%
  pivot_longer(!cluster, names_to = 'var', values_to = 'value') %>%
  group_by(var) %>%
  filter(value <= IQR(value)+quantile(value, 0.75) &
           value >= quantile(value, 0.25)-IQR(value)) %>%
  ggplot() +
  geom_boxplot(aes(y = value, x = cluster, fill = factor(cluster))) +
  scale_fill_discrete(guide = 'none') +
  xlab('') +
  ylab('') +
  facet_wrap(vars(factor(var, levels = c('Delta Q/Initial Q', 'Delta Q', 'Duration',
                                         'Initial N', 'Initial Q', 'Duration Since Last',
                                         'Eastern Total', 'Western Total', 'Eastern 90d Antecedent',
                                         'Tuttle Delta', 'Clinton Delta', 'Milford % of All Outflows',
                                         'Reservoir Outflows/Total EKSRB Precip', 'Western Max Intensity')
                         )), scales = 'free_y', ncol = 3)
  #facet_wrap(vars(factor(var, levels = order_names[order_names %in% names(all_vars[all_vars %in% sig_preds])])), scales = 'free_y')

library(tidyverse)
library(nnet)
library(reshape)


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
         season = factor(season)) %>%
         #spring = factor(ifelse(season == 'Spring', 1, 0))) %>%
  filter(complete.cases(.)) %>%
  mutate(cluster = cls_id$value) %>%
  #remove variables that are correlated > 0.8
  #select(!contains(correlated_predictors)) %>%
  #remove vars of non-clinical importance
  select(!c(FI, HI))
  #select(!c(FI, HI, wateryear, month, season, yday))


#Set Reference Categories
event_summary_cluster$cluster <- relevel(factor(event_summary_cluster$cluster), ref = 8)
event_summary_cluster$wateryear <- relevel(factor(event_summary_cluster$wateryear), ref = '2017')
event_summary_cluster$season <- relevel(event_summary_cluster$season, ref = 'Spring')


# Random testing stuff ----------------------------------------------------
mod0 <- multinom(cluster ~ 1, data = event_summary_cluster)
mod1 <- multinom(cluster ~ 1 + p_10270104_1d_total, data = event_summary_cluster, model = TRUE)
exp(coef(mod1))
coef(mod1)
summary(mod1)$coefficients
summary(mod1)$standard.errors
abs(summary(mod1)$coefficients/summary(mod1)$standard.errors)
exp(10*coef(mod1))
exp(summary(mod1)$coefficients - 1.96*summary(mod1)$standard.errors)
exp(summary(mod1)$coefficients)
exp(summary(mod1)$coefficients + 1.96*summary(mod1)$standard.errors)

(summary(mod1)$coefficients - 1.96*summary(mod1)$standard.errors)
(summary(mod1)$coefficients)
(summary(mod1)$coefficients + 1.96*summary(mod1)$standard.errors)

2*(1-pnorm(abs(summary(mod1)$coefficients/summary(mod1)$standard.errors)))

#visualize probs and logit
cluster <- 1
test <- tibble(
  var = event_summary_cluster$p_10270104_1d_total,
  p = predict(mod1, type = 'probs')[,cluster+1]/(predict(mod1, type = 'probs')[,cluster+1]+predict(mod1, type = 'probs')[,1]),
  g = log(p/(1-p)))
ggplot(data = test) +
  geom_point(aes(x = var, y = p))
ggplot(data = test) +
  geom_point(aes(x = var, y = g)) +
  geom_abline(slope = coef(mod1)[cluster,2], intercept = coef(mod1)[cluster,1])

local_ps <- tibble(
  p1 = predict(mod1, type = 'probs')[,2]/(predict(mod1, type = 'probs')[,2]+predict(mod1, type = 'probs')[,1]),
  p2 = predict(mod1, type = 'probs')[,3]/(predict(mod1, type = 'probs')[,3]+predict(mod1, type = 'probs')[,1]),
  p3 = predict(mod1, type = 'probs')[,4]/(predict(mod1, type = 'probs')[,4]+predict(mod1, type = 'probs')[,1]),
  p4 = predict(mod1, type = 'probs')[,5]/(predict(mod1, type = 'probs')[,5]+predict(mod1, type = 'probs')[,1]),
  p5 = predict(mod1, type = 'probs')[,6]/(predict(mod1, type = 'probs')[,6]+predict(mod1, type = 'probs')[,1]),
  p6 = predict(mod1, type = 'probs')[,7]/(predict(mod1, type = 'probs')[,7]+predict(mod1, type = 'probs')[,1]),
  p7 = predict(mod1, type = 'probs')[,8]/(predict(mod1, type = 'probs')[,8]+predict(mod1, type = 'probs')[,1]),
  p8 = 0.5 #(1-p1 + 1-p2 + 1-p3 + 1-p4 + 1-p5 + 1-p6 + 1-p7)/7
) %>%
  mutate(g1 = log(p1/(1-p1)),
         g2 = log(p2/(1-p2)),
         g3 = log(p3/(1-p3)),
         g4 = log(p4/(1-p4)),
         g5 = log(p5/(1-p5)),
         g6 = log(p6/(1-p6)),
         g7 = log(p7/(1-p7)),
         g8 = log(p8/(1-p8))) %>%
  mutate(p1_cond = exp(g1)/(exp(g1) + exp(g2) + exp(g3) + exp(g4) + exp(g5) + exp(g6) + exp(g7) + exp(g8)),
         p1_global = predict(mod1, type = 'probs')[,2])

 
# Visualize a single predictor (+ season) model ---------------------------

#Function for visualizing a single predictor (+ season) multinomial logistic regression
multi_logit_viz_single <- function(predictor_select, low_quant = 0.1, high_quant = 0.9){

#Get predictor values as vector
predictor_vals <- unlist(as.vector(event_summary_cluster[,predictor_select]))

#Set up new data frame
model_frame <- tibble(
  cluster = event_summary_cluster$cluster,
  #season = event_summary_cluster$season,
  predictor = predictor_vals
)

model <- multinom(cluster ~ predictor, data = model_frame, model = TRUE)

exp(coef(model))
fitted(model)

var_quant_low <- quantile(predictor_vals, probs = low_quant)
var_quant_high <- quantile(predictor_vals, probs = high_quant)
var_quant_steps <- seq(var_quant_low, var_quant_high, by = (var_quant_high - var_quant_low)/100)

synth_data <- tibble(season = rep(c('Winter', 'Spring', 'Summer', 'Fall'), each = 101), predictor = rep(var_quant_steps, 4))

synth_predict <- cbind(synth_data, predict(model, newdata = synth_data, 'probs')) %>%
  mutate(season = factor(season, levels = c('Winter', 'Spring', 'Summer', 'Fall'))) %>%
  pivot_longer(!c(season, predictor), names_to = 'cluster', values_to = 'probability')

plot <- ggplot(data = synth_predict) +
  geom_line(aes(x = predictor, y = probability, color = cluster), linewidth = 1) +
  #scale_color_manual(values = c('#8DA9C4','#60992D','#C33C54','#F4B860')) +
  xlab(predictor_select) #+
  #facet_wrap(vars(season))
  
  return(plot)
}

multi_logit_viz_single('p_10270104_1d_total', low_quant = 0.1, high_quant = 0.9)


# Model with all predictors together --------------------------------------

mod1 <- multinom(cluster ~ 1, data = event_summary_cluster)
mod2 <- multinom(cluster ~ 1 + ., data = event_summary_cluster)

anova(mod1, mod2)

comp <- tibble(
  actual = event_summary_cluster$cluster,
  predict = predict(mod2)) %>%
  mutate(match = actual == predict)
mean(comp$match)


# Model with only significant predictors ----------------------------------

#find significant predictors
predictor_list <- colnames(select(event_summary_cluster, !cluster))

null_mod <- multinom(cluster ~ 1, data = event_summary_cluster)

single_pred_summary <- tibble(
  predictor = predictor_list,
  p = NA
)

for(pred in 1:length(predictor_list)){
  predictor <- predictor_list[pred]
  f <- paste0('cluster ~ ', predictor)
  pred_mod <- multinom(f, data = event_summary_cluster)
  
  single_pred_summary[pred,2] <- anova(null_mod, pred_mod)[2,7]
}
single_pred_summary <- arrange(single_pred_summary, p)


#Model with all significant predictors together
good_preds <- unlist(as.vector(filter(single_pred_summary, p <= 0.05)[,1]))

event_summary_trim <- event_summary_cluster %>%
  select(all_of(good_preds))
colnames(event_summary_trim) <- good_preds
event_summary_trim <- cbind(event_summary_trim, cluster = event_summary_cluster$cluster)

mod1 <- multinom(cluster ~ 1, data = event_summary_cluster)
mod2 <- multinom(cluster ~ 1 + ., data = event_summary_trim)

anova(mod1, mod2)

mod_perform_sigs <- tibble(event = seq(1:nrow(event_summary_cluster))) %>%
  cbind(., fitted(mod2)) %>%
  pivot_longer(!event, names_to = 'cluster', values_to = 'p') %>%
  group_by(event) %>%
  arrange(desc(p), .by_group = TRUE) %>%
  summarise(p1_clust = cluster[1],
            p1_p = p[1],
            p2_clust = cluster[2],
            p2_p = p[2],
            p3_clust = cluster[3],
            p3_p = p[3],
            p4_clust = cluster[4],
            p4_p = p[4],
            p5_clust = cluster[5],
            p5_p = p[5],
            p6_clust = cluster[6],
            p6_p = p[6],
            p7_clust = cluster[7],
            p7_p = p[7],
            p8_clust = cluster[8],
            p8_p = p[8]) %>%
  mutate(actual = event_summary_cluster$cluster) %>%
  select(!event) %>%
  mutate(
    match1 = p1_clust == actual,
    match2 = match1 | p2_clust == actual,
    match3 = match2 | p3_clust == actual,
    match4 = match3 | p4_clust == actual,
    match5 = match4 | p5_clust == actual,
    match6 = match5 | p6_clust == actual,
    match7 = match6 | p7_clust == actual,
    match8 = match7 | p8_clust == actual) %>%
  select(actual, everything(.))

mod_perform_sigs_clust <- tibble(
  cluster = seq(1:8),
  accuracy = c(mean(mod_perform_sigs$match1),
               mean(mod_perform_sigs$match2),
               mean(mod_perform_sigs$match3),
               mean(mod_perform_sigs$match4),
               mean(mod_perform_sigs$match5),
               mean(mod_perform_sigs$match6),
               mean(mod_perform_sigs$match7),
               mean(mod_perform_sigs$match8)),
  clust_n = c(39, 37, 36, 27, 27, 23, 17, 10)) %>%
  mutate(cumul_n = cumsum(clust_n),
         prop_n = cumul_n/sum(clust_n))
mod_perform_sigs_clust$accuracy - mod_perform_sigs_clust$prop_n
ggplot(data = mod_perform_sigs_clust) +
  geom_line(aes(x = cluster, y = accuracy)) +
  geom_line(aes(x = cluster, y = prop_n), color = 'red') +
  ylim(0,1)
    


# Model with best predictors via forward selection ------------

#Manual selection of predictors
# manual_predictors <- c('p_10270104_1d_total', 'initial_N', 'p_10260008_1d_total', 'initial_Q', 'p_10270104_maxint', 'p_10260008_maxint',
#                        'p_10270104_365d_ante', 'p_unint_vol', 'r_m_1d_total', 'r_t_1d_delta', 'delta_rat_Q', 'p_10260008_1d_avgint', 'p_10270104_1d_avgint',
#                        'r_t_1d_delta_rat', 'duration_since_last', 'max_Q', 'r_c_1d_delta_rat', 'reservoir_precip_ratio', 'month', 'r_sum_1d_total',
#                        'runoff_precip_ratio', 'delta_Q', 'p_10260008_180d_ante', 'p_10270104_180d_ante', 'p_10270104_90d_ante', 'r_c_1d_delta',
#                        'r_m_share', 'r_t_1d_total', 'r_t_share', 'reservoir_runoff_ratio', 'total_Q')

#Run forward selection to find best number of predictors
event_summary_sort <- event_summary_cluster %>%
  select(all_of(single_pred_summary$predictor)) %>%
  #select(all_of(manual_predictors)) %>%
  cbind(., cluster = event_summary_cluster$cluster)

forward_summary <- tibble(
  predictor = single_pred_summary$predictor,
  #predictor = manual_predictors,
  dev = NA,
  p = NA
)



for(pred in 1:length(forward_summary$predictor)){
forward_data1 <- as_data_frame(cbind(cluster = event_summary_sort$cluster, event_summary_sort[,1:(pred-1)]))
forward_data2 <- as_data_frame(cbind(cluster = event_summary_sort$cluster, event_summary_sort[,1:pred]))

mod0 <- multinom(cluster ~ 1, data = forward_data1)
mod1 <- multinom(cluster ~ 1 + ., data = forward_data1)
mod2 <- multinom(cluster ~ 1 + ., data = forward_data2)

forward_summary[pred,2] <- deviance(mod2)
forward_summary[pred,3] <- anova(mod1, mod2)[2,7]


if(pred == 1){
  forward_summary[pred,2] <- anova(mod0,mod2)[2,7]
}
print(paste0(pred,'/',length(forward_summary$predictor),' completed!!'))
}
forward_summary <- forward_summary %>%
  mutate(dev_diff = dev - lag(dev))

ggplot(data = forward_summary) +
  geom_line(aes(x = seq(length(dev)), y = dev))


#manual_predictors <- c('initial_N', 'delta_rat_Q', forward_summary$predictor[forward_summary$dev_diff < 0 & !is.na(forward_summary$dev_diff)])

#Model with bestest predictors together, selected at minimum of deviance
#best_preds <- unlist(as.vector(filter(forward_summary, p <= 0.05)[,1]))
best_preds <- forward_summary$predictor[1:17]
manual_preds <- c(forward_summary$predictor[c(1:13, 15, 17)], 'reservoir_precip_ratio', 'reservoir_runoff_ratio', 'max_Q')

event_summary_trimmer <- event_summary_cluster %>%
  select(all_of(best_preds)) %>%
  cbind(cluster = event_summary_cluster$cluster, .)

event_summary_manual <- event_summary_cluster %>%
  select(all_of(manual_preds)) %>%
  cbind(cluster = event_summary_cluster$cluster, .)

mod1 <- multinom(cluster ~ 1, data = event_summary_cluster)
mod2 <- multinom(cluster ~ ., data = event_summary_trimmer)
mod_old <- multinom(cluster ~ ., data = event_summary_manual)

anova(mod4, mod_old)
anova(mod1, mod2)
anova(mod2, mod3)

mod_old_coef <- exp(summary(mod_old)$coefficients)

z <- summary(mod3)$coefficients/summary(mod3)$standard.errors
sig <- (1 - pnorm(abs(z),0,1))*2

mod_perform_fwd <- tibble(event = seq(1:nrow(event_summary_cluster))) %>%
  cbind(., fitted(mod_old)) %>%
  pivot_longer(!event, names_to = 'cluster', values_to = 'p') %>%
  group_by(event) %>%
  arrange(desc(p), .by_group = TRUE) %>%
  summarise(p1_clust = cluster[1],
            p1_p = p[1],
            p2_clust = cluster[2],
            p2_p = p[2],
            p3_clust = cluster[3],
            p3_p = p[3],
            p4_clust = cluster[4],
            p4_p = p[4],
            p5_clust = cluster[5],
            p5_p = p[5],
            p6_clust = cluster[6],
            p6_p = p[6],
            p7_clust = cluster[7],
            p7_p = p[7],
            p8_clust = cluster[8],
            p8_p = p[8]) %>%
  mutate(actual = event_summary_cluster$cluster) %>%
  select(!event) %>%
  mutate(
    match1 = p1_clust == actual,
    match2 = match1 | p2_clust == actual,
    match3 = match2 | p3_clust == actual,
    match4 = match3 | p4_clust == actual,
    match5 = match4 | p5_clust == actual,
    match6 = match5 | p6_clust == actual,
    match7 = match6 | p7_clust == actual,
    match8 = match7 | p8_clust == actual,
    match_n = ifelse(match1, 1, ifelse(match2, 2, ifelse(match3, 3, ifelse(match4, 4, ifelse(match5, 5, ifelse(match6, 6, ifelse(match7, 7, 8)))))))) %>%
  select(actual, everything(.))

mod_perform_fwd_clust <- tibble(
  cluster = seq(1:8),
  accuracy = c(mean(mod_perform_fwd$match1),
               mean(mod_perform_fwd$match2),
               mean(mod_perform_fwd$match3),
               mean(mod_perform_fwd$match4),
               mean(mod_perform_fwd$match5),
               mean(mod_perform_fwd$match6),
               mean(mod_perform_fwd$match7),
               mean(mod_perform_fwd$match8)),
  pn_acc =  c(mean(mod_perform_fwd$p1_p[mod_perform_fwd$p1_clust == mod_perform_fwd$actual]),
              mean(mod_perform_fwd$p2_p[mod_perform_fwd$p2_clust == mod_perform_fwd$actual]),
              mean(mod_perform_fwd$p3_p[mod_perform_fwd$p3_clust == mod_perform_fwd$actual]),
              mean(mod_perform_fwd$p4_p[mod_perform_fwd$p4_clust == mod_perform_fwd$actual]),
              mean(mod_perform_fwd$p5_p[mod_perform_fwd$p5_clust == mod_perform_fwd$actual]),
              mean(mod_perform_fwd$p6_p[mod_perform_fwd$p6_clust == mod_perform_fwd$actual]),
              mean(mod_perform_fwd$p7_p[mod_perform_fwd$p7_clust == mod_perform_fwd$actual]),
              mean(mod_perform_fwd$p8_p[mod_perform_fwd$p8_clust == mod_perform_fwd$actual])),
  clust_n = c(39, 37, 36, 27, 27, 23, 17, 10)) %>%
  mutate(cumul_n = cumsum(clust_n),
         prop_n = cumul_n/sum(clust_n),
         cumul_pn_acc = cumsum(pn_acc))

mod_perform_fwd_clust$accuracy - mod_perform_fwd_clust$prop_n

ggplot(data = mod_perform_fwd_clust) +
  geom_line(aes(x = cluster, y = accuracy)) +
  geom_line(aes(x = cluster, y = prop_n), color = 'red') +
  ylim(0,1) +
  xlab('# of guesses')

mod_perform_fwd %>%
  select(match_n, ends_with('_p')) %>%
  pivot_longer(!match_n, names_to = 'n', values_to = 'p') %>%
  mutate(match = ifelse(as.integer(str_sub(n, 2, 2)) == as.integer(match_n), TRUE, FALSE)) %>%
  ggplot() +
  geom_boxplot(aes(y = p, fill = match)) +
  facet_wrap(vars(n))


# Visualize best model ----------------------------------------------------

#Create visualization using model with all the 'good' predictors
#One line for each cluster, one facet for each predictor
#First, try running a model with every other var held at its median and compare to what the original viz puts out

multi_logit_viz <- function(predictor_select, low_quant = 0.1, high_quant = 0.9, lump_months = FALSE){

  #Get predictor values as vector
  predictor_vals <- unlist(as.vector(event_summary_cluster[,predictor_select]))
  
  #Set up new data frame
  model_frame <- event_summary_trimmer %>%
    select(!predictor_select) %>%
    mutate(predictor = predictor_vals)
  
  model <- multinom(cluster ~ ., data = model_frame)
  
  exp(coef(model))
  fitted(model)
  
  var_quant_low <- quantile(predictor_vals, probs = low_quant)
  var_quant_high <- quantile(predictor_vals, probs = high_quant)
  var_quant_steps <- seq(var_quant_low, var_quant_high, by = (var_quant_high - var_quant_low)/100)
  
  
  ##Synth factors by generating every combination of them
  factor_combs <- expand.grid(levels(event_summary_trimmer$wateryear), levels(event_summary_trimmer$month))
  
  ##Synth numerics by taking the median
  predictor_medians <- apply(select(model_frame, !c(cluster, predictor, wateryear, month)), 2, median)
  predictor_medians <- t(as.data.frame(predictor_medians))
  rownames(predictor_medians) <- NULL
  predictor_medians <- predictor_medians[rep(1, each = 101*108),]
  
  synth_data <- tibble(wateryear = rep(factor_combs$Var1, each = 101), 
                       month = rep(factor_combs$Var2, each = 101),
                       predictor = rep(var_quant_steps, 108)) %>%
    cbind(., predictor_medians)

  synth_predict <- cbind(synth_data, predict(model, newdata = synth_data, 'probs')) %>%
    select(wateryear, month, predictor, as.character(seq(8))) %>%
    pivot_longer(!c(wateryear, month, predictor), names_to = 'cluster', values_to = 'probability') %>%
    group_by(month, predictor, cluster) %>%
    summarize(p = mean(probability)) %>%
    mutate(season = ifelse(month %in% c('Dec', 'Jan', 'Feb'), 'Winter', ifelse(month %in% c('Mar', 'Apr', 'May'), 'Spring', ifelse(month %in% c('Jun', 'Jul', 'Aug'), 'Summer', 'Fall'))))
  
  if(lump_months == TRUE){
    synth_predict <- synth_predict %>%
      group_by(predictor, cluster, season) %>%
      summarize(p = mean(p))
    
    plot <- ggplot(data = synth_predict) +
      geom_line(aes(x = predictor, y = p, color = cluster), linewidth = 1) +
      xlab(predictor_select) +
      facet_wrap(vars(factor(season, levels = c('Winter', 'Spring', 'Summer', 'Fall'))))
    
    return(plot)
  } else {
  
  plot <- ggplot(data = synth_predict) +
    geom_line(aes(x = predictor, y = p, color = cluster), linewidth = 1) +
    xlab(predictor_select) +
    facet_wrap(vars(factor(month, levels = c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'))))
  
  
  return(plot)
  }
}

multi_logit_viz('duration_since_last', 0, 0.9, lump_months = TRUE)

predictor <- event_summary_cluster$max_Q
ggplot(data = event_summary_cluster, aes(x = predictor)) +
  geom_histogram() +
  geom_vline(xintercept = median(predictor))


nullmod <- multinom(cluster ~ 1, data = event_summary_cluster)
singmod <- multinom(cluster ~ 1 + season, data = event_summary_cluster)

anova(nullmod, singmod)

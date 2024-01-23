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
event_summary_cluster$cluster <- relevel(factor(event_summary_cluster$cluster), ref = 8)
# event_summary_cluster$wateryear <- relevel(factor(event_summary_cluster$wateryear), ref = '2017')
event_summary_cluster$season <- relevel(event_summary_cluster$season, ref = 'Spring')


rm(list = ls()[ls() != 'event_summary_cluster'])



# 1: Find all p < 0.25 ----------------------------------------------------

predictor_list <- colnames(select(event_summary_cluster, !c(cluster, FI, HI)))

null_mod <- multinom(cluster ~ 1, data = select(event_summary_cluster, !c(FI, HI)))

single_pred_summary <- tibble(
  predictor = predictor_list,
  p = NA,
  aic = NA
)

for(pred in 1:length(predictor_list)){
  predictor <- predictor_list[pred]
  f <- paste0('cluster ~ ', predictor)
  pred_mod <- multinom(f, data = event_summary_cluster)
  
  single_pred_summary[pred,'p'] <- anova(null_mod, pred_mod)[2,7]
  single_pred_summary[pred, 'aic'] <- pred_mod$AIC
}
single_pred_summary <- arrange(single_pred_summary, p)


# 2/3: Fit all p < 0.25 together, remove individual p > 0.05 and check against orig model ------------------------------------------------------

good_preds <- unlist(as.vector(filter(single_pred_summary, p <= 0.25)[,1]))

event_summary_trim <- event_summary_cluster %>%
  select(all_of(good_preds))
colnames(event_summary_trim) <- good_preds
event_summary_trim <- cbind(event_summary_trim, cluster = event_summary_cluster$cluster)

mod1 <- multinom(cluster ~ 1 + ., data = event_summary_trim)
anova(null_mod, mod1)
mod1_coef <- summary(mod1)$coefficients
mod1_stde <- summary(mod1)$standard.errors
mod1_wald <- abs(mod1_coef/mod1_stde)
mod1_wald_p <- 2*(1 - pnorm(mod1_wald))

#iteratively remove vars and check if model changed
mod2 <- multinom(cluster ~ 1 + ., data = select(event_summary_trim, !c(r_c_1d_total, r_m_1d_total, r_c_1d_delta_rat, 
                                                                       p_10270104_maxint, p_10260015_365d_ante, p_10270104_180d_ante,
                                                                       p_10270104_1d_avgint, r_p_1d_delta_rat, r_p_1d_delta, p_10270104_365d_ante,
                                                                       p_10270104_15d_ante, p_10270104_30d_ante, r_p_share, p_max_y, r_t_share)))
anova(mod1, mod2)

mod2_coef <- summary(mod2)$coefficients
mod2_stde <- summary(mod2)$standard.errors
mod2_lowce <- mod2_coef - 1.96*mod2_stde
mod2_hice <- mod2_coef + 1.96*mod2_stde
mod2_wald <- abs(mod2_coef/mod2_stde)
mod2_wald_p <- 2*(1 - pnorm(mod2_wald))

mod3 <- multinom(cluster ~ 1 + ., data = select(event_summary_trim, !c(r_c_1d_total, r_m_1d_total, r_c_1d_delta_rat,
                                                                       p_10270104_maxint, p_10260015_365d_ante, p_10270104_180d_ante,
                                                                       p_10270104_1d_avgint, r_p_1d_delta_rat, r_p_1d_delta, p_10270104_365d_ante,
                                                                       p_10270104_15d_ante, p_10270104_30d_ante, r_p_share, p_max_y, r_t_share, spring)))
anova(mod2, mod3)
anova(null_mod, mod2)
anova(null_mod, mod3)

mod3_coef <- exp(summary(mod3)$coefficients)
mod3_stde <- summary(mod3)$standard.errors
mod3_wald <- abs(log(mod3_coef)/mod3_stde)
mod3_wald_p <- 2*(1 - pnorm(mod3_wald))


# 4: add originally insignificant vars one by one -----------------------
step2_vars <- colnames(mod3_coef)[-1]
step2_vars <- step2_vars[!str_detect(step2_vars, 'season')]
step2_vars <- c(step2_vars, 'season')

#r_t_1d_delta and duration are significant improvements. Also, swap max_Q for delta_Q
add_var <- select(event_summary_cluster, cluster, step2_vars[step2_vars!='max_Q'], duration, delta_Q, r_t_1d_delta)
mod4 <- multinom(cluster ~ 1 + ., data = add_var)

anova(mod3, mod4)
anova(null_mod, mod4)

mod4_coef <- exp(summary(mod4)$coefficients)
mod4_stde <- summary(mod4)$standard.errors
mod4_wald <- abs(log(mod4_coef)/mod4_stde)
mod4_wald_p <- 2*(1 - pnorm(mod4_wald))



# 5: Check continuous variables for linearity --------------------------------
step4_vars <- colnames(mod4_coef)[-1]
step4_vars <- step4_vars[!str_detect(step4_vars, 'season')]
step4_vars <- c(step4_vars, 'season')

transform_vars <- step4_vars[!str_detect(step4_vars, 'season')]
transform_devs <- tibble(
  var = transform_vars,
  expm2 = NA,
  expm1 = NA,
  expmhalf = NA,
  exphalf = NA,
  exp1 = NA,
  exp2 = NA,
  exp3 = NA,
  loge = NA
)

for(var in seq(1:length(transform_vars))){
  transform_var_data <- event_summary_cluster %>%
    select(transform_vars[var])
  
  transform_var_data[transform_var_data  <= 0] <- 0.01
  
  #Expm2
  transform_frame <- select(event_summary_cluster, cluster, step4_vars[!str_detect(step4_vars, transform_vars[var])]) %>%
    cbind(transform_var_data^(-2))
  transform_devs[var,2] <- deviance(multinom(cluster ~ 1 + ., data = transform_frame))
  
  #Expm1
  transform_frame <- select(event_summary_cluster, cluster, step4_vars[!str_detect(step4_vars, transform_vars[var])]) %>%
    cbind(transform_var_data^(-1))
  transform_devs[var,3] <- deviance(multinom(cluster ~ 1 + ., data = transform_frame))
  
  #Expm2
  transform_frame <- select(event_summary_cluster, cluster, step4_vars[!str_detect(step4_vars, transform_vars[var])]) %>%
    cbind(transform_var_data^(-0.5))
  transform_devs[var,4] <- deviance(multinom(cluster ~ 1 + ., data = transform_frame))
  
  #Expm2
  transform_frame <- select(event_summary_cluster, cluster, step4_vars[!str_detect(step4_vars, transform_vars[var])]) %>%
    cbind(transform_var_data^(0.5))
  transform_devs[var,5] <- deviance(multinom(cluster ~ 1 + ., data = transform_frame))
  
  #Expm2
  transform_frame <- select(event_summary_cluster, cluster, step4_vars[!str_detect(step4_vars, transform_vars[var])]) %>%
    cbind(transform_var_data^(1))
  transform_devs[var,6] <- deviance(multinom(cluster ~ 1 + ., data = transform_frame))
  
  #Expm2
  transform_frame <- select(event_summary_cluster, cluster, step4_vars[!str_detect(step4_vars, transform_vars[var])]) %>%
    cbind(transform_var_data^(2))
  transform_devs[var,7] <- deviance(multinom(cluster ~ 1 + ., data = transform_frame))
  
  #Expm2
  transform_frame <- select(event_summary_cluster, cluster, step4_vars[!str_detect(step4_vars, transform_vars[var])]) %>%
    cbind(transform_var_data^(3))
  transform_devs[var,8] <- deviance(multinom(cluster ~ 1 + ., data = transform_frame))
  
  #loge
  transform_frame <- select(event_summary_cluster, cluster, step4_vars[!str_detect(step4_vars, transform_vars[var])]) %>%
    cbind(log(transform_var_data))
  transform_devs[var,9] <- deviance(multinom(cluster ~ 1 + ., data = transform_frame))
}


transformed <- select(event_summary_cluster, cluster, step4_vars) %>%
  mutate(initial_Q_log = log(initial_Q),
         r_t_1d_delta_log = log(ifelse(r_t_1d_delta <= 0, 0.01, r_t_1d_delta))) %>%
  select(!c(initial_Q, r_t_1d_delta))

mod5 <- multinom(cluster ~ 1 + ., data = transformed)

mod5_coef <- exp(summary(mod5)$coefficients)
mod5_stde <- summary(mod5)$standard.errors
mod5_wald <- abs(log(mod5_coef)/mod5_stde)
mod5_wald_p <- 2*(1 - pnorm(mod5_wald))

#comp <- (coef(mod6) - coef(mod5)[,colnames(coef(mod5)) != 'reservoir_precip_ratio'])/coef(mod5)[,colnames(coef(mod5)) != 'reservoir_precip_ratio']


# 6: Check for interactions -----------------------------------------------
step5_vars <- colnames(transformed)[colnames(transformed) != 'cluster']

var_combs <- combn(step5_vars, 2)

interaction_summary <- tibble(
  interaction = as.character(1:ncol(var_combs)),
  dev = NA)

for(comb in 1:ncol(var_combs)){
  c <- paste0(var_combs[1,comb], '*', var_combs[2,comb])
  f <- paste0('cluster ~ 1 + . + ', var_combs[1,comb], '*', var_combs[2,comb])
  interaction_summary[comb,1] <- c
  interaction_summary[comb,2] <- deviance(multinom(f, data = transformed))
}

interaction_summary <- interaction_summary %>%
  mutate(dev_diff = deviance(mod5) - dev,
         p = 1-pchisq(dev_diff, 7))

sig_interactions <- interaction_summary$interaction[interaction_summary$p < 0.05]


mod6 <- multinom(cluster ~ 1 + . + 
                   delta_rat_Q*season + delta_rat_Q*duration +
                   p_10260015_1d_total*reservoir_precip_ratio +
                   reservoir_precip_ratio*p_10270104_90d_ante, data = transformed)

mod6_coef <- exp(summary(mod6)$coefficients)
mod6_stde <- summary(mod6)$standard.errors
mod6_wald <- abs(log(mod6_coef)/mod6_stde)
mod6_wald_p <- 2*(1 - pnorm(mod6_wald))

# Remove interactions one at a time, see if change is insig
# interactions removed: p_10270104_90d_ante*initial_Q_log, reservoir_precip_ratio*p_10260015_maxint,
# delta_rat_Q*p_10270104_90d_ante, r_m_share*p_10270104_1d_total, reservoir_precip_ratio*p_10270104_90d_ante
event_summary_cluster$cluster <- relevel(factor(event_summary_cluster$cluster), ref = 1)
transformed$cluster <- factor(transformed$cluster, levels = as.character(c(8,1,2,3,4,5,6,7)))
mod7 <- multinom(cluster ~ 1 + . + 
                   delta_rat_Q*season + delta_rat_Q*duration +
                   p_10260015_1d_total*reservoir_precip_ratio, data = transformed)

anova(null_mod, mod7)

mod7_coef <- exp(summary(mod7)$coefficients)
mod7_stde <- summary(mod7)$standard.errors
mod7_wald <- abs(log(mod7_coef)/mod7_stde)
mod7_wald_p <- 2*(1 - pnorm(mod7_wald))

#Significant coefficients
coef_sig_pos_or_neg <- ifelse(mod7_coef > 1, 1, -1) * ifelse(abs(mod7_wald_p) > 0.05, 0, 1)

sum(coef_sig_pos_or_neg^2)

mod8 <- multinom(cluster ~ 1 + . + 
                   delta_rat_Q*season + delta_rat_Q*duration +
                   p_10260015_1d_total*reservoir_precip_ratio, data = select(transformed, !duration_since_last))
anova(mod7, mod8)

mod8_coef <- exp(summary(mod8)$coefficients)
mod8_stde <- summary(mod8)$standard.errors
mod8_wald <- abs(log(mod8_coef)/mod8_stde)
mod8_wald_p <- 2*(1 - pnorm(mod8_wald))

coef_sig_pos_or_neg_mod8 <- ifelse(mod8_coef > 1, 1, -1) * ifelse(abs(mod8_wald_p) > 0.05, 0, 1)

sum(coef_sig_pos_or_neg_mod8^2)

# 7: Assess goodness of fit -----------------------------------------------

#Hosmer-Lemeshow Test
library(generalhoslem)

mod7_fit <- logitgof(event_summary_cluster$cluster, fitted(mod7), g = 10)
(mod7_fit$observed - mod7_fit$expected)/(mod7_fit$expected^(0.5))

unloadNamespace('generalhoslem')
unloadNamespace('MASS')


#Basic Accuracy
mod_perform_sigs <- tibble(event = seq(1:nrow(event_summary_cluster))) %>%
  cbind(., fitted(mod7)) %>%
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
  est_acc = c(mean(mod_perform_sigs$p1_p),
              mean(mod_perform_sigs$p2_p),
              mean(mod_perform_sigs$p3_p),
              mean(mod_perform_sigs$p4_p),
              mean(mod_perform_sigs$p5_p),
              mean(mod_perform_sigs$p6_p),
              mean(mod_perform_sigs$p7_p),
              mean(mod_perform_sigs$p8_p)),
  clust_n = c(39, 37, 36, 27, 27, 23, 17, 10)) %>%
  mutate(cumul_n = cumsum(clust_n),
         prop_n = cumul_n/sum(clust_n),
         cumul_est_acc = cumsum(est_acc))
mod_perform_sigs_clust$accuracy - mod_perform_sigs_clust$prop_n

ggplot(data = mod_perform_sigs_clust) +
  geom_line(aes(x = cluster, y = accuracy)) +
  geom_line(aes(x = cluster, y = cumul_est_acc), color = 'blue') +
  geom_line(aes(x = cluster, y = prop_n), color = 'red') +
  ylim(0,1) +
  xlab('# of guesses')


mod_perform_sigs %>%
  select(match_n, ends_with('_p')) %>%
  pivot_longer(!match_n, names_to = 'n', values_to = 'p') %>%
  mutate(match = ifelse(as.integer(str_sub(n, 2, 2)) == as.integer(match_n), TRUE, FALSE)) %>%
  ggplot() +
  geom_boxplot(aes(y = p, fill = match)) +
  facet_wrap(vars(n))


# Visualize/Summarize -----------------------------------------------------

# Relevel cluster assignments for visualization
event_summary_cluster$cluster <- factor(event_summary_cluster$cluster, levels = as.character(seq(1:8)))
transformed$cluster <- factor(transformed$cluster, levels = as.character(seq(1:8)))


# Scatter showing hits and misses
scatter <- select(event_summary_cluster, FI, HI, cluster) %>%
  cbind(match1 = mod_perform_sigs$match1, match2 = mod_perform_sigs$match2, match3 = mod_perform_sigs$match3)

ggplot(data = scatter) +
  geom_point(aes(x = FI, y = HI, color = cluster, shape = match2), size = 2) +
  scale_shape_manual(values = c(13, 16), guide = 'none') +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  xlim(c(-1,1)) +
  ylim(c(-1,1))

group_by(scatter, cluster) %>%
  summarise(acc1 = mean(match1))
  


# Box plots
#All variables in model
select(transformed, !season) %>%
  pivot_longer(!cluster, names_to = 'var', values_to = 'value') %>%
  group_by(var) %>%
  filter(value <= IQR(value)+quantile(value, 0.75) &
           value >= quantile(value, 0.25)-IQR(value)) %>%
  ggplot() +
  geom_boxplot(aes(y = value, x = cluster, fill = factor(cluster))) +
  scale_fill_discrete(guide = 'none') +
  facet_wrap(vars(var), scales = 'free_y')

#All variables, grouped into clusters
rescale <- transformed %>%
  mutate(id = seq(1:nrow(transformed))) %>%
  select(!c(initial_Q_log, r_t_1d_delta_log, season, cluster)) %>%
  cbind(initial_Q = event_summary_cluster$initial_Q, r_t_1d_delta = event_summary_cluster$r_t_1d_delta) %>%
  pivot_longer(!id, names_to = 'var', values_to = 'value') %>%
  group_by(var) %>%
  filter(value <= IQR(value)+quantile(value, 0.75) &
           value >= quantile(value, 0.25)-IQR(value)) %>%
  pivot_wider(id_cols = 'id', names_from = 'var', values_from = 'value') %>%
  select(!id) %>%
  rescaler(., type = 'range') %>%
  cbind(cluster = transformed$cluster)

rescale %>%
  pivot_longer(!cluster, names_to = 'var', values_to = 'value') %>%
  group_by(var) %>%
  ggplot() +
  geom_boxplot(aes(y = value, x = var, fill = var)) +
  ylab('Characteristic Rescaled to 0-1') +
  facet_wrap(vars(cluster), scales = 'free_y') +
  scale_fill_manual(values = c(rep('#BCE784', 3), rep('#5DD39E', 3), rep('#348AA7', 4), rep('#525174', 4)), guide = 'none') +
  scale_x_discrete(labels = c('Delta Q', 'Delta Q/Initial Q', 'Duration', 'Duration Since Last', 'Initial N', 'Initial Q',
                              'Western Precip Total', 'Western Max Intensity', 'Eastern Precip Total', 'Eastern 90d Antecedent',
                              'Clinton Delta', 'Milford % of Outflows', 'Tuttle Delta', 'Reservoir Precip Ratio')) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

# Odds ratios of interaction terms

# Linear model ------------------------------------------------------------

lin_data <- select(event_summary_cluster, all_of(sig_preds), FI, HI) 
lin_fi <- lm(FI ~ 1 + . + 
               delta_rat_Q*season + delta_rat_Q*duration +
               p_10260015_1d_total*reservoir_precip_ratio, data = select(lin_data, !HI))
lin_fi_perf <- tibble(
  obs = event_summary_cluster$FI,
  est = predict(lin_fi)
)
ggplot(data = lin_fi_perf) +
  geom_point(aes(x = obs, y = est)) +
  xlim(c(-1,1)) +
  ylim(c(-1,1))
cor(lin_fi_perf)


lin_hi <- lm(HI ~ 1 + . + 
               delta_rat_Q*season + delta_rat_Q*duration +
               p_10260015_1d_total*reservoir_precip_ratio, data = select(lin_data, !FI))
lin_hi_perf <- tibble(
  obs = event_summary_cluster$HI,
  est = predict(lin_hi)
)
ggplot(data = lin_hi_perf) +
  geom_point(aes(x = obs, y = est)) +
  xlim(c(-1,1)) +
  ylim(c(-1,1))
cor(lin_hi_perf)

lin_comb <- tibble(
  fi_est = predict(lin_fi),
  hi_est = predict(lin_hi),
  fi_obs = event_summary_cluster$FI,
  hi_obs = event_summary_cluster$HI
) %>%
  mutate(dist = ((fi_est - fi_obs)^2 + (hi_est - hi_obs)^2)^0.5)
mean(lin_comb$dist)

ggplot(data = lin_comb) +
  geom_point(aes(x = fi_est, y = hi_est), color = 'blue') +
  geom_point(aes(x = fi_obs, y = hi_obs)) +
  xlim(c(-1,1)) +
  ylim(c(-1,1))






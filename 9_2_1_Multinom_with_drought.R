library(tidyverse)
library(nnet)

# Load data ---------------------------------------------------------------


# * Basic characteristics

CQ_summary <- read_csv('./DataFiles/hydro_data/event_char/CQ_summary.csv') %>%
  select(FI, HI = HI_mean)

simple_temporal <- read_csv('./DataFiles/hydro_data/event_char/simple_temporal.csv')
simple_size <- read_csv('./DataFiles/hydro_data/event_char/simple_size.csv')

N_yield <- read_csv('./DataFiles/hydro_data/event_char/N_yield.csv') %>%
  cbind(., CQ_summary)

# * Gridded Precip Predictors
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


# * Outflows Predictors 
outflows_inevent_totals <- read_csv('./DataFiles/hydro_data/event_char/outflows_inevent_totals.csv') %>%
  select(contains('1d'))

outflows_deltas <- read_csv('./DataFiles/hydro_data/event_char/outflows_deltas.csv') %>%
  select(contains('1d'))

flow_ratios <- read_csv('./DataFiles/hydro_data/event_char/flow_ratios.csv')


# * Drought Predictors 
spei1 <- read_csv('./DataFiles/hydro_data/event_char/spei1.csv') %>%
  select(contains(c('10260015', '10270104')))

spei3 <- read_csv('./DataFiles/hydro_data/event_char/spei3.csv') %>%
  select(contains(c('10260015', '10270104')))

spei6 <- read_csv('./DataFiles/hydro_data/event_char/spei6.csv') %>%
  select(contains(c('10260015', '10270104')))

spei12 <- read_csv('./DataFiles/hydro_data/event_char/spei12.csv') %>%
  select(contains(c('10260015', '10270104')))

# * Correlated Predictors 

correlated_predictors <- unlist(as.vector(read_csv('./DataFiles/hydro_data/event_char/correlated_predictors.csv'))) 

# * Combine 

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
  spei12
) %>%
  mutate(wateryear = factor(wateryear),
         month = factor(month),
         season = factor(season)) %>%
  select(!contains(correlated_predictors))

# * create clusters 

hclust_data <- event_summary %>%
  filter(complete.cases(.)) %>%
  select(FI, HI)

dist <- daisy(hclust_data, metric = 'euclidean')

cls <- agnes(dist, method = "ward")

k <- 8

cls_id <- cutree(cls, k = k)

# * Add cluster ID back in
event_summary_cluster <- event_summary %>%
  filter(complete.cases(.)) %>%
  mutate(cluster = cls_id)

rm(list = ls()[ls() != 'event_summary_cluster'])

#Set Reference categories
event_summary_cluster$cluster <- relevel(factor(event_summary_cluster$cluster), ref = 8)
# event_summary_cluster$wateryear <- relevel(factor(event_summary_cluster$wateryear), ref = '2017')
event_summary_cluster$season <- relevel(event_summary_cluster$season, ref = 'Spring')



#Remove irrelevants
event_summary_cluster <- event_summary_cluster %>%
  select(!c(yday, month, wateryear))

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


mod_manual_prev <- multinom(cluster ~ 1 + season + d_10260015_spei12 + p_10270104_1d_total +
                         duration_since_last + delta_rat_Q + r_m_share + initial_N +
                         r_t_share + delta_Q+ delta_rat_Q*season + r_c_1d_delta +
                         log10(initial_Q) + season*initial_Q,
                         data = event_summary_cluster)
mod_manual <- multinom(cluster ~ 1 + season + d_10260015_spei12 + p_10270104_1d_total +
                       duration_since_last + delta_rat_Q + r_m_share + initial_N +
                       r_t_share + delta_Q + delta_rat_Q*season + r_c_1d_delta +
                       log10(initial_Q) + season*initial_Q + season*initial_N,
                       data = event_summary_cluster)

anova(mod_manual, mod_manual_prev)
anova(mod_manual, mod)


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
                                                                       p_10270104_15d_ante, p_10270104_30d_ante, r_p_share, p_max_y, r_t_share,
                                                                       d_10260015_spei3, d_10260015_spei6, d_10270104_spei6, d_10270104_spei12,
                                                                       p_10260015_maxint, p_10260015_1d_total, d_10260015_spei1)))
anova(mod1, mod2)
anova(mod2, null_mod)

mod2_coef <- summary(mod2)$coefficients
mod2_stde <- summary(mod2)$standard.errors
mod2_lowce <- mod2_coef - 1.96*mod2_stde
mod2_hice <- mod2_coef + 1.96*mod2_stde
mod2_wald <- abs(mod2_coef/mod2_stde)
mod2_wald_p <- 2*(1 - pnorm(mod2_wald))

mod3 <- multinom(cluster ~ 1 + ., data = select(event_summary_trim, !c(r_c_1d_total, r_m_1d_total, r_c_1d_delta_rat,
                                                                       p_10270104_maxint, p_10260015_365d_ante, p_10270104_180d_ante,
                                                                       p_10270104_1d_avgint, r_p_1d_delta_rat, r_p_1d_delta, p_10270104_365d_ante,
                                                                       p_10270104_15d_ante, p_10270104_30d_ante, r_p_share, p_max_y, r_t_share,
                                                                       d_10260015_spei3, d_10260015_spei6, d_10270104_spei6, d_10270104_spei12,
                                                                       p_10260015_maxint, p_10260015_1d_total, d_10260015_spei1)))
anova(mod2, mod3)
anova(null_mod, mod2)
anova(null_mod, mod3)

mod3_coef <- exp(summary(mod3)$coefficients)
mod3_stde <- summary(mod3)$standard.errors
mod3_wald <- abs(log(mod3_coef)/mod3_stde)
mod3_wald_p <- 2*(1 - pnorm(mod3_wald))

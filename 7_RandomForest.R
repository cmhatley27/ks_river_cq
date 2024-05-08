library(dplyr)
library(tidyverse)
library(tidymodels)
library(partykit)
library(bonsai)
library(ranger)
library(pROC)

source('load_and_combine.R')

# Load data ---------------------------------------------------------------

event_summary <- with_tz(read_csv('./DataFiles/hydro_data/event_summary.csv'), 'America/Chicago') %>%
  mutate(wateryear = factor(wateryear),
         month = factor(month),
         season = factor(season))

event_summary_trim <- with_tz(read_csv('./DataFiles/hydro_data/event_summary_trim.csv'), 'America/Chicago') %>%
  mutate(wateryear = factor(wateryear),
         month = factor(month),
         season = factor(season))

new_test <- cbind(CQ_summary, event_precip_huc8) %>%
  select(-event, -event_number, -HI_sd, -HI_cv)


# Conditional importances for FI ---------------------------------------------------------------

#create data
FI_forest_data <- event_summary_trim %>%
  select(!HI) %>%
  filter(complete.cases(.)) 

#run for multiple seeds and check top variables
nseeds <- 10
npred <- 43
top_predictors_FI <- data.frame(matrix(NA, nrow = nseeds, ncol = npred + 2))
names(top_predictors_FI) <- c('seed', 'oob_mse', seq(npred))

for(seed in seq(nseeds)){

skip_to_next <- FALSE
  
set.seed(9)
data_split <- initial_split(FI_forest_data, prop = 4/5)
train_data <- training(data_split)
test_data <- testing(data_split)

FI_forest_fit <- partykit::cforest(
  FI ~ .,
  data = train_data,
  ntree = 1000,
  control = ctree_control(mincriterion = 0.95))

#extract variable importances
tryCatch(
FI_vi <- partykit::varimp(FI_forest_fit, 
                       conditional = T),
error = function(e){skip_to_next <<- TRUE})

if(skip_to_next) { next }

#collect into data frame
FI_forest_results <- tibble::tibble(predictor = names(FI_vi),
                               conditional_importance = FI_vi) %>%
  arrange(desc(conditional_importance)) %>%
  mutate(predictor = factor(predictor, levels = predictor))

FI_regression_predictions <- data.frame(
  predicted_oob = predict(FI_forest_fit, OOB = TRUE, type = 'response'),
  actual = train_data$FI
) %>%
  mutate(diff_oob = actual - predicted_oob)
oob_mse <- mean(FI_regression_predictions$diff_oob^2)

top_predictors_FI[seed,1] <- seed
top_predictors_FI[seed,2] <- oob_mse
top_predictors_FI[seed,2:npred+1] <- FI_forest_results$predictor[1:npred]
print(paste0(seed, '/', nseeds, ' finished'))

}


ggplot(data = FI_forest_results) +
  geom_col(aes(x = predictor, y = conditional_importance)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle('Variable Importance')


# FI Regression w/ forward selection -----------------------------------------------------------
FI_regression_data <- FI_forest_data %>%
  select(FI, c(FI_forest_results$predictor))

nseeds <- 50
npred <- ncol(FI_regression_data) - 1

FI_regression_npred_results <- data.frame(
  var = names(FI_regression_data[2:ncol(FI_regression_data)]),
  train_mse = vector(length = npred),
  test_mse = vector(length = npred),
  test_mse_stdev = vector(length = npred),
  test_mse_cv = vector(length = npred))

for(predictors in seq(npred)){
  
  FI_regression_stability <- data.frame(
    seed = vector(length = nseeds),
    train_mse = vector(length = nseeds),
    test_mse = vector(length = nseeds))
  
  FI_data_select <- FI_regression_data %>%
    select(1:(14+1))
  
  for(seed in seq(nseeds)){
    set.seed(24)
    data_split <- initial_split(FI_data_select, prop = 4/5)
    train_data <- training(data_split)
    test_data <- testing(data_split)
    
    FI_engine_ranger <-
      rand_forest(trees = 1000) %>%
      set_engine('ranger') %>%
      set_mode('regression')
    
    FI_recipe_ranger <-
      recipe(FI ~ ., data = train_data) %>%
      step_normalize(all_numeric_predictors())
    
    FI_workflow_ranger <-
      workflow() %>%
      add_model(FI_engine_ranger) %>%
      add_recipe(FI_recipe_ranger)
    
    FI_fit_ranger <-
      FI_workflow_ranger %>%
      parsnip::fit(data = train_data)
    
    FI_train_results <- tibble(
      train_prediction = FI_fit_ranger$fit$fit$fit$predictions,
      train_actual = FI_fit_ranger$pre$mold$outcomes$FI,
      error2 = (train_prediction - train_actual)^2)
    
    FI_test_results <-
      tibble::tibble(
        test_prediction = predict(FI_fit_ranger, new_data = test_data)$.pred,
        test_actual = test_data$FI) %>%
      mutate(test_error2 = (test_prediction - test_actual)^2)
    FI_regression_stability$seed[seed] <- seed
    FI_regression_stability$train_mse[seed] <- FI_fit_ranger$fit$fit$fit$prediction.error
    FI_regression_stability$test_mse[seed] <- mean(FI_test_results$test_error2)

    print(paste0(seed,'/',nseeds,' seeds finished'))
  }
  
  FI_regression_npred_results$train_mse[predictors] <- mean(FI_regression_stability$train_mse)
  FI_regression_npred_results$test_mse[predictors] <- mean(FI_regression_stability$test_mse)
  FI_regression_npred_results$test_mse_stdev[predictors] <- sqrt(var(FI_regression_stability$test_mse))
  FI_regression_npred_results$test_mse_cv[predictors] <- sqrt(var(FI_regression_stability$test_mse)) / mean(FI_regression_stability$test_mse)

  print(paste0(predictors,'/',npred,' predictors finished'))
}

cor(FI_test_results$test_actual, FI_test_results$test_prediction)^2

ggplot(data = FI_test_results, aes(x = test_actual, y = test_prediction)) +
  geom_point() +
  xlim(c(-1,1)) +
  ylim(c(-1,1)) +
  geom_smooth(method = 'lm', se = FALSE, linetype = 'dashed', color = 'black') +
  geom_abline(slope = 1, color = 'red') +
  #geom_point(data = FI_train_results, aes(x = train_actual, y = train_prediction), color = 'blue', shape = 1) +
  theme_classic()

# Conditional Importances for HI ------------------------------------------

  #create data
  HI_forest_data <- event_summary_trim %>%
    select(!FI) %>%
    filter(complete.cases(.)) 
  
  #run for multiple seeds and check top variables
  nseeds <- 1
  npred <- 100
  top_predictors_HI <- data.frame(matrix(NA, nrow = nseeds, ncol = npred + 2))
  names(top_predictors_HI) <- c('seed', 'oob_mse', seq(npred))
  
  for(seed in seq(nseeds)){
    
    skip_to_next <- FALSE
    
    set.seed(1)
    data_split <- initial_split(HI_forest_data, prop = 4/5)
    train_data <- training(data_split)
    test_data <- testing(data_split)
    
    HI_forest_fit <- partykit::cforest(
      HI ~ .,
      data = train_data,
      ntree = 1000,
      control = ctree_control(mincriterion = 0.95))
    
    #extract variable importances
    tryCatch(
      HI_vi <- partykit::varimp(HI_forest_fit, 
                                conditional = T),
      error = function(e){skip_to_next <<- TRUE})
    
    if(skip_to_next) { next }
    
    #collect into data frame
    HI_forest_results <- tibble::tibble(predictor = names(HI_vi),
                                        conditional_importance = HI_vi) %>%
      arrange(desc(conditional_importance)) %>%
      mutate(predictor = factor(predictor, levels = predictor))
    
    HI_regression_predictions <- data.frame(
      predicted_oob = predict(HI_forest_fit, OOB = TRUE, type = 'response'),
      actual = train_data$HI
    ) %>%
      mutate(diff_oob = actual - predicted_oob)
    oob_mse <- mean(HI_regression_predictions$diff_oob^2)
    
    top_predictors_HI[seed,1] <- seed
    top_predictors_HI[seed,2] <- oob_mse
    top_predictors_HI[seed,2:npred+1] <- HI_forest_results$predictor[1:npred]
    print(paste0(seed, '/', nseeds, ' finished'))
    
  }
  
  
  ggplot(data = HI_forest_results) +
    geom_col(aes(x = predictor, y = conditional_importance)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle('Variable Importance')


  # HI Regression w/ forward selection -----------------------------------------------------------
  HI_regression_data <- HI_forest_data %>%
    select(HI, c(HI_forest_results$predictor))
  
  nseeds <- 50
  npred <- ncol(HI_regression_data) - 1
  
  HI_regression_npred_results <- data.frame(
    var = names(HI_regression_data[2:ncol(HI_regression_data)]),
    train_mse = vector(length = npred),
    test_mse = vector(length = npred),
    test_mse_stdev = vector(length = npred),
    test_mse_cv = vector(length = npred))
  
  for(predictors in seq(npred)){
    
    HI_regression_stability <- data.frame(
      seed = vector(length = nseeds),
      train_mse = vector(length = nseeds),
      test_mse = vector(length = nseeds))
    
    HI_data_select <- HI_regression_data %>%
      select(1:(predictors+1))
    
    for(seed in seq(nseeds)){
      set.seed(seed)
      data_split <- initial_split(HI_data_select, prop = 4/5)
      train_data <- training(data_split)
      test_data <- testing(data_split)
      
      HI_engine_ranger <-
        rand_forest(trees = 1000) %>%
        set_engine('ranger') %>%
        set_mode('regression')
      
      HI_recipe_ranger <-
        recipe(HI ~ ., data = train_data) %>%
        step_normalize(all_numeric_predictors())
      
      HI_workflow_ranger <-
        workflow() %>%
        add_model(HI_engine_ranger) %>%
        add_recipe(HI_recipe_ranger)
      
      HI_fit_ranger <-
        HI_workflow_ranger %>%
        parsnip::fit(data = train_data)
      
      HI_train_results <- tibble(
        train_prediction = HI_fit_ranger$fit$fit$fit$predictions,
        train_actual = HI_fit_ranger$pre$mold$outcomes$HI,
        error2 = (train_prediction - train_actual)^2)
      
      HI_test_results <-
        tibble::tibble(
          test_prediction = predict(HI_fit_ranger, new_data = test_data)$.pred,
          test_actual = test_data$HI) %>%
        mutate(test_error2 = (test_prediction - test_actual)^2)
      HI_regression_stability$seed[seed] <- seed
      HI_regression_stability$train_mse[seed] <- HI_fit_ranger$fit$fit$fit$prediction.error
      HI_regression_stability$test_mse[seed] <- mean(HI_test_results$test_error2)
      
      print(paste0(seed,'/',nseeds,' seeds finished'))
    }
    
    HI_regression_npred_results$train_mse[predictors] <- mean(HI_regression_stability$train_mse)
    HI_regression_npred_results$test_mse[predictors] <- mean(HI_regression_stability$test_mse)
    HI_regression_npred_results$test_mse_stdev[predictors] <- sqrt(var(HI_regression_stability$test_mse))
    HI_regression_npred_results$test_mse_cv[predictors] <- sqrt(var(HI_regression_stability$test_mse)) / mean(HI_regression_stability$test_mse)
    
    print(paste0(predictors,'/',npred,' predictors finished'))
  }
  
  cor(HI_test_results$test_actual, HI_test_results$test_prediction)^2
  
  ggplot(data = HI_test_results, aes(x = test_actual, y = test_prediction)) +
    geom_point() +
    xlim(c(-1,1)) +
    ylim(c(-1,1)) +
    geom_smooth(method = 'lm', se = FALSE, linetype = 'dashed', color = 'black') +
    geom_abline(slope = 1) #+
    geom_point(data = HI_train_results, aes(x = train_actual, y = train_prediction), color = 'blue', shape = 1)



# 4 way classification ----------------------------------------------------

#* Using all predictors at once --------------------------------------------


class_data <- event_summary_trim %>%
  filter(complete.cases(.)) %>%
  mutate(behavior = ifelse(FI >= 0, 'flushing', 'diluting')) %>%
  mutate(behavior = ifelse(HI >= 0, paste0(behavior, '-cw'), paste0(behavior, '-ccw'))) %>%
  select(behavior, everything()) %>%
  select(!c(FI, HI))

class_counts <- class_data %>%
  count(behavior) %>%
  mutate(frac = n/nrow(class_data))

set.seed(1)
data_split <- initial_split(class_data, prop = 4/5)
train_data <- training(data_split)
test_data <- testing(data_split)

class_engine <-
  rand_forest(trees = 1000) %>%
  set_engine('ranger') %>%
  set_mode('classification')

class_recipe <-
  recipe(behavior ~ ., data = train_data) %>%
  step_normalize(all_numeric_predictors())

class_workflow <-
  workflow() %>%
  add_model(class_engine) %>%
  add_recipe(class_recipe)
  
class_fit <-
  class_workflow %>%
  parsnip::fit(data = train_data)

class_results <-
  tibble::tibble(
    prediction = predict(class_fit, new_data = test_data)$.pred_class,
    actual = test_data$behavior) %>%
  mutate(match = prediction == actual,
         flush_pred = word(prediction, sep = '-'),
         flush_actual = word(actual, sep = '-'),
         flush_match = flush_pred == flush_actual,
         hyst_pred = word(prediction, 2, sep = '-'),
         hyst_actual = word(actual, 2, sep = '-'),
         hyst_match = hyst_pred == hyst_actual)
mean(class_results$match)
mean(class_results$flush_match)
mean(class_results$hyst_match)



#* Forward selection using QUADRANTS -------------------------------------------------------

class_data <- event_summary_trim %>%
  dplyr::select(FI, HI,
                c(FI_forest_results$predictor)) %>%
  filter(complete.cases(.)) %>%
  mutate(behavior = ifelse(FI >= 0, 'flushing', 'diluting')) %>%
  mutate(behavior = ifelse(HI >= 0, paste0(behavior, '-cw'), paste0(behavior, '-ccw'))) %>%
  select(behavior, everything()) %>%
  select(!c(FI, HI))

class_data <- event_summary_cluster %>%
  dplyr::select(behavior = cluster,
                c(vars_obs$var)) %>%
  mutate(behavior = factor(behavior)) %>%
  filter(complete.cases(.))

nseeds <- 50
npred <- ncol(class_data) - 1

class_ranger_npred_acc <- data.frame(
  var = names(class_data[2:ncol(class_data)]),
  train_acc = vector(length = npred),
  test_acc = vector(length = npred),
  test_acc_stdev = vector(length = npred),
  test_acc_cv = vector(length = npred),
  test_auc = vector(length = npred),
  fcw_auc = vector(length = npred),
  fccw_auc = vector(length = npred),
  dcw_auc = vector(length = npred),
  dccw_auc = vector(length = npred)

)

for(predictors in seq(npred)){
  
  class_ranger_stability <- data.frame(
    seed = vector(length = nseeds),
    train_acc = vector(length = nseeds),
    overall_acc = vector(length = nseeds),
    flush_acc = vector(length = nseeds),
    hyst_acc = vector(length = nseeds),
    fcw_auc = vector(length = nseeds),
    fccw_auc = vector(length = nseeds),
    dcw_auc = vector(length = nseeds),
    dccw_auc = vector(length = nseeds),
    mean_auc = vector(length = nseeds))
  
  classes <- c('flushing-cw', 'flushing-ccw', 'diluting-cw', 'diluting-ccw')
  
  class_data_select <- class_data %>%
    select(1:(predictors+1))
  
  for(seed in seq(nseeds)){
    set.seed(seed)
    data_split <- initial_split(class_data_select, prop = 4/5)
    train_data <- training(data_split)
    test_data <- testing(data_split)
    
    class_engine_ranger <-
      rand_forest(trees = 1000) %>%
      set_engine('ranger') %>%
      set_mode('classification')
    
    class_recipe_ranger <-
      recipe(behavior ~ ., data = train_data) %>%
      step_normalize(all_numeric_predictors())
    
    class_workflow_ranger <-
      workflow() %>%
      add_model(class_engine_ranger) %>%
      add_recipe(class_recipe_ranger)
    
    class_fit_ranger <-
      class_workflow_ranger %>%
      parsnip::fit(data = train_data)
    
    class_results_ranger <-
      tibble::tibble(
        prediction = predict(class_fit_ranger, new_data = test_data)$.pred_class,
        actual = test_data$behavior) %>%
      mutate(match = prediction == actual,
             flush_pred = word(prediction, sep = '-'),
             flush_actual = word(actual, sep = '-'),
             flush_match = flush_pred == flush_actual,
             hyst_pred = word(prediction, 2, sep = '-'),
             hyst_actual = word(actual, 2, sep = '-'),
             hyst_match = hyst_pred == hyst_actual)
    class_ranger_stability$seed[seed] <- seed
    class_ranger_stability$train_acc[seed] <- 1 - class_fit_ranger$fit$fit$fit$prediction.error
    class_ranger_stability$overall_acc[seed] <- mean(class_results_ranger$match)
    class_ranger_stability$flush_acc[seed] <- mean(class_results_ranger$flush_match)
    class_ranger_stability$hyst_acc[seed] <- mean(class_results_ranger$hyst_match)
    
    for(class in seq(classes)){
      class_select <- classes[class]
      class_results_selected <- class_results_ranger %>%
        select(prediction, actual) %>%
        mutate(prediction = ifelse(prediction == class_select, 1, 0),
               actual = ifelse(actual == class_select, 1, 0))
      
      class_roc <- roc(data = class_results_selected, response = actual, predictor = prediction, quiet = TRUE)
      class_auc <- as.numeric(auc(class_roc))
      class_ranger_stability[seed,class+5] <- class_auc
    }
    class_ranger_stability$mean_auc[seed] <- mean(as.numeric(class_ranger_stability[seed, 6:9]))
    
    print(paste0(seed,'/',nseeds,' seeds finished'))
  }
  
  class_ranger_npred_acc$train_acc[predictors] <- mean(class_ranger_stability$train_acc)
  class_ranger_npred_acc$test_acc[predictors] <- mean(class_ranger_stability$overall_acc)
  class_ranger_npred_acc$test_acc_stdev[predictors] <- sqrt(var(class_ranger_stability$overall_acc))
  class_ranger_npred_acc$test_acc_cv[predictors] <- sqrt(var(class_ranger_stability$overall_acc)) / mean(class_ranger_stability$overall_acc)
  class_ranger_npred_acc$test_auc[predictors] <- mean(class_ranger_stability$mean_auc)
  class_ranger_npred_acc[predictors,7:10] <- as.numeric(colMeans(class_ranger_stability[,6:9]))
  
  
  print(paste0(predictors,'/',npred,' predictors finished'))
}

ggplot(data = class_ranger_npred_acc, aes(x = seq(nrow(class_ranger_npred_acc)))) +
  geom_path(aes(y = 1 - train_acc), color = 'red', linewidth = 1) +
  geom_point(aes(y = 1 - train_acc), color = 'red', size = 2) +
  geom_path(aes(y = 1 - test_acc), color = 'blue', linewidth = 1) +
  geom_point(aes(y = 1 - test_acc), color = 'blue', size = 2) +
  geom_vline(xintercept = 14, linetype = 'dashed') +
  ylim(c(0.45,0.65)) +
  xlab('# of Predictors') +
  ylab('Prediction Error') +
  theme_classic()
ggsave('./Figures/RF_npred_performace.tiff', device = 'tiff', height = 5, width = 8, units = 'in', compression = 'lzw', dpi = 700)

ggplot(data = class_ranger_npred_acc, aes(x = seq(nrow(class_ranger_npred_acc)))) +
  geom_path(aes(y = fcw_auc), color = 'red', linewidth = 1) +
  geom_point(aes(y = fcw_auc), color = 'red', size = 2) +
  geom_vline(xintercept = 14, linetype = 'dashed') +
  #ylim(c(0.48,0.68)) +
  xlab('# of Predictors') +
  ylab('AUC') +
  theme_classic()

ggplot(data = class_ranger_stability) +
  geom_histogram(aes(x = mean_auc), bins = 10)

#* Forward selection using CLUSTERS -------------------------------------------------------
#* 
class_data <- event_summary %>%
  dplyr::select(behavior = cluster,
                everything()) %>%
  select(!c(starts_with(c('FI', 'HI')), event)) %>%
  mutate(behavior = factor(behavior)) %>%
  filter(!is.na(behavior))

nseeds <- 2
npred <- ncol(class_data) - 1

class_ranger_npred_acc <- data.frame(
  var = names(class_data[2:ncol(class_data)]),
  train_acc = vector(length = npred),
  test_acc = vector(length = npred),
  test_acc_stdev = vector(length = npred),
  test_acc_cv = vector(length = npred)
)

for(predictors in seq(npred)){
  
  class_ranger_stability <- data.frame(
    seed = vector(length = nseeds),
    train_acc = vector(length = nseeds),
    overall_acc = vector(length = nseeds))
  
  class_data_select <- class_data %>%
    select(1:(predictors+1))
  
  for(seed in seq(nseeds)){
    set.seed(seed)
    data_split <- initial_split(class_data_select, prop = 4/5)
    train_data <- training(data_split)
    test_data <- testing(data_split)
    
    class_engine_ranger <-
      rand_forest(trees = 1000) %>%
      set_engine('ranger') %>%
      set_mode('classification')
    
    class_recipe_ranger <-
      recipe(behavior ~ ., data = train_data) %>%
      step_normalize(all_numeric_predictors())
    
    class_workflow_ranger <-
      workflow() %>%
      add_model(class_engine_ranger) %>%
      add_recipe(class_recipe_ranger)
    
    class_fit_ranger <-
      class_workflow_ranger %>%
      parsnip::fit(data = train_data)
    
    class_results_ranger <-
      tibble::tibble(
        prediction = predict(class_fit_ranger, new_data = test_data)$.pred_class,
        actual = test_data$behavior) %>%
      mutate(match = prediction == actual)
    class_ranger_stability$seed[seed] <- seed
    class_ranger_stability$train_acc[seed] <- 1 - class_fit_ranger$fit$fit$fit$prediction.error
    class_ranger_stability$overall_acc[seed] <- mean(class_results_ranger$match)
    
    print(paste0(seed,'/',nseeds,' seeds finished'))
  }
  
  class_ranger_npred_acc$train_acc[predictors] <- mean(class_ranger_stability$train_acc)
  class_ranger_npred_acc$test_acc[predictors] <- mean(class_ranger_stability$overall_acc)
  class_ranger_npred_acc$test_acc_stdev[predictors] <- sqrt(var(class_ranger_stability$overall_acc))
  class_ranger_npred_acc$test_acc_cv[predictors] <- sqrt(var(class_ranger_stability$overall_acc)) / mean(class_ranger_stability$overall_acc)
  
  print(paste0(predictors,'/',npred,' predictors finished'))
}

ggplot(data = class_ranger_npred_acc, aes(x = seq(nrow(class_ranger_npred_acc)))) +
  geom_path(aes(y = 1 - train_acc), color = 'red', linewidth = 1) +
  geom_point(aes(y = 1 - train_acc), color = 'red', size = 2) +
  geom_path(aes(y = 1 - test_acc), color = 'blue', linewidth = 1) +
  geom_point(aes(y = 1 - test_acc), color = 'blue', size = 2) +
  geom_vline(xintercept = 14, linetype = 'dashed') +
  ylim(c(0,1)) +
  xlab('# of Predictors') +
  ylab('Prediction Error') +
  theme_classic()
ggsave('./Figures/RF_npred_performace.tiff', device = 'tiff', height = 5, width = 8, units = 'in', compression = 'lzw', dpi = 700)



# Classification with partykit --------------------------------------------
# Compare with party #
set.seed(6969)
data_split <- initial_split(class_data, prop = 4/5)
train_data <- training(data_split)
test_data <- testing(data_split)

class_engine_party <-
  rand_forest(trees = 1000) %>%
  set_engine('partykit') %>%
  set_mode('classification')

class_recipe_party <-
  recipe(behavior ~ ., data = train_data) %>%
  step_normalize(all_numeric_predictors())

class_workflow_party <-
  workflow() %>%
  add_model(class_engine_party) %>%
  add_recipe(class_recipe_party)

class_fit_party <-
  class_workflow_party %>%
  parsnip::fit(data = train_data)

class_results_party <-
  tibble::tibble(
    prediction = predict(class_fit_party, new_data = test_data)$.pred_class,
    actual = test_data$behavior) %>%
  mutate(match = prediction == actual,
         flush_pred = word(prediction, sep = '-'),
         flush_actual = word(actual, sep = '-'),
         flush_match = flush_pred == flush_actual,
         hyst_pred = word(prediction, 2, sep = '-'),
         hyst_actual = word(actual, 2, sep = '-'),
         hyst_match = hyst_pred == hyst_actual)
mean(class_results_party$match)
mean(class_results_party$flush_match)
mean(class_results_party$hyst_match)

#check stability
nseeds <- 10
class_party_stability <- data.frame(
  seed = vector(length = nseeds),
  overall_acc = vector(length = nseeds),
  flush_acc = vector(length = nseeds),
  hyst_acc = vector(length = nseeds)
) 
for(seed in seq(nseeds)){
  set.seed(seed)
  data_split <- initial_split(class_data, prop = 4/5)
  train_data <- training(data_split)
  test_data <- testing(data_split)
  
  class_engine_party <-
    rand_forest(trees = 1000) %>%
    set_engine('partykit') %>%
    set_mode('classification')
  
  class_recipe_party <-
    recipe(behavior ~ ., data = train_data) %>%
    step_normalize(all_numeric_predictors())
  
  class_workflow_party <-
    workflow() %>%
    add_model(class_engine_party) %>%
    add_recipe(class_recipe_party)
  
  class_fit_party <-
    class_workflow_party %>%
    parsnip::fit(data = train_data)
  
  class_results_party <-
    tibble::tibble(
      prediction = predict(class_fit_party, new_data = test_data)$.pred_class,
      actual = test_data$behavior) %>%
    mutate(match = prediction == actual,
           flush_pred = word(prediction, sep = '-'),
           flush_actual = word(actual, sep = '-'),
           flush_match = flush_pred == flush_actual,
           hyst_pred = word(prediction, 2, sep = '-'),
           hyst_actual = word(actual, 2, sep = '-'),
           hyst_match = hyst_pred == hyst_actual)
  class_party_stability$seed[seed] <- seed
  class_party_stability$overall_acc[seed] <- mean(class_results_party$match)
  class_party_stability$flush_acc[seed] <- mean(class_results_party$flush_match)
  class_party_stability$hyst_acc[seed] <- mean(class_results_party$hyst_match)
  print(paste0(seed,'/',nseeds,' seeds finished'))
}

mean(class_party_stability$overall_acc)

ggplot(data = FI_results, aes(x = actual, y = prediction)) +
  geom_point() +
  ylim(c(-1,1))

#tuning
folds <- vfold_cv(train_data, v = 5)
rf_grid <- 
  grid_latin_hypercube(
    min_n(), 
    mtry(range = c(4, 9)), 
    trees(), 
    size = 80)

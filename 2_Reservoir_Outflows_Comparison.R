library(tidyverse)
library(lubridate)
library(dataRetrieval)
library(hydroEvents)
library(rnoaa)
library(ggpubr)
library(FluMoDL)

# Load data ---------------------------------------------------------------

clinton_gage_daily <- read_csv('./DataFiles/hydro_data/reservoir_outflows/clinton_gage_daily.csv')
milford_gage_daily <- read_csv('./DataFiles/hydro_data/reservoir_outflows/milford_gage_daily.csv')
perry_gage_daily <- read_csv('./DataFiles/hydro_data/reservoir_outflows/perry_gage_daily.csv')
tuttle_gage_daily <- read_csv('./DataFiles/hydro_data/reservoir_outflows/tuttle_gage_daily.csv')
clinton_dam_daily <- read_csv('./DataFiles/hydro_data/reservoir_outflows/clinton_dam_daily.csv')
milford_dam_daily <- read_csv('./DataFiles/hydro_data/reservoir_outflows/milford_dam_daily.csv')
perry_dam_daily <- read_csv('./DataFiles/hydro_data/reservoir_outflows/perry_dam_daily.csv')
tuttle_dam_daily <- read_csv('./DataFiles/hydro_data/reservoir_outflows/tuttle_dam_daily.csv')
DeSoto_daily <- read_csv('./DataFiles/hydro_data/DeSoto_gage/DeSoto_daily.csv')
lawrence_precip_daily <- read_csv('./DataFiles/hydro_data/precip/lawrence_precip_daily.csv')
all_hydro_daily <- read_csv('./DataFiles/hydro_data/all_hydro_daily.csv')
all_hydro_daily_linterp <- read_csv('./DataFiles/hydro_data/all_hydro_daily_linterp.csv')


# Create USACE v USGS comparison data frames ------------------------------

clinton_comparison <- clinton_dam_daily %>%
  left_join(clinton_gage_daily, by = 'date') %>%
  select(date, wateryear, damQ = outflow, gageQ = Q) %>%
  filter(wateryear >= 2021)
milford_comparison <- milford_dam_daily %>%
  left_join(milford_gage_daily, by = 'date') %>%
  select(date, wateryear, damQ = outflow, gageQ = Q) %>%
  filter(wateryear >= 2021)
perry_comparison <- perry_dam_daily %>%
  left_join(perry_gage_daily, by = 'date') %>%
  select(date, wateryear, damQ = outflow, gageQ = Q) %>%
  filter(wateryear >= 2021)
tuttle_comparison <- tuttle_dam_daily %>%
  left_join(tuttle_gage_daily, by = 'date') %>%
  select(date, wateryear, damQ = outflow, gageQ = Q) %>%
  filter(wateryear >= 2021)


# USACE v USGS statistics -------------------------------------------------

clinton_r2 <- cor(clinton_comparison$damQ, clinton_comparison$gageQ)^2
milford_r2 <- cor(milford_comparison$damQ, milford_comparison$gageQ)^2
perry_r2 <- cor(perry_comparison$damQ, perry_comparison$gageQ)^2
tuttle_r2 <- cor(tuttle_comparison$damQ, tuttle_comparison$gageQ)^2


# Time series of flows ----------------------------------------------------

startyear <- 2021
endyear <- 2021

ggplot() +
  #Reservoir Outflows from USACE dam info
  geom_line(data = clinton_comparison, aes(x = date, y = damQ), color = 'red') +
  geom_line(data = milford_comparison, aes(x = date, y = damQ), color = 'gold') +
  geom_line(data = perry_comparison, aes(x = date, y = damQ), color = 'blue') +
  geom_line(data = tuttle_comparison, aes(x = date, y = damQ), color = 'purple') +
  #Reservoir Outflows from USGS gages downstream of dams
  geom_line(data = clinton_comparison, aes(x = date, y = gageQ), color = 'red', linetype = 'dashed') +
  geom_line(data = milford_comparison, aes(x = date, y = gageQ), color = 'gold', linetype = 'dashed') +
  geom_line(data = perry_comparison, aes(x = date, y = gageQ), color = 'blue', linetype = 'dashed') +
  geom_line(data = tuttle_comparison, aes(x = date, y = gageQ), color = 'purple', linetype = 'dashed') +
  theme_classic() +
  labs(title = 'USACE v USGS Outflows',
       x = NULL) +
  scale_x_date(limits = as_date(c(str_c(startyear-1, '-10-01'), str_c(endyear, '-09-30'))),
               date_breaks = str_c(endyear-startyear+1,' months'),
               date_labels = '%b %y') +
  coord_cartesian(ylim = c(0,250))


# Time series of percent difference of gage Q from dam Q ------------------

startyear <- 2021
endyear <- 2021

ggplot() +
  geom_line(data = clinton_comparison, aes(x = date, y = (gageQ-damQ)/damQ), color = 'red')+
  geom_line(data = milford_comparison, aes(x = date, y = (gageQ-damQ)/damQ), color = 'gold')+
  geom_line(data = perry_comparison, aes(x = date, y = (gageQ-damQ)/damQ), color = 'blue')+
  geom_line(data = tuttle_comparison, aes(x = date, y = (gageQ-damQ)/damQ), color = 'purple')+
  theme_classic() +
  labs(title = 'USACE v USGS Outflows',
       x = NULL) +
  scale_x_date(limits = as_date(c(str_c(startyear-1, '-10-01'), str_c(endyear, '-09-30'))),
               date_breaks = str_c(endyear-startyear+1,' months'),
               date_labels = '%b %y') +
  scale_y_continuous(limits = c(0,25))


# Scatter of absolute flows -----------------------------------------------

ggplot() +
  geom_point(data = clinton_comparison, aes(x = damQ, y = gageQ), color = 'red') +
  geom_point(data = milford_comparison, aes(x = damQ, y = gageQ), color = 'gold') +
  geom_point(data = perry_comparison, aes(x = damQ, y = gageQ), color = 'blue') +
  geom_point(data = tuttle_comparison, aes(x = damQ, y = gageQ), color = 'purple') +
  geom_abline(slope = 1, intercept = 0) +
  theme_classic()


# Scatter of normalized flows ---------------------------------------------

ggplot() +
  geom_point(data = clinton_comparison, aes(x = damQ/max(damQ), y = gageQ/max(gageQ)), color = 'red') +
  geom_point(data = milford_comparison, aes(x = damQ/max(damQ), y = gageQ/max(gageQ)), color = 'gold') +
  geom_point(data = perry_comparison, aes(x = damQ/max(damQ), y = gageQ/max(gageQ)), color = 'blue') +
  geom_point(data = tuttle_comparison, aes(x = damQ/max(damQ), y = gageQ/max(gageQ)), color = 'purple') +
  geom_abline(slope = 1, intercept = 0) +
  theme_classic()

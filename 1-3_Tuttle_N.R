library(tidyverse)
library(dataRetrieval)

all_hydro_15min <- with_tz(read_csv('./DataFiles/hydro_data/all_hydro_15min.csv'), tz = 'America/Chicago')

params <- readNWISpCode('all')
N_params <- params %>%
  filter(parameter_group_nm == 'Nutrient') %>%
  filter(str_detect(srsname, regex('Nit', ignore_case = TRUE)))

tuttle_N_info <- whatNWISdata(siteNumber = '06887000', parameterCd = c(N_params$parameter_cd)) %>%
  left_join(N_params, by = join_by('parm_cd' == 'parameter_cd'))

tuttle_N <- readNWISqw('06887000', '00618', startDate = '1990-01-01', tz = 'America/Chicago') %>%
  select(dateTime = startDateTime, N = result_va)

ggplot(data = all_hydro_15min) +
  geom_line(aes(x = dateTime, y = tuttleQ)) +
  geom_line(aes(x = dateTime, y = N*50), color = 'red') +
  geom_point(data = tuttle_N, aes(x = dateTime, y = N*50), color = 'darkgreen') +
  scale_x_datetime(limits = as_datetime(c('2019-07-01', '2020-07-01')))






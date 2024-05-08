library(tidyverse)
source('load_and_combine.R')
source('Theme+Settings.R')
min <- read_csv(file.path('DataFiles','hydro_data','all_hydro_15min.csv'))
daily <- read_csv(file.path('DataFiles','hydro_data','all_hydro_daily.csv'))
event_summary <- mutate(event_summary, month = factor(month, levels = c('Jan', 'Feb', 'Mar', 
                                                                        'Apr', 'May', 'Jun',
                                                                        'Jul', 'Aug', 'Sep',
                                                                        'Oct', 'Nov', 'Dec')))

ggplot(daily, aes(x = N)) +
  geom_histogram() +
  geom_histogram(data = event_summary, aes(x = max_N), fill = 'red')

ggplot(event_summary, aes(x = max_N)) +
  geom_histogram()

table(cut(event_summary$max_N, 10))

summary(min$N)
quantile(min$N, seq(0,1,by=0.05),na.rm = TRUE)
quantile(event_summary$max_N, seq(0,1,by=0.05),na.rm = TRUE)


ggplot(data = min, aes(x = riverQ, y = N)) +
  geom_point()


#number of high N events per year
high_n <- filter(event_summary, max_N >= 2)
ggplot(high_n) +
  geom_bar(aes(x = wateryear))

ggplot(event_summary) +
  geom_bar(aes(x = wateryear), color = 'black', fill = 'grey80') +
  geom_bar(data = high_n, aes(x = wateryear), fill = 'maroon', color = 'black')

ggplot(event_summary) +
  geom_bar(aes(x = season), color = 'black', fill = 'grey80') +
  geom_bar(data = high_n, aes(x = season), fill = 'maroon', color = 'black') +
  facet_wrap(vars(wateryear))
#event stats per year
year_sum <- group_by(event_summary, wateryear) %>%
  summarise(across(c(initial_Q, total_Q, delta_rat_Q, initial_N, p_10270104_1d_total), ~mean(.x, na.rm = TRUE)))
vars_sel <- c('initial_Q', 'd_10260010_spei12')
select(event_summary, c(wateryear, all_of(vars_sel))) %>%
  pivot_longer(!wateryear, names_to = 'var', values_to = 'val') %>%
  ggplot(.) +
  geom_boxplot(aes(x = wateryear, y = val)) +
  facet_wrap(vars(var), scales = 'free_y')
ggplot(event_summary) +
  geom_boxplot(aes(x = wateryear, y = initial_Q))



#water year stats
years <- group_by(daily, wateryear) %>%
  summarise(mean_Q = mean(riverQ, na.rm = TRUE),
            total_Q = sum(riverQ, na.rm = TRUE))
ggplot(years) +
  geom_col(aes(x = factor(wateryear), y = mean_Q))


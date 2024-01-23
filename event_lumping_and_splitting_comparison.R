library(tidyverse)
library(ggpubr)


# Load data ---------------------------------------------------------------
all_hydro_15min <- read_csv('./DataFiles/hydro_data/all_hydro_15min.csv')

events_1 <- read_csv('./DataFiles/hydro_data/event_summary.csv') %>%
  select(FI, HI = HI_mean, start_dateTime, end_dateTime) %>%
  mutate(set = 'set1')

events_2 <- read_csv('./DataFiles/hydro_data/event_summary_compound1.csv')%>%
  select(FI, HI = HI_mean, start_dateTime, end_dateTime) %>%
  mutate(set = 'set3')

events_3 <- read_csv('./DataFiles/hydro_data/event_summary_compound2.csv')%>%
  select(FI, HI = HI_mean, start_dateTime, end_dateTime) %>%
  mutate(set = 'set2')

events_all <- rbind(events_1, events_2, events_3) %>%
  mutate(wateryear = if_else(month(start_dateTime) > 9, year(start_dateTime) + 1, year(start_dateTime)))


# Plot Indices --------------------------------------------------------------------
ggplot(data = events_all, aes(x = FI, y = HI, color = set)) +
  geom_point()

ggplot(data = events_all, aes(x = FI, fill = set)) +
  geom_histogram() +
  facet_wrap(vars(set), scales = 'free_y')

ggplot(data = events_all, aes(x = HI, fill = set)) +
  geom_histogram() +
  facet_wrap(vars(set), scales = 'free_y')

ggplot(data = events_all, aes(x = FI, y = HI, color = set)) +
  geom_point()

pivot_longer(events_all, c(FI, HI), names_to = 'index', values_to = 'val') %>%
  ggplot(aes(x = index, y = val, fill = set)) +
  geom_boxplot()
       
combs <- combn(unique(events_all$set), 2)
for(comb in 1:ncol(combs)){
  print(paste0('FI ',combs[1,comb],'-',combs[2,comb],': ',
  ks.test(events_all$FI[events_all$set == combs[1,comb]], 
          events_all$FI[events_all$set == combs[2,comb]])$p.value))
  print(paste0('HI ',combs[1,comb],'-',combs[2,comb],': ',
  ks.test(events_all$HI[events_all$set == combs[1,comb]], 
          events_all$HI[events_all$set == combs[2,comb]])$p.value))
}



# Plot Time Series --------------------------------------------------------
s_year <- 2019
e_year <- 2019

set1 <- ggplot() +
  geom_line(data = subset(all_hydro_15min, wateryear %in% seq(s_year, e_year)), 
            aes(x = dateTime, y = riverQ)) +
  geom_rect(data = subset(events_all, wateryear %in% seq(s_year, e_year) & set == 'set1'),
            aes(xmin = start_dateTime, xmax = end_dateTime), color = 'red',
            ymin = -1000, ymax = 5000, fill = 'grey', alpha = 0.5) +
  ylim(c(0, max(subset(all_hydro_15min, wateryear %in% seq(s_year, e_year))$riverQ, na.rm = TRUE))) +
  theme_classic()

set2 <- ggplot() +
  geom_line(data = subset(all_hydro_15min, wateryear %in% seq(s_year, e_year)), 
            aes(x = dateTime, y = riverQ)) +
  geom_rect(data = subset(events_all, wateryear %in% seq(s_year, e_year) & set == 'set2'),
            aes(xmin = start_dateTime, xmax = end_dateTime), color = 'green',
            ymin = -1000, ymax = 5000, fill = 'grey', alpha = 0.5) +
  ylim(c(0, max(subset(all_hydro_15min, wateryear %in% seq(s_year, e_year))$riverQ, na.rm = TRUE))) +
  theme_classic()

set3 <- ggplot() +
  geom_line(data = subset(all_hydro_15min, wateryear %in% seq(s_year, e_year)), 
            aes(x = dateTime, y = riverQ)) +
  geom_rect(data = subset(events_all, wateryear %in% seq(s_year, e_year) & set == 'set3'),
            aes(xmin = start_dateTime, xmax = end_dateTime), color = 'blue',
            ymin = -1000, ymax = 5000, fill = 'grey', alpha = 0.5) +
  ylim(c(0, max(subset(all_hydro_15min, wateryear %in% seq(s_year, e_year))$riverQ, na.rm = TRUE))) +
  theme_classic()

ggarrange(set1, set2, set3, nrow = 3)


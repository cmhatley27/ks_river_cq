# Load --------------------------------------------------------------------
library(tidyverse)
library(patchwork)
library(conover.test)
source('8_combine_all_event_data.R')
source('Theme+Settings.R')

plot_dat <- event_summary %>%
  mutate(cluster = factor(cluster),
         res_end = cut(reservoir_runoff_ratio, breaks = c(0,0.5,1.5), labels = c('Precip','Reservoir')))

res_means <- plot_dat %>%
  group_by(res_end) %>%
  summarise(mean_fi = mean(FI_n),
            mean_hi = mean(HI_n))


# scatter -----------------------------------------------------------------
scatter <- ggplot(plot_dat, aes(x = FI_n, y = HI_n, fill = season, shape = res_end)) +
  geom_hline(yintercept = 0, color = 'grey') +
  geom_vline(xintercept = 0, color = 'grey') +
  geom_point(size = 2, stroke = 0.25) +
  # geom_point(data = res_means, aes(x = mean_fi, y = mean_hi, color = res_end), shape = 15, size = 3) +
  # geom_point(x = mean(plot_dat$FI_n), y = mean(plot_dat$HI_n), color = 'red', shape = 15, size = 3) + #global mean
  scale_y_continuous(position = 'right', limits = c(-1,1)) +
  scale_fill_manual(name = 'Season', values = c('#8DA9C4','#60992D','#C33C54','#F4B860'), guide = 'none') +
  scale_shape_manual(values = c(21,24), guide = 'none') +
  xlab('Flushing Index') +
  ylab('Hysteresis Index') +
  xlim(-1,1) +
  theme(axis.title = element_text(face = 'bold', size = 9))
  # theme(text = element_text(size = 18)) #For poster
scatter
scatter_width <- 4
scatter_height <- 3.5
ggsave(file.path('Figures','final','index scatter','index_scatter.png'), height = scatter_height, width = scatter_width, units = 'in')


# revised box plots -------------------------------------------------------
subset_width <- 1
season_width <- 1.5
#FI subset
ggplot(plot_dat, aes(fill = res_end, x = FI_n, y=fct_rev(res_end))) +
  geom_vline(xintercept = 0, color = 'grey') +
  geom_boxplot(size = 0.25) +
  scale_fill_manual(limits = c('Reservoir','Precip'),
                    values = c('#d73027','#4575b4'), guide = 'none') +
  scale_y_discrete(position = 'right', name = '', labels = c(-0.5, 0.5)) +
  scale_x_continuous(position = 'top', name = NULL, labels = NULL, limits = c(-1,1))
ggsave(file.path('Figures','final','index scatter','fi_subsets.png'), height = subset_width, width = scatter_width, units = 'in')

#FI seasons
ggplot(plot_dat, aes(fill = fct_rev(season), x = FI_n, y=fct_rev(res_end))) +
  geom_vline(xintercept = 0, color = 'grey') +
  geom_boxplot(size = 0.25) +
  scale_fill_manual(values = season_palette, guide = 'none') +
  scale_y_discrete(position = 'right', name = '', labels = c(-0.5, 0.5)) +
  scale_x_continuous(position = 'top', name = NULL, labels = NULL, limits = c(-1,1))
ggsave(file.path('Figures','final','index scatter','fi_seasons.png'), height = season_width, width = scatter_width, units = 'in')


#HI subset
ggplot(plot_dat, aes(fill = res_end, y = HI_n, x=res_end)) +
  geom_hline(yintercept = 0, color = 'grey') +
  geom_boxplot(size = 0.25) +
  scale_fill_manual(limits = c('Reservoir','Precip'),
                    values = c('#d73027','#4575b4'), guide = 'none') +
  scale_x_discrete(position = 'bottom', name = '', labels = c(-0.5,0.5)) +
  scale_y_continuous(position = 'left', name = NULL, labels = NULL, limits = c(-1,1))
ggsave(file.path('Figures','final','index scatter','hi_subsets.png'), height = scatter_height, width = subset_width, units = 'in')

#HI seasons
ggplot(plot_dat, aes(fill = season, y = HI_n, x=res_end)) +
  geom_hline(yintercept = 0, color = 'grey') +
  geom_boxplot(size = 0.25) +
  scale_fill_manual(values = season_palette, guide = 'none') +
  scale_x_discrete(position = 'bottom', name = '', labels = c(-0.5,0.5)) +
  scale_y_continuous(position = 'left', name = NULL, labels = NULL, limits = c(-1,1))
ggsave(file.path('Figures','final','index scatter','hi_seasons.png'), height = scatter_height, width = season_width, units = 'in')



# inset loops -------------------------------------------------------------
hydro_dat <- read_csv('./DataFiles/hydro_data/event_data/event_data.csv')
selected_events <- c(20,139,21,127)
# selected_events <- c(20,139,167,255) #events from methods fig
event_sel <- 3
event_dat <- filter(hydro_dat, event == selected_events[event_sel])

ggplot(event_dat, aes(x = Q_norm, y = N_smooth_norm, color = dateTime)) +
  geom_path() +
  scale_color_gradient2(low = 'black',mid = 'blue', high = 'red', midpoint = event_dat$dateTime[event_dat$Q_norm == 1],guide = 'none') +
  scale_x_continuous(name = NULL, breaks = NULL) +
  scale_y_continuous(name = NULL, breaks = NULL) +
  theme(axis.title.y = element_text(size = 6, angle = 0, vjust = 0.5),
        axis.title.x = element_text(size = 6),
        plot.margin = unit(rep(0,4), 'inches'))
ggsave(file.path('Figures','final','index scatter',paste0('inset_q',event_sel,'.png')), height = 0.45, width = 0.45, units = 'in')
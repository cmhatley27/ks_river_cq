library(tidyverse)
source('Theme+Settings.R')
source('8_combine_all_event_data.R')

characteristics_sel <- c('initial_Q','duration_since_last',
                         'max_Q','delta_Q','delta_rat_Q','duration',
                         'reservoir_runoff_ratio','r_t_share','r_m_share','r_p_share')

cor_dat <- event_summary %>%
  mutate(res_end = cut(reservoir_runoff_ratio, breaks = c(0,0.5,1.5), labels = c('Precip','Reservoir'))) %>%
  select(c(FI = FI_n, HI = HI_n, res_end, all_of(characteristics_sel))) %>%
  pivot_longer(c(FI,HI), names_to = 'target', values_to = 'target_val')

plot_dat <- cor_dat %>%
  group_by(target, res_end) %>%
  summarise(across(all_of(characteristics_sel), ~cor(.x, target_val, method = 'spearman'))) %>%
  pivot_longer(!c(target, res_end), names_to = 'var', values_to = 'cor') %>%
  mutate(var = factor(var, levels = characteristics_sel),
         cor = round(cor, digits = 3))


# make color sheets -------------------------------------------------------
target_vars <- c('FI','HI')
for(target_sel in 1:2){
  ggplot(subset(plot_dat, target == target_vars[target_sel]), aes(x = var, y = fct_rev(res_end), fill = cor)) +
    geom_tile(width = 2, height = 1, color = 'black') +
    scale_fill_gradient2(limits = c(-0.32, 0.32), low = '#a50026', mid = '#f3faaf', high = '#006837', guide = 'none') +
    theme_void()
  ggsave(file.path('Figures','final','correlation table',paste0(target_vars[target_sel],'_colors.png')),
         height = 1, width = 6, units = 'in', dpi = 500)
}

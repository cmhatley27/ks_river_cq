library(tidyverse)
library(smapr)

all_hydro_daily <- read_csv('./DataFiles/hydro_data/all_hydro_daily.csv')

set_smap_credentials('cmhatley', 'Rfh91895!Rfh91895!', overwrite = TRUE)

find_smap('SPL4SMAU', dates = '2016-05-27', version = 7)

smap <- download_smap(find_smap('SPL4SMAU', dates = '2016-05-27', version = 7))

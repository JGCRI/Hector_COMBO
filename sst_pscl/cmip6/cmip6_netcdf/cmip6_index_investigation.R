# checking cmip6 index

library(dplyr)
library(tidyr)
library(tibble)


files <- as.tibble(read.csv('cmip6/cmip6_netcdf/cmip6_archive_index.csv', stringsAsFactors = F))

files %>%
  filter(variable == 'tos' | variable == 'tas')  %>%
  select(model, variable, experiment, domain) %>%  
  filter(experiment == "ssp245" | experiment == 'historical') %>%
  distinct %>%
  filter(grepl('mon', .$domain)) %>%
  arrange(model) %>%
  as.data.frame()



files %>%
  filter(variable == 'tos')  %>% 
  select(model, experiment) %>%
  filter(experiment == 'ssp245' | experiment == 'historical') %>%
  distinct %>%
  arrange(model) %>%
  as.data.frame()

files %>%
  filter(variable == 'tas')  %>% 
  select(model, experiment)%>%
  filter(experiment == 'ssp245' | experiment == 'historical' | experiment == 'esm-hist') %>%
  distinct

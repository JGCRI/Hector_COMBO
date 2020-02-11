# checking cmip6 index

library(dplyr)
library(tidyr)
library(tibble)


files <- as.tibble(read.csv('cmip6/cmip6_netcdf/cmip6_archive_index.csv', stringsAsFactors = F))

files %>%
  filter(grepl('.nc', file),
         variable == 'tos' | variable == 'tas')  %>%
  select(model, variable, experiment, domain) %>%  
  filter(experiment == "ssp245" | experiment == 'historical') %>%
  distinct %>%
  filter(grepl('mon', .$domain)) %>%
  arrange(model) %>%
  filter((variable == 'tas' & experiment == 'ssp245') | variable == 'tos') ->
  files2 


files2 %>%
  as.data.frame() 


# arrange so that just get a list of models with completed data.
files2 %>%
  unite(united, c('variable', 'experiment'), sep = '~') %>% 
  mutate(val = 1) %>%
  select(-domain) %>%
  distinct %>%
  spread(united, val) %>%
  na.omit %>%
  select(model)

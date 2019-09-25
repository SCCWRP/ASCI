library(tidyverse)

fls <- list.files('C:/Users/Marcus/Desktop/quants_24Sept2019/', full.names = T)

# all data, refcalmeans and min/max
dat <- enframe(fls) %>% 
  group_by(value) %>% 
  mutate(
    rawdat = purrr::map(value, read.csv, stringsAsFactors = F)
  ) %>% 
  ungroup %>% 
  select(-name) %>% 
  mutate(
    value = basename(value), 
    value = gsub('\\_', '.', value),
    value = gsub('reaser|\\.csv$', '', value), 
    value = gsub('\\.quants', '', value)
  )

# refcalmeans only
refcalmeans <- dat %>% 
  filter(grepl('refcalmean', value)) %>% 
  unnest(rawdat) %>% 
  mutate(
    value = gsub('\\.refcalmean$', '', value),
    X = case_when(
      !grepl('\\_raw$', X) ~ paste(X, 'mod', sep = '_'), 
      T ~ X
    )
  ) %>% 
  rename(
    Metric = X,
    RefCalMean = x,
    Assemblage = value
  )

# minmax only
minmax <- dat %>% 
  filter(!grepl('refcalmean', value)) %>% 
  mutate(rawdat = purrr::map(rawdat, gather, key = 'Metric', value = 'minmax')) %>%
  unnest(rawdat) %>% 
  filter(!grepl('^X$', Metric)) %>% 
  separate(minmax, c('Min', 'Max'), sep = ' ') %>% 
  separate(value, c('Assemblage', 'StressResponse', 'type'), sep = '\\.') %>% 
  mutate(
    Min = as.numeric(Min), 
    Max = as.numeric(Max)
  ) %>% 
  unite('Metric', Metric, type)
  
# combine refcalmeans and minmax
omni.ref <- full_join(minmax, refcalmeans, by = c('Metric', 'Assemblage')) %>% 
  select(Metric, Min, Max, StressResponse, RefCalMean, Assemblage)
  
# save
mmilkup$omni.ref <- omni.ref

save(mmilkup, file = 'data/mmilkup.RData', compress = 'xz', version = 2)

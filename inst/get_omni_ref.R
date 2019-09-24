fls <- list.files('C:/Users/Marcus.SCCWRP2K/Desktop/quants_24Sept2019/quants_24Sept2019/', full.names = T)

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
    value = gsub('reaser|\\.csv$', '', value)
  )

refcalmeans <- dat %>% 
  filter(grepl('refcalmean', value)) %>% 
  unnest(rawdat) %>% 
  mutate(value = gsub('\\.refcalmean$', '', value)) %>% 
  rename(
    Metric = X,
    RefCalMean = x,
    Assemblage = value
  )

minmax <- dat %>% 
  filter(!grepl('refcalmean', value)) %>% 
  mutate(
    rawdat = purrr::map(rawdat, gather)
  ) %>% 
  unnest(rawdat)
    {
      
      browser()
    })
  )
  unnest(rawdat)
  
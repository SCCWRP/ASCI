library(dplyr)
library(RPostgreSQL)
library(dbplyr)

# connection
con <- DBI::dbConnect(
  RPostgreSQL::PostgreSQL(),
  dbname = 'smc',
  host = '192.168.1.17',
  user = rstudioapi::askForPassword('Database user'),
  password = rstudioapi::askForPassword("Database password")
)

# gis data
datgiscon <- tbl(con, 'tblgismetrics') %>% 
  collect()

# algae data
datalgcon <- tbl(con, 'tmp_algae') %>%
  collect()

# intersect(datgiscon$stationcode, datalgcon$stationcode) %>% sample(3) %>% sort %>% dput
togrb <- c("404M07357", "801M16916", "909M24937")

demo_algae_tax <- datalgcon %>%
  filter(stationcode %in% togrb) %>%
  select(
    StationCode = stationcode,
    SampleDate = sampledate,
    Replicate = replicate,
    SampleTypeCode = sampletypecode,
    BAResult = baresult,
    Result = result,
    FinalID = finalid
  ) %>% 
  mutate(
    Result = as.numeric(Result)
  )

demo_station <- datgiscon %>%
  filter(stationcode %in% togrb) %>%
  select(
    StationCode = stationcode,
    CondQR50 = condqr50,
    SITE_ELEV = site_elev,
    TEMP_00_09 = temp_00_09,
    KFCT_AVE = kfct_ave,
    AtmCa = atmca,
    PPT_00_09 = ppt_00_09,
    MAX_ELEV = max_elev,
    CaO_Mean = cao_mean,
    MgO_Mean = mgo_mean,
    S_Mean = s_mean,
    UCS_Mean = ucs_mean,
    LPREM_mean =lprem_mean,
    AtmMg = atmmg,
    AtmSO4 = atmso4,
    MINP_WS = minp_ws,
    MEANP_WS = meanp_ws,
    SumAve_P = sumave_p,
    TMAX_WS = tmax_ws,
    XWD_WS = xwd_ws,
    MAXWD_WS = maxwd_ws,
    LST32AVE = lst32ave,
    BDH_AVE = bdh_ave,
    PRMH_AVE = prmh_ave,
    PSA6C = psa6c
  ) %>%
  mutate(XerMtn = NA)

rmv <- chkinp(demo_algae_tax, demo_station, getval  = T)
demo_algae_tax <- demo_algae_tax %>% 
  filter(!FinalID %in% rmv)

save(demo_algae_tax, file = 'data/demo_algae_tax.RData', compress = 'xz', version = 2)
save(demo_station, file = 'data/demo_station.RData', compress = 'xz', version = 2)


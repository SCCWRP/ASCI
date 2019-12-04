#' Score samples using the ASCI tool
#'
#' @param taxa \code{data.frame} for input taxonomy data
#' @param station \code{data.frame} for input station data
#' 
#' @details 
#' One index for three taxonomy types are scored, MMI for diatoms, soft-bodied algae, and hybrid. 
#' If only soft-bodied algae or diatoms are present, only the respective ASCI and metrics are returned.
#' This function outputs the results of \code{\link{mmifun}} functions in a user-friendly format.
#' 
#' @return 
#' A dataframe with all metrics calculated for each provided taxa
#'
#' @export
#' 
#' @importFrom dplyr bind_rows mutate select case_when mutate_all group_by ungroup inner_join summarize full_join
#' @importFrom magrittr "%>%"
#' @importFrom tidyr gather spread unnest unite
#' @import purrr
#' @import tibble
#' 
#' @seealso  \code{\link{mmifun}}
#' 
#' @examples 
#' # calculate all
#' ASCI(demo_algae_tax, demo_station)
#' 
#' # works if either soft-bodied or diatoms are missing
#' # remove diatoms from station sample 909M24937
#' tmp <- subset(demo_algae_tax, !(StationCode == '909M24937' & SampleTypeCode == 'Integrated'))
#' ASCI(tmp, demo_station)
#' 
#' # works if either soft-bodied or diatoms are missing
#' # remove soft-bodied from station sample 801M16916
#' tmp <- subset(demo_algae_tax, !(StationCode == '801M16916' & SampleTypeCode != 'Integrated'))
#' ASCI(tmp, demo_station)
ASCI <- function(taxa, station){

  # run all other checks, get output if passed
  dat <- chkinp(taxa, station)
  txrmv <- dat$txrmv
  dat <- dat$taxa
  
  # calculate GIS from stations
  station <- calcgis(station)
  
  # mmi
  mmind <- mmifun(dat, station)
  
  ##
  # main output (scores)
  mmiscr <- mmind %>% 
    map(~ .x$MMI_scores) %>% 
    map(gather, 'Metric', 'Value', -SampleID) %>% 
    enframe('taxa') %>% 
    unnest(cols = value) %>% 
    mutate_all(~ replace(., is.na(.), NA))
  
  ##
  # supplementary info
  Supp1_mmi <- mmind %>% 
    map(~ .x$MMI_supp) %>% 
    map(gather, 'Metric', 'Value', -SampleID) %>% 
    enframe('taxa') %>% 
    unnest(cols = value) %>% 
    mutate_all(~ replace(., is.na(.), NA))
  
  # get diatom valve counts
  extra1 <- dat %>% 
    filter(SampleTypeCode == 'Integrated') %>% 
    group_by(SampleID) %>%
    summarize(
      D_ValveCount = sum(BAResult, na.rm = T)
    )
  
  # get sampletype concatenated column, soft-bodied entity and biovolume count
  extra2 <- dat %>% 
    group_by(SampleID) %>% 
    summarize(
      SampleType = paste0(unique(SampleTypeCode), collapse = ', '),
      S_EntityCount = sum(BAResult, na.rm = T),
      S_Biovolume = sum(Result, na.rm = T)
    ) %>% 
    full_join(extra1, by = 'SampleID')

  # combine all
  out <- rbind(mmiscr, Supp1_mmi) %>% 
    mutate(
    taxa = case_when(
      taxa == 'diatoms' ~ 'D',
      taxa == 'sba' ~ 'S',
      TRUE ~ 'H'
      )
    ) %>% 
    unite('Met', c('taxa', 'Metric'), sep = '_') %>% 
    group_by(SampleID, Met) %>% 
    mutate(grouped_id = dplyr::row_number()) %>% 
    spread(Met, Value) %>% 
    select(-grouped_id) %>% 
    ungroup() 
  out1 <- extra2 %>% 
    inner_join(out, by = 'SampleID') %>% 
    filter(SampleID != 1)
  
  # get original stationcode, date, and replicate
  # add unrecognized taxa
  out1 <- getids(out1, concatenate = FALSE) %>% 
    left_join(txrmv, by = 'SampleID')
  
  # unholy column selection
  colsel <- c("SampleID", "StationCode", "SampleDate", "Replicate", "SampleType", 
              
              "D_ValveCount", "S_EntityCount", "S_Biovolume", "D_NumberTaxa", 
              "S_NumberTaxa", "H_NumberTaxa", "UnrecognizedTaxa", "D_ASCI", "S_ASCI", "H_ASCI", 
              
              # "D_cnt.spp.IndicatorClass_TP_low_raw", "D_cnt.spp.IndicatorClass_TP_low_raw_score", 
              # "D_prop.spp.Saprobic.BM_raw", "D_prop.spp.Saprobic.BM_raw_score", 
              # "D_prop.spp.SPIspecies4_mod", "D_prop.spp.SPIspecies4_mod_score", 
             #  "D_Salinity.BF.richness_mod", "D_Salinity.BF.richness_mod_score", 
             #  "D_pcnt.attributed.IndicatorClass_TP_Low", "D_pcnt.attributed.Salinity.BF", 
             #  "D_pcnt.attributed.Saprobic.BM", "D_pcnt.attributed.SPIspecies4", 
             #  "H_OxyRed.DO_30.richness_mod", "H_OxyRed.DO_30.richness_mod_score", "H_prop.spp.BCG4_mod",
             #  "H_prop.spp.BCG4_mod_score",  "H_prop.spp.IndicatorClass_DOC_high_raw","H_prop.spp.IndicatorClass_DOC_high_raw_score", 
             #  "H_Salinity.BF.richness_mod", "H_Salinity.BF.richness_mod_score",  
             #  "H_pcnt.attributed.OxyRed.DO_30", "H_pcnt.attributed.BCG4", "H_pcnt.attributed.IndicatorClass_DOC_high",
             #  "H_pcnt.attributed.Salinity.BF",
             #  "S_prop.spp.BCG45_raw", "S_prop.spp.BCG45_raw_score", "S_prop.spp.Green_raw", 
             #  "S_prop.spp.Green_raw_score", "S_cnt.spp.IndicatorClass_DOC_high_raw", 
             #  "S_cnt.spp.IndicatorClass_DOC_high_raw_score", "S_pcnt.attributed.BCG45", "S_pcnt.attributed.Green", 
             #  "S_pcnt.attributed.IndicatorClass_DOC_high"
              NULL
              )

  out1 <- out1[, colsel]
  
  
  return(out1)
  
}
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
      S_EntityCount = sum(BAResult[which(SampleTypeCode != 'Integrated')], na.rm = T),
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
  #out1 <- getids(out1, concatenate = FALSE) %>% 
  out1 <- out1 %>%
    left_join(txrmv, by = 'SampleID') %>%
    #dplyr::select(-c("StationCode","SampleDate","Replicate")) %>%
    left_join(
      dat %>% dplyr::distinct(
        SampleID,StationCode,SampleDate,Replicate
      ),
      by = "SampleID"
    )
  # In the above code we also grab stationcode sampledate and replicate from original input data
  # getids function was sometimes not grabbing them correctly
  
  # Besides, as Rafi has suggested in the past, extracting fields by parsing a string may not
  # be the best way to go. He suggested doing it this way and 
  # I think it may be a good approach
    
  
  
  
  
  # Tack on a column that warns them if there is only diatom or soft body at a site
  warnings_column <- dat %>%
    # if Result is 0 it counts as it not being there
    dplyr::filter(((Result != 0) & is.na(BAResult)) | ((BAResult != 0) & is.na(Result))) %>% 
    dplyr::group_by(SampleID) %>%
    dplyr::mutate(
      dia_or_sba = ifelse(SampleTypeCode == "Integrated", "dia", "sba")
    ) %>% 
    dplyr::summarise(
      Comments = dplyr::case_when(
        ("dia" %in% unique(dia_or_sba)) & (!("sba" %in% unique(dia_or_sba))) ~ "Warning - Only Diatom data present",
        ("sba" %in% unique(dia_or_sba)) & (!("dia" %in% unique(dia_or_sba))) ~ "Warning - Only Soft Body data present",
        TRUE ~ ""
      )
    )
  
  beginning_cols <- c('SampleID', 'StationCode','SampleDate','Replicate',
                      'SampleType','D_ValveCount','S_EntityCount','S_Biovolume','D_NumberTaxa',
                      'S_NumberTaxa','H_NumberTaxa','UnrecognizedTaxa',
                      'D_ASCI','S_ASCI','H_ASCI')
  
  # Here we tack on that warnings column
  out1 <- out1 %>% dplyr::left_join(warnings_column, by = "SampleID") %>%
    dplyr::select(
      c(beginning_cols, names(out1)[which(!names(out1) %in% beginning_cols)], "Comments")
    ) 
  
  names(out1) <- gsub("pcnt","pct",names(out1))
  names(out1) <- gsub("prop","prp",names(out1))
  names(out1) <- gsub("attributed","att",names(out1))
  names(out1) <- gsub("score","scr",names(out1))
  names(out1) <- gsub("IndicatorClass","IC",names(out1))
  names(out1) <- gsub("OxyRed","OxRd",names(out1))
  names(out1) <- gsub("OxyReq","OxRq",names(out1))
  names(out1) <- gsub("DO_100orDO_75","DO100_75",names(out1))
  names(out1) <- gsub("\\.","_",names(out1))

  
  return(out1)
  
}

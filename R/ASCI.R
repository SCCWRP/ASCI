#' Score samples using the ASCI tool
#'
#' @param taxa \code{data.frame} for input taxonomy data
#' @param station \code{data.frame} for input station data
#' @param tax chr string indicating output to return from a specific taxa, must one to many of \code{'diatoms'}, \code{'sba'}, or \code{'hybrid'}, defaults to all
#' @param ... additional arguments passed to other funcions
#' 
#' @details 
#' One index for three taxonomy types are scored, MMI for diatoms, soft-bodied algae, and hybrid. 
#' This function outputs the reulsts of \code{\link{mmifun}} functions in a user-friendly format.
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
#' results <- ASCI(demo_algae_tax)
#' 
ASCI <- function(taxa, station, tax = c('diatoms', 'sba', 'hybrid'), ...){
  
  # check tax argument
  if(any(!tax %in% c('diatoms', 'sba', 'hybrid')))
    stop('tax must match diatoms, sba, and/or hybrid')
  
  # run all other checks, get output if passed
  dat <- chkinp(taxa)
  
  # mmi
  mmind <- mmifun(dat, station, ...)
  
  ##
  # main output (scores)
  mmiscr <- mmind %>% 
    map(~ .x$MMI_scores) %>% 
    map(gather, 'Metric', 'Value', -SampleID) %>% 
    enframe('taxa') %>% 
    unnest %>% 
    mutate_all(~ replace(., . < 0 | is.na(.), NA))
  
  ##
  # supplementary info
  Supp1_mmi <- mmind %>% 
    map(~ .x$MMI_supp) %>% 
    map(gather, 'Metric', 'Value', -SampleID) %>% 
    enframe('taxa') %>% 
    unnest %>% 
    mutate_all(~ replace(., . < 0 | is.na(.), NA))
  
  # subset taxa if needed
  if(length(tax) < 3){
    
    mmiscr <- mmiscr %>% 
      filter(taxa %in% tax) %>% 
      mutate_all(~ replace(., . < 0 | is.na(.), NA))
    
    Supp1_mmi <- Supp1_mmi %>% 
      filter(taxa %in% tax) %>% 
      mutate_all(~ replace(., . < 0 | is.na(.), NA))
    
  }
  
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
      SampleType = paste0(unique(SampleTypeCode), collapse = '|'),
      S_EntityCount = sum(BAResult, na.rm = T),
      S_Biovolume = sum(Result, na.rm = T),
      UnrecognizedTaxa = paste0(setdiff(FinalID, STE$FinalID), collapse = '|')
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
  out1 <- getids(out1, concatenate = FALSE)
  
  # reorder columns
  out1 <- out1 %>% 
    select(SampleID, StationCode, SampleDate, Replicate, SampleType, D_ValveCount, S_EntityCount, S_Biovolume, 
           D_NumberTaxa, S_NumberTaxa, H_NumberTaxa, UnrecognizedTaxa, D_ASCI, S_ASCI, H_ASCI, dplyr::contains('D_'), 
           dplyr::contains('S_'), dplyr::contains('H_'), dplyr::everything()
    )
  
  return(out1)
  
}


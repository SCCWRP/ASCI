#' @title Check input taxonomy and site data
#'
#' @description Check input taxonomy and site for required information
#' 
#' @param taxa \code{data.frame} for input taxonomy data
#' @param station \code{data.frame} for input station data
#' @param getval logical to return a vector of values not satisfied by checks, useful for data prep
#'
#' @return A two element list of the original data (named \code{taxa}) and removed taxa by \code{SampleID} (named \code{txrmv})
#' if all checks are met.  The original data also includes a new column for \code{SampleID}.  An 
#' error message is returned if the datasetsdo not meet requirements or a vector of values that caused the error if \code{getval = TRUE}.  
#' Site data will include only those sites in the taxonomic data.
#' 
#' @details 
#' The following are checked:
#' \itemize{
#' \item Required columns in taxonomy data: StationCode, SampleDate, Replicate, SampleTypeCode, BAResult, Result, FinalID
#' \item Taxonomic names are present in the \code{\link{STE}} reference file
#' \item Sites include both diatom and soft-bodied algae data (warning if not)
#' \item No missing abundance values for diatoms (for rarification)
#' \item One of \code{CondQR50} or all predictors for the conductivity model in the station data
#' \item One of \code{XerMtn} or \code{PSA6C} in the station data
#' \item Additional required columns for the station data: StationCode, CondQR50, SITE_ELEV, TEMP_00_09, KFCT_AVE, AtmCa, PPT_00_09, MAX_ELEV
#' \item No missing data in additional required columns for stationdata
#' }
#' 
#' @export
#' 
#' @seealso \code{\link{calcgis}}
#'
#' @importFrom dplyr arrange filter group_by mutate select summarise
#' @importFrom magrittr "%>%"
#' @importFrom tidyr gather
#' @importFrom stringr str_squish
#' 
#' @examples
#' # all checks passed, data returned with SampleID
#' chkinp(demo_algae_tax, demo_station)
#' 
#' # errors
#' \dontrun{
#' # missing columns in taxa data
#' tmp <- demo_algae_tax[, 1, drop = FALSE]
#' chkinp(tmp, demo_station)
#' chkinp(tmp, demo_station, getval = TRUE)
#' 
#' # incorrect taxonomy
#' tmp <- demo_algae_tax
#' tmp[1, 'FinalID'] <- 'asdf'
#' chkinp(tmp, demo_station)
#' chkinp(tmp, demo_station, getval = TRUE)
#' 
#' # missing diatom data at sites, returns only a warning
#' tmp <- merge(demo_algae_tax, STE, all.x = T) %>%
#'   filter(Class %in% 'Bacillariophyceae')
#' chkinp(tmp, demo_station)
#' 
#' 
#' # missing abundance data for diatoms
#' tmp <- demo_algae_tax
#' tmp$BAResult <- NA
#' chkinp(tmp, demo_station)
#' chkinp(tmp, demo_station, getval = TRUE)
#' 
#' # stations not shared between taxa and station
#' tmp <- demo_station[-1, ]
#' chkinp(demo_algae_tax, tmp)
#' 
#' # missing both of XerMtn and PSA6C in station
#' tmp <- demo_station[, !names(demo_station) %in% c('XerMtn', 'PSA6C')]
#' chkinp(demo_algae_tax, tmp)
#' 
#' # missing CondQR50 and incomplete predictor fields
#' tmp <- demo_station[, !names(demo_station) %in% c('CondQR50', 'TMAX_WS')]
#' chkinp(demo_algae_tax, tmp)
#' 
#' # missing remaining station predictors
#' tmp <- demo_station[, !names(demo_station) %in% c('AtmCa')]
#' chkinp(demo_algae_tax, tmp)
#' 
#' # missing data in remaining station predictors
#' tmp <- demo_station
#' tmp$AtmCa[2] <- NA
#' chkinp(demo_algae_tax, tmp)
#' }

chkinp <- function(taxa, station, getval = FALSE){
  
  # Replace -88's with NA
  # doesnt work for character columns
  # taxa <- taxa %>% mutate_all(~dplyr::na_if(.,-88))
  
  # trying a more 'old school' approach
  taxa[taxa == -88] <- NA_real_
  
  # Clean up the STE
  STE <- STE %>% mutate(
    FinalID = str_trim(FinalID),
    FinalIDassigned = str_trim(FinalIDassigned)
  )
  

  ##
  # check if required columns are present in taxa
  taxcols <- c('StationCode', 'SampleDate', 'Replicate',
               'SampleTypeCode', 'BAResult', 'Result', 'FinalID')
  chk <- taxcols %in% names(taxa)
  if(any(!chk)){
    
    chk <- taxcols[!chk]
    if(getval) return(chk)
    
    msg <- paste(chk, collapse = ', ') %>% 
      paste('Required columns not found taxa:', .)
    stop(msg, call. = FALSE)
  }
  
  taxa <- taxa %>%
    left_join(
      taxa_crosswalk %>% 
        rename(FinalID = NewName)
      ,
      by = 'FinalID'
    ) %>% 
    mutate(
      FinalID = case_when(
        is.na(OldName_in_ASCI_Calculator) ~ FinalID,
        T ~ OldName_in_ASCI_Calculator
      )
    ) %>%
    select(-OldName_in_ASCI_Calculator)
  
  
  
  # Reassure that all columns are the correct datatype
  taxa <- taxa %>%
    dplyr::mutate(
      StationCode = str_squish(as.character(StationCode)),
      # Nothing in the code necessarily demands that the SampleDate needs to be a Date field (to my knowledge)
      # For that reason, we will not coerce it to a date, in case they put different date formats or something like that
      #SampleDate = as.POSIXct(SampleDate),
      Replicate = as.integer(as.character(Replicate)),
      Replicate = replace(Replicate, which(is.na(Replicate)), -88),
      SampleTypeCode = str_squish(as.character(SampleTypeCode)),
      BAResult = as.integer(as.character(BAResult)),
      Result = as.numeric(as.character(Result)),
      FinalID = str_squish(as.character(FinalID)),
      SampleID = paste(StationCode, SampleDate, Replicate, sep = "_")
    ) %>% 
    dplyr::select(SampleID, StationCode, SampleDate, Replicate, SampleTypeCode, BAResult, Result, FinalID)
  
  
  
  ##
  # add id values after columns are checked
  #taxa <- getids(taxa)  
  
  # check if They put integrated SampleTypeCode for a FinalID of a Soft Body Algae organism
  # Standard to distinguish diatoms and soft body technically is the phylum - Bacillariophyta
  # Algae with phylum "Bacillariophyta" are diatoms
  # If the Phylum is not Bacillariophyta, then it is soft body
  tmp <- taxa %>% 
    left_join(STE %>% select(FinalID, Phylum), by = 'FinalID') %>% 
    filter((Phylum != 'Bacillariophyta') & (SampleTypeCode == 'Integrated'))
  
  chk <- (nrow(tmp) > 0)
  if(chk){
    offendingTaxa <- tmp$FinalID %>% unique() %>% paste( collapse = ', ')
    
    msg <- paste(
      'Error in data submitted. ASCI scores not calculated. Soft algae taxa present in the Integrated fraction: \n', 
      offendingTaxa,
      'Please correct and resubmit'
    )
    
    stop(msg)
  }

  # check if They put Microalgae, Macroalgae, or Qualitative SampleTypeCode for a Diatom
  # Standard to distinguish diatoms and soft body technically is the phylum - Bacillariophyta
  # Algae with phylum "Bacillariophyta" are diatoms
  # If the Phylum is not Bacillariophyta, then it is soft body
  tmp <- taxa %>% 
    left_join(STE %>% select(FinalID, Phylum), by = 'FinalID') %>% 
    filter((Phylum == 'Bacillariophyta') & (SampleTypeCode %in% c('Microalgae', 'Macroalgae', 'Qualitative') ))
  
  chk <- (nrow(tmp) > 0)
  if(chk){
    offendingTaxa <- tmp$FinalID %>% unique() %>% paste( collapse = ', ')
    
    msg <- paste(
      'Error in data submitted. ASCI scores not calculated. Diatom taxa present in the Microalgae, Macroalgae, or Qualitative fraction: \n', 
      offendingTaxa,
      'Please correct and resubmit'
    )
    
    stop(msg)
  }
  
  
  
  ##
  # check if sites have both diatom and sba data
  tmp <-taxa %>% 
    merge(STE, all.x = TRUE) %>% 
    mutate(diaind = ifelse(SampleTypeCode %in% 'Integrated', 'dia', 'sba')) %>% 
    group_by(SampleID) %>%
    summarise(n = length(unique(diaind)))
  chk <- tmp$n < 2
  if(any(chk)){
    
    chk <- tmp[chk, ] %>% 
      .$SampleID
    msg <- paste(chk, collapse = ', ') %>% 
      paste('Only diatom or soft-bodied algae present at sites:', .)
    warning(msg)
  }
  
  ## 
  # check if abundance diatom data available in taxonomy
  chk <- taxa %>% 
    merge(STE, all.x = T) %>% 
    filter(Class %in% 'Bacillariophyceae') %>% 
    group_by(SampleID) %>% 
    summarise(chk = any(is.na(BAResult)))
  if(any(chk$chk)){

    chk <- chk$SampleID[chk$chk]
    if(getval) return(chk)
    
    msg <- paste(chk, collapse = ', ') %>% 
      paste('Missing abundance data for diatoms', .)
    #stop(msg, call. = FALSE)
    
  }

  # check if stationcode in stations match those in taxa
  chk <- setdiff(taxa$StationCode, station$StationCode)
  if(length(chk) > 0){
    
    if(getval) return(chk)
      
    msg <- paste(chk, collapse = ', ') %>% 
      paste('Station names not shared between taxa and station:', .)
    stop(msg, call. = FALSE)
    
  }

  
  # check for duplicates in stations data
  chk <- station %>% 
    group_by(StationCode) %>% 
    filter(n() > 1) %>%
    distinct(StationCode)
  if (nrow(chk) > 0) {
    stop(
      paste(
        "Duplicate records in station data for the following sites:\n",
        paste(chk, collapse = ',\n')
      )
    )
  }
  
    
  ##
  # check station columns
  
  condvals <- c("CaO_Mean", "MgO_Mean", "S_Mean", "UCS_Mean", "LPREM_mean", 
                "AtmCa", "AtmMg", "AtmSO4", "MINP_WS", "MEANP_WS", "SumAve_P", 
                "TMAX_WS", "XWD_WS", "MAXWD_WS", "LST32AVE", "BDH_AVE", "KFCT_AVE", 
                "PRMH_AVE")
  regvals <- c("StationCode", 
               "AREA_SQKM",
                "KFCT_AVE",
                #"XerMtn",
                "PPT_00_09",
                #"CondQR50",
                "TMAX_WS",
                "MAX_ELEV",
                "SITE_ELEV",
                "AtmCa")

  # must have one of or both XerMtn and PSA6C
  chk <- c('XerMtn', 'PSA6C') %in% names(station)
  if(sum(chk) == 0){
    
    stop('Station data must include one of XerMtn or PSA6C', call. = FALSE)
    
  }

  # must have one of or both CondQR50 and condvals  
  chk <- !'CondQR50' %in% names(station) & sum(names(station) %in% condvals) < 18
  if(chk){

    msg <- paste(condvals, collapse = ', ') %>% 
      paste('Station data must include CondQR50 or all of the following:', .)
    stop(msg, call. = FALSE)
    
  }
  
  # must have all of regvals
  chk <- setdiff(regvals, names(station))
  if(length(chk) != 0){
    
    msg <- paste(chk, collapse = ', ') %>% 
      paste('Station data missing the following:', .)
    stop(msg, call. = FALSE)    
    
  }
  
  # for regvals, no missing data
  # message if NA values in conductivity predictors
  chk <- station[, names(station) %in% regvals] %>% 
    apply(2, anyNA) %>% 
    which %>% 
    names
  if(length(chk) > 0){
    
    msg <- chk %>% 
      paste(collapse = ', ') %>% 
      paste('Missing values in GIS predictors:', .)
    stop(msg, call. = FALSE)
    
  }

  # Station columns that need to be character:
  # StationCode, PSA6, PSA6C, County
  # These seem to be read in as factors sometimes
  station <- station %>% dplyr::mutate_if(is.factor, as.character)
  
  
  ##
  # find and remove unidentified taxonomy
  chk <- setdiff(taxa$FinalID, STE$FinalID)
  #if(length(chk) > 0){
  if (TRUE) {  # in the case that length(chk) == 0, txrmv is not getting defined. Thus causing it to break when it references it later on
    # removed FinalID by SampleID
    txrmv <- taxa %>% 
      filter(FinalID %in% chk) %>% 
      select(SampleID, FinalID) %>% 
      unique %>% 
      group_by(SampleID) %>% 
      summarise(UnrecognizedTaxa = paste(FinalID, collapse = ', ') )
    
    # purge unidentified from taxa
    taxa <- taxa %>% 
      filter(!FinalID %in% chk)
    
  }
  
  ##
  # return list of taxa and removed FinalID by SampleID
  if (exists("txrmv") ==T) { 
  out <- list(
    taxa = taxa,
    txrmv = txrmv
  ) 
  
  return(out) } 
  
  if (exists("txrmv") ==F) { 
  out <- list(
    taxa = taxa
  ) 
  
  return(out) } 
  
  
  
}



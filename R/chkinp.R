#' @title Check input taxonomy and site data
#'
#' @description Check input taxonomy and site data for required information
#' 
#' @param taxain \code{data.frame} for input taxonomy data
#' @param getval logical to return a vector of values not satisfied by checks, useful for data prep
#'
#' @return The original data with only relevant columns are returned if all checks are met, including a new column for \code{SampleID} (see \code{\link{getids}}).  An error message is returned if the datasets do not meet requirements or a vector of values that caused the error if \code{getval = TRUE}.  Site data will include only those sites in the taxonomic data.
#' 
#' @details 
#' The following are checked:
#' \itemize{
#' \item Required columns in taxonomy data: StationCode, SampleDate, Replicate, SampleTypeCode, BAResult, Result, FinalID
#' \item Required columns in site data: StationCode, SampleDate, Replicate, AREA_SQKM, AtmCa, CondQR50, DayOfYear, KFCT_AVE, LogWSA, LST32AVE, MAX_ELEV, MEANP_WS, MINP_WS, New_Lat, New_Long, PPT_00_09, SITE_ELEV, TMAX_WS, XWD_WS 
#' \item SampleID in taxonomic data are present in site data
#' \item Taxonomic names are present in the \code{\link{STE}} reference file
#' \item Sites include both diatom and soft-bodied algae data
#' \item No missing abundance values for diatoms (for rarification)
#' \item No missing values for environmental predictor variables
#' }
#' 
#' @export
#' 
#' @seealso \code{\link{getids}}
#'
#' @importFrom dplyr arrange filter group_by mutate select summarise
#' @importFrom magrittr "%>%"
#' @importFrom tidyr gather
#' 
#' @examples
#' # all checks passed, data returned with SampleID
#' chkinp(demo_algae_tax)
#' 
#' # errors
#' \dontrun{
#' # missing columns in taxa data
#' tmp <- demo_algae_tax[, 1, drop = FALSE]
#' chkinp(tmp)
#' chkinp(tmp, getval = TRUE)
#' 
#' # incorrect taxonomy
#' tmp <- demo_algae_tax
#' tmp[1, 'FinalID'] <- 'asdf'
#' chkinp(tmp)
#' chkinp(tmp, getval = TRUE)
#' 
#' # missing diatom data at sites
#' tmp <- merge(demo_algae_tax, STE, all.x = T) %>%
#'   filter(!Class %in% 'Bacillariophyceae')
#' chkinp(tmp)
#' chkinp(tmp, getval = TRUE)
#' 
#' # missing abundance data for diatoms
#' tmp <- demo_algae_tax
#' tmp$BAResult <- NA
#' chkinp(tmp)
#' chkinp(tmp, getval = TRUE)

chkinp <- function(taxain, getval = FALSE){
  
  ##
  # check if required columns are present in taxain
  taxcols <- c('StationCode', 'SampleDate', 'Replicate',
               'SampleTypeCode', 'BAResult', 'Result', 'FinalID')
  chk <- taxcols %in% names(taxain)
  if(any(!chk)){
    
    chk <- taxcols[!chk]
    if(getval) return(chk)
    
    msg <- paste(chk, collapse = ', ') %>% 
      paste('Required columns not found taxain:', .)
    stop(msg, call. = FALSE)
  }
  
  ##
  # add id values after columns are checked
  taxain <- getids(taxain)  


  ##
  # check taxonomy names
  chk <- setdiff(taxain$FinalID, STE$FinalID)
  if(length(chk) > 0){

    if(getval) return(chk)

    msg <- paste(chk, collapse = ', ') %>%
      paste('Unrecognized taxa:', .)
    stop(msg, call. = FALSE)

    }

  ##
  # check if sites have both diatom and sba data
  tmp <-taxain %>% 
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
  chk <- taxain %>% 
    merge(STE, all.x = T) %>% 
    filter(Class %in% 'Bacillariophyceae') %>% 
    group_by(SampleID) %>% 
    summarise(chk = any(is.na(BAResult)))
  if(any(chk$chk)){

    chk <- chk$SampleID[chk$chk]
    if(getval) return(chk)
    
    msg <- paste(chk, collapse = ', ') %>% 
      paste('Missing abundance data for diatoms', .)
    stop(msg, .call = FALSE)
    
  }

  
  ##
  # return if all checks met
  out <- taxain
  
  return(out)
  
}



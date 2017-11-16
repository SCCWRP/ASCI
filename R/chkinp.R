#' @title Check input taxonomy and site data
#'
#' @description Check input taxonomy and site data for required information
#' 
#' @param taxain \code{data.frame} for input taxonomy data
#' @param sitein \code{data.frame} for input site data
#' @param getval logical to return a vector of values not satisfied by checks, useful for data prep
#'
#' @return The original data are returned if all checks are met, including a new column for \code{SampleID} (see \code{\link{getids}}).  An error message is returned if the datasets do not meet requirements or a vector of values that caused the error if \code{getval = TRUE}.
#' 
#' @details 
#' The following are checked:
#' \itemize{
#' \item Required columns in taxonomy data: StationCode, SampleDate, Replicate,SampleTypeCode, BAResult, Result, FinalID
#' \item SampleID in taxonomic data are present in site data
#' \item Taxonomic names are present in the \code{\link{STE}} reference file
#' \item No missing abundance values for diatoms
#' }
#' 
#' @export
#' 
#' @seealso \code{\link{getids}}
#'
#' @importFrom dplyr group_by summarise
#' @importFrom magrittr "%>%"
#' 
#' @examples
#' # all checks passed, data returned with SampleID
#' chkinp(demo_algae_tax, demo_algae_sitedata)
#' 
#' # errors
#' \dontrun{
#' # missing columns
#' tmp <- demo_algae_tax[, -c(1, 2)]
#' chkinp(tmp, demo_algae_sitedata)
#' chkinp(tmp, demo_algae_sitedata, getval = TRUE)
#' 
#' # site data not found for taxonomic site
#' tmp <- demo_algae_sitedata[-1, ]
#' chkinp(demo_algae_tax, tmp)
#' chkinp(demo_algae_tax, tmp, getval = TRUE)
#' 
#' # incorrect taxonomy
#' tmp <- demo_algae_tax
#' tmp[1, 'FinalID'] <- 'asdf'
#' chkinp(tmp, demo_algae_sitedata)
#' chkinp(tmp, demo_algae_sitedata, getval = TRUE)
#' 
#' # missing abundance data
#' tmp <- demo_algae_tax
#' tmp$BAResult <- NA
#' chkinp(tmp, demo_algae_sitedata)
#' chkinp(tmp, demo_algae_sitedata, getval = TRUE)
#' }
chkinp <- function(taxain, sitein, getval = FALSE){
  
  ##
  # check if required columns are present
  cols <- c('StationCode', 'SampleDate', 'Replicate','SampleTypeCode', 'BAResult', 'Result', 'FinalID')
  chk <- cols %in% names(taxain)
  if(any(!chk)){
    
    chk <- cols[!chk]
    if(getval) return(chk)
    
    msg <- paste(chk, collapse = ', ') %>% 
      paste('Required columns not found:', .)
    stop(msg, call. = FALSE)
    
  }
  
  ##
  # add id values after columns are checked
  if(!'SampleID' %in% names(taxain))
    taxain <- getids(taxain)
  if(!'SampleID' %in% names(sitein))
    sitein <- getids(sitein)
  sitein <- sitein %>% 
    filter(SampleID %in% taxain$SampleID)
  
  ##
  # check if all sites in taxain are in sitein
  chk <- unique(taxain$SampleID) %in% sitein$SampleID
  if(any(!chk)){

    chk <- unique(taxain$SampleID)[!chk]
    if(getval) return(chk)

    msg <- paste(chk, collapse = ', ') %>%
      paste('SampleID in taxonomic data not found in site data: ', .)
    stop(msg, call. = FALSE)

  }

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
  out <- list(
    taxain = taxain, 
    sitein = sitein
  )
  
  return(out)
  
}



#' @title Check input taxonomy and site data
#'
#' @description Check input taxonomy and site data for required information
#' 
#' @param taxain \code{data.frame} for input taxonomy data
#' @param sitein \code{data.frame} for input site data
#' @param getval logical to return a vector of values not satisfied by checks, useful for data prep
#'
#' @return The original data with only relevant columns are returned if all checks are met, including a new column for \code{SampleID} (see \code{\link{getids}}).  An error message is returned if the datasets do not meet requirements or a vector of values that caused the error if \code{getval = TRUE}.
#' 
#' @details 
#' The following are checked:
#' \itemize{
#' \item Required columns in taxonomy data: StationCode, SampleDate, Replicate, SampleTypeCode, BAResult, Result, FinalID
#' \item Required columns in site data: StationCode, SampleDate, Replicate, AREA_SQKM, AtmCa, CondQR50, DayOfYear, KFCT_AVE, LogWSA, LST32AVE, MAX_ELEV, MEANP_WS, MINP_WS, New_Lat, New_Long, PPT_00_09, SITE_ELEV, TMAX_WS, XWD_WS 
#' \item SampleID in taxonomic data are present in site data
#' \item Taxonomic names are present in the \code{\link{STE}} reference file
#' \item Sites with soft-bodied algae also include diatom data
#' \item No missing abundance values for diatoms
#' \item No missing values for environmental predictor variables
#' }
#' 
#' @export
#' 
#' @seealso \code{\link{getids}}
#'
#' @importFrom dplyr group_by select summarise
#' @importFrom magrittr "%>%"
#' @importFrom tidyr gather
#' 
#' @examples
#' # all checks passed, data returned with SampleID
#' chkinp(demo_algae_tax, demo_algae_sitedata)
#' 
#' # errors
#' \dontrun{
#' # missing columns in taxa data
#' tmp <- demo_algae_tax[, 1, drop = FALSE]
#' chkinp(tmp, demo_algae_sitedata)
#' chkinp(tmp, demo_algae_sitedata, getval = TRUE)
#' 
#' # missing columns in site data
#' tmp <- demo_algae_sitedata[, 1, drop = FALSE]
#' chkinp(demo_algae_tax, tmp)
#' chkinp(demo_algae_tax, tmp, getval = TRUE)
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
#' # missing diatom data at a site
#' tmp <- merge(demo_algae_tax, STE, all.x = T) %>%
#'   filter(!Class %in% 'Bacillariophyceae')
#' chkinp(tmp, demo_algae_sitedata)
#' chkinp(tmp, demo_algae_sitedata, getval = TRUE)
#' 
#' # missing abundance data for diatoms
#' tmp <- demo_algae_tax
#' tmp$BAResult <- NA
#' chkinp(tmp, demo_algae_sitedata)
#' chkinp(tmp, demo_algae_sitedata, getval = TRUE)
#' 
#' # missing environmental data
#' tmp <- demo_algae_sitedata
#' tmp$SITE_ELEV <- NA
#' tmp$AREA_SQKM <- NA
#' chkinp(demo_algae_tax, tmp)
#' chkinp(demo_algae_tax, tmp, getval = TRUE)
#' }
chkinp <- function(taxain, sitein, getval = FALSE){
  
  ##
  # check if required columns are present in taxain
  taxcols <- c('StationCode', 'SampleDate', 'Replicate','SampleTypeCode', 'BAResult', 'Result', 'FinalID')
  chk <- taxcols %in% names(taxain)
  if(any(!chk)){
    
    chk <- taxcols[!chk]
    if(getval) return(chk)
    
    msg <- paste(chk, collapse = ', ') %>% 
      paste('Required columns not found taxain:', .)
    stop(msg, call. = FALSE)
    
  }
  
  ## 
  # check if required columns are present in sitein
  sitcols <- c(
    row.names(diatom_rf_oe$importance), 
    row.names(sba_rf_oe$importance),
    row.names(hybrid_rf_oe$importance)
    ) %>% 
    unique %>% 
    sort %>% 
    c('StationCode', 'SampleDate', 'Replicate', .)
  chk <- sitcols %in% names(sitein)
  if(any(!chk)){
    
    chk <- sitcols[!chk]
    if(getval) return(chk)
    
    msg <- paste(chk, collapse = ', ') %>% 
      paste('Required columns not found in sitein:', .)
    stop(msg, call. = FALSE)
    
  }
  
  ##
  # add id values after columns are checked
  taxain <- getids(taxain)    
  sitein <- getids(sitein)
  # sitein <- sitein %>% 
  #   filter(SampleID %in% taxain$SampleID)
  
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
  # check if diatom data are available at all sites
  taxain <- merge(taxain, STE, all.x = T)
  dia_sit <- subset(taxain, Class %in% 'Bacillariophyceae') %>% 
    .$SampleID %>% 
    unique
  oth_sit <-subset(taxain, !Class  %in% 'Bacillariophyceae') %>% 
    .$SampleID %>% 
    unique
  chk <- oth_sit[!oth_sit %in% dia_sit]
  if(length(chk) > 0){
    
    if(getval) return(chk)
    
    msg <- paste(chk, collapse = ', ') %>% 
      paste('No diatoms at sites:', .)
    stop(msg, .call = FALSE)
    
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
  # check if predictor data available for oe rf models
  chk <- sitein[, c('SampleID', sitcols)] %>% 
    gather('var', 'val', -SampleID) %>% 
    filter(is.na(val))
  if(nrow(chk) > 0){
    
    if(getval) return(chk)
    
    msg <- unique(chk$SampleID)
    for(i in seq_along(msg)){

      vrs <- chk %>% 
        filter(SampleID %in% msg[i]) %>% 
        .$var
      msg[i] <- paste(vrs, collapse = ', ') %>% 
        paste0(msg[i], ':\t', .)
      
    }

    msg <- paste(msg, collapse = '\n\n') %>% 
      paste('\n\nMissing environmental data:\n\n', .)
    stop(msg, .call = FALSE)
    
  }

  ##
  # return if all checks met
  out <- list(
    taxain = taxain[, c('SampleID', taxcols)], 
    sitein = sitein[, c('SampleID', sitcols)]
  )
  
  return(out)
  
}



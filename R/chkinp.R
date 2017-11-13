#' @title Check input taxonomy and site data
#'
#' @description Check input taxonomy and site data for required information
#' 
#' @param taxain \code{data.frame} for input taxonomy data
#' @param sitein \code{data.frame} for input site data
#'
#' @return \code{NULL} if all checks are met, otherwise an informative error message is returned
#' 
#' @details 
#' The following is checked:
#' \itemize{
#' \item Required columns in taxonomy data: StationCode, SampleDate, Replicate,SampleTypeCode, BAResult, Result, FinalID
#' }
#' 
#' @export
#'
#' @importFrom magrittr "%>%"
#' 
#' @examples
#' chkinp(demo_algae_tax)
chkinp <- function(taxain, sitein){
  
  # check if required columns are present
  cols <- c('StationCode', 'SampleDate', 'Replicate','SampleTypeCode', 'BAResult', 'Result', 'FinalID')
  chk <- cols %in% names(taxain)
  if(!any(chk)){
    
    msg <- cols[!cols] %>% 
      paste(collapse = ', ')
    
    stop('Columns ', msg, ' not found', call. = FALSE)
    
  }
  
  return(NULL)
  
}



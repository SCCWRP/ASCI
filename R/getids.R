#' @title Assign unique sample id
#' 
#' @description Assing unique sample id to rows of taxonomy and site data
#'
#' @param datin \code{data.frame} of taxonomy or site data
#'
#' @details Assigns unique sample id based on a concatenation of StationCode, SampleDate, and Replicate
#'
#' @return The original input with a new 'SampleID' column
#' 
#' @export
#' 
#' @importFrom magrittr "%>%"
#' @importFrom tidyr unite
#' 
#' @examples
#' getids(demo_algae_sitedata)
#' getids(demo_algae_tax)
getids <- function(datin){
  
  datin <- datin %>% 
    unite('SampleID', StationCode, SampleDate, Replicate, sep = '_', remove = FALSE)

  return(datin)
  
}
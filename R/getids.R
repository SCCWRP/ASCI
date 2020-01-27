#' @title Assign unique sample id
#' 
#' @description Assing unique sample id to rows of taxonomy and site data
#'
#' @param datin \code{data.frame} of taxonomy or site data
#' @param concatenate logical indicating if \code{SampleID} is created or parsed 
#' @details Assigns unique sample id based on a concatenation of StationCode, SampleDate, and Replicate, used internally in \code{\link{chkinp}} and \code{\link{ASCI}}
#'
#' @return The original input with a new \code{SampleID} column if \code{concatenate = TRUE}, or the original input with columns for \code{StationCode}, \code{Date}, and \code{Replicate} if \code{concatenate = FALSE}
#' 
#' @export
#' 
#' @importFrom magrittr "%>%"
#' @importFrom tidyr unite
#' 
#' @seealso \code{\link{chkinp}}, \code{\link{calcgis}}
#' 
#' @examples
#' getids(demo_algae_tax)
#' 
getids <- function(datin, concatenate = T){
  
  if(!'SampleID' %in% names(datin) & concatenate)
  
    datin <- datin %>% 
      unite('SampleID', StationCode, SampleDate, Replicate, sep = '_', remove = FALSE)

  if('SampleID' %in% names(datin) & !concatenate)

    datin <- datin %>%
      mutate(
        col = SampleID,
        Replicate = gsub('.*_([0-9]+)$', '\\1', col),
        Replicate = as.numeric(Replicate),
        col = gsub('_[0-9]+$', '', col),
        SampleDate = gsub('.*_(.*)$', '\\1', col), 
        StationCode = gsub(paste(SampleDate, collapse = '_'), '', col), 
        StationCode = gsub('_$', '', StationCode)
      ) %>% 
      dplyr::select(-col)
        
  return(datin)
  
}
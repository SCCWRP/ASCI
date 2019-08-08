#' Score samples using the ASCI tool
#'
#' @param taxain \code{data.frame} for input taxonomy data
#' @param tax chr string indicating output to return from a specific taxa, must one to many of \code{'diatoms'}, \code{'sba'}, or \code{'hybrid'}, defaults to all
#' @param ... additional arguments passed to other funcions
#' 
#' @details 
#' One index for three taxonomy types are scored, MMI for diatoms, soft-bodied algae, and hybrid. 
#' This function outputs the reulsts of \code{\link{mmifun}} functions in a user-friendly format.
#' 
#' @return 
#' A \code{\link{asci}} object with specific methods.  See the examples for accessing.
#'
#' @export
#' 
#' @importFrom dplyr bind_rows mutate select
#' @importFrom magrittr "%>%"
#' @importFrom tidyr gather spread unnest
#' @import purrr
#' @import tibble
#' 
#' @seealso  \code{\link{mmifun}}
#' 
#' @examples 
#' results <- ASCI(demo_algae_tax)
#' scores(results)
#' Supp1_mmi(results)
ASCI <- function(taxain, tax = c('diatoms', 'sba', 'hybrid'), ...){
  
  # check tax argument
  if(any(!tax %in% c('diatoms', 'sba', 'hybrid')))
    stop('tax must match diatoms, sba, and/or hybrid')
  
  # run all other checks, get output if passed
  dat <- chkinp(taxain)
  
  # mmi
  mmind <- mmifun(dat, ...)
  
  ##
  # main output (scores)
  mmiscr <- mmind %>% 
    map(~ .x$MMI_scores) %>% 
    map(gather, 'met', 'val', -SampleID) %>% 
    enframe('taxa') %>% 
    unnest %>% 
    spread(met, val)
  
  ##
  # supplementary info
  Supp1_mmi <- mmind %>% 
    map(~ .x$MMI_supp) %>% 
    map(gather, 'Metric', 'Value', -SampleID) %>% 
    enframe('taxa') %>% 
    unnest
  
  # subset taxa if needed
  if(length(tax) < 3){
    
    mmiscr <- mmiscr %>% 
      filter(taxa %in% tax) %>% 
      mutate_all(~ replace(., . < 0 | is.na(.), NA))
    
    Supp1_mmi <- Supp1_mmi %>% 
      filter(taxa %in% tax) %>% 
      mutate_all(~ replace(., . < 0 | is.na(.), NA))
    
  }
  
  ##
  # create asci class output
  out <- asci(scores = mmiscr, Supp1_mmi = Supp1_mmi, taxa = tax)
  
  return(out)
  
}
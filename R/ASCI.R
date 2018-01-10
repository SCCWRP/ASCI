#' Score samples using the ASCI tool
#'
#' @param taxain \code{data.frame} for input taxonomy data
#' @param sitein \code{data.frame} for input site data
#' @param tax chr string indicating output to return from a specific taxa, must one to many of \code{'diatoms'}, \code{'sba'}, or \code{'hybrid'}, defaults to all
#' @param ... additional arguments passed to other funcions, e.g., \code{\link{rfpred}}
#' 
#' @details 
#' Three taxonomy types are scored for the pMMI: diatoms, soft-bodied algae, and hybrid. This function combines output from the \code{\link{pmmifun}} function in a user-friendly format.
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
#' @seealso  \code{\link{pmmifun}}
#' 
#' @examples 
#' results <- ASCI(demo_algae_tax, demo_algae_sitedata)
#' scores(results)
#' Supp1_mmi(results)
ASCI <- function(taxain, sitein, tax = c('diatoms', 'sba', 'hybrid'), ...){
  
  ## sanity checks
  
  # check tax argument
  if(any(!tax %in% c('diatoms', 'sba', 'hybrid')))
    stop('tax must match diatoms, sba, and/or hybrid')
  
  # run all other checks, get output if passed
  dat <- chkinp(taxain, sitein)
  
  ##
  # individual output
  
  # pmmi
  pmmind <- pmmifun(dat$taxa, dat$site, ...)
    
  ##
  # main output (scores)

  # pmmi
  pmmiscr <- pmmind %>% 
    map(~ .x$MMI_scores) %>% 
    map(gather, 'met', 'val', -SampleID) %>% 
    enframe('taxa') %>% 
    unnest

  # combine scrs, long format
  scr <- pmmiscr %>% 
    spread(met, val) %>% 
    select(taxa, SampleID, MMI, MMI_Percentile)
  
  ##
  # supplementary info

  # pmmi
  Supp1_mmi <- pmmind %>% 
    map(~ .x$MMI_supp) %>% 
    map(gather, 'Metric', 'Value', -SampleID) %>% 
    enframe('taxa') %>% 
    unnest
  
  # subset taxa if needed
  if(length(tax) < 3){
    
    scr <- scr %>% 
      filter(taxa %in% tax)
    
    Supp1_mmi <- Supp1_mmi %>% 
      filter(taxa %in% tax)
    
  }
    
  ##
  # create asci class output
  out <- asci(scores = scr, Supp1_mmi = Supp1_mmi, taxa = tax)
  
  return(out)
  
}
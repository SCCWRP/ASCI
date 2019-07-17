#' Title
#'
#' @param taxain 
#' @param sitein 
#' @param tax 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
ASCI_1 <- function(taxain, sitein, tax = c('diatoms', 'sba', 'hybrid'), ...){
  
  ## sanity checks
  
  # check tax argument
  if(any(!tax %in% c('diatoms', 'sba', 'hybrid')))
    stop('tax must match diatoms, sba, and/or hybrid')
  
  # run all other checks, get output if passed
  dat <- chkinp(taxain, sitein)
  
  ##
  # individual output
  
  # oe
  # oeind <- oefun(dat$taxa, dat$site, ...)
  
  # mmi
  mmind <- mmifun(dat$taxa, dat$site, ...)
  
  ##
  # main output (scores)
  
  # oe
  # oescr <- oeind %>% 
  #   map(function(x){
  #     
  #     x$OE_scores %>% 
  #       select(SampleID, O, E, OoverE, OoverE_Percentile) %>% 
  #       gather('met', 'val', -SampleID)
  #     
  #   }) %>% 
  #   enframe('taxa') %>% 
  #   unnest
  # 
  # pmmi
  mmiscr <- mmind %>% 
    map(~ .x$MMI_scores) %>% 
    map(gather, 'met', 'val', -SampleID) %>% 
    enframe('taxa') %>% 
    unnest
  
  # combine scrs, long format
  # scr <- bind_rows(oescr, pmmiscr) %>% 
  #   spread(met, val) %>% 
  #   select(taxa, SampleID, MMI, MMI_Percentile, O, E, OoverE, OoverE_Percentile)
  # 
  ##
  # supplementary info
  
  # pmmi
  Supp1_mmi <- mmind %>% 
    map(~ .x$MMI_supp) %>% 
    map(gather, 'Metric', 'Value', -SampleID) %>% 
    enframe('taxa') %>% 
    unnest
  
  # oe capture probs
  # Supp1_OE <- oeind %>% 
  #   map(~ .x$Capture_Probs) %>% 
  #   enframe('taxa') %>% 
  #   unnest
  
  # oe group occurrence probs
  # Supp2_OE <- oeind %>% 
  #   map(~ .x$Group_Occurrence_Probs) %>% 
  #   map(gather, 'pGroup', 'Prob', -SampleID) %>% 
  #   enframe('taxa') %>% 
  #   unnest
  
  # oe null
  # null_OE <- oeind %>% 
  #   map(function(x){
  #     
  #     x$OE_scores %>% 
  #       select(SampleID, Onull, Enull, OoverE.null) %>% 
  #       gather('met', 'val', -SampleID)
  #     
  #   }) %>% 
  #   enframe('taxa') %>% 
  #   unnest
  
  # subset taxa if needed
  if(length(tax) < 3){
    
    scr <- scr %>% 
      filter(taxa %in% tax)
    
    Supp1_mmi <- Supp1_mmi %>% 
      filter(taxa %in% tax)
    
    # Supp1_OE <- Supp1_OE %>% 
    #   filter(taxa %in% tax)
    # 
    # Supp2_OE <- Supp2_OE %>% 
    #   filter(taxa %in% tax)
    # 
    # null_OE <- null_OE %>% 
    #   filter(taxa %in% tax)
    
  }
  
  ##
  # create asci class output
  out <- asci_1(scores = scr, Supp_mmi = Supp1_mmi, taxa = tax)
  
  return(out)
  
}
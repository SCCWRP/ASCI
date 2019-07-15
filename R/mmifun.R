#' Run the ASCI MMI 
#' 
#' Run the ASCI MMI index for diatoms, soft-bodied algae, and hybrids
#' 
#' @param taxain \code{data.frame} for input taxonomy data
#' 
#' @details 
#' Three index scores are calculated and returned as a named list
#' 
#' @return 
#' A list with three elements named as \code{d.results.scored}, \code{sba.results.scored}, and \code{hybrid.results.scored}. Each element is a \code{data.frame} with metric scores by site.
#'
#' @export
#' 
#' @importFrom reshape2 acast melt
#' @importFrom dplyr bind_rows filter mutate group_by select summarise ungroup
#' @importFrom tidyr gather spread nest unnest separate unite
#' @import purrr
#' @import tibble
#' 
#' @examples 
#' taxain <- getids(demo_algae_tax)
#' mmifun(taxain)
mmifun <- function(taxain){
  
  options(gsubfn.engine = "R")
  
  # Step 1. Import taxonomy data -----------------------------------------------------------
  bugs <- taxain 
  
  # Step 2. Get diatom, sba, hybrid --------------------------------------------------------
  bugs<- merge(bugs, STE, all.x = T) %>% 
    filter(
      SampleTypeCode != "Qualitative",
      !is.na(FinalIDassigned),
      !is.null(FinalIDassigned))
  
  # subset into assemblages # Nov 2 right now no zero scores getting thrown out -- add a flag 
  bugs.d <- bugs %>% 
    filter(
      Phylum == 'Bacillariophyta',
      BAResult != 0
    )
  bugs.sba <- bugs %>% 
    filter(
      Phylum != 'Bacillariophyta',
      BAResult != 0
    )
  bugs <- bugs %>% 
    mutate(ComboResult = as.numeric(pmax(BAResult, BAResult, na.rm = T))) %>% 
    filter(ComboResult != 0)
  
  # Step 3. Convert to species abd matrix at Species level  -----------------------------------------------------------
  bugs.d.m <- as.data.frame(acast(bugs.d, 
                                  SampleID ~ FinalIDassigned, 
                                  value.var = "BAResult", 
                                  fun.aggregate=sum))
  bugs.sba.m <- as.data.frame(acast(bugs.sba, 
                                    SampleID ~ FinalIDassigned, 
                                    value.var = "Result", 
                                    fun.aggregate=sum))
  bugs.hybrid.m <- as.data.frame(acast(bugs, 
                                       SampleID ~ FinalIDassigned, 
                                       value.var = "ComboResult", 
                                       fun.aggregate=sum))
  
  # calculate metrics
  # still need stations in Sussy's code, have to confirm -- mmi_calcmetrics might need it
  d.metrics <- mmi_calcmetrics('diatoms', bugs.d.m, stations)
  sba.metrics <- mmi_calcmetrics('sba', bugs.sba.m, stations)
  hybrid.metrics <- mmi_calcmetrics('hybrid', bugs.hybrid.m, stations)
  
  
  # Load winning metrics -----------------------------------------------------------
  
  d.win <- mmilkup$d.win 
  sba.win <- mmilkup$sba.win
  hybrid.win <- mmilkup$hybrid.win
  
  d.results <- d.metrics %>%
    select(d.win) %>%
    filter()
  
  
  
  
}
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
mmifun <- function(taxain, sitein){
  
  options(gsubfn.engine = "R")
  
  # Step 1. Import taxonomy data -----------------------------------------------------------
  bugs <- taxain 
  
  # Step 2. Get diatom, sba, hybrid --------------------------------------------------------
  bugs <- merge(bugs, STE, all.x = T) %>% 
    filter(
      SampleTypeCode != "Qualitative",
      !is.na(FinalIDassigned),
      !is.null(FinalIDassigned))
  
  # subset into assemblages 
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
                                  fun.aggregate=sum)) %>% 
    rownames_to_column()
  bugs.sba.m <- as.data.frame(acast(bugs.sba, 
                                    SampleID ~ FinalIDassigned, 
                                    value.var = "Result", 
                                    fun.aggregate=sum)) %>% 
    rownames_to_column()
  bugs.hybrid.m <- as.data.frame(acast(bugs, 
                                       SampleID ~ FinalIDassigned, 
                                       value.var = "ComboResult", 
                                       fun.aggregate=sum)) %>% 
    rownames_to_column()
  
  stations <- sitein
  
  # calculate metrics
  d.metrics <- mmi_calcmetrics('diatoms', bugs.d.m, stations)
  sba.metrics <- mmi_calcmetrics('sba', bugs.sba.m, stations)
  hybrid.metrics <- mmi_calcmetrics('hybrid', bugs.hybrid.m, stations)
  
  
  # Load winning metrics -----------------------------------------------------------
  
  d.win <- mmilkup$d.win 
  sba.win <- mmilkup$sba.win
  hybrid.win <- mmilkup$hybrid.win
  
  d.results <- d.metrics %>%
    rownames_to_column() %>% 
    select(d.win) %>%
    filter(rowname %in% bugs.d.m$rowname) %>% 
    column_to_rownames()
    
  sba.results <- sba.metrics %>% 
    rownames_to_column() %>% 
    select(sba.win) %>% 
    filter(rowname %in% bugs.sba.m$rowname) %>% 
    column_to_rownames()
  
  
  hybrid.results <- hybrid.metrics %>% 
    rownames_to_column() %>% 
    select(hybrid.win) %>% 
    filter(rowname %in% bugs.hybrid.m$rowname) %>% 
    column_to_rownames()
  
  d.scored <- score_metric(d.results, bugs.d.m, mmilkup$d.inc, inc = T) %>% 
    left_join(score_metric(bugs.d.m, mmilkup$d.dec, inc = F), by = 'rowname') %>% 
    column_to_rownames() %>% 
    replace(. > 1, 1) %>% 
    replace(. < 0, 0) %>% 
    select(sort(colnames()))
  
  sba.scored <- score_metric(sba.results, bugs.sba.m, mmilkup$sba.inc, inc = T) %>% 
   left_join(score_metric(bugs.sba.m, mmilkup$sba.dec, inc = F), by = 'rowname') %>% 
    column_to_rownames() %>% 
    replace(. > 1, 1) %>% 
    replace(. < 0, 0) %>% 
    select(sort(colnames()))
  
  hybrid.scored <- score_metric(hybrid.results, bugs.hybrid.m, mmilkup$hybrid.inc, inc = T) %>% 
   left_join(score_metric(bugs.hybrid.m, mmilkup$hybrid.dec, inc = F), by = 'rowname') %>% 
    column_to_rownames() %>% 
    replace(. > 1, 1) %>% 
    replace(. < 0, 0) %>% 
    select(sort(colnames()))
  
  d.rf.mean <- mmilkup$d.rf.mean %>% 
    select(colnames(d.scored))
  
  sba.rf.mean <- mmilkup$sba.rf.mean %>% 
    select(colnames(sba.scored))
  
  hybrid.rf.mean <- mmilkup$hybrid.rf.mean %>% 
    select(colnames(hybrid.scored))
  
  d.scored.scaled <- sweep(d.scored, MARGIN = 2, FUN = "/",
                           STATS = colMeans(d.rf.mean, na.rm = T)) 
  sba.scored.scaled <- sweep(sba.scored, MARGIN = 2, FUN = "/",
                           STATS = colMeans(sba.rf.mean, na.rm = T))
  hybrid.scored.scaled <- sweep(hybrid.scored, MARGIN = 2, FUN="/",
                                STATS = colMeans(hybrid.rf.mean, na.rm = T))
  
}
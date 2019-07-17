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
#' sitein <- getids(demo_algae_sitedata)
#' mmifun(taxain, sitein)
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
  browser()
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
  
  omni.ref <- mmilkup$omni.ref
  
  d.scored <- score_metric(taxa = 'diatoms', bugs.d.m, d.results, omni.ref) %>% 
    replace(. > 1, 1) %>% 
    replace(. < 0, 0) %>% 
    select(sort(colnames()))
  
  sba.scored <- score_metric(taxa = 'sba', bugs.sba.m, sba.results, omni.ref) %>% 
    replace(. > 1, 1) %>% 
    replace(. < 0, 0) %>% 
    select(sort(colnames()))
  
  hybrid.scored <- score_metric(taxa = 'hybrid', bugs.hybrid.m, hybrid.results, omni.ref) %>% 
    replace(. > 1, 1) %>% 
    replace(. < 0, 0) %>% 
    select(sort(colnames()))
  
  d.rf.mean <- omni.ref %>%
    filter(Assemblage == 'diatoms',
           Metrics %in% colnames(d.scored))
    
  sba.rf.mean <- omni.ref %>%
    filter(Assemblage == 'sba',
           Metrics %in% colnames(sba.scored))
  
  hybrid.rf.mean <- omni.ref %>%
    filter(Assemblage == 'hybrid',
           Metrics %in% colnames(hybrid.scored))
  
  d.scored.scaled <- d.scored %>% 
    column_to_rownames() %>% 
    sweep(., MARGIN = 2, FUN = "/",
          STATS = colMeans(d.rf.mean, na.rm = T)) 
  
  sba.scored.scaled <- sba.scored %>% 
    column_to_rownames() %>% 
    sweep(., MARGIN = 2, FUN = "/",
          STATS = colMeans(sba.rf.mean, na.rm = T))
  
  hybrid.scored.scaled <- hybrid.scored %>% 
    column_to_rownames() %>% 
    sweep(., MARGIN = 2, FUN="/",
          STATS = colMeans(hybrid.rf.mean, na.rm = T))
  
  
  # put all results in long format
  out <- list(
    diatoms_obs = d.results, 
    diatoms_scr = d.results.scored,
    sba_obs = sba.results,
    sba_scr = sba.results.scored,
    hybrid_obs = hybrid.results, 
    hybrid_scr = hybrid.results.scored
  ) %>% 
    enframe %>% 
    mutate(
      value = map(value, rownames_to_column, 'SampleID'),
      value = map(value, gather, 'met', 'val', -SampleID)
    ) %>% 
    unnest %>% 
    separate(name, c('taxa', 'results'), sep = '_') 
  
  # metric housekeeping
  # still need ceiling?
  out <- out %>% 
    mutate(
      val = ifelse(results == 'scr', pmin(val, 1), val), # ceiling at 1
      val = ifelse(results == 'scr', pmax(val, 0), val) # floor at 0
    )
  
  # get mmi total score
  mmiout <- out %>% 
    filter(results %in% 'scr') %>% 
    group_by(taxa, SampleID) %>% 
    summarise(
      MMI = mean(val)
    ) %>% 
    mutate(MMI_Percentile = pnorm(MMI, mean(MMI), sd(MMI))) %>% 
    ungroup %>% 
    split(.$taxa) %>% 
    map(select, -taxa)
  
  # make out a list
  out <- out %>% 
    unite('met', met, results, sep = '_') %>% 
    mutate(
      met = gsub('_obs$', '', met),
      met = gsub('_scr$', '_score', met)
    ) %>% 
    split(.$taxa) %>% 
    map(select, -taxa) %>% 
    map(spread, met, val)
  browser()
  # list of lists for input to ASCI
  out <- list(
    diatoms = list(mmiout$diatoms, out$diatoms),
    sba = list(mmiout$sba, out$sba),
    hybrid = list(mmiout$hybrid, out$hybrid)
  )
  
  # assign names to list elements
  out <- out %>% 
    map(function(x){
      names(x) <- c('MMI_scores', 'MMI_supp')
      return(x)
    }
    )
  
  
}
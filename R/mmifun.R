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
  bugs <- merge(bugs, STE[,c("FinalID", "FinalIDassigned", "Genus", "Phylum", "Class")], all.x = T) %>% 
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
      Result != 0
    )
  bugs <- bugs %>% 
    mutate(ComboResult = as.numeric(pmax(BAResult, Result, na.rm = T))) %>% 
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
    select(SampleID, colnames(d.win)) %>%
    filter(SampleID %in% rownames(bugs.d.m)) %>% 
    column_to_rownames('SampleID') 
    
  sba.results <- sba.metrics %>% 
    select(SampleID, colnames(sba.win)) %>%
    filter(SampleID %in% rownames(bugs.sba.m)) %>% 
    column_to_rownames('SampleID')  
  
  hybrid.results <- hybrid.metrics %>% 
    select(SampleID, colnames(hybrid.win)) %>%
    filter(SampleID %in% rownames(bugs.hybrid.m)) %>% 
    column_to_rownames('SampleID')
  
  
  omni.ref <- mmilkup$omni.ref
  
  d.scored <- score_metric(taxa = 'diatoms', bugs.d.m, d.results, omni.ref) %>% 
    select(-rowname) %>% 
    replace(. > 1, 1) %>% 
    replace(. < 0, 0) %>% 
    select(sort(colnames(.))) 
  
  sba.scored <- score_metric(taxa = 'sba', bugs.sba.m, sba.results, omni.ref) %>% 
    select(-rowname) %>% 
    replace(. > 1, 1) %>% 
    replace(. < 0, 0) %>% 
    select(sort(colnames(.)))
  
  hybrid.scored <- score_metric(taxa = 'hybrid', bugs.hybrid.m, hybrid.results, omni.ref) %>%
    select(-rowname) %>% 
    replace(. > 1, 1) %>% 
    replace(. < 0, 0) %>% 
    select(sort(colnames(.)))
  
  d.rf.mean <- omni.ref %>%
    filter(Assemblage == 'diatoms',
           Metric %in% colnames(d.scored)) %>% 
    select(RefCalMean)
    
  sba.rf.mean <- omni.ref %>%
    filter(Assemblage == 'sba',
           Metric %in% colnames(sba.scored)) %>% 
    select(RefCalMean)
  
  hybrid.rf.mean <- omni.ref %>%
    filter(Assemblage == 'hybrid',
           Metric %in% colnames(hybrid.scored)) %>% 
    select(RefCalMean)
  
  d.scored.scaled <- d.scored %>% 
    # column_to_rownames() %>% 
    sweep(., MARGIN = 2, FUN = "/",
          STATS = colMeans(d.rf.mean, na.rm = T)) 
  
  sba.scored.scaled <- sba.scored %>% 
    # column_to_rownames() %>% 
    sweep(., MARGIN = 2, FUN = "/",
          STATS = colMeans(sba.rf.mean, na.rm = T))
  
  hybrid.scored.scaled <- hybrid.scored %>% 
    # column_to_rownames() %>% 
    sweep(., MARGIN = 2, FUN="/",
          STATS = colMeans(hybrid.rf.mean, na.rm = T))
  
  
  # put all results in long format
  out <- list(
    diatoms_obs = d.results, 
    diatoms_scr = d.scored.scaled,
    sba_obs = sba.results,
    sba_scr = sba.scored.scaled,
    hybrid_obs = hybrid.results, 
    hybrid_scr = hybrid.scored.scaled
  ) %>% 
    enframe %>% 
    mutate(
      value = map(value, rownames_to_column, 'SampleID'),
      value = map(value, gather, 'met', 'val', -SampleID)
    ) %>% 
    unnest %>% 
    separate(name, c('taxa', 'results'), sep = '_') 
  
  # metric housekeeping
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


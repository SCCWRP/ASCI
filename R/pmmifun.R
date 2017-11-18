#' Run the ASCI pMMI 
#' 
#' Run the ASCI pMMI index for diatoms, soft-bodied algae, and hybrids
#' 
#' @param taxain \code{data.frame} for input taxonomy data
#' @param sitein \code{data.frame} for input site data
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
#' pmmifun(taxain, sitein)
pmmifun <- function(taxain, sitein){

  options(gsubfn.engine = "R")
  
  # Step 1. Import taxonomy data -----------------------------------------------------------
  bugs <- taxain
  
  # Step 2. Import stations data -----------------------------------------------------------
  stations <- sitein
  
  # Step 3. get diatom, sba, hybrid --------------------------------------------------------
  bugs<- merge(bugs, STE, all.x = T)
  bugs.d<-subset(bugs, Class=="Bacillariophyceae")
  bugs.sba<-subset(bugs, Class!="Bacillariophyceae")
  bugs$ComboResult<-as.numeric(pmax(bugs$BAResult,bugs$Result, na.rm=T))
  
  # Step 4. Rarify diatom data -----------------------------------------------------------
  bugs.d.sub<-rarify(inbug=bugs.d, sample.ID="SampleID", abund="BAResult", subsiz=500)
  
  # Step 5. Convert to species abd matrix at Species level  -----------------------------------------------------------
  bugs.d.m<-as.data.frame(acast(bugs.d.sub, SampleID~FinalIDassigned, value.var="BAResult", fun.aggregate=sum))
  bugs.sba.m<-as.data.frame(acast(bugs.sba, SampleID~FinalIDassigned, value.var="Result", fun.aggregate=sum))
  bugs.hybrid.m<-as.data.frame(acast(bugs, SampleID~FinalIDassigned, value.var="ComboResult", fun.aggregate=sum))

  # calculate metrics
  d.metrics<-pmmi_calcmetrics('diatoms', bugs.d.m, stations)
  sba.metrics<-pmmi_calcmetrics('sba', bugs.sba.m, stations)
  hybrid.metrics<-pmmi_calcmetrics('hybrid', bugs.hybrid.m, stations)
  
  # Load winning metrics -----------------------------------------------------------
  d.win <- pmmilkup$d.win 
  sba.win <- pmmilkup$sba.win
  hybrid.win <- pmmilkup$hybrid.win
  d.win<-colnames(d.win[,-(length(names(d.win)))])
  sba.win<-colnames(sba.win[,-(length(names(sba.win)))])
  hybrid.win<-colnames(hybrid.win[,-(length(names(hybrid.win)))])
  
  # subset winning metrics for new sites
  d.results <- data.frame(d.metrics[,d.win], row.names = d.metrics$SampleID)
  sba.results <- data.frame(sba.metrics[,sba.win], row.names = sba.metrics$SampleID)
  hybrid.results <- data.frame(hybrid.metrics[,hybrid.win], row.names = hybrid.metrics$SampleID)
  
  # score results ------------------------------------------------------
  
  # hard-coded values are in pmmilkup$quants
  # for increasers, min is 5th of ref and max is 95th of str
  # for decreasers, min is 5th of str and max is 95th of ref
  
  d.results.scored<-data.frame(row.names = row.names(d.results))
  sba.results.scored<-data.frame(row.names = row.names(sba.results))
  hybrid.results.scored<-data.frame(row.names = row.names(hybrid.results))
  
  # increase (obs - max) / ( min - max)
  d.results.scored$prop.spp.Salinity.BF <- ((d.results$prop.spp.Salinity.BF - 12.9216849) / (-0.12408199 - 12.9216849)) / 1.2768156
  d.results.scored$prop.spp.HighMotility <- ((d.results$prop.spp.HighMotility - 0.4322735) / (0 - 0.4322735)) / 1.0160970
  d.results.scored$prop.ind.most.tol <- ((d.results$prop.ind.most.tol - 0.4798261) / (0 - 0.4798261)) / 1.0169511
  sba.results.scored$cnt.spp.IndicatorClass_TP_high <- ((sba.results$cnt.spp.IndicatorClass_TP_high - 5) / (0 - 5)) / 1.1299149
  sba.results.scored$prop.spp.IndicatorClass_DOC_high <- ((sba.results$prop.spp.IndicatorClass_DOC_high - 0.7736111) / (0.03125000 - 0.7736111)) / 1.1658549
  sba.results.scored$prop.spp.Green <- ((sba.results$prop.spp.Green - 0.7000000) / (0 - 0.7000000)) / 1.1019174
  hybrid.results.scored$prop.spp.IndicatorClass_DOC_high <- ((hybrid.results$prop.spp.IndicatorClass_DOC_high - 0.2500000) / (0.01445135 - 0.2500000)) / 1.5067045
  
  #decrease (obs - min) / (max - min)
  d.results.scored$prop.spp.BCG3 <- ((d.results$prop.spp.BCG3 - 0) / (0.4670833 - 0)) / 0.7808209
  sba.results.scored$prop.spp.BCG3 <- ((sba.results$prop.spp.BCG3 - 0) / (0.5000000 - 0)) / 0.7387326
  hybrid.results.scored$prop.spp.Trophic.I <- ((hybrid.results$prop.spp.Trophic.I - 0) / (0.1764706 - 0)) / 0.8838483
  hybrid.results.scored$prop.spp.ZHR <- ((hybrid.results$prop.spp.ZHR - 0) / (0.1878378 - 0)) / 1.2028715
  hybrid.results.scored$prop.spp.BCG3 <- ((hybrid.results$prop.spp.BCG3 - 0.05579365) / (0.4196944 - 0.05579365)) / 0.7299234
  
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
  out <- out %>% 
    mutate(
      val = ifelse(results == 'scr', pmin(val, 1), val), # ceiling at 1
      val = ifelse(results == 'scr', pmax(val, 0), val) # floor at 0
      )
  
  # get pmmi total score
  pmmiout <- out %>% 
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
    diatoms = list(pmmiout$diatoms, out$diatoms),
    sba = list(pmmiout$sba, out$sba),
    hybrid = list(pmmiout$hybrid, out$hybrid)
  )
  
  # assign names to list elements
  out <- out %>% 
    map(function(x){
      names(x) <- c('MMI_scores', 'MMI_supp')
      return(x)
      }
    )
  
  return(out)
  
}
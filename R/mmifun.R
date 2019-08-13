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
  bugs <- merge(bugs, STE[,c("FinalID", "FinalIDassigned", "Genus", "Phylum", "Class")], all.x = T) %>% 
    filter(
      SampleTypeCode != "Qualitative",
      !is.na(FinalIDassigned), 
      !is.null(FinalIDassigned))
  
  # subset into assemblages 
  
  chkmt <- function(df) {
    if(nrow(df) == 0) {
      df[1,] = rep(-88)
    }
    return(df)
  }
  
  bugs.d <- bugs %>% 
    filter(
      SampleTypeCode == 'Integrated'
    )
  bugs.d <- chkmt(bugs.d)
  
  
  bugs.sba <- bugs %>% 
    filter(
      SampleTypeCode != 'Integrated'
    )
  bugs.sba <- chkmt(bugs.sba)
  
  if(bugs.d$FinalID[1] == -88 | bugs.sba$FinalID[1] == -88) {
    bugs <- bugs[1,]
    bugs[1, ] <- rep(-88)
    bugs <- bugs %>% 
      mutate(ComboResult = as.numeric(pmax(BAResult, Result, na.rm = T)))
  } else {
    smpid <- intersect(bugs.d$SampleID, bugs.sba$SampleID)
    bugs <- bugs %>% 
      filter(SampleID %in% smpid) %>% 
      group_by(SampleID) %>% 
      mutate(ComboResult = as.numeric(pmax(BAResult, Result, na.rm = T))) %>% 
      filter(ComboResult != 0) %>% 
      ungroup()
  }
  
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
  
  stations <- taxain %>% 
    select(StationCode, SampleDate, Replicate, SampleID) %>% 
    unique()

  # calculate metrics
  d.metrics <- mmi_calcmetrics('diatoms', bugs.d.m, stations)
  sba.metrics <- mmi_calcmetrics('sba', bugs.sba.m, stations)
  hybrid.metrics <- mmi_calcmetrics('hybrid', bugs.hybrid.m, stations)
  
  
  # Load winning metrics -----------------------------------------------------------
  d.win <- c('cnt.spp.BCG3', 'prop.Cyclotella', 
             'prop.Surirella', 'prop.spp.OxyReq.DO_10',
             'richness', 'Cyclotella.richness', 
             'Surirella.richness', 'OxyReq.DO_10.richness')
  sba.win <- c('cnt.spp.BCG5', 'cnt.spp.IndicatorClass_Cu_high', 
               'cnt.spp.IndicatorClass_DOC_high', 'cnt.spp.IndicatorClass_TP_high',
               'richness')
  hybrid.win <- c('cnt.ind.most.tol', 'cnt.spp.IndicatorClass_Cu_high',
                  'prop.spp.OrgN.NHHONF', 'prop.Cyclotella',
                  'richness', 'Cyclotella.richness', 
                  'OrgN.NHHONF.richness')
  
  d.results <- d.metrics %>%
    select(SampleID, d.win) %>%
    filter(SampleID %in% rownames(bugs.d.m)) %>%
    mutate(
      pcnt.attributed.BCG3 = cnt.spp.BCG3/richness,
      pcnt.attributed.Cyclotella = Cyclotella.richness/richness,
      pcnt.attributed.Surirella = Surirella.richness/richness,
      pcnt.attributed.OxyReg.DO_10 = OxyReq.DO_10.richness/richness
    ) %>% 
    select(-c('Cyclotella.richness', 
              'Surirella.richness', 'OxyReq.DO_10.richness')) %>% 
    rename(NumberTaxa = richness) %>% 
    column_to_rownames('SampleID')
  d.results <- chkmt(d.results)
    
  sba.results <- sba.metrics %>% 
    select(SampleID, sba.win) %>%
    filter(SampleID %in% rownames(bugs.sba.m)) %>% 
    mutate(
      pcnt.attributed.BCG5 = cnt.spp.BCG5/richness,
      pcnt.attributed.HiCu = cnt.spp.IndicatorClass_Cu_high/richness,
      pcnt.attributed.HiDOC = cnt.spp.IndicatorClass_DOC_high/richness,
      pcnt.attributed.HiTP.DO_10 = cnt.spp.IndicatorClass_TP_high/richness
    ) %>% 
    rename(NumberTaxa = richness) %>% 
    column_to_rownames('SampleID')
  sba.results <- chkmt(sba.results)
  
  hybrid.results <- hybrid.metrics %>% 
    select(SampleID, hybrid.win) %>%
    filter(SampleID %in% rownames(bugs.hybrid.m)) %>% 
    mutate(
      pcnt.attributed.HiTolerance = cnt.ind.most.tol/richness,
      pcnt.attributed.Cyclotella = Cyclotella.richness/richness,
      pcnt.attributed.HiCu = cnt.spp.IndicatorClass_Cu_high/richness,
      pcnt.attributed.NHHONF = OrgN.NHHONF.richness/richness
    ) %>% 
    select(-c('Cyclotella.richness', 
              'OrgN.NHHONF.richness')) %>% 
    rename(NumberTaxa = richness) %>% 
    column_to_rownames('SampleID')
  hybrid.results <- chkmt(hybrid.results)
  
  
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
    arrange(Metric) %>% 
    column_to_rownames('Metric') %>% 
    select(RefCalMean) %>% 
    t()
    
  sba.rf.mean <- omni.ref %>%
    filter(Assemblage == 'sba',
           Metric %in% colnames(sba.scored)) %>% 
    arrange(Metric) %>% 
    column_to_rownames('Metric') %>% 
    select(RefCalMean) %>% 
    t()
  
  hybrid.rf.mean <- omni.ref %>%
    filter(Assemblage == 'hybrid',
           Metric %in% colnames(hybrid.scored)) %>% 
    arrange(Metric) %>% 
    column_to_rownames('Metric') %>% 
    select(RefCalMean) %>% 
    t()

  d.scored.scaled <- d.scored %>% 
    # column_to_rownames() %>% 
    sweep(., MARGIN = 2, FUN = "/",
          STATS = colMeans(d.rf.mean, na.rm = T)) 
  if(bugs.d[1,1] == -88){
    d.scored.scaled <- d.scored.scaled[1,] 
    d.scored.scaled[1,] <- rep(-88)
  }
  
  sba.scored.scaled <- sba.scored %>% 
    sweep(., MARGIN = 2, FUN = "/",
          STATS = colMeans(sba.rf.mean, na.rm = T))
  if(bugs.sba[1,1] == -88){
    sba.scored.scaled <- sba.scored.scaled[1,]
    sba.scored.scaled[1,] <- rep(-88)
  }

  hybrid.scored.scaled <- hybrid.scored %>% 
    sweep(., MARGIN = 2, FUN="/",
          STATS = colMeans(hybrid.rf.mean, na.rm = T))
  if(bugs[1,1] == -88){
    hybrid.scored.scaled <- hybrid.scored.scaled[1, ]
    hybrid.scored.scaled[1, ] <- rep(-88)
  }
  
  
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
  
  
  # get mmi total score
  mmiout <- out %>% 
    filter(results %in% 'scr') %>% 
    group_by(taxa, SampleID) %>% 
    summarise(
      ASCI = mean(val)
    ) %>% 
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


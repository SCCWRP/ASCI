#' Run the ASCI MMI 
#' 
#' Run the ASCI MMI index for diatoms, soft-bodied algae, and hybrids
#' 
#' @param taxa \code{data.frame} for input taxonomy data
#' @param station \code{data.frame} for input station data
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
#' @importFrom dplyr bind_rows contains filter mutate group_by rename_at select select_at summarise vars ungroup
#' @importFrom tidyr gather spread nest unnest separate unite
#' @import purrr
#' @import tibble
#' 
#' @examples 
#' # check input taxonomy data
#' dat <- chkinp(demo_algae_tax, demo_station)
#' dat <- dat$taxa
#' 
#' # calculate GIS from stations
#' station <- calcgis(demo_station)
#' 
#' # calc metrics
#' out <- mmifun(dat, station)
#' out
mmifun <- function(taxa, station){
  
  options(gsubfn.engine = "R")
  
  # Import taxonomy data -----------------------------------------------------------
  bugs <- taxa 
  
  # Get diatom, sba, hybrid --------------------------------------------------------
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
  
  # diatoms only in integrated
  bugs.d <- bugs %>% 
    filter(
      SampleTypeCode == 'Integrated'
    )
  bugs.d <- chkmt(bugs.d)
  
  # soft-bodied not in integrated
  bugs.sba <- bugs %>% 
    filter(
      SampleTypeCode != 'Integrated'
    )
  bugs.sba <- chkmt(bugs.sba)
  
  # create hybrid, but first see if both exist, if not create dummy data frame
  if(bugs.d$FinalID[1] == -88 | bugs.sba$FinalID[1] == -88) {
    bugs <- bugs[1,]
    bugs[1, ] <- rep(-88)
    bugs <- bugs %>% 
      mutate(ComboResult = as.numeric(pmax(BAResult, Result, na.rm = T)))
  } else { # otherwise subset both
    smpid <- intersect(bugs.d$SampleID, bugs.sba$SampleID)
    bugs <- bugs %>% 
      filter(SampleID %in% smpid) %>% 
      group_by(SampleID) %>% 
      mutate(ComboResult = as.numeric(pmax(BAResult, Result, na.rm = T))) %>% 
      filter(ComboResult != 0) %>% 
      ungroup()
  }
  
  # Convert to species abd matrix at Species level  -----------------------------------------------------------
  bugs.d.m <- as.data.frame(acast(bugs.d, 
                                  SampleID ~ FinalIDassigned, 
                                  value.var = "BAResult", 
                                  fun.aggregate=sum, na.rm = T))
  bugs.sba.m <- as.data.frame(acast(bugs.sba, 
                                    SampleID ~ FinalIDassigned, 
                                    value.var = "Result", 
                                    fun.aggregate=sum, na.rm = T))
  bugs.hybrid.m <- as.data.frame(acast(bugs, 
                                       SampleID ~ FinalIDassigned, 
                                       value.var = "ComboResult", 
                                       fun.aggregate=sum, na.rm = T))
  
  stations <- taxa %>% 
    select(StationCode, SampleDate, Replicate, SampleID) %>% 
    unique()
  
  # calculate metrics
  d.metrics <- mmi_calcmetrics('diatoms', bugs.d.m, stations)
  sba.metrics <- mmi_calcmetrics('sba', bugs.sba.m, stations)
  hybrid.metrics <- mmi_calcmetrics('hybrid', bugs.hybrid.m, stations)
  
  # Setup GIS predictors by station id --------------------------------------
  
  stationid <- taxa %>% 
    select(SampleID, StationCode) %>% 
    unique %>% 
    full_join(station, by ='StationCode')
  
  # Load winning metrics -----------------------------------------------------------
  
  # metric names for mmi_calcmetrics
  # d.win <- c('prop.spp.SPIspecies4', 'Salinity.BF.richness', 
  #            'prop.spp.Saprobic.BM', 'cnt.spp.IndicatorClass_TP_low',
  #            'richness', 'cnt.spp.SPIspecies4', 'Saprobic.BM.richness')
  # sba.win <- c('prop.spp.Green', 'cnt.spp.IndicatorClass_DOC_high', 
  #              'prop.spp.BCG45', 'richness', 'cnt.spp.Green', 'cnt.spp.BCG45')
  # hybrid.win <- c('prop.spp.BCG4', 'Salinity.BF.richness', 
  #                 'prop.spp.IndicatorClass_DOC_high', 'OxyRed.DO_30.richness', 
  #                 'richness', 'cnt.spp.BCG4', 'cnt.spp.IndicatorClass_DOC_high')
  
  d.win<-c("prop.spp.BCG12", "prop.spp.OxyReq.DO_100orDO_75", "prop.spp.Salinity.BF", "prop.spp.Trophic.E", "richness")
  sba.win<-c("cnt.spp.IndicatorClass_DOC_high" , "prop.spp.BCG45", "prop.spp.Green", "richness")
  hybrid.win<-c("OxyRed.DO_30.richness","prop.spp.BCG4","prop.spp.IndicatorClass_DOC_high", "Salinity.BF.richness", "richness")
  
  
  # names with suffix
  # d.win.suf <- c('prop.spp.SPIspecies4_mod', 'Salinity.BF.richness_mod', 
  #                'prop.spp.Saprobic.BM_raw', 'cnt.spp.IndicatorClass_TP_low_raw')
  # sba.win.suf <- c('prop.spp.Green_raw', 'cnt.spp.IndicatorClass_DOC_high_raw', 
  #                  'prop.spp.BCG45_raw')
  # hybrid.win.suf <- c('prop.spp.BCG4_mod', 'Salinity.BF.richness_mod', 
  #                     'prop.spp.IndicatorClass_DOC_high_raw', 'OxyRed.DO_30.richness_mod')
  
  d.win.suf<-c("prop.spp.BCG12_mod", "prop.spp.OxyReq.DO_100orDO_75_raw", "prop.spp.Salinity.BF_mod", "prop.spp.Trophic.E_mod")
  sba.win.suf<-c("cnt.spp.IndicatorClass_DOC_high_raw" , "prop.spp.BCG45_raw", "prop.spp.Green_raw")
  hybrid.win.suf<-c("OxyRed.DO_30.richness_mod","prop.spp.BCG4_mod","prop.spp.IndicatorClass_DOC_high_raw", "Salinity.BF.richness_mod")

  
  # Calculated observed and predicted metrics -------------------------------
  
  ##
  # diatoms
  
  # get observe diatom metrics and percent attributed
  d.results <- d.metrics %>%
    select(SampleID, d.win) %>%
    filter(SampleID %in% rownames(bugs.d.m)) %>%
    mutate(
      pcnt.attributed.prop.spp.BCG12 = prop.spp.BCG12/richness,
      pcnt.attributed.prop.spp.OxyReq.DO_100orDO_75 = prop.spp.OxyReq.DO_100orDO_75/richness,
      pcnt.attributed.prop.spp.Salinity.BF = prop.spp.Salinity.BF/richness,
      pcnt.attributed.prop.spp.Trophic.E = prop.spp.Trophic.E/richness
    ) # %>% 
   #  select(-c('prop.spp.BCG12','prop.spp.Salinity.BF','prop.spp.Trophic.E_mod')) # mystery line 
  names(d.results) <- paste0(names(d.results), '_raw') 
  d.results <- d.results %>% 
    rename(
      NumberTaxa = richness_raw, 
      SampleID = SampleID_raw
    )
  
  # predicted diatom metrics
  d.predmet <- stationid %>% 
    mutate(
      prop.spp.BCG12_pred = predict(rfmods$diatoms.prop.spp.BCG12, newdata = .[, c("CondQR50","MAX_ELEV")]), 
      prop.spp.Salinity.BF_pred = predict(rfmods$diatoms.prop.spp.Salinity.BF, newdata = .[, c("XerMtn","KFCT_AVE","CondQR50","LST32AVE","AtmCa","SITE_ELEV")]), 
      prop.spp.Trophic.E_pred = predict(rfmods$diatoms.prop.spp.Trophic.E, newdata = .[, c("SITE_ELEV", "KFCT_AVE")]) 
    ) %>% 
    select(SampleID, prop.spp.SPI.species4_pred, Salinity.BF.richness_pred, prop.spp.Trophic.E_pred)
  
  # join with observed, take residuals for raw/pred metrics
  d.results <- d.results %>% 
    left_join(d.predmet, by = 'SampleID') %>%
    mutate(
      prop.spp.BCG12_mod = prop.spp.BCG12_raw - prop.spp.BCG12_pred, 
      prop.spp.Salinity.BF_mod = prop.spp.Salinity.BF_raw - prop.spp.Salinity.BF_pred, 
      prop.spp.Trophic.E_mod = prop.spp.Trophic.E_raw - prop.spp.Trophic.E_pred
    ) %>% 
    column_to_rownames('SampleID')
  
  # final selection
  colsel <- names(d.results) %in% d.win.suf | grepl('^NumberTaxa|^pcnt\\.attributed', names(d.results))
  d.results <- d.results[, colsel] %>% 
    rename_at(vars(contains('pcnt.attributed')), function(x) gsub('\\_raw$', '', x))
  d.results <- chkmt(d.results)
  
  ##
  # sba, no predicted metrics
  
  # get observed soft-bodied metrics and percent attributed  
  sba.results <- sba.metrics %>% 
    select(SampleID, sba.win) %>%
    filter(SampleID %in% rownames(bugs.sba.m)) %>% 
    mutate(
      pcnt.attributed.cnt.spp.IndicatorClass_DOC_high = cnt.spp.IndicatorClass_DOC_high/richness,
      pcnt.attributed.prop.spp.BCG45 = prop.spp.BCG45/richness,
      pcnt.attributed.prop.spp.Green = prop.spp.Green/richness
    ) # %>% 
   #  select(-c('foo'))  # mystery line 
  names(sba.results) <- paste0(names(sba.results), '_raw') 
  sba.results <- sba.results %>% 
    rename(
      NumberTaxa = richness_raw, 
      SampleID = SampleID_raw
    ) %>% 
    column_to_rownames('SampleID')
  
  # final selection
  colsel <- names(sba.results) %in% sba.win.suf | grepl('^NumberTaxa|^pcnt\\.attributed', names(sba.results))
  sba.results <- sba.results[, colsel] %>% 
    rename_at(vars(contains('pcnt.attributed')), function(x) gsub('\\_raw$', '', x))
  sba.results <- chkmt(sba.results)
  
  ##
  # hybrid
  
  #   hybrid.win.suf<-c("OxyRed.DO_30.richness_mod","prop.spp.BCG4_mod","prop.spp.IndicatorClass_DOC_high_raw", "Salinity.BF.richness_mod")
  
  # get observed hybrid results and percent attributed
  hybrid.results <- hybrid.metrics %>% 
    select(SampleID, hybrid.win) %>%
    filter(SampleID %in% rownames(bugs.hybrid.m)) %>% 
    mutate(
      pcnt.attributed.OxyRed.DO_30.richness = OxyRed.DO_30.richness/richness,
      pcnt.attributed.prop.spp.BCG4 = prop.spp.BCG4/richness,
      pcnt.attributed.prop.spp.IndicatorClass_DOC_high = prop.spp.IndicatorClass_DOC_high/richness,
      pcnt.attributed.Salinity.BF.richness = Salinity.BF.richness/richness
    ) # %>% 
   # select(-c('OxyRed.DO_30.richness', 'prop.spp.BCG4', 'Salinity.BF.richness')) # mystery line 
  names(hybrid.results) <- paste0(names(hybrid.results), '_raw') 
  hybrid.results <- hybrid.results %>% 
    rename(
      NumberTaxa = richness_raw, 
      SampleID = SampleID_raw
    )
  
  # predicted hybrid metrics
  hybrid.predmet <- stationid %>% 
    mutate(
      OxyRed.DO_30.richness_pred = predict(rfmods$hybrid.OxyRed.DO_30.richness, newdata = .[, c("AtmCa","PPT_00_09")]), 
      prop.spp.BCG4_pred = predict(rfmods$hybrid.prop.spp.BCG4, newdata = .[, c("MAX_ELEV","CondQR50")]), 
      Salinity.BF.richness_pred = predict(rfmods$hybrid.Salinity.BF.richness, newdata = .[, c("XerMtn", "KFCT_AVE")]) 
    ) %>% 
    select(SampleID, OxyRed.DO_30.richness_pred, prop.spp.BCG4_pred, Salinity.BF.richness_pred)
  
  # join with observed, take residuals for raw/pred metrics
  hybrid.results <- hybrid.results %>% 
    left_join(hybrid.predmet, by = 'SampleID') %>%
    mutate(
      OxyRed.DO_30.richness_mod = OxyRed.DO_30.richness_raw - OxyRed.DO_30.richness_pred, 
      prop.spp.BCG4_mod = prop.spp.BCG4_raw - prop.spp.BCG4_pred, 
      Salinity.BF.richness_mod = Salinity.BF.richness_raw - Salinity.BF.richness_pred
    ) %>% 
    column_to_rownames('SampleID') %>% 
    select_at(vars(-contains('pred')))
  
  # final selection
  colsel <- names(hybrid.results) %in% hybrid.win.suf | grepl('^NumberTaxa|^pcnt\\.attributed', names(hybrid.results))
  hybrid.results <- hybrid.results[, colsel] %>% 
    rename_at(vars(contains('pcnt.attributed')), function(x) gsub('\\_raw$', '', x))
  hybrid.results <- chkmt(hybrid.results)
  
  # Score metrics -----------------------------------------------------------
  
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
    unnest(cols = value) %>% 
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

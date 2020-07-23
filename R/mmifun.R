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
#' @importFrom vegan diversity specnumber 
#'
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
      df[1,] = rep(NA)
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
  
  # ASCI should always calculate hybrid even if they are missing assemblage 
  smpid <- bugs$SampleID
  bugs <- bugs %>% 
    filter(SampleID %in% smpid) %>% 
    group_by(SampleID) %>% 
    mutate(ComboResult = as.numeric(pmax(BAResult, Result, na.rm = T))) %>% 
    filter(ComboResult != 0) %>% 
    ungroup()
  
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
  d.mmi <- mmi_calcmetrics('diatoms', bugs.d.m, stations)
  sba.mmi <- mmi_calcmetrics('sba', bugs.sba.m, stations)
  hybrid.mmi <- mmi_calcmetrics('hybrid', bugs.hybrid.m, stations)
  
  # recent change to allow access to the stations combined variable that gets generated in the mmi_calmetrics function
  # it returns a list now
  # 7/1/2020
  d.metrics <- d.mmi$metrics
  d.stations_combined <- d.mmi$stations_combined
  
  sba.metrics <- sba.mmi$metrics
  sba.stations_combined <- sba.mmi$stations_combined
  
  hybrid.metrics <- hybrid.mmi$metrics
  hybrid.stations_combined <- hybrid.mmi$stations_combined
  
  
  
  # Setup GIS predictors by station id --------------------------------------
  
  stationid <- taxa %>% 
    select(SampleID, StationCode) %>% 
    unique %>% 
    full_join(station, by ='StationCode')
  
  # Load winning metrics -----------------------------------------------------------
  
  # metric names for mmi_calcmetrics
  
  d.win<-c("cnt.spp.most.tol",
           "EpiRho.richness",
           "prop.spp.IndicatorClass_TN_low",
           "prop.spp.Planktonic",
           "prop.spp.Trophic.E",
           "Salinity.BF.richness",
           "shannon",
           "simpson",
           "richness")
  sba.win<-c("prop.spp.IndicatorClass_DOC_high",
             "prop.spp.IndicatorClass_NonRef",
             "prop.spp.IndicatorClass_TP_high",
             "prop.spp.ZHR",
             "shannon",
             "simpson",
             "richness")
  hybrid.win<-c("cnt.spp.IndicatorClass_TP_high",
                "cnt.spp.most.tol",
                "EpiRho.richness",
                "OxyRed.DO_30.richness",
                "prop.spp.Planktonic",
                "prop.spp.Trophic.E",
                "prop.spp.ZHR",
                "Salinity.BF.richness",
                "shannon",
                "simpson",
                "richness")
  
  # names with suffix mod or raw
  
  d.win.suf<-c("cnt.spp.most.tol_mod",
               "EpiRho.richness_mod",
               "prop.spp.IndicatorClass_TN_low_mod",
               "prop.spp.Planktonic_mod",
               "prop.spp.Trophic.E_mod",
               "Salinity.BF.richness_mod")
  sba.win.suf<-c("prop.spp.IndicatorClass_DOC_high_raw",
                 "prop.spp.IndicatorClass_NonRef_raw",
                 "prop.spp.IndicatorClass_TP_high_raw",
                 "prop.spp.ZHR_raw")
  hybrid.win.suf<-c("cnt.spp.IndicatorClass_TP_high_mod",
                    "cnt.spp.most.tol_mod",
                    "EpiRho.richness_mod",
                    "OxyRed.DO_30.richness_mod",
                    "prop.spp.Planktonic_mod",
                    "prop.spp.Trophic.E_mod",
                    "prop.spp.ZHR_raw",
                    "Salinity.BF.richness_mod")
  
  # Calculated observed and predicted metrics -------------------------------
  ##
  # diatoms
  
  # get observe diatom metrics and percent attributed
  d.results <- d.metrics %>%
    select(SampleID, d.win) %>%
    filter(SampleID %in% rownames(bugs.d.m)) %>%
    mutate(
      pcnt.attributed.cnt.spp.most.tol =  cnt.spp.most.tol/100,
      pcnt.attributed.EpiRho.richness =  EpiRho.richness/100,
      pcnt.attributed.prop.spp.IndicatorClass_TN_low =  prop.spp.IndicatorClass_TN_low/100,
      pcnt.attributed.prop.spp.Planktonic =  prop.spp.Planktonic/100,
      pcnt.attributed.prop.spp.Trophic.E =  prop.spp.Trophic.E/100,
      pcnt.attributed.Salinity.BF.richness =  Salinity.BF.richness/100
      
    )  # %>% 
  
  names(d.results) <- paste0(names(d.results), '_raw') 
  
  d.results <- d.results %>% 
    dplyr::rename(
      NumberTaxa = richness_raw,
      SampleID = SampleID_raw
    )
  
  # predicted diatom metrics
  d.predmet <- stationid %>% 
    mutate(
      cnt.spp.most.tol = predict(rfmods$diatoms.cnt.spp.most.tol, newdata = .[,c("XerMtn","PPT_00_09")]),
      EpiRho.richness = predict(rfmods$diatoms.EpiRho.richness, newdata = .[,c("AREA_SQKM","TMAX_WS")]),
      prop.spp.IndicatorClass_TN_low = predict(rfmods$diatoms.prop.spp.IndicatorClass_TN_low,newdata = .[,c("CondQR50","MAX_ELEV")]),
      prop.spp.Planktonic = predict(rfmods$diatoms.prop.spp.Planktonic,newdata = .[,c("CondQR50","SITE_ELEV")]),
      prop.spp.Trophic.E = predict(rfmods$diatoms.prop.spp.Trophic.E,newdata = .[,c("KFCT_AVE","CondQR50")]),
      Salinity.BF.richness = predict(rfmods$diatoms.Salinity.BF.richness,newdata = .[,c("XerMtn","KFCT_AVE","CondQR50")])
    ) %>% 
    select(SampleID, 
           cnt.spp.most.tol, EpiRho.richness, prop.spp.IndicatorClass_TN_low, 
           prop.spp.Planktonic, prop.spp.Trophic.E, Salinity.BF.richness) 
  
  # join with observed, take residuals for raw/pred metrics
  d.results <- d.results %>% 
    left_join(d.predmet, by = 'SampleID') %>%
    mutate(
      cnt.spp.most.tol_mod = cnt.spp.most.tol_raw - cnt.spp.most.tol,
      EpiRho.richness_mod = EpiRho.richness_raw - EpiRho.richness,
      prop.spp.IndicatorClass_TN_low_mod = prop.spp.IndicatorClass_TN_low_raw - prop.spp.IndicatorClass_TN_low,
      prop.spp.Planktonic_mod = prop.spp.Planktonic_raw - prop.spp.Planktonic,
      prop.spp.Trophic.E_mod = prop.spp.Trophic.E_raw - prop.spp.Trophic.E,
      Salinity.BF.richness_mod = Salinity.BF.richness_raw - Salinity.BF.richness
    ) %>% 
    column_to_rownames('SampleID')
  
  # final selection
  colsel <- names(d.results) %in% d.win.suf | grepl('^NumberTaxa|^pcnt\\.attributed', names(d.results)) | grepl('pred',
                                                                                                                names(d.results))
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
      pcnt.attributed.prop.spp.IndicatorClass_DOC_high = prop.spp.IndicatorClass_DOC_high/100,
      pcnt.attributed.prop.spp.IndicatorClass_NonRef = prop.spp.IndicatorClass_NonRef/100,
      pcnt.attributed.prop.spp.IndicatorClass_TP_high = prop.spp.IndicatorClass_TP_high/100, 
      pcnt.attributed.prop.spp.ZHR_raw = prop.spp.ZHR/100 # prop.spp.ZHR_raw wasn't a column
    )  # %>% 
  # select(-c('foo'))  # mystery line 
  
  names(sba.results) <- paste0(names(sba.results), '_raw') 
  
  
  sba.results <- sba.results %>% 
    rename(
      NumberTaxa = richness_raw, 
      #NumberTaxa = 10, 
      SampleID = SampleID_raw
    ) %>% 
    column_to_rownames('SampleID')
  
  # final selection
  colsel <- names(sba.results) %in% sba.win.suf | grepl('^NumberTaxa|^pcnt\\.attributed', names(sba.results)) | grepl('pred',
                                                                                                                      names(sba.results))
  
  sba.results <- sba.results[, colsel] %>% 
    rename_at(vars(contains('pcnt.attributed')), function(x) gsub('\\_raw$', '', x))
  sba.results <- chkmt(sba.results)
  
  ##
  # hybrid
  
  hybrid.win.suf<-c("cnt.spp.IndicatorClass_TP_high_mod",
                    "cnt.spp.most.tol_mod",
                    "EpiRho.richness_mod",
                    "OxyRed.DO_30.richness_mod",
                    "prop.spp.Planktonic_mod",
                    "prop.spp.Trophic.E_mod",
                    "prop.spp.ZHR_raw",
                    "Salinity.BF.richness_mod")
  
  # get observed hybrid results and percent attributed
  hybrid.results <- hybrid.metrics %>% 
    select(SampleID, hybrid.win) %>%
    filter(SampleID %in% rownames(bugs.hybrid.m)) %>% 
    mutate(
      pcnt.attributed.cnt.spp.IndicatorClass_TP_high = cnt.spp.IndicatorClass_TP_high/100,
      pcnt.attributed.cnt.spp.most.tol = cnt.spp.most.tol/100,
      pcnt.attributed.EpiRho.richness = EpiRho.richness/100,
      pcnt.attributed.OxyRed.DO_30.richness = OxyRed.DO_30.richness/100,
      pcnt.attributed.prop.spp.Planktonic = prop.spp.Planktonic/100,
      pcnt.attributed.prop.spp.Trophic.E = prop.spp.Trophic.E/100,
      pcnt.attributed.prop.spp.ZHR_raw = prop.spp.ZHR/100, # prop.spp.ZHR_raw wasn't a column
      pcnt.attributed.Salinity.BF.richness = Salinity.BF.richness/100
      
    ) #  %>% 
  # select(-c('OxyRed.DO_30.richness', 'prop.spp.BCG4', 'Salinity.BF.richness')) # mystery line 
  names(hybrid.results) <- paste0(names(hybrid.results), '_raw') 
  hybrid.results <- hybrid.results %>% 
    rename(
      NumberTaxa = richness_raw,
      #NumberTaxa = 10,
      SampleID = SampleID_raw
    )
  
  # predicted hybrid metrics
  
  hybrid.predmet <- stationid %>% 
    mutate(
      # hybrid.cnt.spp.IndicatorClass_TP_high is supposed to be a randomForest model object thing
      # However, right now it is saying that it is NULL........
      cnt.spp.IndicatorClass_TP_high = predict(rfmods$hybrid.cnt.spp.IndicatorClass_TP_high, newdata = .[, c("PPT_00_09", "KFCT_AVE")]), 
      cnt.spp.most.tol = predict(rfmods$hybrid.cnt.spp.most.tol, newdata = .[, c("CondQR50", "XerMtn")]), 
      EpiRho.richness = predict(rfmods$hybrid.EpiRho.richness, newdata = .[, c("AREA_SQKM", "TMAX_WS")]), 
      OxyRed.DO_30.richness = predict(rfmods$hybrid.OxyRed.DO_30.richness, newdata = .[, c("AtmCa", "PPT_00_09")]), 
      prop.spp.Planktonic = predict(rfmods$hybrid.prop.spp.Planktonic, newdata = .[, c("CondQR50", "SITE_ELEV")]), 
      prop.spp.Trophic.E = predict(rfmods$hybrid.prop.spp.Trophic.E, newdata = .[, c("CondQR50", "KFCT_AVE")]), 
      Salinity.BF.richness = predict(rfmods$hybrid.Salinity.BF.richness, newdata = .[, c("XerMtn", "KFCT_AVE")]) 
      
    ) %>% 
    select(SampleID, cnt.spp.IndicatorClass_TP_high,cnt.spp.most.tol,EpiRho.richness, 
           OxyRed.DO_30.richness,prop.spp.Planktonic,prop.spp.Trophic.E,Salinity.BF.richness )
  
  # join with observed, take residuals for raw/pred metrics
  hybrid.results <- hybrid.results %>% 
    left_join(hybrid.predmet, by = 'SampleID') %>%
    mutate(
      cnt.spp.IndicatorClass_TP_high_mod = cnt.spp.IndicatorClass_TP_high_raw - cnt.spp.IndicatorClass_TP_high, 
      cnt.spp.most.tol_mod = cnt.spp.most.tol_raw - cnt.spp.most.tol, 
      EpiRho.richness_mod = EpiRho.richness_raw - EpiRho.richness, 
      OxyRed.DO_30.richness_mod = OxyRed.DO_30.richness_raw - OxyRed.DO_30.richness, 
      prop.spp.Planktonic_mod = prop.spp.Planktonic_raw - prop.spp.Planktonic, 
      prop.spp.Trophic.E_mod = prop.spp.Trophic.E_raw - prop.spp.Trophic.E, 
      Salinity.BF.richness_mod = Salinity.BF.richness_raw - Salinity.BF.richness
      
    ) %>% 
    column_to_rownames('SampleID') # %>% 
  # select_at(vars(-contains('pred')))
  
  # final selection
  colsel <- names(hybrid.results) %in% hybrid.win.suf | grepl('^NumberTaxa|^pcnt\\.attributed', names(hybrid.results)) | grepl('pred', names(hybrid.results))
  hybrid.results <- hybrid.results[, colsel] %>% 
    rename_at(vars(contains('pcnt.attributed')), function(x) gsub('\\_raw$', '', x))
  hybrid.results <- chkmt(hybrid.results)
  
  # Score metrics [make from 0 to 1] -----------------------------------------------------------
  
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
  
  # Find average. Scale index by dividing by refcalmean ----------------------------------

  # average 
  d.scored$MMI <- rowMeans(d.scored, na.rm = T)
  sba.scored$MMI <- rowMeans(sba.scored, na.rm = T)
  hybrid.scored$MMI <- rowMeans(hybrid.scored, na.rm = T) 
  
  # refcalmeans
  d.rf.mean <- 0.752705813
  sba.rf.mean <- 0.70994176
  hybrid.rf.mean <- 0.73063118
  
  # ASCI has to be MMI divided by Ref Cal Mean
  d.scored$ASCI <- d.scored$MMI / d.rf.mean
  sba.scored$ASCI <- sba.scored$MMI / sba.rf.mean
  hybrid.scored$ASCI <- hybrid.scored$MMI / hybrid.rf.mean
  
  d.scored.scaled<-d.scored
  sba.scored.scaled<-sba.scored
  hybrid.scored.scaled<-hybrid.scored
  
  # These were problematic
  # d.predmet and hybrid.predmet
  # They needed to have SampleID to be the rownames
  d.predmet <- d.predmet %>% column_to_rownames("SampleID")
  hybrid.predmet <- hybrid.predmet %>% column_to_rownames("SampleID")

  # put all results in long format
 
  out <- list(
    diatoms_obs = d.results, 
    diatoms_pred = d.predmet,
    diatoms_scr = d.scored.scaled %>% select(-ASCI),
    sba_obs = sba.results,
    sba_scr = sba.scored.scaled %>% select(-ASCI),
    hybrid_obs = hybrid.results, 
    hybrid_pred = hybrid.predmet,
    hybrid_scr = hybrid.scored.scaled %>% select(-ASCI)
  ) %>% 
    enframe %>% 
    mutate(
      value = map(value, rownames_to_column, 'SampleID'),
      value = map(value, gather, 'met', 'val', -SampleID)
    ) %>% 
    unnest(cols = value) %>% 
    separate(name, c('taxa', 'results'), sep = '_') 
  
  # get mmi total score
  # new way of calculating asci, replaced by lines 382-384
  # mmiout <- out %>% 
  #   filter(results %in% 'scr') %>% 
  #   group_by(taxa, SampleID) %>% 
  #   summarise(
  #     ASCI = mean(val)
  #   ) %>% 
  #   ungroup %>% 
  #   split(.$taxa) %>% 
  #   map(select, -taxa)


  # Robert 7/23/2020
  # Looks like this is where the ASCI scores were getting duplicated
  # Notice that above, sba.scored.scaled gets added to out
  # it also got the ASCI column of it added here
  # all columns from "out" get those suffixes put on them, namely "_scr"
  # I prevented ASCI from getting included up there
  mmiout <- list(
    diatoms = d.scored.scaled %>% 
      rownames_to_column('SampleID') %>% 
      select(SampleID, ASCI),
    
    hybrid = hybrid.scored.scaled %>% 
      rownames_to_column('SampleID') %>% 
      select(SampleID, ASCI),
    
    sba = sba.scored.scaled %>% 
      rownames_to_column('SampleID') %>% 
      select(SampleID, ASCI)
  )
  
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
  
  
  # out <- list(
  #   diatoms = list(out$diatoms),
  #   sba = list(out$sba),
  #   hybrid = list(out$hybrid)
  # )
  
  
  # assign names to list elements
  out <- out %>% 
    map(function(x){
      names(x) <- c('MMI_scores', 'MMI_supp')
      return(x)
    }
    )
}

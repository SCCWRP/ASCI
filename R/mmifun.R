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
mmifun <- function(dat, station){
  
  options(gsubfn.engine = "R")
  
  # Import taxonomy data -----------------------------------------------------------
  algae <- dat
  
  # Get diatom, sba, hybrid --------------------------------------------------------
  algae <- merge(
      algae, 
      STE[,c("FinalID", "FinalIDassigned", "Genus", "Phylum", "Class")], 
      all.x = T
    ) %>% 
    filter(
      SampleTypeCode != "Qualitative",
      !is.na(FinalIDassigned),
      !is.null(FinalIDassigned)
    )
  
  # subset into assemblages 
  
  # check empty, check if the dataframe is empty or not
  # if it is, make a row of NA's
  chkmt <- function(df) {
    if(nrow(df) == 0) {
      df[1,] = rep(NA)
    }
    return(df)
  }
  
  # diatoms only in integrated
  algae.d <- algae %>% 
    filter(
      SampleTypeCode == 'Integrated'
    )
  algae.d <- chkmt(algae.d)
  
  # soft-bodied not in integrated
  algae.sba <- algae %>% 
    filter(
      SampleTypeCode != 'Integrated'
    )
  algae.sba <- chkmt(algae.sba)
  
  # ASCI should always calculate hybrid even if they are missing assemblage 
  algae <- algae %>%  
    group_by(SampleID) %>% 
    mutate(
      ComboResult = as.numeric(pmax(BAResult, Result, na.rm = T))
    ) %>% 
    filter(ComboResult != 0) %>% 
    ungroup()
  
  # Convert to species abd matrix at Species level  -----------------------------------------------------------
  algae.d.m <- as.data.frame(
    acast(
      algae.d, 
      SampleID ~ FinalIDassigned, 
      value.var = "BAResult", 
      fun.aggregate=sum, na.rm = T)
    )
  algae.sba.m <- as.data.frame(
    acast(
      algae.sba, 
      SampleID ~ FinalIDassigned, 
      value.var = "Result", 
      fun.aggregate=sum, na.rm = T)
    )
  algae.hybrid.m <- as.data.frame(
    acast(
      algae, 
      SampleID ~ FinalIDassigned, 
      value.var = "ComboResult", 
      fun.aggregate=sum, na.rm = T)
    )
  
  
  # calculate raw metrics
  stations_for_calcmetrics <- dat %>% 
    select(StationCode, SampleDate, Replicate, SampleID) %>% 
    unique()
  
  d.metrics <- mmi_calcmetrics('diatoms', algae.d.m, stations_for_calcmetrics)
  sba.metrics <- mmi_calcmetrics('sba', algae.sba.m, stations_for_calcmetrics)
  hybrid.metrics <- mmi_calcmetrics('hybrid', algae.hybrid.m, stations_for_calcmetrics)

  
  # Setup GIS predictors by station id --------------------------------------
  gis_predictors <- dat %>% 
    select(SampleID, StationCode) %>% 
    unique %>% 
    full_join(station, by ='StationCode')
  
  # Load winning metrics -----------------------------------------------------------
  
  # metric names for mmi_calcmetrics
  d.win<-c(
    "cnt.spp.most.tol",
    "EpiRho.richness",
    "prop.spp.IndicatorClass_TN_low",
    "prop.spp.Planktonic",
    "prop.spp.Trophic.E",
    "Salinity.BF.richness",
    "shannon",
    "simpson",
    "richness"
  )
  sba.win<-c(
    "prop.spp.IndicatorClass_DOC_high",
    "prop.spp.IndicatorClass_NonRef",
    "prop.spp.IndicatorClass_TP_high",
    "prop.spp.ZHR",
    "shannon",
    "simpson",
    "richness"
  )
  hybrid.win<-c(
    "cnt.spp.IndicatorClass_TP_high",
    "cnt.spp.most.tol",
    "EpiRho.richness",
    "OxyRed.DO_30.richness",
    "prop.spp.Planktonic",
    "prop.spp.Trophic.E",
    "prop.spp.ZHR",
    "Salinity.BF.richness",
    "shannon",
    "simpson",
    "richness"
  )
  
  # Calculated observed and predicted metrics -------------------------------
  ##
  # diatoms
  
  # get observe diatom metrics and percent attributed
  d.results <- d.metrics %>%
    select(SampleID, all_of(d.win)) %>%
    filter(SampleID %in% rownames(algae.d.m)) %>%
    # stick the raw suffix on all the column names, since that is what they are
    rename_at(
      vars(-c(SampleID, richness)),
      function(x) {
        paste0(x,"_raw")
      }
    ) %>%
    # percent attributed is simply raw / 100
    mutate(
      cnt.spp.most.tol_pct_att = cnt.spp.most.tol_raw/100,
      EpiRho.richness_pct_att = EpiRho.richness_raw/100,
      prop.spp.IndicatorClass_TN_low_pct_att = prop.spp.IndicatorClass_TN_low_raw/100,
      prop.spp.Planktonic_pct_att =  prop.spp.Planktonic_raw/100,
      prop.spp.Trophic.E_pct_att =  prop.spp.Trophic.E_raw/100,
      Salinity.BF.richness_pct_att =  Salinity.BF.richness_raw/100
      
    ) %>% 
    rename(
      NumberTaxa = richness
    )
  
  # predicted diatom metrics
  d.predmet <- gis_predictors %>% 
    mutate(
      cnt.spp.most.tol_pred = predict(rfmods$diatoms.cnt.spp.most.tol, newdata = .[,c("XerMtn","PPT_00_09")]),
      EpiRho.richness_pred = predict(rfmods$diatoms.EpiRho.richness, newdata = .[,c("AREA_SQKM","TMAX_WS")]),
      prop.spp.IndicatorClass_TN_low_pred = predict(rfmods$diatoms.prop.spp.IndicatorClass_TN_low,newdata = .[,c("CondQR50","MAX_ELEV")]),
      prop.spp.Planktonic_pred = predict(rfmods$diatoms.prop.spp.Planktonic,newdata = .[,c("CondQR50","SITE_ELEV")]),
      prop.spp.Trophic.E_pred = predict(rfmods$diatoms.prop.spp.Trophic.E,newdata = .[,c("KFCT_AVE","CondQR50")]),
      Salinity.BF.richness_pred = predict(rfmods$diatoms.Salinity.BF.richness,newdata = .[,c("XerMtn","KFCT_AVE","CondQR50")])
    ) %>% 
    select(SampleID, 
           cnt.spp.most.tol_pred, EpiRho.richness_pred, prop.spp.IndicatorClass_TN_low_pred, 
           prop.spp.Planktonic_pred, prop.spp.Trophic.E_pred, Salinity.BF.richness_pred) 
  
  # join with observed, take residuals for raw/pred metrics
  d.results <- d.results %>% 
    left_join(d.predmet, by = 'SampleID') %>%
    column_to_rownames('SampleID')
  
  # Later code depends on rownames being the sampleID
  # above code doesn't like the rownames to be sampleID's
  d.predmet <- d.predmet %>% column_to_rownames("SampleID")
  
  d.results <- chkmt(d.results)
  
  ##
  # sba, no predicted metrics
  
  # get observed soft-bodied metrics and percent attributed  
  
  sba.results <- sba.metrics %>% 
    select(SampleID, all_of(sba.win)) %>%
    filter(SampleID %in% rownames(algae.sba.m)) %>% 
    # stick the raw suffix on all the column names, since that is what they are
    rename_at(
      vars(-c(SampleID, richness)),
      function(x) {
        paste0(x,"_raw")
      }
    ) %>%
    # pct_att (percent attributed is simply raw over 100)
    mutate(
      prop.spp.IndicatorClass_DOC_high_pct_att = prop.spp.IndicatorClass_DOC_high_raw/100,
      prop.spp.IndicatorClass_NonRef_pct_att = prop.spp.IndicatorClass_NonRef_raw/100,
      prop.spp.IndicatorClass_TP_high_pct_att = prop.spp.IndicatorClass_TP_high_raw/100, 
      prop.spp.ZHR_raw_pct_att = prop.spp.ZHR_raw/100 # prop.spp.ZHR_raw wasn't a column
    ) %>% 
    rename(
      NumberTaxa = richness
    ) %>% 
    column_to_rownames('SampleID')
  
  sba.results <- chkmt(sba.results)
  
  ##
  
  
  # get observed hybrid results and percent attributed
  hybrid.results <- hybrid.metrics %>% 
    select(SampleID, all_of(hybrid.win)) %>%
    filter(SampleID %in% rownames(algae.hybrid.m)) %>% 
    # stick the raw suffix on all the column names, since that is what they are
    rename_at(
      vars(-c(SampleID, richness)),
      function(x) {
        paste0(x,"_raw")
      }
    ) %>%
    # pct_att (percent attributed is simply raw over 100)
    mutate(
      cnt.spp.IndicatorClass_TP_high_pct_att = cnt.spp.IndicatorClass_TP_high_raw/100,
      cnt.spp.most.tol_pct_att = cnt.spp.most.tol_raw/100,
      EpiRho.richness_pct_att = EpiRho.richness_raw/100,
      OxyRed.DO_30.richness_pct_att = OxyRed.DO_30.richness_raw/100,
      prop.spp.Planktonic_pct_att = prop.spp.Planktonic_raw/100,
      prop.spp.Trophic.E_pct_att = prop.spp.Trophic.E_raw/100,
      prop.spp.ZHR_raw_pct_att = prop.spp.ZHR_raw/100,
      Salinity.BF.richness_pct_att = Salinity.BF.richness_raw/100
      
    ) %>% 
    rename(
      NumberTaxa = richness
    )
  
  # predicted hybrid metrics
  
  hybrid.predmet <- gis_predictors %>% 
    mutate(
      cnt.spp.IndicatorClass_TP_high_pred = predict(rfmods$hybrid.cnt.spp.IndicatorClass_TP_high, newdata = .[, c("PPT_00_09", "KFCT_AVE")]), 
      cnt.spp.most.tol_pred = predict(rfmods$hybrid.cnt.spp.most.tol, newdata = .[, c("CondQR50", "XerMtn")]), 
      EpiRho.richness_pred = predict(rfmods$hybrid.EpiRho.richness, newdata = .[, c("AREA_SQKM", "TMAX_WS")]), 
      OxyRed.DO_30.richness_pred = predict(rfmods$hybrid.OxyRed.DO_30.richness, newdata = .[, c("AtmCa", "PPT_00_09")]), 
      prop.spp.Planktonic_pred = predict(rfmods$hybrid.prop.spp.Planktonic, newdata = .[, c("CondQR50", "SITE_ELEV")]), 
      prop.spp.Trophic.E_pred = predict(rfmods$hybrid.prop.spp.Trophic.E, newdata = .[, c("CondQR50", "KFCT_AVE")]), 
      Salinity.BF.richness_pred = predict(rfmods$hybrid.Salinity.BF.richness, newdata = .[, c("XerMtn", "KFCT_AVE")]) 
    ) %>% 
    select(
      SampleID, cnt.spp.IndicatorClass_TP_high_pred, cnt.spp.most.tol_pred, EpiRho.richness_pred, 
      OxyRed.DO_30.richness_pred, prop.spp.Planktonic_pred, prop.spp.Trophic.E_pred, Salinity.BF.richness_pred
    )
  
  # join with observed, take residuals for raw/pred metrics
  hybrid.results <- hybrid.results %>% 
    left_join(hybrid.predmet, by = 'SampleID') %>%
    column_to_rownames('SampleID') 
  
  # Later code depends on rownames being the sampleID
  # above code doesn't like the rownames to be sampleID's
  hybrid.predmet <- hybrid.predmet %>% column_to_rownames("SampleID")

  hybrid.results <- chkmt(hybrid.results)
  
  
  # Score metrics [make from 0 to 1] -----------------------------------------------------------
  omni.ref <- mmilkup$omni.ref
  d.scored <- score_metric(taxa = 'diatoms', algae.d.m, d.results, omni.ref) %>% 
    select(-rowname) %>% 
    replace(. > 1, 1) %>% 
    replace(. < 0, 0) %>% 
    select(sort(colnames(.)))
  names(d.scored) <- gsub("_raw","_scr",names(d.scored))

  
  sba.scored <- score_metric(taxa = 'sba', algae.sba.m, sba.results, omni.ref) %>% 
    select(-rowname) %>% 
    replace(. > 1, 1) %>%
    replace(. < 0, 0) %>% 
    select(sort(colnames(.)))
  names(sba.scored) <- gsub("_raw","_scr",names(sba.scored))
  
  
  hybrid.scored <- score_metric(taxa = 'hybrid', algae.hybrid.m, hybrid.results, omni.ref) %>%
    select(-rowname) %>% 
    replace(. > 1, 1) %>% 
    replace(. < 0, 0) %>% 
    select(sort(colnames(.)))
  names(hybrid.scored) <- gsub("_raw","_scr",names(hybrid.scored))
  
  # Find average. Scale index by dividing by refcalmean ----------------------------------


  # Here we calculate ASCI
  # Take the mean of each "score metric" and divide by the refcal mean
  # refcalmeans
  d.rf.mean <- 0.752705813
  sba.rf.mean <- 0.70994176
  hybrid.rf.mean <- 0.73063118
  
  # average 
  d.scored$ASCI <- rowMeans(d.scored, na.rm = T) / d.rf.mean
  sba.scored$ASCI <- rowMeans(sba.scored, na.rm = T) / sba.rf.mean
  hybrid.scored$ASCI <- rowMeans(hybrid.scored, na.rm = T) / hybrid.rf.mean


  # put all results in long format
 
  out <- list(
    diatoms_obs = d.results, 
    #diatoms_pred = d.predmet,
    diatoms_scr = d.scored %>% select(-ASCI),
    sba_obs = sba.results,
    sba_scr = sba.scored %>% select(-ASCI),
    hybrid_obs = hybrid.results, 
    #hybrid_pred = hybrid.predmet,
    hybrid_scr = hybrid.scored %>% select(-ASCI)
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
    diatoms = d.scored %>% 
      rownames_to_column('SampleID') %>% 
      select(SampleID, ASCI),
    
    hybrid = hybrid.scored %>% 
      rownames_to_column('SampleID') %>% 
      select(SampleID, ASCI),
    
    sba = sba.scored %>% 
      rownames_to_column('SampleID') %>% 
      select(SampleID, ASCI)
  )
  

  out <- out %>% 
    # i think the suffixes are already in place
    # unite('met', met, results, sep = '_') %>% 
    # mutate(
    #   met = gsub('_obs$', '', met),
    #   met = gsub('_scr$', '_score', met)
    # ) %>% 
    split(.$taxa) %>% 
    map(select, -c(taxa,results)) %>% 
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

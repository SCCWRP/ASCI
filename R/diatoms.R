#' Calculate ASCI metrics for diatoms
#' 
#' @param algae \code{data.frame} input taxonomy data, with traits and indicators attached.
#' @param gismetrics \code{data.frame} the dataframe with the gis predictors 
#' 
#' @return a dataframe in wide format with the ASCI metrics for diatoms
#' 
#' @import dplyr
#' 
#' @export
diatoms <- function(algae, gismetrics) {
  
  # First things first - Keep only integrated sampletypecodes
  # The scores appear to come out closer to the expected results when we DONT include this step
  algae <- algae %>%
    filter(
      SampleTypeCode == 'Integrated',
      Phylum == 'Bacillariophyta'
    ) %>%
    select(
      -c(SampleTypeCode, BAResult, Result)
    ) %>%
    group_by(
      SampleID
    ) %>%
    distinct(
      FinalIDassigned, .keep_all = TRUE
    ) %>%
    ungroup()

  # Now we get the metrics
  d.metrics <- algae %>%
    # groupby includes StationCode so that the column is preserved, to later join with gismetrics
    group_by(StationCode, SampleID) %>%
    summarize(
      # ---- Get the raw and percent attributed metrics ----
      # Raw metric is, well, the actual calculated metric value
      # Percent attributed is the percentage of species that actually had data available in the column used 
      #   to calculate that metric
      
      # Count of most tolerant species
      cnt.spp.most.tol_raw = sum(na.omit(designation=='Stressed')),
      cnt.spp.most.tol_pct_att = length(na.omit(designation)) / length(designation),
      
      # Count of genus Epithemia and Rhopalodia
      EpiRho.richness_raw = sum(na.omit(Genus=='Epithemia'|Genus=='Rhopalodia')),
      EpiRho.richness_pct_att = length(na.omit(Genus)) / length(Genus),
      
      # Proportion of species with a value of "low" in the IndicatorClass_TN column
      prop.spp.IndicatorClass_TN_low_raw = sum(na.omit(IndicatorClass_TN=='low'))/length(na.omit(IndicatorClass_TN)),
      prop.spp.IndicatorClass_TN_low_pct_att = length(na.omit(IndicatorClass_TN)) / length(IndicatorClass_TN),
      
      # Proportion of Species with the Habitat = P
      prop.spp.Planktonic_raw = sum(na.omit(Habitat=='P'))/length(na.omit(Habitat)),
      prop.spp.Planktonic_pct_att = length(na.omit(Habitat)) / length(Habitat),
      
      # Proportion of species with the TrophicState value equal to "E"
      prop.spp.Trophic.E_raw = sum(na.omit(TrophicState=='E'))/length(na.omit(TrophicState)),
      prop.spp.Trophic.E_pct_att = length(na.omit(TrophicState)) / length(TrophicState),
      
      # How many species have a Salinity value of "BF"
      Salinity.BF.richness_raw = sum(na.omit(Salinity=='BF')),
      Salinity.BF.richness_pct_att = length(na.omit(Salinity)) / length(Salinity)
    ) %>%
    ungroup() %>%
    # join with the gismetrics dataframe by the StationCode column
    # gismetrics dataframe has no SampleID column
    full_join(gismetrics, by = 'StationCode') %>%
    group_by(StationCode, SampleID) %>%
    # Next we get the Predicted values and the Mod values (Raw - Pred)
    mutate(
      # ---- Get Predicted Metrics ----
      # Count of most tolerant species
      # GIS Predictors: XerMtn and PPT_00_09
      cnt.spp.most.tol_pred = ifelse(
        is.na(cnt.spp.most.tol_raw),
        NA_real_,
        predict(
          rfmods$diatoms.cnt.spp.most.tol, c(XerMtn,PPT_00_09)
        )
      ),
      cnt.spp.most.tol_mod = cnt.spp.most.tol_raw - cnt.spp.most.tol_pred,
      
      # Count of genus Epithemia and Rhopalodia
      # GIS Predictors: AREA_SQKM and TMAX_WS
      EpiRho.richness_pred = ifelse(
        is.na(EpiRho.richness_raw),
        NA_real_,
        predict(
          rfmods$diatoms.EpiRho.richness, c(AREA_SQKM, TMAX_WS)
        )
      ),
      EpiRho.richness_mod = EpiRho.richness_raw - EpiRho.richness_pred,
      
      # Proportion of species with a value of "low" in the IndicatorClass_TN column
      # GIS Predictors: CondQR50 and MAX_ELEV
      prop.spp.IndicatorClass_TN_low_pred = ifelse(
        is.na(prop.spp.IndicatorClass_TN_low_raw),
        NA_real_,
        predict(
          rfmods$diatoms.prop.spp.IndicatorClass_TN_low, c(CondQR50, MAX_ELEV)
        )
      ),
      prop.spp.IndicatorClass_TN_low_mod = prop.spp.IndicatorClass_TN_low_raw - prop.spp.IndicatorClass_TN_low_pred,
      
      # Proportion of Species with the Habitat = P
      # GIS Predictors: CondQR50 and SITE_ELEV
      prop.spp.Planktonic_pred = ifelse(
        is.na(prop.spp.Planktonic_raw),
        NA_real_,
        predict(
          rfmods$diatoms.prop.spp.Planktonic, c(CondQR50, SITE_ELEV)
        )
      ),
      prop.spp.Planktonic_mod = prop.spp.Planktonic_raw - prop.spp.Planktonic_pred,
      
      # Proportion of species with the TrophicState value equal to "E"
      # GIS Predictors: KFCT_AVE and CondQR50
      prop.spp.Trophic.E_pred = ifelse(
        is.na(prop.spp.Trophic.E_raw),
        NA_real_,
        predict(
          rfmods$diatoms.prop.spp.Trophic.E, c(KFCT_AVE,CondQR50)
        )
      ),
      prop.spp.Trophic.E_mod = prop.spp.Trophic.E_raw - prop.spp.Trophic.E_pred,
      
      # How many species have a Salinity value of "BF"
      # GIS Predictors: XerMtn, KFCT_AVE and CondQR50
      Salinity.BF.richness_pred = ifelse(
        is.na(Salinity.BF.richness_raw),
        NA_real_,
        predict(
          rfmods$diatoms.Salinity.BF.richness, c(XerMtn,KFCT_AVE,CondQR50)
        )
      ),
      Salinity.BF.richness_mod = Salinity.BF.richness_raw - Salinity.BF.richness_pred
    ) %>%
    ungroup()  %>%
    select(-c(names(gismetrics), StationCode))
  

  # if rows have all NA values, drop the rows.
  d.metrics <- d.metrics[
    rowSums(is.na(d.metrics)) != ncol(d.metrics)
  ,]

    
  # Last but not least we return the metrics
  return(d.metrics)
    
}

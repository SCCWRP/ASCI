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
hybrid <- function(algae, gismetrics) {
  hybrid.metrics <- algae %>%
    select(
      -c(SampleTypeCode, BAResult, Result)
    ) %>%
    distinct(
      FinalIDAssigned, .keep_all = TRUE
    ) %>%
    # groupby includes StationCode so that the column is preserved, to later join with gismetrics
    group_by(StationCode, SampleID) %>%
    summarize(
      # ---- Raw and Pct Att ----
      # Get the raw and percent attributed metrics
      # Raw metric is, well, the actual calculated metric value
      # Percent attributed is the percentage of species that actually had data available in the column used 
      #   to calculate that metric
      
      # Count of species with a value of "high" in the IndicatorClass_TP column
      cnt.spp.IndicatorClass_TP_high_raw = sum(na.omit(IndicatorClass_TP=='high')),
      cnt.spp.IndicatorClass_TP_high_pct_att = length(na.omit(IndicatorClass_TP)) / length(IndicatorClass_TP),
      
      # Count of most tolerant species
      cnt.spp.most.tol_raw = sum(na.omit(designation=='Stressed')),
      cnt.spp.most.tol_pct_att = length(na.omit(designation)) / length(designation),
      
      # Count of genus Epithemia and Rhopalodia
      EpiRho.richness_raw = sum(na.omit(Genus=='Epithemia'|Genus=='Rhopalodia')),
      EpiRho.richness_pct_att = length(na.omit(Genus)) / length(Genus),
          
      # How many species have an OxygenRequirement of DO_30?
      # Should the Metric be called "OxyReq"? Not OxyRed...?
      OxyRed.DO_30.richness_raw = sum(na.omit(OxygenRequirements=='DO_30')),
      OxyRed.DO_30.richness_pct_att = length(na.omit(OxygenRequirements)) / length(OxygenRequirements),

      # Proportion of Species with the Habitat = P
      prop.spp.Planktonic_raw = sum(na.omit(Habitat=='P'))/length(na.omit(Habitat)),
      prop.spp.Planktonic_pct_att = length(na.omit(Habitat)) / length(Habitat),
      
      # Proportion of species with the TrophicState value equal to "E"
      prop.spp.Trophic.E_raw = sum(na.omit(TrophicState=='E'))/length(na.omit(TrophicState)),
      prop.spp.Trophic.E_pct_att = length(na.omit(TrophicState)) / length(TrophicState),
            
      # Proportion of species with a ZHR value of "yes"
      prop.spp.ZHR_raw = sum(na.omit(ZHR=='yes'))/length(na.omit(ZHR)),
      prop.spp.ZHR_pct_att = length(na.omit(ZHR)) / length(ZHR),

      # How many species have a Salinity value of "BF"
      Salinity.BF.richness_raw = sum(na.omit(Salinity=='BF')),
      Salinity.BF.richness_pct_att = length(na.omit(Salinity)) / length(Salinity)
    ) %>%
    ungroup() %>%
    # join with the gismetrics dataframe by the StationCode column
    # gismetrics dataframe has no SampleID column
    inner_join(gismetrics, by = 'StationCode') %>%
    group_by(SampleID) %>%
    # Next we get the Predicted values and the Mod values (Raw - Pred)
    mutate(
      # ---- Get Predicted Metrics ----
      # Count of species with a value of "high" in the IndicatorClass_TN column
      # GIS Predictors: "PPT_00_09", "KFCT_AVE"
      cnt.spp.IndicatorClass_TP_high_pred = predict(
        rfmods$hybrid.cnt.spp.IndicatorClass_TP_high, c(PPT_00_09, KFCT_AVE)
      ),
      cnt.spp.IndicatorClass_TP_high_mod = cnt.spp.IndicatorClass_TP_high_raw - cnt.spp.IndicatorClass_TP_high_pred,
      
      # Count of most tolerant species
      # GIS Predictors: "CondQR50", "XerMtn"
      cnt.spp.most.tol_pred = predict(
        rfmods$hybrid.cnt.spp.most.tol, c(CondQR50, XerMtn)
      ),
      cnt.spp.most.tol_mod = cnt.spp.most.tol_raw - cnt.spp.most.tol_pred,
      
      # Count of genus Epithemia and Rhopalodia
      # GIS Predictors: "AREA_SQKM", "TMAX_WS"
      EpiRho.richness_pred = predict(
        rfmods$hybrid.EpiRho.richness, c(AREA_SQKM, TMAX_WS)
      ),
      EpiRho.richness_mod = EpiRho.richness_raw - EpiRho.richness_pred,
      
      # Count of speices with OxygenRequirements DO_30
      # GIS Predictors: "AtmCa", "PPT_00_09"
      OxyRed.DO_30.richness_pred = predict(
        rfmods$hybrid.OxyRed.DO_30.richness, c(AtmCa, PPT_00_09)
      ),
      OxyRed.DO_30.richness_mod = OxyRed.DO_30.richness_raw - OxyRed.DO_30.richness_pred,

      # Proportion of Species with the Habitat = P
      # GIS Predictors: "CondQR50", "SITE_ELEV"
      prop.spp.Planktonic_pred = predict(
        rfmods$hybrid.prop.spp.Planktonic, c(CondQR50, SITE_ELEV)
      ),
      prop.spp.Planktonic_mod = prop.spp.Planktonic_raw - prop.spp.Planktonic_pred,
      
      # Proportion of species with the TrophicState value equal to "E"
      # GIS Predictors: "CondQR50", "KFCT_AVE"
      prop.spp.Trophic.E_pred = predict(
        rfmods$hybrid.prop.spp.Trophic.E, c(CondQR50, KFCT_AVE)
      ),
      prop.spp.Trophic.E_mod = prop.spp.Trophic.E_raw - prop.spp.Trophic.E_pred,
      
      # How many species have a Salinity value of "BF"
      # GIS Predictors: "XerMtn", "KFCT_AVE"
      Salinity.BF.richness_pred = predict(
        rfmods$hybrid.Salinity.BF.richness, c(XerMtn, KFCT_AVE)
      ),
      Salinity.BF.richness_mod = Salinity.BF.richness_raw - Salinity.BF.richness_pred
    ) %>%
    # we don't want those GIS Metric names as columns in our final output
    select(-names(gismetrics)) %>%
    ungroup()
  
  
  # Last but not least we return the metrics
  return(hybrid.metrics)
  
}
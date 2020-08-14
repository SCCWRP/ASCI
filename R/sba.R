#' Calculate ASCI metrics for soft body algae
#' 
#' @param algae \code{data.frame} input taxonomy data, with traits and indicators attached
#' 
#' @return a dataframe in wide format with the ASCI metrics for diatoms
#' 
#' @import dplyr
#' 
#' @export
sba <- function(algae) {  
  # First things first - Keep only non integrated sampletypecodes
  algae <- algae %>%
    filter(
      SampleTypeCode != 'Integrated',
      Phylum != 'Bacillariophyta'
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
  # No predictive metrics for soft body algae
  # no GIS data needed.
  # Makes life easy here
  sba.metrics <- algae %>%
    group_by(SampleID) %>%
    summarize(
      # Get the raw and percent attributed metrics
      # Raw metric is, well, the actual calculated metric value
      # Percent attributed is the percentage of species that actually had data available in the column used 
      #   to calculate that metric
      
      # Proportion of species with a value of "high" in the IndicatorClass_DOC column
      prop.spp.IndicatorClass_DOC_high_raw = sum(na.omit(IndicatorClass_DOC=='high'))/length(na.omit(IndicatorClass_DOC)),
      prop.spp.IndicatorClass_DOC_high_pct_att = length(na.omit(IndicatorClass_DOC)) / length(IndicatorClass_DOC),
      
      # Proportion of species with a value of "NRF" in the IndicatorClass_Ref column
      prop.spp.IndicatorClass_NonRef_raw = sum(na.omit(IndicatorClass_Ref=='NRF'))/length(na.omit(IndicatorClass_Ref)),
      prop.spp.IndicatorClass_NonRef_pct_att = length(na.omit(IndicatorClass_Ref)) / length(IndicatorClass_Ref),
      
      # Proportion of species with a value of "high" in the IndicatorClass_TP column
      prop.spp.IndicatorClass_TP_high_raw = sum(na.omit(IndicatorClass_TP=='high'))/length(na.omit(IndicatorClass_TP)),
      prop.spp.IndicatorClass_TP_high_pct_att = length(na.omit(IndicatorClass_TP)) / length(IndicatorClass_TP),
      
      # Proportion of Species with ZHR == 'yes'
      prop.spp.ZHR_raw = sum(na.omit(ZHR=='yes'))/length(na.omit(ZHR)),
      prop.spp.ZHR_pct_att = length(na.omit(ZHR)) / length(ZHR),
      
    ) %>%
    ungroup() 
  
  
  # Last but not least we return the metrics
  return(sba.metrics)
  
}
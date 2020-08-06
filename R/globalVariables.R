globalVariables(
  c(
    'STE', 'Class', 'acast', 'mmilkup', 'SampleTypeCode', 'Phylum', 'Result', 'ComboResult',
    'melt', 'value', 'specnumber', 'ddply', 'summarize', 'stations', 'diversity',
    'SampleID', 'met', 'name', 'pnorm', 
    'results', 'sd', 'taxa', 'val', '.', 'SampleID_raw', 'mets', 'rfmods', 
    'MMI', 'MMI_Percentile',
    'Replicate', 'SampleDate', 'StationCode', 'BAResult', 'RefCalMean',
    'diaind', 'variable', 'FinalIDassigned', 'StationID', 'ave', 'cls', 'data', 'fst', 'grp', 
    'ind', 'lm', 'mmi', 'prc_amg', 'prc_wth', 'psa', 'res_tst', 'rfmod', 'scr.x', 
    'scr.y', 'typ', 'metest', 'Assemblage', 'rowname', 'Metric', 'FinalID', 'Met', 'Value',
    'grouped_id', 'scores', 'D_ASCI', 'D_NumberTaxa',
    'D_ValveCount', 'H_ASCI', 'H_NumberTaxa', 'S_ASCI', 'S_Biovolume', 'S_EntityCount',
    'S_NumberTaxa', 'SampleType', 'UnrecognizedTaxa', 
    
    # traits 
    'NitrogenUptakeMetabolism',
    'OxygenRequirements', 
    'Salinity', 
    'Saprobity', 
    'TrophicState', 
    'Genus', 
    'Motility', 
    'designation', 
    'Salinity2', 
    'ZHR', 
    'CRUS', 
    'Green', 'BCG', 'BCG12', 'BCG45', 'CRUS', 'IndicatorClass_TP', 'IndicatorClass_Cu', 
    'IndicatorClass_DOC', 'IndicatorClass_Ref', 'IndicatorClass_TN', 'Heterocystous', 
    'NitrogenUptakeMetabolism2','Habitat', 
    
    # metrics 
    "cnt.spp.most.tol", "cnt.spp.most.tol_raw", "cnt.spp.most.tol_pred",
    "EpiRho.richness", "EpiRho.richness_raw", "EpiRho.richness_pred",
    "prop.spp.IndicatorClass_TN_low", "prop.spp.IndicatorClass_TN_low_raw", "prop.spp.IndicatorClass_TN_low_pred",
    "Salinity.BF.richness", "Salinity.BF.richness_raw","Salinity.BF.richness_pred",
    "cnt.spp.IndicatorClass_TP_high", "cnt.spp.IndicatorClass_TP_high_raw", "cnt.spp.IndicatorClass_TP_high_pred",
    "OxyRed.DO_30.richness", "OxyRed.DO_30.richness_raw","OxyRed.DO_30.richness_pred",
    "prop.spp.Planktonic", "prop.spp.Planktonic_raw","prop.spp.Planktonic_pred",
    "prop.spp.Trophic.E", "prop.spp.Trophic.E_raw", "prop.spp.Trophic.E_pred",
    "prop.spp.ZHR_raw",
    "Salinity.BF.richness", "Salinity.BF.richness_raw", "Salinity.BF.richness_pred",
    "prop.spp.IndicatorClass_DOC_high_raw",
    "prop.spp.IndicatorClass_NonRef_raw",
    "prop.spp.IndicatorClass_TP_high_raw",
    "prop.spp.ZHR_raw", 
    "richness",

    # GIS Predictors
    "AREA_SQKM","AtmCa","CondQR50","SITE_ELEV","TMAX_WS","XerMtn",
    "KFCT_AVE","PPT_00_09","MAX_ELEV",

    # Other various column names that show up in strings of dplyr pipes
    # and Travis wants to see them show up in this file
    "Score", "Metric", "Metric_Type", "Used for ASCI"

  )
)

#' @importFrom stats ave lm na.omit pnorm predict runif sd t.test var
NULL

#' @importFrom utils data tail
NULL

#' @importFrom methods .S4methods new
NULL

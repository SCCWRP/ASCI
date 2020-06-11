library(tidyverse)

#Many fields require manual filling. These are left as blank ("") for manual entry in Excel.
swampify_ASCI<-function(x){
  
  non_analytename_cols <- c('StationCode', 'SampleID', 'SampleDate', 
  'ProjectCode', 'EventCode', 'ProtocolCode', 'AgencyCode', 'SampleComments', 
  'LocationCode', 'GeometryShape', 'CollectionTime', 'CollectionMethodCode', 
  'Replicate', 'HabitatCollectionComments', 'MatrixName', 'MethodName', 'AnalyteName',
  'FractionName', 'UnitName', 'VariableResult', 'Result', 'ResQualCode', 
  'QACode', 'ComplianceCode', 'BatchVerificationCode', 'CollectionDeviceName', 
  'HabitatResultComments', 'Comments', 'version_number')
  
  
  #Core, Suppl1_mmi, and Suppl1_grps
  xdf<- x %>%
    select(-SampleType, -UnrecognizedTaxa) %>%
    # names() %>%dput()
    pivot_longer(
      cols = names(x)[which(!(names(x) %in% c(non_analytename_cols, 'SampleType', 'UnrecognizedTaxa') ))],
      names_to = "AnalyteName",
      values_to = "Result"
    ) %>%
    mutate(
      AnalyteName = case_when(
        AnalyteName == "D_ASCI"~"ASCI_D",
        AnalyteName == "S_ASCI"~"ASCI_S",
        AnalyteName == "H_ASCI"~"ASCI_H",
        T~paste0("ASCI_", AnalyteName)
      ),
      SampleDate=SampleDate,
      ProjectCode="",
      EventCode="BA",
      ProtocolCode="", #Or always "SWAMP_2016_WS"?
      AgencyCode="",
      SampleComments="",
      LocationCode="X",
      GeometryShape="Point",
      CollectionTime="",
      CollectionMethodCode="",
      Replicate=Replicate,
      HabitatCollectionComments="",
      MatrixName="benthic",
      MethodName=paste0(
        "ASCI_software_v",
        strsplit(packageVersion("ASCI") %>% as.character() ,split=".")[[1]][1],".x"),
      FractionName="None",
      UnitName="none",
      VariableResult="",
      ResQualCode="=",
      QACode="None",
      ComplianceCode="Pend",
      BatchVerificationCode="NR",
      CollectionDeviceName="", 
      HabitatResultComments=Comments,
           
    ) %>%
  mutate(
    UnitName = case_when(
      grepl("ASCI_D|ASCI_S|ASCI_H|_scr", AnalyteName) ~ "score",
      grepl('pct|pcnt|prp|prop', AnalyteName) ~ "%",
      grepl('count|cnt|richness', tolower(AnalyteName)) ~ "count",
      AnalyteName %in% c("ASCI_S_Biovolume") ~ "mg/m2",
      T~"error"
    )
  )
  
  xdf
  }

###EXAMPLE
#generate ASCI results
library(ASCI)
#example(ASCI)
results<-ASCI(tmp, demo_station)


results_swampified<-  swampify_ASCI(results) 

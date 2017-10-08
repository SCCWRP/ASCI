#' Calculate ASCI metrics
#' 
#' @param taxa chr string indicating taxa to use to calculate metrics
#' @param tax_dat data.frame of taxonomy data for diatoms, soft-bodied algae, or hybrid
#'
#' @export
#' 
#' @details \code{tax_dat} is converted to \code{taxonomy_pa} for presence/absence
#' 
#' @return a data.frame of scores for the relevant taxa
#' 
#' @examples 
#' pmmi_calcmetrics('diatoms')
pmmi_calcmetrics <- function(taxa = c('diatoms', 'sba', 'hybrid'), tax_dat){
  
  taxa <- match.arg(taxa)
  
  # convert taxonomy data to presence/absence
  taxonomy_pa <- as.data.frame(ifelse(tax_dat > 0, 1, 0))
  
  indicators.count=read.csv('lookups/indic.with.count.csv',header=TRUE,strip.white=TRUE,check.names=FALSE)
  indicators=indicators.count
  
  ##################################################################################################################################################################
  # data prep
  ##################################################################################################################################################################
  
  ###General Data Prep For Evaluating Metrics
  taxonomy_pa$SampleID<-row.names(taxonomy_pa)
  taxonomy_pa=taxonomy_pa[order(taxonomy_pa$SampleID),] #ordering by SampleID
  taxonomy_pa_melt=melt(taxonomy_pa,id=c('SampleID')) #melting data into long format
  taxonomy_pa_melt<-as.data.frame(droplevels(subset(taxonomy_pa_melt, value!=0)))
  #taxonomy_pa_melt=sqldf("select SampleID, variable as FinalIDassigned from taxonomy_pa_melt where value=1 order by SampleID") #selecting only those rows with non-zero values
  names(taxonomy_pa_melt)[2] <- "FinalIDassigned"
  
  ###Stations Data Prep
  #the next set of lines is the combining of tables to create one gaint table that is to be used in the metrics calculations below
  stations_combined1=sqldf("select * from (stations left join taxonomy_pa_melt using(SampleID)) left join traits using(FinalIDassigned)")
  stations_combined2=sqldf("select * from stations_combined1 left join indicators using(FinalIDassigned)")
  
  stations_combined3=sqldf("select * from stations_combined2
                           left join
                           taxonomy_pa_melt
                           on stations_combined2.SampleID=taxonomy_pa_melt.SampleID
                           and stations_combined2.FinalIDassigned=taxonomy_pa_melt.FinalIDassigned") # i edited the .variable to .FinalIDassigned
  
  #the next two lines create a genus variable by splitting the full taxa name
  stations_list=strsplit(stations_combined3$FinalIDassigned," ")
  stations_combined3$Genus=lapply(stations_list,FUN=function(x) x[1])
  
  ################################################################################
  # Metric calculations --------
  ################################################################################
  
  ### Evaluating Specialty Metrics
  z<-length(colnames(taxonomy_pa))
  shannon=diversity(taxonomy_pa[,-z],index='shannon') #Computing Shannon
  simpson=diversity(taxonomy_pa[,-z],index='simpson') #Computing Simpson
  richness=specnumber(taxonomy_pa[,-z]) #Computing Richness
  specialty_metrics=data.frame(SampleID=taxonomy_pa[,z],shannon,simpson,richness) #Combining the three into one table with SampleIDs attached
  
  ### Metric Calculations ------------------------------------------------------------------------------------------------------------
  metrics=ddply(.data = stations_combined3, .variables = ~SampleID,.fun = summarize,
  
                OrgN.NAHON.richness=sum(na.omit(NitrogenUptakeMetabolism=='NAHON')), #count of NAHON - species
                OrgN.NALON.richness=sum(na.omit(NitrogenUptakeMetabolism=='NALON')), #count of NALON - species
                OrgN.NHHONF.richness=sum(na.omit(NitrogenUptakeMetabolism=='NHHONF')), #count of NHHONF - species
                OrgN.NHHONForNHHONO.richness=sum(na.omit(NitrogenUptakeMetabolism=='NHHONF'|NitrogenUptakeMetabolism=='NHHONO')), #count of NHHONF and NHHONO - species
                OrgN.NHHONO.richness=sum(na.omit(NitrogenUptakeMetabolism=='NHHONO')), #count of NHHONO - species
                OxyReq.DO_100.richness=sum(na.omit(OxygenRequirements=='DO_100')), #count of DO_100 - species
                OxyReq.DO_100orDO_75.richness=sum(na.omit(OxygenRequirements=='DO_100'|OxygenRequirements=='DO_75')), #count of DO_100 and DO_75 - species
                OxyReq.DO_75.richness=sum(na.omit(OxygenRequirements=='DO_75')), #count of DO_75 - species
                OxyReq.DO_50.richness=sum(na.omit(OxygenRequirements=='DO_50')), #count of DO_50 - species
                OxyRed.DO_30.richness=sum(na.omit(OxygenRequirements=='DO_30')), #count of DO_30 - species
                OxyReq.DO_30andDO_10.richness=sum(na.omit(OxygenRequirements=='DO_30'|OxygenRequirements=='DO_10')), #count of DO_30 and DO_10 - species
                OxyReq.DO_10.richness=sum(na.omit(OxygenRequirements=='DO_10')), #count of DO_10 - species
                Salinity.B.richness=sum(na.omit(Salinity=='B')), #count of B - species
                Salinity.BF.richness=sum(na.omit(Salinity=='BF')), #count of BF - species
                Salinity.F.richness=sum(na.omit(Salinity=='F')), #count of F - species
                Salinity.FB.richness=sum(na.omit(Salinity=='FB')), #count of FB - species
                Saprobic.AM.richness=sum(na.omit(Saprobity=='AM')), #count of AM - species
                Saprobic.AMorAMPS.richness=sum(na.omit(Saprobity=='AM'|Saprobity=='AMPS')), #count of AM and AMPS - species
                Saprobic.AMPS.richness=sum(na.omit(Saprobity=='AMPS')), #count of AMPS - species
                Saprobic.BM.richness=sum(na.omit(Saprobity=='BM')), #count of BM - species
                Saprobic.OS.richness=sum(na.omit(Saprobity=='OS')), #count of OS - species
                Saprobic.OSandPS.richness=sum(na.omit(Saprobity=='OS'|Saprobity=='PS')), #count of OS and PS - species
                Saprobic.PS.richness=sum(na.omit(Saprobity=='PS')), #count of PS - species
                Trophic.E.richness=sum(na.omit(TrophicState=='E')), #count of E - species
                Trophic.EorI.richness=sum(na.omit(TrophicState=='E'|TrophicState=='I')), #count of E and I - species
                Trophic.I.richness=sum(na.omit(TrophicState=='I')), #count of I - species
                Trophic.M.richness=sum(na.omit(TrophicState=='M')), #count of M - species
                Trophic.ME.richness=sum(na.omit(TrophicState=='ME')), #count of ME - species
                Trophic.O.richness=sum(na.omit(TrophicState=='O')), #count of O - species
                Trophic.OorOMorPH.richness=sum(na.omit(TrophicState=='O'|TrophicState=='OM'|TrophicState=='PH')), #count of O, OM, and PH - species
                Trophic.OM.richness=sum(na.omit(TrophicState=='OM')), #count of OM
                Trophic.PH.richness=sum(na.omit(TrophicState=='PH')), #count of PH
                prop.spp.OrgN.NAHON=sum(na.omit(NitrogenUptakeMetabolism=='NAHON'))/length(NitrogenUptakeMetabolism), #proportion of NAHON - species
                prop.spp.OrgN.NALON=sum(na.omit(NitrogenUptakeMetabolism=='NALON'))/length(NitrogenUptakeMetabolism), #proportion of NALON - species
                prop.spp.OrgN.NHHONF=sum(na.omit(NitrogenUptakeMetabolism=='NHHONF'))/length(NitrogenUptakeMetabolism), #proportion of NHHONF - species
                prop.spp.OrgN.NHHONForNHHONO=sum(na.omit(NitrogenUptakeMetabolism=='NHHONF'|NitrogenUptakeMetabolism=='NHHONO'))/length(NitrogenUptakeMetabolism), #proportion of NHHONF and NHHONO - species
                prop.spp.OrgN.NHHONO=sum(na.omit(NitrogenUptakeMetabolism=='NHHONO'))/length(NitrogenUptakeMetabolism), #proportion of NHHONO - species
                prop.spp.OxyReq.DO_100=sum(na.omit(OxygenRequirements=='DO_100'))/length(OxygenRequirements), #proportion of DO_100 - species
                prop.spp.OxyReq.DO_100orDO_75=sum(na.omit(OxygenRequirements=='DO_100'|OxygenRequirements=='DO_75'))/length(OxygenRequirements), #proportion of DO_100 and DO_75 - species
                prop.spp.OxyReq.DO_75=sum(na.omit(OxygenRequirements=='DO_75'))/length(OxygenRequirements), #proportion of DO_75 - species
                prop.spp.OxyReq.DO_50=sum(na.omit(OxygenRequirements=='DO_50'))/length(OxygenRequirements), #proportion of DO_50 - species
                prop.spp.OxyReq.DO_atleast50=sum(na.omit(OxygenRequirements==c('DO_50','DO_75','DO_100')))/length(OxygenRequirements), #proportion of at least DO_50 - species
                prop.spp.OxyRed.DO_30=sum(na.omit(OxygenRequirements=='DO_30'))/length(OxygenRequirements), #proportion of DO_30 - species
                prop.spp.OxyReq.DO_30orDO_10=sum(na.omit(OxygenRequirements=='DO_30'|OxygenRequirements=='DO_10'))/length(OxygenRequirements), #proportion of DO_30 and DO_10 - species
                prop.spp.OxyReq.DO_10=sum(na.omit(OxygenRequirements=='DO_10'))/length(OxygenRequirements), #proportion of DO_10 - species
                prop.spp.Salinity.B=sum(na.omit(Salinity=='B'))/length(Salinity), #proportion of B - species
                prop.spp.Salinity.BF=sum(na.omit(Salinity=='BF'))/length(Salinity), #proportion of BF - species
                prop.spp.Salinity.F=sum(na.omit(Salinity=='F'))/length(Salinity), #proportion of F - species
                prop.spp.Salinity.FB=sum(na.omit(Salinity=='FB'))/length(Salinity), #proportion of FB - species
                prop.spp.Saprobic.AM=sum(na.omit(Saprobity=='AM'))/length(Saprobity), #proportion of AM - species
                prop.spp.Saprobic.AMorAMPS=sum(na.omit(Saprobity=='AM'|Saprobity=='AMPS'))/length(Saprobity), #proportion of AM and AMPS - species
                prop.spp.Saprobic.AMPS=sum(na.omit(Saprobity=='AMPS'))/length(Saprobity), #proportion of AMPS - species
                prop.spp.Saprobic.BM=sum(na.omit(Saprobity=='BM'))/length(Saprobity), #proportion of BM - species
                prop.spp.Saprobic.OS=sum(na.omit(Saprobity=='OS'))/length(Saprobity), #proportion of OS - species
                prop.spp.Saprobic.OSorPS=sum(na.omit(Saprobity=='OS'|Saprobity=='PS'))/length(Saprobity), #proportion of OS and PS - species
                prop.spp.Saprobic.PS=sum(na.omit(Saprobity=='PS'))/length(Saprobity), #proportion of PS - species
                prop.spp.Trophic.E=sum(na.omit(TrophicState=='E'))/length(TrophicState), #proportion of E - species
                prop.spp.Trophic.EorI=sum(na.omit(TrophicState=='E'|TrophicState=='I'))/length(TrophicState), #proportion of E and I - species
                prop.spp.Trophic.I=sum(na.omit(TrophicState=='I'))/length(TrophicState), #proportion of I - species
                prop.spp.Trophic.M=sum(na.omit(TrophicState=='M'))/length(TrophicState), #proportion of M - species
                prop.spp.Trophic.ME=sum(na.omit(TrophicState=='ME'))/length(TrophicState), #proportion of ME - species
                prop.spp.Trophic.O=sum(na.omit(TrophicState=='O'))/length(TrophicState), #proportion of O - species
                prop.spp.Trophic.OorOMorPH=sum(na.omit(TrophicState=='O'|TrophicState=='OM'|TrophicState=='PH'))/length(TrophicState), #proportion of O, OM, and PH - species
                prop.spp.Trophic.OM=sum(na.omit(TrophicState=='OM'))/length(TrophicState), #proportion of OM
                prop.spp.Trophic.PH=sum(na.omit(TrophicState=='PH'))/length(TrophicState), #proportion of PH
                Achnanthes.richness=sum(na.omit(Genus=='Achnanthes')), #count of genus Achnanthes
                Amphora.richness=sum(na.omit(Genus=='Amphora')), #count of genus Amphora
                Cocconeis.richness=sum(na.omit(Genus=='Cocconeis')), #count of genus Cocconeis
                Cyclotella.richness=sum(na.omit(Genus=='Cyclotella')), #count of genus Cyclotella
                Cymbella.richness=sum(na.omit(Genus=='Cymbella')), #count of genus Cymbella
                Epithemia.richness=sum(na.omit(Genus=='Epithemia')), #count of genus Epithemia
                Eunotia.richness=sum(na.omit(Genus=='Eunotia')), #count of genus Eunotia
                Fragilaria.richness=sum(na.omit(Genus=='Fragilaria')), #count of genus Fragilaria
                Frustulia.richness=sum(na.omit(Genus=='Frustulia')), #count of genus Frustulia
                Gomphonema.richness=sum(na.omit(Genus=='Gomphonema')), #count of genus Gomphonema
                Navicula.richness=sum(na.omit(Genus=='Navicula')), #count of genus Navicula
                Nitzschia.richness=sum(na.omit(Genus=='Nitzschia')), #count of genus Nitzschia
                Rhoicosphenia.richness=sum(na.omit(Genus=='Rhoicosphenia')), #count of genus Rhoicosphenia
                Rhopalodia.richness=sum(na.omit(Genus=='Rhopalodia')), #count of genus Rhopalodia
                Surirella.richness=sum(na.omit(Genus=='Surirella')), #count of genus Surirella
                Synedra.richness=sum(na.omit(Genus=='Synedra')), #count of genus Synedra
                EpiRho.richness=sum(na.omit(Genus=='Epithemia'|Genus=='Rhopalodia')), #count of genus Epithemia and Rhopalodia
                prop.Achnanthes=sum(na.omit(Genus=='Achnanthes'))/length(Genus), #proportion of genus Achnanthes
                prop.Amphora=sum(na.omit(Genus=='Amphora'))/length(Genus), #proportion of genus Amphora
                prop.Cocconeis=sum(na.omit(Genus=='Cocconeis'))/length(Genus), #proportion of genus Cocconeis
                prop.Cyclotella=sum(na.omit(Genus=='Cyclotella'))/length(Genus), #proportion of genus Cyclotella
                prop.Cymbella=sum(na.omit(Genus=='Cymbella'))/length(Genus), #proportion of genus Cymbella
                prop.Epithemia=sum(na.omit(Genus=='Epithemia'))/length(Genus), #proportion of genus Epithemia
                prop.Eunotia=sum(na.omit(Genus=='Eunotia'))/length(Genus), #proportion of genus Eunotia
                prop.Fragilaria=sum(na.omit(Genus=='Fragilaria'))/length(Genus), #proportion of genus Fragilaria
                prop.Frustulia=sum(na.omit(Genus=='Frustulia'))/length(Genus), #proportion of genus Frustulia
                prop.Gomphonema=sum(na.omit(Genus=='Gomphonema'))/length(Genus), #proportion of genus Gomphonema
                prop.Navicula=sum(na.omit(Genus=='Navicula'))/length(Genus), #proportion of genus Navicula
                prop.Nitzschia=sum(na.omit(Genus=='Nitzschia'))/length(Genus), #proportion of genus Nitzschia
                prop.Rhoicosphenia=sum(na.omit(Genus=='Rhoicosphenia'))/length(Genus), #proportion of genus Rhoicosphenia
                prop.Rhopalodia=sum(na.omit(Genus=='Rhopalodia'))/length(Genus), #proportion of genus Rhopalodia
                prop.Surirella=sum(na.omit(Genus=='Surirella'))/length(Genus), #proportion of genus Surirella
                prop.Synedra=sum(na.omit(Genus=='Synedra'))/length(Genus), #proportion of genus Synedra
                prop.AchOverAchPlusNav=sum(na.omit(Genus=='Achnanthes'))/sum(na.omit(Genus=='Achnanthes'|Genus=='Navicula')), #proportion of genus Achnanthes divided by genus Achnanthes or Navicula
                prop.CymOverCymPlusNav=sum(na.omit(Genus=='Cymbella'))/sum(na.omit(Genus=='Cymbella'|Genus=='Navicula')), #proportion of genus Cymbella divided by genus Cymbella or Navicula
                prop.EpiOverEpiPlusRho=sum(na.omit(Genus=='Epithemia'))/sum(na.omit(Genus=='Epithemia'|Genus=='Rhopalodia')), #proportion of genus Epithemia divided by genus Epithemia or Rhopalodia
                cnt.spp.HighMotility=sum(na.omit(Motility=='HM')), #count of high motility - species
                cnt.spp.ModMotility=sum(na.omit(Motility=='MM')), #count of moderate motility - species
                cnt.spp.MinMotility=sum(na.omit(Motility=='LM')), #count of low motility - species
                cnt.spp.NonMotile=sum(na.omit(Motility=='NM')), #count of no motility - species
                prop.spp.HighMotility=sum(na.omit(Motility=='HM'))/length(Motility), #proportion of high motility - species
                prop.spp.ModMotility=sum(na.omit(Motility=='MM'))/length(Motility), #proportion of moderate motility - species
                prop.spp.MinMotility=sum(na.omit(Motility=='LM'))/length(Motility), #proportion of low motility - species
                prop.spp.NonMotile=sum(na.omit(Motility=='NM'))/length(Motility), #proportion of no motility - species
                prop.NNS=sum(na.omit(Genus=='Navicula'|Genus=='Nitzschia'|Genus=='Surirella'))/length(Genus), #proportion of NNS
                prop.ind.least.tol=sum(na.omit(value[designation=='Reference']))/sum(na.omit(value)), #proportion of individuals that are least tolerant
                prop.ind.most.tol=sum(na.omit(value[designation=='Stressed']))/sum(na.omit(value)), #proportion of individuals that are most tolerant
                prop.spp.least.tol=sum(na.omit(designation=='Reference'))/length(designation), #proportion of species that are least tolerant
                prop.spp.most.tol=sum(na.omit(designation=='Stressed'))/length(designation), #proportion of species that are most tolerant
                cnt.spp.BCG1=sum(na.omit(BCG=='1')), #count of BCG 1 - species
                prop.spp.BCG1=sum(na.omit(BCG=='1'))/length(BCG), #proportion of BCG 1 - species
                cnt.spp.BCG2=sum(na.omit(BCG=='2')), #count of BCG 2 - species
                prop.spp.BCG2=sum(na.omit(BCG=='2'))/length(BCG), #proportion of BCG 2 - species
                cnt.spp.BCG3=sum(na.omit(BCG=='3')), #count of BCG 3 - species
                prop.spp.BCG3=sum(na.omit(BCG=='3'))/length(BCG), #proportion of BCG 3 - species
                cnt.spp.BCG4=sum(na.omit(BCG=='4')), #count of BCG 4 - species
                prop.spp.BCG4=sum(na.omit(BCG=='4'))/length(BCG), #proportion of BCG 4 - species
                cnt.spp.BCG5=sum(na.omit(BCG=='5')), #count of BCG 5 - species
                prop.spp.BCG5=sum(na.omit(BCG=='5'))/length(BCG), #proportion of BCG 5 - species
                cnt.spp.BCG6=sum(na.omit(BCG=='6')), #count of BCG 6 - species
                prop.spp.BCG6=sum(na.omit(BCG=='6'))/length(BCG), #proportion of BCG 6 - species
                cnt.spp.BCG12=sum(na.omit(BCG12=='yes')), #count of BCG 6 - species
                prop.spp.BCG12=sum(na.omit(BCG12=='yes'))/length(BCG), #proportion of BCG 6 - species
                cnt.spp.BCG45=sum(na.omit(BCG45=='yes')), #count of BCG 6 - species
                prop.spp.BCG45=sum(na.omit(BCG45=='yes'))/length(BCG), #proportion of BCG 6 - species
                cnt.spp.IndicatorClass_TP_high=sum(na.omit(IndicatorClass_TP=='high')), #from betty
                prop.spp.IndicatorClass_TP_high=sum(na.omit(IndicatorClass_TP=='high'))/length(IndicatorClass_TP),
                cnt.spp.IndicatorClass_TP_low=sum(na.omit(IndicatorClass_TP=='low')), #from betty
                prop.spp.IndicatorClass_TP_low=sum(na.omit(IndicatorClass_TP=='low'))/length(IndicatorClass_TP),
                cnt.spp.IndicatorClass_Cu_high=sum(na.omit(IndicatorClass_Cu=='high')),
                prop.spp.IndicatorClass_Cu_high=sum(na.omit(IndicatorClass_Cu=='high'))/length(IndicatorClass_Cu),
                cnt.spp.IndicatorClass_DOC_high=sum(na.omit(IndicatorClass_DOC=='high')),
                prop.spp.IndicatorClass_DOC_high=sum(na.omit(IndicatorClass_DOC=='high'))/length(IndicatorClass_DOC),
                cnt.spp.IndicatorClass_Ref=sum(na.omit(IndicatorClass_Ref=='RF')),
                prop.spp.IndicatorClass_Ref=sum(na.omit(IndicatorClass_Ref=='RF'))/length(IndicatorClass_Ref),
                cnt.spp.IndicatorClass_TN_high=sum(na.omit(IndicatorClass_TN=='high')),
                prop.spp.IndicatorClass_TN_high=sum(na.omit(IndicatorClass_TN=='high'))/length(IndicatorClass_TN),
                cnt.spp.Heterocy=sum(na.omit(Heterocystous=='Yes')),
                prop.spp.Heterocy=sum(na.omit(Heterocystous=='Yes'))/length(Heterocystous),
                cnt.spp.NHeterotroph=sum(na.omit(NitrogenUptakeMetabolism2=='Heterotroph')),
                prop.spp.NHeterotroph=sum(na.omit(NitrogenUptakeMetabolism2=='Heterotroph'))/length(NitrogenUptakeMetabolism2),
                cnt.spp.Halo=sum(na.omit(Salinity2=='Halo')),
                prop.spp.Halo=sum(na.omit(Salinity2=='Halo'))/length(Salinity2),
                cnt.spp.ZHR=sum(na.omit(ZHR=='Yes')), #from betty
                prop.spp.ZHR=sum(na.omit(ZHR=='Yes'))/length(IndicatorClass_TP),
                cnt.spp.CRUS=sum(na.omit(CRUS=='Yes')), #from betty
                prop.spp.CRUS=sum(na.omit(CRUS=='Yes'))/length(IndicatorClass_TP),
                cnt.spp.Green=sum(na.omit(Green=='Yes')), #from betty
                prop.spp.Green=sum(na.omit(Green=='Yes'))/length(IndicatorClass_TP)
  
  )
  
  # Combining With Specialty Metrics and Outputting
  if (taxa == 'sba') {
    out=sqldf("select * from metrics left join specialty_metrics using(SampleID)") #Linking metrics and specialty_metrics tables
  
    } #Saving completer metrics dataframe as CSV
  if (taxa == 'hybrid') {
    out=sqldf("select * from metrics left join specialty_metrics using(SampleID)") #Linking metrics and specialty_metrics tables
  
    } #Saving completer metrics dataframe as CSV
  
  # now for the one modeled metric
  if (taxa == 'diatoms') {
       #   for (i in d.win.metrics.modeled) {
       #   filename<-paste0("~/Documents/R/ASCI/pMMI.all/diatoms.models.out/", "diatoms.", i, "RF.model.Rdata")
       #   load(filename)
       #   assign(paste0(i, ".d.RF"), rf.out.top)
       #   print(paste(i, "is modeled"))}
    filename<-paste0("data/", "diatoms.", "prop.spp.Salinity.BF", "RF.model.Rdata")
    load(filename)
    assign(paste0("prop.spp.Salinity.BF", ".d.RF"), rf.out.top)
    required.predictors <- row.names(prop.spp.Salinity.BF.d.RF$importance)
    required.predictors
    missingpredictors <- setdiff(required.predictors, colnames(stations))
    if (length(missingpredictors) > 0 ) {print(paste("missing predictors", missingpredictors)) }
    stations.predictors<-stations[,required.predictors]
    row.names(stations.predictors) <- stations$SampleID
    predicted<-predict(prop.spp.Salinity.BF.d.RF, stations.predictors)
    observed<-metrics[,"prop.spp.Salinity.BF"]
    resid<-observed-predicted
    metrics$prop.spp.Salinity.BF<-resid
  
    out = sqldf("select * from metrics left join specialty_metrics using(SampleID)") #Linking metrics and specialty_metrics tables
    
    }
  
  return(out)

}

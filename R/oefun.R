#' Run the ASCI O/E 
#' 
#' Run the ASCI O/E index for diatoms, soft-bodied algae, and hybrids
#' 
#' @param taxain chr string for path to input taxonomy data
#' @param sitein chr string for path to input site data
#' 
#' @details 
#' Three index scores are calculated and returned as a named list
#' 
#' @return 
#' A list with three elements named as \code{Scores.diatoms}, \code{Scores.sba}, and \code{Scores.hybrid}.  Each element includes an additional list of three data.frames names \code{OE.scores}, \code{Capture.Probs}, and \code{Group.Occurrence.Probs}.
#' 
#' @export
#' 
#' @examples 
#' oefun(demo_algae_tax, demo_algae_sitedata)
oefun <- function(taxain, sitein){

  # Step 1. Import taxonomy data -----------------------------------------------------------
  bugs<- taxain
  reqfields<- c("StationCode", "SampleDate", "Replicate","SampleTypeCode", "BAResult", "Result", "FinalID")
  missingtaxafields<-setdiff(reqfields, colnames(bugs))
  if( length(missingtaxafields) >0 ) { print(paste("Missing fields", missingtaxafields))}
  bugs$SampleID <- paste(bugs$StationCode, bugs$SampleDate, bugs$Replicate, sep="_")
  
  # Step 2. Import stations data -----------------------------------------------------------
  stations <- sitein
  stations$SampleID <- paste(stations$StationCode, stations$SampleDate, stations$Replicate, sep="_")
  row.names(stations) <- stations$SampleID
  missingsites<-setdiff(bugs$StationCode, stations$StationCode)
  if(length(missingsites) > 0 ) {print(paste("Missing station codes", missingsites)) }
  stations<-subset(stations, stations$StationCode %in% bugs$StationCode)
  
  # Step 3. Make ASCI-readable taxa names -----------------------------------------------------------
  STE<-read.csv("lookups/algae_STE.csv", stringsAsFactors = F)
  #bugs$FinalID2<-toupper(bugs$FinalID)
  #STE$FinalID2<-toupper(STE$FinalID)
  unrecognizedtaxa <- setdiff(bugs$FinalID, STE$FinalID)
  if (length(unrecognizedtaxa) > 0 ) { print(paste("Unrecognized taxa", unrecognizedtaxa))}  
  bugs<- merge(bugs, STE[,c("FinalID", "FinalIDassigned", "Genus", "Phylum", "Class")], all.x = T) # non matches get purged for now  #this is now case sensitive, could change
  bugs.d<-subset(bugs, Class=="Bacillariophyceae")
  bugs.sba<-subset(bugs, Class!="Bacillariophyceae")
  bugs$ComboResult<-as.numeric(pmax(bugs$BAResult,bugs$Result, na.rm=T))
    
  # Step 4. Rarify diatom data -----------------------------------------------------------
  bugs.d.sub<-rarify(inbug=bugs.d, sample.ID="SampleID", abund="BAResult", subsiz=500) 
  
  # Step 5. Convert to species abd matrix at Genus level  -----------------------------------------------------------
  bugs.d.m<-as.data.frame(acast(bugs.d.sub, SampleID~Genus, value.var="BAResult", fun.aggregate=sum))
  bugs.sba.m<-as.data.frame(acast(bugs.sba, SampleID~Genus, value.var="Result", fun.aggregate=sum))
  bugs.hybrid.m<-as.data.frame(acast(bugs, SampleID~Genus, value.var="ComboResult", fun.aggregate=sum))
  
  # Make presence/absence 
  bugs.d.m <- ifelse(bugs.d.m > 0,1,0)
  bugs.sba.m <- ifelse(bugs.sba.m > 0,1,0)
  bugs.hybrid.m <- ifelse(bugs.hybrid.m > 0,1,0)

  # Step 5. Calculate O/E for diatoms -----------------------------------------------------------
  required.d.oe.predictors <- row.names(diatom.rf.oe$importance)
  missingpredictors <- setdiff(required.d.oe.predictors, colnames(stations))
  if (length(missingpredictors) > 0 ) {print(paste("missing predictors", missingpredictors)) }
  
  stations.d.oe.predictors<-stations[,required.d.oe.predictors]
  row.names(stations.d.oe.predictors) <- stations$SampleID
  
  calib.bugs.d.tax.refcal<-read.csv("lookups/diatoms.bugs.rc.csv", stringsAsFactors = F, row.names=1)
  calib.stations.d.refcal.BG<- read.csv("lookups/diatoms.stations.rc.csv", stringsAsFactors = F, row.names=1)
  
  Scores.diatoms <- rfpred(
    bugcal.pa = calib.bugs.d.tax.refcal,
    grps.final = calib.stations.d.refcal.BG$BG,
    preds.final = required.d.oe.predictors,
    ranfor.mod = diatom.rf.oe, 
    prednew = stations.d.oe.predictors,
    bugnew = bugs.d.m
    )

  # Step 5. Calculate O/E for sba ----------------------------------------------------------------------------------------------------------------------
  required.sba.oe.predictors <- row.names(sba.rf.oe$importance)
  required.sba.oe.predictors
  missingpredictors <- setdiff(required.sba.oe.predictors, colnames(stations))
  if (length(missingpredictors) > 0 ) {print(paste("missing predictors", missingpredictors)) }
  
  stations.sba.oe.predictors<-stations[,required.sba.oe.predictors]
  row.names(stations.sba.oe.predictors) <- stations$SampleID
  
  calib.bugs.sba.tax.refcal<-read.csv("lookups/sba.bugs.rc.csv", stringsAsFactors = F, row.names=1)
  calib.stations.sba.refcal.BG<- read.csv("lookups/sba.stations.rc.csv", stringsAsFactors = F, row.names=1)
  
  Scores.sba <- rfpred(
    bugcal.pa = calib.bugs.sba.tax.refcal,
    grps.final = calib.stations.sba.refcal.BG$BG,
    preds.final = required.sba.oe.predictors,
    ranfor.mod = sba.rf.oe, 
    prednew = stations.sba.oe.predictors,
    bugnew = bugs.sba.m
    )

  # Step 5. Calculate O/E for hybrid ----------------------------------------------------------------------------------------------------------------------
  required.hybrid.oe.predictors <- row.names(hybrid.rf.oe$importance)
  required.hybrid.oe.predictors
  missingpredictors <- setdiff(required.hybrid.oe.predictors, colnames(stations))
  if (length(missingpredictors) > 0 ) {print(paste("missing predictors", missingpredictors)) }
  
  stations.hybrid.oe.predictors<-stations[,required.hybrid.oe.predictors]
  row.names(stations.hybrid.oe.predictors) <- stations$SampleID
  
  calib.bugs.hybrid.tax.refcal<-read.csv("lookups/hybrid.bugs.rc.csv", stringsAsFactors = F, row.names=1)
  calib.stations.hybrid.refcal.BG<- read.csv("lookups/hybrid.stations.rc.csv", stringsAsFactors = F, row.names=1)

  Scores.hybrid <- rfpred(
    bugcal.pa = calib.bugs.hybrid.tax.refcal,
    grps.final = calib.stations.hybrid.refcal.BG$BG,
    preds.final = required.hybrid.oe.predictors,
    ranfor.mod = hybrid.rf.oe, 
    prednew = stations.hybrid.oe.predictors,
    bugnew = bugs.hybrid.m
    )

  out <- list(
    Scores.diatoms = Scores.diatoms,
    Scores.sba = Scores.sba,
    Scores.hybrid = Scores.hybrid
  )
  
  return(out)

}







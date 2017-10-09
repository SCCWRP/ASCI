#' Run the ASCI pMMI 
#' 
#' Run the ASCI pMMI index for diatoms, soft-bodied algae, and hybrids
#' 
#' @param taxain chr string for path to input taxonomy data
#' @param sitein chr string for path to input site data
#' 
#' @details 
#' Three index scores are calculated and returned as a named list
#' 
#' @return 
#' A list with three elements named as \code{d.results.scored}, \code{sba.results.scored}, and \code{hybrid.results.scored}. Each element is a \code{data.frame} with metric scores by site.
#'
#' @export
#' 
#' @examples 
#' pmmifun(demo_algae_tax, demo_algae_sitedata)
pmmifun <- function(taxain, sitein){

  options(gsubfn.engine = "R")
  
  # Step 1. Import taxonomy data -----------------------------------------------------------
  bugs <- taxain
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
  unrecognizedtaxa <- setdiff(bugs$FinalID, STE$FinalID)
  if (length(unrecognizedtaxa) > 0 ) { print(paste("Unrecognized taxa", unrecognizedtaxa))}
  bugs<- merge(bugs, STE[,c("FinalID", "FinalIDassigned", "Genus", "Phylum", "Class")], all.x = T) # non matches get purged for now  #this is now case sensitive, could change
  bugs.d<-subset(bugs, Class=="Bacillariophyceae")
  bugs.sba<-subset(bugs, Class!="Bacillariophyceae")
  bugs$ComboResult<-as.numeric(pmax(bugs$BAResult,bugs$Result, na.rm=T))
  
  # Step 4. Rarify diatom data -----------------------------------------------------------
  bugs.d.sub<-rarify(inbug=bugs.d, sample.ID="SampleID", abund="BAResult", subsiz=500)
  
  # Step 5. Convert to species abd matrix at Species level  -----------------------------------------------------------
  bugs.d.m<-as.data.frame(acast(bugs.d.sub, SampleID~FinalIDassigned, value.var="BAResult", fun.aggregate=sum))
  bugs.sba.m<-as.data.frame(acast(bugs.sba, SampleID~FinalIDassigned, value.var="Result", fun.aggregate=sum))
  bugs.hybrid.m<-as.data.frame(acast(bugs, SampleID~FinalIDassigned, value.var="ComboResult", fun.aggregate=sum))
  
  # Import traits table
  traits <- pmmilkup$traits

  # calculate metrics using runpMMI.calcmetrics.R
  d.metrics<-pmmi_calcmetrics('diatoms', bugs.d.m)
  sba.metrics<-pmmi_calcmetrics('sba', bugs.sba.m)
  hybrid.metrics<-pmmi_calcmetrics('hybrid', bugs.hybrid.m)
  
  # Load winning metrics -----------------------------------------------------------
  d.win <- pmmilkup$d.win 
  sba.win <- pmmilkup$sba.win
  hybrid.win <- pmmilkup$hybrid.win
  d.win<-colnames(d.win[,-(length(names(d.win)))])
  sba.win<-colnames(sba.win[,-(length(names(sba.win)))])
  hybrid.win<-colnames(hybrid.win[,-(length(names(hybrid.win)))])
  
  # subset winning metrics for new sites
  d.results<-d.metrics[,d.win]
  sba.results<-sba.metrics[,sba.win]
  hybrid.results<- hybrid.metrics[,hybrid.win]
  
  # score results ------------------------------------------------------
  
  quants <- pmmilkup$quants
  # for increasers, min is 5th of ref and max is 95th of str
  # for decreasers, min is 5th of str and max is 95th of ref
  
  d.results.scored<-data.frame(row.names(d.results))
  sba.results.scored<-data.frame(row.names(sba.results))
  hybrid.results.scored<-data.frame(row.names(hybrid.results))
  
  # increase (obs - max) / ( min - max)
  d.results.scored$prop.spp.Salinity.BF <- ((d.results$prop.spp.Salinity.BF - 12.9216849) / (-0.12408199 - 12.9216849)) / 1.2768156
  d.results.scored$prop.spp.HighMotility <- ((d.results$prop.spp.HighMotility - 0.4322735) / (0 - 0.4322735)) / 1.0160970
  d.results.scored$prop.ind.most.tol <- ((d.results$prop.ind.most.tol - 0.4798261) / (0 - 0.4798261)) / 1.0169511
  sba.results.scored$cnt.spp.IndicatorClass_TP_high <- ((sba.results$cnt.spp.IndicatorClass_TP_high - 5) / (0 - 5)) / 1.1299149
  sba.results.scored$prop.spp.IndicatorClass_DOC_high <- ((sba.results$prop.spp.IndicatorClass_DOC_high - 0.7736111) / (0.03125000 - 0.7736111)) / 1.1658549
  sba.results.scored$prop.spp.Green <- ((sba.results$prop.spp.Green - 0.7000000) / (0 - 0.7000000)) / 1.1019174
  hybrid.results.scored$prop.spp.IndicatorClass_DOC_high <- ((hybrid.results$prop.spp.IndicatorClass_DOC_high - 0.2500000) / (0.01445135 - 0.2500000)) / 1.5067045
  
  #decrease (obs - min) / (max - min)
  d.results.scored$prop.spp.BCG3 <- ((d.results$prop.spp.BCG3 - 0) / (0.4670833 - 0)) / 0.7808209
  sba.results.scored$prop.spp.BCG3 <- ((sba.results$prop.spp.BCG3 - 0) / (0.5000000 - 0)) / 0.7387326
  hybrid.results.scored$prop.spp.Trophic.I <- ((hybrid.results$prop.spp.Trophic.I - 0) / (0.1764706 - 0)) / 0.8838483
  hybrid.results.scored$prop.spp.ZHR <- ((hybrid.results$prop.spp.ZHR - 0) / (0.1878378 - 0)) / 1.2028715
  hybrid.results.scored$prop.spp.BCG3 <- ((hybrid.results$prop.spp.BCG3 - 0.05579365) / (0.4196944 - 0.05579365)) / 0.7299234
  
  row.names(d.results.scored)<-d.results.scored[,1]; d.results.scored<-d.results.scored[,2:5]
  row.names(sba.results.scored)<-sba.results.scored[,1]; sba.results.scored<-sba.results.scored[,2:5]
  row.names(hybrid.results.scored)<-hybrid.results.scored[,1]; hybrid.results.scored<-hybrid.results.scored[,2:5]
  
  d.results.scored[d.results.scored>1]<-1
  d.results.scored[d.results.scored<0]<-0
  sba.results.scored[sba.results.scored>1]<-1
  sba.results.scored[sba.results.scored<0]<-0
  hybrid.results.scored[hybrid.results.scored>1]<-1
  hybrid.results.scored[hybrid.results.scored<0]<-0
  
  # compile results ------------------------------------------------------
  
  d.results.scored$diatom.pMMI<-rowMeans(d.results.scored)
  sba.results.scored$sba.pMMI<-rowMeans(sba.results.scored)
  hybrid.results.scored$hybri.pMMI<-rowMeans(hybrid.results.scored)

  # output
  out <- list(
    d.results.scored = d.results.scored, 
    sba.results.scored = sba.results.scored,
    hybrid.results.scored = hybrid.results.scored
  )
  
  return(out)
  
}
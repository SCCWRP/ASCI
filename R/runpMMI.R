# # Run the ASCI pMMI 
# 
# #### Load packages 
# library("sampling")
# library("rms")
# library("randomForest")
# library(sqldf)
# library(vegan)
# library(indicspecies)
# library(reshape2)
# library(plyr) # WARNING: problems arise if you have dplyr already loaded, so maybe start from scratch
# options(gsubfn.engine = "R")
# source("R/OE.load.and.source.R")
# 
# # Step 1. Import taxonomy data -----------------------------------------------------------
# bugs<- read.csv("demo_algae_tax.csv", stringsAsFactors = F) 
# reqfields<- c("StationCode", "SampleDate", "Replicate","SampleTypeCode", "BAResult", "Result", "FinalID")
# missingtaxafields<-setdiff(reqfields, colnames(bugs))
# if( length(missingtaxafields) >0 ) { print(paste("Missing fields", missingtaxafields))}
# bugs$SampleID <- paste(bugs$StationCode, bugs$SampleDate, bugs$Replicate, sep="_")
# 
# # Step 2. Import stations data -----------------------------------------------------------
# stations<-read.csv("demo_algae_sitedata.csv", stringsAsFactors = F)
# stations$SampleID <- paste(stations$StationCode, stations$SampleDate, stations$Replicate, sep="_")
# row.names(stations) <- stations$SampleID
# missingsites<-setdiff(bugs$StationCode, stations$StationCode)
# if(length(missingsites) > 0 ) {print(paste("Missing station codes", missingsites)) }
# stations<-subset(stations, stations$StationCode %in% bugs$StationCode)
# 
# # Step 3. Make ASCI-readable taxa names -----------------------------------------------------------
# STE<-read.csv("lookups/algae_STE.csv", stringsAsFactors = F)
# #bugs$FinalID2<-toupper(bugs$FinalID)
# #STE$FinalID2<-toupper(STE$FinalID)
# unrecognizedtaxa <- setdiff(bugs$FinalID, STE$FinalID)
# if (length(unrecognizedtaxa) > 0 ) { print(paste("Unrecognized taxa", unrecognizedtaxa))}  
# bugs<- merge(bugs, STE[,c("FinalID", "FinalIDassigned", "Genus", "Phylum", "Class")], all.x = T) # non matches get purged for now  #this is now case sensitive, could change
# bugs.d<-subset(bugs, Class=="Bacillariophyceae")
# bugs.sba<-subset(bugs, Class!="Bacillariophyceae")
# bugs$ComboResult<-as.numeric(pmax(bugs$BAResult,bugs$Result, na.rm=T))
#   
# # Step 4. Rarify diatom data -----------------------------------------------------------
# bugs.d.sub<-rarify(inbug=bugs.d, sample.ID="SampleID", abund="BAResult", subsiz=500) 
# 
# # Step 5. Convert to species abd matrix at Species level  -----------------------------------------------------------
# bugs.d.m<-as.data.frame(acast(bugs.d.sub, SampleID~FinalIDassigned, value.var="BAResult", fun.aggregate=sum))
# bugs.sba.m<-as.data.frame(acast(bugs.sba, SampleID~FinalIDassigned, value.var="Result", fun.aggregate=sum))
# bugs.hybrid.m<-as.data.frame(acast(bugs, SampleID~FinalIDassigned, value.var="ComboResult", fun.aggregate=sum))
# 
# # Convert to presence/absence -----------------------------------------------------------
# bugs.d.m<-as.data.frame(ifelse(bugs.d.m>0,1,0))
# bugs.sba.m<-as.data.frame(ifelse(bugs.sba.m>0,1,0))
# bugs.hybrid.m<-as.data.frame(ifelse(bugs.hybrid.m>0,1,0))
# 
# # Import traits table 
# traits<-read.csv('lookups/combotraits.fromaaron4.csv',header=TRUE,strip.white=TRUE,check.names=FALSE)
# 
# # calculate metrics using runpMMI.calcmetrics.R
# diatoms = T 
# sba = F
# hybrid = F 
# source("R/runpMMI.calcmetrics.R")
# diatoms = F
# sba = T
# hybrid = F 
# source("R/runpMMI.calcmetrics.R")
# diatoms = F
# sba = F
# hybrid = T 
# source("R/runpMMI.calcmetrics.R")
# 
# # Import calculated metrics (that you just generated)
# d.metrics<-read.csv("diatom.metrics.csv", row.names=1, stringsAsFactors = F)
# sba.metrics<-read.csv("sba.metrics.csv", row.names=1, stringsAsFactors = F)
# hybrid.metrics<-read.csv("hybrid.metrics.csv", row.names=1, stringsAsFactors = F)
# 
# # Load winning metrics -----------------------------------------------------------
# d.win<-read.csv("lookups/diatoms.combined.win.metrics.scaled.csv", row.names=1, stringsAsFactors = F)
# sba.win<-read.csv("lookups/sba.combined.win.metrics.scaled.csv", row.names=1, stringsAsFactors = F)
# hybrid.win<-read.csv("lookups/hybrid.combined.win.metrics.scaled.csv", row.names=1, stringsAsFactors = F)
# d.win<-colnames(d.win[,-(length(names(d.win)))])
# sba.win<-colnames(sba.win[,-(length(names(sba.win)))])
# hybrid.win<-colnames(hybrid.win[,-(length(names(hybrid.win)))])
# 
# # subset winning metrics for new sites 
# d.results<-d.metrics[,d.win]
# sba.results<-sba.metrics[,sba.win]
# hybrid.results<- hybrid.metrics[,hybrid.win]
# 
# # score results ------------------------------------------------------
# 
# quants<-read.csv("lookups/quants/quants.csv", stringsAsFactors = F)
# # for increasers, min is 5th of ref and max is 95th of str
# # for decreasers, min is 5th of str and max is 95th of ref 
# 
# d.results.scored<-data.frame(row.names(d.results))
# sba.results.scored<-data.frame(row.names(sba.results))
# hybrid.results.scored<-data.frame(row.names(hybrid.results))
# 
# # increase (obs - max) / ( min - max) 
# d.results.scored$prop.spp.Salinity.BF <- ((d.results$prop.spp.Salinity.BF - 12.9216849) / (-0.12408199 - 12.9216849)) / 1.2768156
# d.results.scored$prop.spp.HighMotility <- ((d.results$prop.spp.HighMotility - 0.4322735) / (0 - 0.4322735)) / 1.0160970
# d.results.scored$prop.ind.most.tol <- ((d.results$prop.ind.most.tol - 0.4798261) / (0 - 0.4798261)) / 1.0169511
# sba.results.scored$cnt.spp.IndicatorClass_TP_high <- ((sba.results$cnt.spp.IndicatorClass_TP_high - 5) / (0 - 5)) / 1.1299149
# sba.results.scored$prop.spp.IndicatorClass_DOC_high <- ((sba.results$prop.spp.IndicatorClass_DOC_high - 0.7736111) / (0.03125000 - 0.7736111)) / 1.1658549
# sba.results.scored$prop.spp.Green <- ((sba.results$prop.spp.Green - 0.7000000) / (0 - 0.7000000)) / 1.1019174
# hybrid.results.scored$prop.spp.IndicatorClass_DOC_high <- ((hybrid.results$prop.spp.IndicatorClass_DOC_high - 0.2500000) / (0.01445135 - 0.2500000)) / 1.5067045
# 
# #decrease (obs - min) / (max - min)
# d.results.scored$prop.spp.BCG3 <- ((d.results$prop.spp.BCG3 - 0) / (0.4670833 - 0)) / 0.7808209
# sba.results.scored$prop.spp.BCG3 <- ((sba.results$prop.spp.BCG3 - 0) / (0.5000000 - 0)) / 0.7387326
# hybrid.results.scored$prop.spp.Trophic.I <- ((hybrid.results$prop.spp.Trophic.I - 0) / (0.1764706 - 0)) / 0.8838483
# hybrid.results.scored$prop.spp.ZHR <- ((hybrid.results$prop.spp.ZHR - 0) / (0.1878378 - 0)) / 1.2028715
# hybrid.results.scored$prop.spp.BCG3 <- ((hybrid.results$prop.spp.BCG3 - 0.05579365) / (0.4196944 - 0.05579365)) / 0.7299234
# 
# row.names(d.results.scored)<-d.results.scored[,1]; d.results.scored<-d.results.scored[,2:5]
# row.names(sba.results.scored)<-sba.results.scored[,1]; sba.results.scored<-sba.results.scored[,2:5]
# row.names(hybrid.results.scored)<-hybrid.results.scored[,1]; hybrid.results.scored<-hybrid.results.scored[,2:5]
# 
# d.results.scored[d.results.scored>1]<-1
# d.results.scored[d.results.scored<0]<-0
# sba.results.scored[sba.results.scored>1]<-1
# sba.results.scored[sba.results.scored<0]<-0
# hybrid.results.scored[hybrid.results.scored>1]<-1
# hybrid.results.scored[hybrid.results.scored<0]<-0
# 
# # compile results ------------------------------------------------------
# 
# d.results.scored$diatom.pMMI<-rowMeans(d.results.scored)
# sba.results.scored$sba.pMMI<-rowMeans(sba.results.scored)
# hybrid.results.scored$hybri.pMMI<-rowMeans(hybrid.results.scored)
# 
# 
# # output results ------------------------------------------------------
# 
# write.csv(d.results.scored, "diatom.results.scored.csv")
# write.csv(sba.results.scored, "sba.results.scored.csv")
# write.csv(hybrid.results.scored, "hybrid.results.scored.csv")
# 
# 
# 

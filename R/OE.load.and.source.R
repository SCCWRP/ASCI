###################################################################################
#  LOAD AND SOURCE # 
###################################################################################
library(reshape2)
library(ggplot2)
library("cluster")
library(randomForest)
library(vegan)

source("R/dapply.R")
source("R/assess.one.sample.4.1.r")

jaccfun<-function(siti,sitj) {
  shared<-sum((siti>0)&(sitj>0));
  uniquei<-sum((siti>0)&(sitj==0));
  uniquej<-sum((siti==0)&(sitj>0));
  1-(shared/(shared+uniquei+uniquej)); #return Jaccard dissimilarity;
} #end of function;

sornfun<-function(siti,sitj) {
  shared<-sum((siti>0)&(sitj>0));
  uniquei<-sum((siti>0)&(sitj==0));
  uniquej<-sum((siti==0)&(sitj>0));
  1-(2*shared/(2*shared+uniquei+uniquej)); #return Sorenson dissimilarity;
} #end of function;

bcfun<-function(siti,sitj) {
  bcnum<-sum(abs(siti-sitj));
  bcdenom<-sum(siti+sitj);
  ifelse(bcdenom>0, (bcnum/bcdenom),0); #return BC dissimilarity;
} #end of function;

rep.sam.sd<-function(occprb.cal,Pc) {
  cal.cut<-(occprb.cal>=Pc); #TRUE/FALSE matrix denoting which taxa are above predicted P cutoff;
  #use occurrence probabilities only of taxa above the cutoff (cal.cut='TRUE');
  E.cal<-apply(occprb.cal*cal.cut,1,sum); #vector of predicted E ;
  # numerator of site-specific replicate sampling var. Result is a site vector
  RS.cal<-apply(occprb.cal*cal.cut,1,function(x)sum(x*(1-x)));
  SDRS<-sqrt(mean(RS.cal/(E.cal^2))); #replicate sampling SD is sqrt(mean(site-specific replicate sampling variances));
  print(' ',quote=F)
  print(' Replicate sampling SD of O/E: ',quote=F)
  print(SDRS,digits=4);  }; #end of function;


## Assess individual sites;
#This function assesses a single site or sample from a new (test) data set to which
# model.predict.v4.1() has already been applied.
# It compares observed occurrences with the model-predicted probabilities of occurrence for all taxa;

#Input parameters are:
#  case -- A selected site or sample ID, for which a prediction has already been made using model.predict.v4.1(). ;
# result.prd -- Output of model.predict.v4.1() for new samples that include the chosen case;
# bugnew  -- Sample-by-taxa matrix (abundance or presence/absence) of new samples that was submitted to model.predict.v4.1().
# Pc -- Cutoff for capture probabilties for inclusion of taxa in O/E;

#The function outputs a data frame with one row for each taxon, and the following columns:
 # observed presence(1) or absence(0);
 # predicted capture probability;
# Big.diff = "Yes", if there is a big difference (>=0.5 in magnitude) between observed and predicted;
# In.OtoE = "Yes" if the taxon would be included in the O/E calculation for this sample, given the stated value of Pc;

#NOTE that rows (taxa) in this data frame are sorted by the magnitude of (observed-predicted),
   # as suggested in Van Sickle, J. (2008), JNABS 27:227-235;

#The function prints prints the data frame, followed by O,E and O/E for this case;

#Example usage:
# site1.result<-assess.one.sample(case="01020CSR",result.prd=OE.assess.test, bugnew=bug.test, Pc=0.5);

     assess.one.sample.4.1<-function(case, result.prd, bugnew, Pc){
     #first, need to restructure bugnew to match taxa of the model predictions;
      #convert bug matrix to P/A (1/0);
      temp.pa<-bugnew;
      temp.pa[temp.pa>0]<-1;
      #reshape bugnew.pa columns (taxa) to match those in bugcal.pa, and be in same order;
      # New bug data might have fewer or more columns, and it might have extraneous columns;
       # create a new empty site x taxa matrix with bugcal.pa columns and bugnew rows, fill it with zeros;
       nsite.new<-dim(temp.pa)[[1]];
       ntaxa<-dim(result.prd$Capture.Probs)[[2]];
       bugnew.pa<-matrix(rep(0,nsite.new*ntaxa),nrow=nsite.new,dimnames=list(rownames(temp.pa),dimnames(result.prd$Capture.Probs)[[2]]));
      #loop through columns of new matrix and fill with columns of the original test data matrix;
       col.match<-match(dimnames(bugnew.pa)[[2]],dimnames(temp.pa)[[2]]);
       for(kcol in 1:ntaxa) if(!is.na(col.match[kcol]))bugnew.pa[,kcol]<-temp.pa[,col.match[kcol]];
     #Restructure of new bug data matrix is complete;
    ############# Now work on the special case;
       observed<-as.vector(t(bugnew.pa[case,]));
       predicted<-as.vector(t(result.prd$Capture.Probs[case,]));
       ttt<-data.frame(cbind(observed ,predicted),Big.diff=ifelse(abs(observed-predicted)>=0.5,"Yes","  "),
           In.OtoE=ifelse(predicted>=Pc,"Yes","  "), row.names=dimnames(bugnew.pa)[[2]]);
        ttt<-ttt[order(abs(observed-predicted),decreasing=TRUE),];
      print(paste(c("Observed occurrences and predicted probabilities for case ", case),collapse=""),quote=F);
      print(ttt,quote=F);
      print("",quote=F);
      print(result.prd$OE.scores[case,],quote=F)
      return(ttt)};



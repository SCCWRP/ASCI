#' Predict O/E scores for new sites
#'  
#' Predict O/E scores for a new set of sites based on a Random Forest predictive model
#' 
#' @param bugcal.pa \code{data.frame} of taxonomy at calibration reference sites
#' @param grps.final \code{data.frame} of environmetnal data at calibration reference sites
#' @param preds.final chr string of predictors used in \code{ranfor.mod}
#' @param ranfor.mod \code{\link[randomForest]{randomForest}} model used for prediction
#' @param prednew \code{data.frame} containing predictor variables (columns) at all new sites/samples (rows) at which predictions are desired (e.g., 'test' sites)
#' @param bugnew \code{data.frame} of taxonomy presence/absence data at all new sites/samples (rows) at which predictions are desired (e.g., 'test' sites)
#' @param Pc numeric value for capture probability cutoff
#' @param Cal.OOB logical value default \code{FALSE}, set to \code{TRUE} only if you are predicting for the model calibration data and desire out-of-bag predictions 
#' @param saverfe logical value indicating if RFE scores are saved as a text file in the working directory
#'
#' @author J. Van Sickle, USEPA, \email{VanSickle.John@epa.gov}
#'
#' @note Version 4.2 - Added option for out-of-bag predictions on CAL data. 12/30/10
#' 
#' @export
#' 
#' @return function output is a list containing three elements
#' \itemize{ 
#' \item \code{OE.scores} A data frame for all samples, containing O, E, O/E and BC from the predictive and null models, as well as outlier flags
#' \item \code{Capture.Probs} Matrix of model-predicted capture (occurrence) probabilties for all taxa in all samples
#' \item \code{Group.Occurrnce.Probs} Matrix of predicted probabilities of occurrence for each sample in each calibration-site group
#' }
rfpred <- function(bugcal.pa, grps.final, preds.final, ranfor.mod, prednew, bugnew, Pc = 0.5, Cal.OOB = FALSE, saverfe = FALSE) {

  #first convert bug matrix to P/A (1/0)
  temp.pa<-bugnew
  temp.pa[temp.pa>0]<-1
  rm(bugnew)

  #1. - initial definitions
  names(grps.final)<-row.names(bugcal.pa)
  nsite.cal<-length(grps.final) #number of calibration sites
  npreds<-length(preds.final) #number of predictor variables
  grpsiz<-table(grps.final) #tabulate sites per group
  ngrps<-length(grpsiz)  #number of groups

  #2. Alignment of new predictor and bug data with model data
  #2a) Align the rows (samples) of the new bug data to the new predictor data
  temp.pa<-temp.pa[row.names(prednew),]
  #2b)reshape bugnew columns (taxa) to match those in bugcal.pa, and be in same order
  # New bug data might have fewer or more columns (taxa) than calibration bug data
  # create a new empty site x taxa matrix with bugcal.pa columns and bugnew rows, fill it with zeros
  nsite.new<-dim(temp.pa)[[1]]
  ntaxa<-dim(bugcal.pa)[[2]]
  bugnew.pa<-matrix(rep(0,nsite.new*ntaxa),nrow=nsite.new,dimnames=list(rownames(temp.pa),colnames(bugcal.pa)))
  #loop through columns of new matrix and fill with columns of the original test data matrix
  col.match<-match(dimnames(bugnew.pa)[[2]],dimnames(temp.pa)[[2]])
  for(kcol in 1:ntaxa) if(!is.na(col.match[kcol]))bugnew.pa[,kcol]<-temp.pa[,col.match[kcol]]
  
  ## STEP 3. -- Use RF to predict the group (cluster) membership for all new sites. 
  # Does not use RIVPACS assumption of weighting the membership probabilities by Calibration group size, as a prior probability
  # Also, RF predictions do not have an outlier test, unlike DFA predictions
  # Predicted probs are outputted as a matrix, sites are rows, columns are groups
  
  #If Cal.OOB is true, do OOB predictions, appropriate ONLY for CAL data
  # If it is false, do a new prediction
  if(Cal.OOB==TRUE) grpprobs<-ranfor.mod$votes else grpprobs<-predict(ranfor.mod,newdata=prednew[,preds.final],type='prob')

  ############
  #STEP 4 -- Compute predicted occurrence probability for each modeled taxon at each new sample
  # "modeled OTU's" consist of all taxa that were found at >=1 calibration sample
  #To do this, first calculate the occurrence freqs of all modeled taxa in the Calibration sample groups
  grpocc<-apply(bugcal.pa,2,function(x){tapply(x,grps.final,function(y){sum(y)/length(y)})})
  
  #finally, compute the matrix of predicted occurrence (capture) probabilities, for all new samples and all modeled taxa
  #This is the matrix-algebra form of the RIVPACS combining formula (e.g., Clarke et al. 2003, Eq. 4)
  site.pred.dfa<-grpprobs%*%grpocc

  #######################

  # STEP 5. Compute O, E, O/E and BC for all samples. 
  # Also compute O/E and BC for the null model

  #5.1 loop over all samples. Compute and store  O, predicted E, predicted BC for each sample. 
  #temporary data frame to hold nonnull results for all samples. 
  nsit.new<-dim(prednew)[[1]]
  OE.stats<-data.frame(OBS=rep(NA,nsit.new), E.prd=rep(NA,nsit.new),BC.prd=rep(NA,nsit.new),row.names=row.names(prednew))
  for(i in 1:nsit.new) {
     #i<-1
     cur.prd<-site.pred.dfa[i,] #vector of predicted probs for current sample
     spdyn<-names(cur.prd)[cur.prd>=Pc]  #subset of taxa with Pi>=Pcutoff
     cur.prd<-cur.prd[spdyn] #vector of Pi for subset of included taxa
     cur.obs<-bugnew.pa[i,spdyn] #vector of observed P/A for those taxa
     OE.stats$OBS[i]<-sum(cur.obs) #observed richness (O)
     OE.stats$E.prd[i]<-sum(cur.prd) #Expected richness (E)
     OE.stats$BC.prd[i]<-sum(abs(cur.obs-cur.prd))/ (OE.stats$OBS[i]+OE.stats$E.prd[i]) #BC value
           } #finish sample loop

  #5.2 - Compute Expected richness (E) and BC for null model using taxa >= Pc.
  # Note that the set of taxa included in the null model is fixed for all samples
  pnull<-apply(bugcal.pa,2,sum)/dim(bugcal.pa)[[1]]  #null model predicted occurrnece probabilities, all taxa
  nulltax<-names(pnull[pnull>=Pc]) #subset of taxa with Pnull >= Pc
  Enull<-sum(pnull[nulltax])
  Obsnull<-apply(bugnew.pa[,nulltax],1,sum) #vector of Observed richness, new samples, under null model
  BC.null<-apply(bugnew.pa[,nulltax],1,function(x)sum(abs(x-pnull[nulltax])))/(Obsnull+Enull) #vector of null-model BC
  
  #5.3 - Final data frame contains values of O, E, O/E, Onull, Enull, Onull/Enull, BC.prd and BC.null, for all samples
  #Also includes outlier flags
  
  OE.final<-data.frame(O=OE.stats$OBS,E=OE.stats$E.prd,
                        OoverE=OE.stats$OBS/OE.stats$E.prd, 
     Onull=Obsnull,Enull=rep(Enull,length(Obsnull)),OoverE.null=Obsnull/Enull,
     BC= OE.stats$BC.prd,BC.null=BC.null,
            row.names=row.names(bugnew.pa))

  #print some summary statistics of O/E to wrap up
  c1<-mean(OE.final$OoverE, na.rm = T) 
  c2<-mean(OE.final$OoverE.null, na.rm=T)
  s1<-sqrt(var(OE.final$OoverE, na.rm = T)) 
  s2<-sqrt(var(OE.final$OoverE.null, na.rm = T))
  outz<-paste(c1, s1, c2, s2, digits=3)
  saveas2<-paste0("out_RFscores", ".txt")
  saveas3<-paste(saveas2, outz)
  if(saverfe) capture.output(print(saveas3), file="RFEscores.txt", append=T)
  
  #function output is a list containing OE.final, matrix of predicted capture probs, and predicted group membership probs
  out <- list(OE.scores=OE.final,Capture.Probs=site.pred.dfa,Group.Occurrence.Probs=grpprobs)
  
  return(out)
  
  }




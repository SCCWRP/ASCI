#' Rarify a macroinvertebrate sample
#' 
#' Rarify (subsample) a macroinvertebrate sample down to a fixed count
#'
#' @param inbug \code{data.frame} for input data
#' @param sample.ID chr string of column name with the sample ID
#' @param abund chr string of column name with the abundance data
#' @param subsiz numeric value indicating size of the subsample
#'
#' @author John Van Sickle, USEPA, \email{VanSickle.John@epa.gov}
#' 
#' @note Version 1.0 06/10/05
#' 
#' @return A rarified \code{data.frame}.
rarify <- function(inbug, sample.ID, abund, subsiz){

  # setup parameters
  outbug<-inbug
  sampid<-unique(inbug[,sample.ID])
  nsamp<-length(sampid)

  #zero out all abundances in output data set
  outbug[,abund]<-0
  
  #loop over samples, rarify each one in turn
  for(i in 1:nsamp){ 

    #extract current sample
    isamp<-sampid[i]

    onesamp<-inbug[inbug[,sample.ID]==isamp,]
    onesamp<-data.frame(onesamp,row.id=seq(1,dim(onesamp)[[1]])) #add sequence numbers as a new column
    #expand the sample into a vector of individuals
    samp.expand<-rep(x=onesamp$row.id,times=onesamp[,abund])
    nbug<-length(samp.expand) #number of bugs in sample
    #vector of uniform random numbers
    ranvec<-runif(n=nbug)
    #sort the expanded sample randomly
    samp.ex2<-samp.expand[order(ranvec)]
    #keep only the first piece of ranvec, of the desired fised count size
    #if there are fewer bugs than the fixed count size, keep them all
    if(nbug>subsiz){subsamp<-samp.ex2[1:subsiz]} else{subsamp<-samp.ex2}
    #tabulate bugs in subsample
    subcnt<-table(subsamp)
    #define new subsample frame and fill it with new reduced counts
    newsamp<-onesamp
    newsamp[,abund]<-0
    newsamp[match(newsamp$row.id,names(subcnt),nomatch=0)>0,abund]<-as.vector(subcnt)
    outbug[outbug[,sample.ID]==isamp,abund]<-newsamp[,abund]
    
    }

  return(outbug)
  
  }

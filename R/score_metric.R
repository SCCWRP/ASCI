
#' Create scored metrics
#'
#' @param taxa chr string indicating taxa to use to calculate metrics
#' @param bugs.m species abd matrix at Species level
#' @param results.metric metrics calculated from \code{mmi_calmetrics} combined with win metrics
#' @param omni.ref a data frame in \code{mmilkup} indicates directions of the metrics
#'
#' @return a data.frame scored
#' @export 
#'
#' @examples
#' \dontrun{
#' score_metric(taxa = 'diatoms', bugs.d.m, d.results, omni.ref)
#' }

score_metric <- function(taxa, bugs.m, results.metric, omni.ref){
  scored <- data.frame(rowname = row.names(bugs.m))
  foo <- which(colnames(results.metric) %in% omni.ref$Metric)
  cols <- names(results.metric)[foo]
  
  for (i in cols){      # i<-"prop.Cyclotella"
    foo <- omni.ref %>% 
      filter(Metric == i, 
             Assemblage == taxa)
    min <- foo$Min[1]
    max <- foo$Max[1]
    observed<-results.metric[i]
    if(foo$StressResponse == 'inc'){
      a <-(observed-max)
      b <-(min-max)
    } else {
      a <-(observed-min)
      b <-(max-min)
    }
    c <- as.data.frame(a/b)
    names(c)[1] <- i
    scored <- cbind(scored,c)
  }
  return(scored)
}
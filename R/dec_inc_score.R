
#' Title
#'
#' @param taxa 
#' @param bugs.m 
#' @param results.metric 
#' @param omni.ref 
#'
#' @return
#' @export
#'
#' @examples
score_metric <- function(taxa, bugs.m, results.metric, omni.ref){
  scored <- data.frame(rowname = row.names(bugs.m))
  foo <- which(colnames(results.metric) %in% omni.ref$Metric)
  cols <- names(win.metric)[foo]
  
  for (i in cols){      # i<-"prop.Cyclotella"
    foo <- omni.ref %>% 
      filter(Metric == i, 
             Assemblage == taxa)
    min <- foo$Min[1]
    max <- foo$Max[1]
    observed<-win.metric[[i]]
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
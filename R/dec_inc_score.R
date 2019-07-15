#,
#,
#,
#,

score_metric <- function(win, species_matrix, dir_df, inc = T){
  scored <- data.frame(row.names(species_matrix))
  foo <- which(colnames(win) %in% rownames(dir_df))
  cols <- names(win)[foo]
  
  for (i in cols){      # i<-"prop.Cyclotella"
    foo <- subset(dir_df, X==i)
    min <- foo$min[1]
    max <- foo$max[1]
    observed<-win[[i]]
    if(inc = T){
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
# caret source 

library(randomForest)
library(caret)
# library(foreach)
# registerDoSEQ()

################## 
# functions and cutoffs #
################## 
ctrl <- rfeControl(functions = rfFuncs, method = "cv",verbose = FALSE, returnResamp = "all")

#this is Rafi's update 
ctrl$functions$selectSize<-
  function (x, metric, tol = tolerance, maximize)
  {
    if (!maximize) {
      best <- min(x[which(x$Variables<=10), metric])
      perf <- (x[, metric] - best)/best * 100
      flag <- perf <= tol
    }
    else {
      best <- max(x[which(x$Variables<=10), metric])
      perf <- (best - x[, metric])/best * 100
      flag <- perf <= tol
    }
    min(x[flag, "Variables"])
  }

tolerance<- 1

set.seed(106)
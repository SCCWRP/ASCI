######
#' Random forest models for metrics
#'
#' List of random forest models for predictive metrics
#'  
#' @format A \code{\link[base]{list}} with six elements, one for each predictive metric and a final element for conductivity estimates:
#' \describe{
#'   \item{\code{d.prop.spp.SPIspecies4}}{randomForest}
#'   \item{\code{d.Salinity.BF.richness}}{randomForest}
#'   \item{\code{h.OxyRed.DO_30.richness}}{randomForest}
#'   \item{\code{h.prop.spp.BCG4}}{randomForest}
#'   \item{\code{h.Salinity.BF.richness}}{randomForest}
#'   \item{\code{temp.cond.qrf}}{quantregForest}
#' }
#'
#' @examples 
#' data(rfmods)
"rfmods"
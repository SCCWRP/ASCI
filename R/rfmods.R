######
#' Random forest models for metrics
#'
#' List of random forest models for predictive metrics
#'  
#'  
#' @format A \code{\link[base]{list}} with six elements, one for each predictive metric and a final element for conductivity estimates:
#' \describe{
#'   \item{\code{diatoms.prop.spp.BCG12}}{randomForest}
#'   \item{\code{diatoms.prop.spp.Salinity.BF}}{randomForest}
#'   \item{\code{diatoms.prop.spp.Trophic.E}}{randomForest}
#'   \item{\code{hybrid.OxyRed.DO_30.richness}}{randomForest}
#'   \item{\code{hybrid.prop.spp.BCG4}}{randomForest}
#'   \item{\code{hybrid.Salinity.BF.richness}}{randomForest}
#'   \item{\code{cond.qrf}}{quantregForest}
#' }
#'
#' @examples 
#' data(rfmods)
"rfmods"
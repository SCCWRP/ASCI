######
#' Site classifications and PSA regions
#'
#' Site classifications and PSA regions
#'  
#' @format A \code{\link{data.frame}} of site classifications and regions.
#'
#' @details The object is used internally in \code{\link{perf}}. The \code{SampleID} follows the format created by \code{\link{getids}}.  Site classications are:
#' \itemize{
#' \item{int}{intermediate}
#' \item{notrecent}{duplicate site visit}
#' \item{rc}{reference calibration}
#' \item{rv}{reference validation}
#' \item{str}{stressed}
#' }
#' PSA regions are those defined for California.
#' 
#' @examples 
#' data(sitcat)
"sitcat"
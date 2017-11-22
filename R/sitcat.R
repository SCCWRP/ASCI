######
#' Site classifications and PSA regions
#'
#' Site classifications and PSA regions
#'  
#' @format A \code{\link{data.frame}} of site classifications and regions.
#'
#' @details The object is used internally in \code{\link{perf}}. The \code{SampleID} follows the format created by \code{\link{getids}}.  Site classications are:
#' \itemize{
#' \item{\code{int}} {intermediate}
#' \item{\code{notrecent}} {duplicate site visit}
#' \item{\code{rc}} {reference calibration}
#' \item{\code{rv}} {reference validation}
#' \item{\code{str}} {stressed}
#' }
#' PSA regions are those defined for California.
#' 
#' @examples 
#' data(sitcat)
"sitcat"
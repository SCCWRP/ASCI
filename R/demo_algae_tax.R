######
#' Example data for algae taxonomy
#'
#' An example dataset for algae taxonomy. Data are in long format, one row per sample. Sites match those in \code{\link{demo_algae_sitedata}}.
#'  
#' @format A \code{\link[base]{data.frame}} with 105 rows and 8 variables:
#' \describe{
#'   \item{\code{StationCode}}{chr}
#'   \item{\code{SampleDate}}{chr}
#'   \item{\code{Replicate}}{int}
#'   \item{\code{CollectionMethodCode}}{chr}
#'   \item{\code{SampleTypeCode}}{chr}
#'   \item{\code{BAResult}}{num}
#'   \item{\code{Result}}{num}
#'   \item{\code{FinalID}}{chr}
#' }
#'
#' @examples 
#' data(demo_algae_tax)
"demo_algae_tax"
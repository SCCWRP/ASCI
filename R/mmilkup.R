######
#' Lookup tables for MMI calculations
#'
#' Lookup tables for MMI calculations used in \code{\link{mmifun}} and \code{\link{mmi_calcmetrics}}
#'  
#' @format A \code{\link[base]{list}} with six elements named \code{traits}, \code{d.win}, \code{sba.win}, \code{hybrid.win}, \code{omni.ref}, and \code{indicators}.
#'
#' @details Each element is a \code{\link[base]{data.frame}} where all are used in \code{\link{mmifun}}, except \code{indicators} which is used in \code{\link{mmi_calcmetrics}}.
#' 
#' @seealso \code{\link{mmifun}}, \code{\link{mmi_calcmetrics}}
#' 
#' @examples 
#' data(mmilkup)
"mmilkup"
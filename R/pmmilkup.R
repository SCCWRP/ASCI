######
#' Lookup tables for pMMI calculations
#'
#' Lookup tables for pMMI calculations used in \code{\link{pmmifun}} and \code{\link{pmmi_calcmetrics}}
#'  
#' @format A \code{\link[base]{list}} with six elements named \code{traits}, \code{d.win}, \code{sba.win}, \code{hybrid.win}, \code{quants}, and \code{indicators}.
#'
#' @details Each element is a \code{\link[base]{data.frame}} where all are used in \code{\link{pmmifun}}, except \code{indicators} which is used in \code{\link{pmmi_calcmetrics}}.
#' 
#' @seealso \code{\link{pmmifun}}, \code{\link{pmmi_calcmetrics}}
#' 
#' @examples 
#' data(pmmilkup)
"pmmilkup"
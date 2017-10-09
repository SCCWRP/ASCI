######
#' Lookup tables for O/E calculations
#'
#' Lookup tables for O/E calculations used in \code{\link{oefun}}
#'  
#' @format A \code{\link[base]{list}} with six elements named \code{diatoms.bugs.rc}, \code{diatoms.stations.rc}, \code{sba.bugs.rc}, \code{sba.stations.rc}, \code{hybrid.bugs.rc}, and \code{hybrid.stations.rc}.
#'
#' @details Each element is a \code{\link[base]{data.frame}} of taxonomic ('bugs') or site ('stations') information for diatoms, soft-bodied algae, or hyrid taxa at reference calibration sites.
#' 
#' @seealso \code{\link{oefun}}
#' 
#' @examples 
#' data(oelkup)
"oelkup"
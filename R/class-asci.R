#' @title An S4 class to represent an \code{asci} object
#' 
#' @description The S4 \code{asci} object includes several methods shown below.
#' 
#' @slot scores \code{data.frame} of ASCI scores
#' @slot Supp1_mmi \code{data.frame} of supplemental MMI metric scores
#' @slot taxa character string of taxa returned with initial function call
#' 
#' @examples 
#' results <- ASCI(demo_algae_tax, demo_algae_sitedata)
#' scores(results)
#' Supp1_mmi(results)
#' Supp1_OE(results)
#' Supp2_OE(results)
asci <- setClass('asci', 
         slots = list(
           scores = 'data.frame',
           Supp1_mmi = 'data.frame', 
           taxa = 'character'
           )
)


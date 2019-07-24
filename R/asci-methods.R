#' @include class-asci.R
NULL

#' @rdname asci-class
setMethod('show', 
          signature = 'asci', 
          definition = function(object){
            
            # defined methods from below
            meths <- .S4methods(class = 'asci')
            meths <-  attr(meths, 'info')$generic
            meths <- meths[!meths %in% 'show']

            # unique samples 
            n <- scores(object)
            n <- unique(n$SampleID)
            n <- length(n)
              
            # on screen
            cat('An object of class', class(object), '\n')
            cat('Scores calculated for', paste(object@taxa, collapse = ', ') , 'indices for', n, 'unique samples\n')
            cat('Use these functions for access:', paste(meths, collapse = ', '), '\n')
            invisible(NULL)
            
          })

#' @param object \code{asci} object created with \code{\link{ASCI}}
#' @rdname asci-class
#' 
#' @export
setGeneric('scores', function(object) standardGeneric('scores'))

#' @rdname asci-class
setMethod('scores', 'asci', function(object) object@scores)
          
#' @rdname asci-class
#' @export
setGeneric('Supp1_mmi', function(object) standardGeneric('Supp1_mmi'))

#' @rdname asci-class
setMethod('Supp1_mmi', 'asci', function(object) object@Supp1_mmi)

#' @title Get performance measures for ASCI
#' 
#' @description Get performance measures for ASCI using statewide results
#'
#' @param object \code{\link{asci}} object for statewide results
#'
#' @return A \code{data.frame} of performance measures by index, index group (diatoms, soft-bodied algae, hybrid), predicted/null, and reference calibration/validation
#' 
#' @export
#'
#' @details This functions estimates seven performance measures that evaluate index accuracy, precision, and responsiveness for calibration and validation datasets. The following are estimated:
#' \itemize{
#' \item{\code{ave}} {Average scores at reference sites}
#' \item{\code{fst}} {F statistic of reference scores between ecoregions}
#' \item{\code{vrs_nat}} {Variance of index scores explained by natural gradients at reference sites}
#' \item{\code{prc_amg}} {Precision among sites as standard deviation of scores at reference sites}
#' \item{\code{prc_wth}} {Precision within sites as standard deviation of within-site residuals for repeat site visits at reference sites}
#' \item{\code{res_tst}} {t-statistic as a measure of precision between scores at reference and stressed sites}
#' \item{\code{res_var}} {Variance of index scores explained by human-activity gradients at all sites as a measure of responsiveness}
#' }
#' This function requires site classifications as reference (calibration and validation), stressed, and repeat visits.  Regions for each site are also required.
#' 
#' @importFrom dplyr filter group_by inner_join left_join mutate rename
#' @importFrom magrittr "%>%"
#' @importFrom purrr map
#' @importFrom tidyr gather nest separate unnest
#' 
#' @seealso \code{\link{sitcat}}
#' 
#' #' @examples
#' #' perf(allscr)
#' setGeneric('perf', function(object) standardGeneric('perf'))
#' 
#' #' @rdname perf
#' setMethod('perf', 'asci', function(object){
#'   
#'   # get predictive index scores from object
#'   scr <- object %>% 
#'     scores %>% 
#'     dplyr::select(SampleID, taxa, OoverE, MMI) %>% 
#'     rename(
#'       oe = OoverE, 
#'       mmi = MMI,
#'       grp = taxa
#'     ) %>% 
#'     gather('ind', 'scr', oe, mmi) %>% 
#'     mutate(typ = 'Predictive')
#'   
#'   # get null scores from results
#'   nll <- object %>% 
#'     .@null_OE %>% 
#'     rename(
#'       grp = taxa,
#'       ind = met, 
#'       scr = val
#'     ) %>% 
#'     dplyr::select(SampleID, grp, ind, scr) %>% 
#'     filter(ind %in% 'OoverE.null') %>% 
#'     mutate(
#'       ind = 'oe',
#'       typ = 'Null'
#'     ) 
#'   
#'   # rf mods
#'   rfmods <- list(
#'     diatoms_oe_Predictive = diatom_rf_oe,
#'     hybrid_oe_Predictive = hybrid_rf_oe,
#'     sba_oe_Predictive = sba_rf_oe,
#'     diatoms_mmi_Predictive = rf_out_top
#'   ) %>% 
#'     enframe('val', 'rfmod') %>% 
#'     separate(val, c('grp', 'ind', 'typ'))
#'   
#'   # combine predictive and null scores
#'   # add site classification, regions (remove CV)
#'   # add rf mods
#'   toprf <- rbind(scr, nll) %>% 
#'     left_join(sitcat, by = 'SampleID') %>% 
#'     na.omit %>% 
#'     filter(!psa %in% 'Central_Valley') %>% 
#'     filter(cls %in% c('rc', 'rv', 'str', 'notrecent')) %>% 
#'     separate(SampleID, c('StationID', 'Replicate'), sep = '_[0-9]+/[0-9]+/[0-9]+_', remove = FALSE) %>% 
#'     group_by(grp, ind, typ) %>% 
#'     nest %>% 
#'     left_join(rfmods, by = c('grp', 'ind', 'typ'))
#'   
#'   # get performance by grp, ind, typ combo
#'   prf <- toprf %>% 
#'     mutate(prf = pmap(list(data, rfmod), function(data, rfmod){
#'       
#'       # nest by cal, val
#'       ref <- data %>% 
#'         filter(cls %in% c('rc', 'rv')) %>% 
#'         group_by(cls) %>% 
#'         nest
#'       
#'       # get stressed
#'       str <- data %>% 
#'         filter(cls %in% 'str')
#'       
#'       # get reps
#'       ntr <- data %>% 
#'         filter(cls %in% 'notrecent')
#'       
#'       # variance with natural gradients
#'       vrs_nat <- tail(rfmod$rsq, 1)
#'       if(is.null(vrs_nat))
#'         vrs_nat <- NA
#'       
#'       # responsiveness as var
#'       res_var <- NA
#'       
#'       # evaluate by calibration/validation
#'       refprf <- ref %>% 
#'         mutate(
#'           evl = map(data, function(scrs){
#'             
#'             # average scores
#'             ave <- scrs %>% 
#'               .$scr %>% 
#'               mean(na.rm = TRUE)
#'             
#'             # f stat
#'             fst <- scrs %>% 
#'               lm(scr ~ psa, .) %>% 
#'               summary %>% 
#'               .$fstatistic %>% 
#'               .[1]   
#'             
#'             # precision among sites
#'             prc_amg <- scrs %>% 
#'               .$scr %>% 
#'               sd(na.rm = TRUE)
#'             
#'             # precision within sites
#'             prc_wth <- scrs %>% 
#'               inner_join(ntr, by = 'StationID') %>% 
#'               dplyr::select(StationID, scr.x, scr.y) %>% 
#'               gather('scr', 'val', scr.x, scr.y) %>% 
#'               group_by(StationID) %>% 
#'               mutate(val = val - mean(val)) %>% 
#'               .$val %>% 
#'               sd(na.rm = TRUE)
#'             
#'             # responsiveness as t.test
#'             res_tst <- scrs %>% 
#'               .$scr %>% 
#'               t.test(., str$scr) %>% 
#'               .$statistic
#'             
#'             # combine for output
#'             out <- list(
#'               ave = ave,
#'               fst = fst,
#'               prc_amg = prc_amg,
#'               prc_wth = prc_wth, 
#'               res_tst = res_tst
#'             ) %>% 
#'               data.frame
#'             
#'             return(out)
#'             
#'           })
#'         ) %>% 
#'         dplyr::select(-data) %>% 
#'         unnest
#'       
#'       # add res_var, vrs_nat
#'       refprf$res_var <- res_var
#'       refprf$vrs_nat <- vrs_nat
#'       
#'       # order columns
#'       refprf <- refprf %>% 
#'         # dplyr::select(cls, ave, fst, vrs_nat, prc_amg, prc_wth, res_tst, res_var)
#'         dplyr::select(cls, ave, fst, prc_amg, prc_wth, res_tst)
#'       
#'       return(refprf)
#'       
#'       })
#' 
#'     ) %>% 
#'     dplyr::select(-data, -rfmod) %>%
#'     unnest
#'   
#'   return(prf)
#'   
#' })
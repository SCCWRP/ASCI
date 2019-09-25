#' @title Calculate missing station data
#'
#' @description Calculate missing station data for \code{XerMtn} and \code{CondQR50}
#' 
#' @param station \code{data.frame} for input station data
#'
#' @return The original station data with calculated fields where applicable.
#' 
#' @details 
#' \code{XerMtn} is calculated from \code{PSA6C} if the latter is present with no \code{NA} value. If \code{XerMtn} is 
#' present with \code{NA} values, \code{PSA6C} is used if present with no \code{NA} values. \code{XerMtn = 1} if
#' \code{PSA6C} is \code{'SN', 'NC'}, otherwise \code{XerMtn = 0} if \code{PSA6C} is \code{'CH', 'CV', 'DM', 'SC'}.
#' 
#' \code{CondQR50} is calculated from the \code{\link[quantregForest]{quantregForest}} object in \code{\link{rfmods}}
#' if all predictors are present in stations and no missing values are in the predictors.  Similar to \code{XerMtn}, 
#' \code{NA} values will be predicted per row only if the column already exists.  The required predictors are \code{'CaO_Mean', 
#' 'MgO_Mean', 'S_Mean', 'UCS_Mean', 'LPREM_mean', 'AtmCa', 'AtmMg', 'AtmSO4', 'MINP_WS', 'MEANP_WS', 'SumAve_P',
#' 'TMAX_WS', 'XWD_WS', 'MAXWD_WS', 'LST32AVE', 'BDH_AVE', 'KFCT_AVE', 'PRMH_AVE'}.
#' 
#' @export
#' 
#' @seealso \code{\link{chkinp}}, \code{\link{getids}}
#'
#' @importFrom dplyr mutate
#' @importFrom magrittr "%>%"
#' @import quantregForest
#' 
#' @examples
#' 
#' # this calculates CondQR50 and XerMtn
#' calcgis(demo_station)
#' 
#' \dontrun{
#' # get XerMtn from PSA6c
#' tmp <- demo_station[, !names(demo_station) %in% 'XerMtn']
#' calcgis(tmp)
#' 
#' # error, cannot get XerMtn if PSA6c not found
#' tmp <- demo_station[, !names(demo_station) %in% c('XerMtn', 'PSA6C')]
#' calcgis(tmp)
#' 
#' # get conductivity
#' tmp <- demo_station
#' calcgis(tmp)
#' 
#' # get conductivity for only NA 
#' tmp <- demo_station
#' tmp$CondQR50[1] <- 200
#' calcgis(tmp)
#' 
#' # error, cannot calculate conductivity if missing predictors
#' tmp <- demo_station[, !names(demo_station) %in% c('TMAX_WS', 'AtmSO4')]
#' calcgis(tmp)
#' 
#' # error, cannot calculate conductivity if missing values in predictors
#' tmp <- demo_station
#' tmp$MINP_WS[2] <- NA
#' tmp$AtmSO4[3] <- NA
#' calcgis(tmp)
#' }
calcgis <- function(station){
  
  nms <- names(station)
  
  ##
  # XerMTN
  
  # if XerMtn not found, calculate from PSA6c if found
  if(!'XerMtn' %in% nms){
    
    if(!'PSA6C' %in% nms)
      stop('Cannot calculate XerMtn if PSA6c is missing', call. = FALSE)
    
    if(any(is.na(station$PSA6C)))
      stop('Cannot calculate XerMtn if any missing values in PSA6c', call. = FALSE)
    
    station <- station %>% 
      mutate(
        XerMtn = case_when(
          PSA6C %in% c('CH', 'CV', 'DM', 'SC') ~ 0, 
          PSA6C %in% c('SN', 'NC') ~ 1
        )
      )
    
  }
  
  # if XerMtn is found, replace any NA using PSA if found
  if('XerMtn' %in% nms){
    
    if(any(is.na(station$XerMtn))){
      
      if(!'PSA6C' %in% nms)
        stop('Cannot calculate XerMtn if PSA6c is missing', call. = FALSE)
    
      if(any(is.na(station$PSA6C)))
        stop('Cannot calculate XerMtn if any missing values in PSA6c', call. = FALSE)
      
      station <- station %>% 
        mutate(
          XerMtn = case_when(
            is.na(XerMtn) & PSA6C %in% c('CH', 'CV', 'DM', 'SC') ~ 0, 
            is.na(XerMtn) & PSA6C %in% c('SN', 'NC') ~ 1
          )
        )

    }
    
  }
  
  ##
  # CondQR50

  # predictors for conductivity model
  condvals <- c("CaO_Mean", "MgO_Mean", "S_Mean", "UCS_Mean", "LPREM_mean",
                "AtmCa", "AtmMg", "AtmSO4", "MINP_WS", "MEANP_WS", "SumAve_P",
                "TMAX_WS", "XWD_WS", "MAXWD_WS", "LST32AVE", "BDH_AVE", "KFCT_AVE",
                "PRMH_AVE")
  
  # message if missing conductivity predictors
  msgmissprd <- setdiff(condvals, nms) %>% 
    paste(collapse = ', ') %>% 
    paste('Missing the following conductivity predictors in station data:', .)
 
  # message if NA values in conductivity predictors
  msgnaprd <- station[, nms %in% condvals] %>% 
    apply(2, anyNA) %>% 
    which %>% 
    names %>% 
    paste(collapse = ', ') %>% 
    paste('Missing values in conductivity predictors:', .)
  
  # conductivity model
  condmod <- rfmods$temp.cond.qrf

  # if CondQR50 not found, calculate from predictors if found
  if(!'CondQR50' %in% nms){

    # must have all predictors
    if(sum(nms %in% condvals) < 18)
      stop(msgmissprd, call. = FALSE)
  
    # must have no NA values in predictors
    if(anyNA(station[, condvals]))
      stop(msgnaprd, call. = FALSE)

    # get predictions
    station <- station %>%
      mutate(
        CondQR50 = predict(condmod, newdata = .[, condvals], what = 0.5)
      )
    
  }
  
  # if CondQR50 is found with NA, calculate from predictors if found
  if(anyNA(station$CondQR50)){
    
    # must have all predictors
    if(sum(nms %in% condvals) < 18)
      stop(msgmissprd, call. = FALSE)
    
    # must have no NA values in predictors
    if(anyNA(station[, condvals]))
      stop(msgnaprd, call. = FALSE)
    
    # get predictions
    station <- station %>%
      mutate(
        CondQR50 = case_when(
          is.na(CondQR50) ~ as.numeric(predict(condmod, newdata = .[, condvals], what = 0.5)), 
          T ~ as.numeric(CondQR50)
        )
      )
    
  }

  ##
  # return calculated info
  out <- station
  
  return(out)
  
}



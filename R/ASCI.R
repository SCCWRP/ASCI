#' Score samples using the ASCI tool
#'
#' @param taxa \code{data.frame} for input taxonomy data
#' @param stations \code{data.frame} for input station data
#' 
#' @details 
#' One index for three taxonomy types are scored, MMI for diatoms, soft-bodied algae, and hybrid. 
#' If only soft-bodied algae or diatoms are present, only the respective ASCI and metrics are returned.
#' This function outputs the results of \code{\link{mmifun}} functions in a user-friendly format.
#' 
#' @return 
#' A dataframe with all metrics calculated for each provided taxa
#'
#' @export
#' 
#' @importFrom dplyr bind_rows mutate select case_when mutate_all group_by ungroup inner_join summarize full_join
#' @importFrom magrittr "%>%"
#' @importFrom tidyr gather spread unnest unite replace_na
#' @import purrr
#' @import tibble
#' 
#' @seealso  \code{\link{mmifun}}
#' 
#' @examples 
#' # calculate all
#' ASCI(demo_algae_tax, demo_station)
#' 
#' # works if either soft-bodied or diatoms are missing
#' # remove diatoms from station sample 909M24937
#' tmp <- subset(demo_algae_tax, !(StationCode == '909M24937' & SampleTypeCode == 'Integrated'))
#' ASCI(tmp, demo_station)
#' 
#' # works if either soft-bodied or diatoms are missing
#' # remove soft-bodied from station sample 801M16916
#' tmp <- subset(demo_algae_tax, !(StationCode == '801M16916' & SampleTypeCode != 'Integrated'))
#' ASCI(tmp, demo_station)
ASCI <- function(taxa, stations){
  # Renamed the station argument to gismetrics, because it is actually GIS metric data on the stations
  # for which ASCI is being processed
  # The station data comes from the table called tblgismetrics in the SMC database

  # run all other checks, get output if passed
  dat <- chkinp(taxa, stations)
  txrmv <- dat$txrmv
  dat <- dat$taxa
  
  # calculate GIS from stations
  # calcgis calculates required gismetrics which are not provided, from the ones that were provided
  # IF possible... If not enough data was provided it throws an error, or at least it is supposed to...
  gismetrics <- calcgis(stations)
  
  # Prepare the Algae data to have metrics ran on it
  # We need to tack on the traits and indicators dataframes, which has certain important info on the species
  # That info is used to calculate metrics which are later used for ASCI
  algae <- dat %>%
    left_join(STE[,c("FinalID","FinalIDassigned")], by = "FinalID") %>%
    left_join(mmilkup$traits, by = 'FinalIDassigned') %>%
    left_join(mmilkup$indicators, by = 'FinalIDassigned') %>%
    filter(
      # I saw the sampletypecode Qualitative filtered out in previous versions of the ASCI calculator
      SampleTypeCode != 'Qualitative',
      # The na.replace function makes it so that if Result is NA, it gets counted as "Not being equal to Zero"
      # As opposed to just an NA
      # The reason we want to get rid of Result and BAResult being Zero is that BAResult or Result of Zero is the same
      # as the organism not being present
      replace_na(Result != 0, T), 
      replace_na(BAResult != 0, T)
    )
  
  # Pass the algae and the gis data to each subsequent function
  # the sba function wont need the gis data since there are no predictive metrics for SBA data
  diatom.scores <- diatoms(algae, gismetrics)
  sba.scores = sba(algae)
  hybrid.scores = hybrid(algae, gismetrics)
  
  # Here we bring them all together
  combined.scores <- diatom.scores %>%
    inner_join(
      sba.scores,
      by = "SampleID"
    ) %>%
    inner_join(
      hybrid.scores,
      by = "SampleID"
    )

  # Here we get the what I like to call the supplementary information
  supplementary_info <- algae %>%
    group_by(SampleID) %>%
    summarize(
      SampleType = paste0(unique(SampleTypeCode), collapse = ', '),
      D_NumberTaxa = sum(!is.na(BAResult[which(SampleTypeCode == 'Integrated')])), # May be incorrect
      H_NumberTaxa = length(FinalID), # May be incorrect
      S_NumberTaxa = H_NumberTaxa - D_NumberTaxa, # May be incorrect
      D_ValveCount = sum(BAResult[which(SampleTypeCode == 'Integrated')], na.rm = T),
      # This line of code seems wrong to me
      # Since it is Soft Body, wouldn't it be Result? Rather than BAResult
      S_EntityCount = sum(BAResult[which(SampleTypeCode != 'Integrated')], na.rm = T),
      S_Biovolume = sum(Result, na.rm = T),
      Comments = dplyr::case_when(
        !("Integrated" %in% SampleTypeCode) ~ "Warning - Only Soft Body data present",
        # if it gets to this next step, it has already determined that Integrated is one of the sampletypecodes
        # if the length of the unique sampletypecodes is 1, then that means that Integrated is the ONLY sampletypecode
        (length(unique(SampleTypeCode)) == 1) ~ "Warning - Only Diatom data present",
        # if it failed to meet the above criteria, its fine and we don't need to give a warning
        TRUE ~ ""
      ),
      version_number = as.character(packageVersion('ASCI'))
    ) %>%
    ungroup()

  # --- Preparing the final output ----
  # The next few lines of code are going to be purely for presentation, and organizing final output
  beginning_cols <- c('SampleID', 'StationCode','SampleDate','Replicate',
                      'SampleType','D_ValveCount','S_EntityCount','S_Biovolume','D_NumberTaxa',
                      'S_NumberTaxa','H_NumberTaxa','UnrecognizedTaxa',
                      'D_ASCI','S_ASCI','H_ASCI')

  ending_cols <- c('Comments','version_number')

  out <- combined.scores %>%
    inner_join(
      # join it with original data to tack on StationCode, SampleDate and Replicate, 
      # which disappeared because we were only grouping based on SampleID
      dat %>% 
        select(SampleID,StationCode,SampleDate,Replicate) %>% 
        unique(),
      by = 'SampleID'
    ) %>%
    inner_join(
      supplementary_info,
      by = 'SampleID'
    ) %>%
    # Tack on the Unrecognized Taxa gathered by chkinp at the beginning
    # left join because the txrmv only has rows if a sample had unrecognized taxa
    left_join(
      txrmv,
      by = 'SampleID'
    ) %>%
    select(
      # Final ordering of the columns
      all_of(beginning_cols),
      starts_with("D_"),
      starts_with("S_"),
      starts_with("H_"),
      all_of(ending_cols)
    )

  return(out)

}

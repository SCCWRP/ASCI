#' Create scored metrics
#'
#' @param metrics mmi metrics dataframe, having metrics to be scored
#' @param assemblage chr string indicating the assemblage type - diatoms, sba or hybrid
#'
#' @return a data.frame scored
#' @export
#' @importFrom tidyr pivot_longer pivot_wider
#'
#' @examples
#' \dontrun{
#' score(taxa = 'diatoms', bugs.d.m, d.results, omni.ref)
#' }

score <- function(metrics, assemblage){
  if (!assemblage %in% c('diatoms','sba','hybrid')){
    stop('Error in score_metric.R - assemblage keyword arg must be diatoms, sba or hybrid')
  }
  
  # we make the assumption that the dataframe of metrics comes in a wide 
  # For this reason, we will pivot it out to long format, to be able to join with the onmi.ref dataframe
  scored.metrics <- metrics %>% 
    pivot_longer(
      -SampleID, 
      names_to = 'Metric',
      values_to = 'Value'
    ) %>%
    # mmilkup$omni.ref is a dataframe with a column called "Metric" which is the metric name
    # and it has other columns with info about that metric that lets us get the scored version of each metric
    # doing an inner join ensures that we only keep the "raw" and "mod" version of each metric,
    #   since those are the only ones that get scored
    inner_join(
      mmilkup$omni.ref %>% filter(Assemblage == assemblage),
      by = 'Metric'
    ) %>%
    mutate(
      # This is how the scores were calculated in the old version of the package
      Score = case_when(
        StressResponse == 'inc' ~ (Value - Max) / (Min - Max),
        StressResponse == 'dec' ~ (Value - Min) / (Max - Min),
        TRUE ~ NA_real_
      )
    ) %>%
    select(
      SampleID, Metric, Score, RefCalMean
    )
  
  # Scored Metric dataframe is in a very convenient form to get the ASCI scores
  # We will just go ahead and calculate ASCI scores
  ASCI.scores <- scored.metrics %>%
    group_by(SampleID) %>%
    summarize(
      # Now, RefCalMean is a constant value for each assemblage type
      # ASCI is mean of scored metrics divided by the RefCalMean (a constant)
      # but remember, in the dataframe, the RefCalMean column had that same constant value for each row
      # dividing by the mean of the RefCalMean column will be the same as dividing by the constant value
      ASCI = mean(Score, na.rm = T) / mean(RefCalMean, na.rm = T)
    )
  
  # now for the final scored metrics output, we will want it to be wide.
  # we kept it in long format for the sake of easy calculation of ASCI, as you can see above
  # now we will pivot it wider because that's how we will need it to be from here on out
  scored.metrics <- scored.metrics %>%
    select(-RefCalMean) %>%
    # put the metric names back as column names
    pivot_wider(
      names_from = Metric, 
      values_from = Score
    ) %>%
    # We need to add on the scr suffix
    # keep in mind for eahc metric there is a mod score and a raw score
    rename_at(
      vars(
        names(.)[which(grepl("_mod|_raw",names(.)))]
      ),
      function(x) {
        paste0(x, "_scr")
      }
    )
  
  # Here we will bring the final output together and order the columns
  final.output <- metrics %>%
    inner_join(
      scored.metrics,
      by = 'SampleID'
    ) %>%
    inner_join(
      ASCI.scores,
      by = 'SampleID'
    ) %>%
    select(
      # Select SampleID, ASCI, and the rest of the columns, ordered.
      'SampleID',
      'ASCI',
      order(names(.))[which(!names(.) %in% c('SampleID','ASCI'))]
    ) %>%
    rename_at(
      vars(-SampleID),
      function(x) {
        # take the first letter of the assemblage keyword arg, and that will be the prefix
        #   for all the final column names (except SampleID of course)
        # That toupper function expression will return D for diatoms, S for sba and H for hybrid
        # assemblage keyword are can be one of "diatoms","sba", or "hybrid"
        paste0(
          toupper(substr(assemblage,1,1)), "_", x
        ) 
      }
    )
  
  return(final.output)
  
}
#' Re format ASCI output for SMC database
#'
#' @param df.wide The wide data that gets output from the ASCI function
#'
#' @return a data.frame long - the same ASCI data in long format to go to the SMC database tables
#' 
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom dplyr case_when select arrange mutate
#' 
#' @export
smcify <- function(df.wide) {
df.long <- df.wide %>%
  pivot_longer(
    cols = starts_with(c("D_","S_","H_")), 
    names_to = c("Metric"),
    values_to = c("Result")
  ) %>% 
  mutate(
    Assemblage = case_when(
      grepl("D_", Metric) ~ "Diatom",
      grepl("S_", Metric) ~ "SBA",
      grepl("H_", Metric) ~ "Hybrid",
      # There are different types of NA's in R
      # NA_character must be specified because this column is a character vector
      TRUE ~ NA_character_ 
    ) 
  ) %>%
  mutate(
    Metric_Type = case_when(
      grepl("_pct_att", Metric) ~ "Pct Att",
      grepl("_raw", Metric) ~ "Raw",
      grepl("_pred", Metric) ~ "Pred",
      grepl("_scr",Metric) ~ "Score",  
      TRUE ~ NA_character_
    )
  ) %>% 
  mutate(
    # Now we alter the Metric column
    # replace the first two characters (D_, S_, H_) with an empty string
    Metric = gsub("^.{0,2}","",Metric),
    
    # Now get rid of the suffix attached to the metric name 
    Metric = gsub("_pct_att|_pred|_raw|_scr","", Metric)
  ) %>%
  select(
    "SampleID","StationCode","SampleDate","Replicate",
    "SampleType","UnrecognizedTaxa","Assemblage",
    "Metric","Metric_Type","Result",
    "Comments","version_number"
  ) %>% 
  arrange(SampleID)
  
  return(df.long)
  
}
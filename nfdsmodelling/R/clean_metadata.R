#' Clean isolate metadata
#'
#' Harmonises serotype names, keeps one country, and derives age in years,
#' age group, collection type and vaccine type (PCV7 / PCV13 / NVT).
#'
#' @param meta Raw metadata, as read from metadata.csv.
#' @param serotype_renames Named character vector, old label = new label.
#'   Renamed isolates are treated as the new serotype from here on.
#' @param keep_country Country to keep.
#' @param child_cutoff Isolates with an age below this (in years) are "child".
#' @return One row per isolate.
#' @importFrom dplyr %>% mutate filter na_if coalesce case_when if_else
#'   add_count select
#' @importFrom readr parse_number
#' @export
clean_metadata <- function(meta, serotype_renames,
                           keep_country = "South Africa", child_cutoff = 6) {
  meta$serotype <- meta$In_Silico_Serotype

  hit <- meta$serotype %in% names(serotype_renames)
  meta$serotype[hit] <- unname(serotype_renames[meta$serotype[hit]])

  meta_clean <- meta %>%
    filter(Country == keep_country) %>%
    mutate(
      # Convert "." and other non-numeric values to NA
      Age_Year_num  = parse_number(na_if(as.character(Age_Year),  ".")),
      Age_Month_num = parse_number(na_if(as.character(Age_Month), ".")),
      Age_Day_num   = parse_number(na_if(as.character(Age_Day),   ".")),

      # Prioritize years, then months, then days
      age_years = coalesce(Age_Year_num, Age_Month_num / 12, Age_Day_num / 365.25),

      year = Col_Year,
      gpsc = as.character(new_GPSC),

      age_group = case_when(
        is.na(age_years)         ~ "unknown",
        age_years < child_cutoff ~ "child",
        TRUE                     ~ "adult"
      ),

      collection_type = if_else(Source == "NP", "carriage", "IPD")
    )

  pcv7_types  <- VACCINE_TYPES$PCV7
  pcv13_extra <- setdiff(VACCINE_TYPES$PCV13, VACCINE_TYPES$PCV7)  # the 6 serotypes PCV13 added

  meta_clean %>%
    mutate(vaccine_type = case_when(serotype %in% pcv7_types  ~ "PCV7",
                                    serotype %in% pcv13_extra ~ "PCV13",
                                    TRUE                      ~ "NVT")) %>%
    # add_count(serotype, gpsc, name = "n_sero_gpsc") %>%
    # add_count(serotype,       name = "n_sero") %>%
    select(
      Sample_accession, year, serotype, gpsc, vaccine_type,
      Age_Year, Age_Month, Age_Day,
      Age_Year_num, Age_Month_num, Age_Day_num,
      age_years, age_group, collection_type, Source
    )
}

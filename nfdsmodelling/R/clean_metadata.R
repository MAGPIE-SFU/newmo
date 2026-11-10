#' Clean isolate metadata
#'
#' Harmonises serotype names, keeps one country, and derives age in years,
#' age group, collection type and vaccine type (PCV7 / PCV13 / NVT).
#'
#' @param meta Raw metadata, as read from metadata.csv.
#' @param keep_country Country to keep.
#' @param child_cutoff Isolates with an age below this (in years) are "child".
#' @return One row per isolate.
#' @importFrom dplyr %>% mutate recode filter na_if coalesce case_when if_else
#'   add_count select
#' @importFrom readr parse_number
#' @export
clean_metadata <- function(meta, keep_country = "South Africa", child_cutoff = 6) {
  meta$serotype <- meta$In_Silico_Serotype

  meta <- meta %>%
    mutate(serotype = recode(serotype,
                             "11C"     = "11B",    # Nick has 11B only; group 11B/11C under that name
                             "10X"     = "10B",
                             "23B1"    = "23B",    # 23B1 is a 23B variant
                             "33A/33F" = "33F",
                             "15B/15C" = "15B/C",  # Nick's spelling
                             "20A"     = "20",     # 20 was split into 20A/20B after Nick's data
                             "20B"     = "20",
                             "35D"     = "35B"     # 35D isolates were typed as 35B until recently
    ))

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

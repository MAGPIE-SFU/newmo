# clean the data
# 29 July 2026
# 1. match the invasiveness data (disease progression rates)
#    to the serotypes in SA data. some of them don't match correctly.
#    see OL notes for rationale on the substitutions. Note: if i used the
#    rate for a different sero I created a col to indicate which one so the
#    original In_silico sero is maintained. new cols match the cols in nicks
#    dataset `invasiveness` and upper/lower. 
# 2. Gets the collection months from GPS/belman and adds them to the metadata.
# 3. Computes an age for each observation based on Age_Year, Age_Month, Age_Day
#    i also made an age indicator for adult/child
#    - `age_at_collection_years` (real) is the fraction of a year they are old
#       at the month of collection. ie, ceiling of that number is their "age"
#       unless they are <1yo.
#    - `child_under6` (logical) under 6 at their collection date
# 4. two new files are created
#    1. `meta_age_rates_clean.csv` the original data plus #1-3 above.
#    2. `clean_kids.csv` only the children for years around the carriage study.
#        see the data visualization script for more info on clean_kids.csv

library(dplyr)
library(readr)
library(stringr)


progression_rates_raw <- read_csv(
  "data/nick-inv-sts-children.csv",
  show_col_types = FALSE
)

south_africa_metadata_raw <- read_csv(
  "data/metadata.csv",
  show_col_types = FALSE
)



# ------------------------------
# clean nicks progression rates
# ------------------------------

# Clean and match childhood pneumo disease-progression rates
# to the South African pneumo metadata.
# note that nick has sperate ones for adults so we would have to
# get those and match them if we want to include adults.
# right now everyone gets the child rates



clean_serotype <- function(x) {
  x |>
    as.character() |>
    str_to_upper() |>
    str_trim() |>
    str_replace_all("\\s+", "") |>
    str_replace_all("_", "/") |>
    str_replace_all("-", "/") |>

    # Standardize grouped serotype labels
    str_replace("^15B/15C$", "15B/C") |>
    str_replace("^15C/15B$", "15B/C") |>
    str_replace("^15C/B$", "15B/C") |>

    # Remove possible prefixes
    str_remove("^SEROTYPE") |>
    str_remove("^ST")
}


progression_rates_clean <- progression_rates_raw |>
  transmute(
    progression_serotype_original = Serotype,
    serotype_progression = clean_serotype(Serotype),
    invasiveness,
    invasiveness_lower,
    invasiveness_upper
  ) |>
  distinct()


# Identify South African serotypes without progression-rate estimates
sa_serotypes <- south_africa_metadata_raw |>
  transmute(
    serotype_original = In_Silico_Serotype,
    serotype_clean = clean_serotype(In_Silico_Serotype)
  ) |>
  distinct()

progression_serotypes <- progression_rates_clean |>
  distinct(serotype_progression)

missing_progression_rates <- sa_serotypes |>
  anti_join(
    progression_serotypes,
    by = c("serotype_clean" = "serotype_progression")
  ) |>
  arrange(serotype_clean)

print(missing_progression_rates, n = Inf)


# These serotypes are excluded because no suitable progression-rate
# estimate was selected.
serotypes_to_exclude <- c(
  "10X",
  "11C",
  "25F",
  "UNTYPABLE",
  "SWISS/NT"
)

# These South African labels use the progression rate from a related
# serotype.
serotype_progression_crosswalk <- c(
  "6G"      = "6B",
  "35D"     = "35B",
  "33D"     = "33A",
  "24"      = "24F",
  "23B1"    = "23B",
  "20A"     = "20",
  "33A/33F" = "33F"
)


# ---------------------------
# assign matching serotypes
# ---------------------------

south_africa_metadata_clean <- south_africa_metadata_raw |>
  mutate(
    serotype_original = as.character(In_Silico_Serotype),
    serotype_clean = clean_serotype(In_Silico_Serotype)
  ) |>

  # Remove serotypes without an accepted progression-rate match
  filter(!serotype_clean %in% serotypes_to_exclude) |>

  # Assign the serotype label used for progression-rate matching
  mutate(
    serotype_progression = recode(
      serotype_clean,
      !!!serotype_progression_crosswalk,
      .default = serotype_clean
    )
  )


# Join progression rates onto the South African metadata
south_africa_with_progression <- south_africa_metadata_clean |>
  left_join(
    progression_rates_clean |>
      select(
        serotype_progression,
        invasiveness,
        invasiveness_lower,
        invasiveness_upper
      ),
    by = "serotype_progression"
  ) |>
  mutate(
    Source = factor(Source),
    Region = factor(Region),
    In_Silico_Serotype = factor(In_Silico_Serotype)
  )


# Validation checks

# Overall matching success
progression_match_summary <- south_africa_with_progression |>
  summarise(
    n_observations = n(),
    n_with_rate = sum(!is.na(invasiveness)),
    n_without_rate = sum(is.na(invasiveness)),
    proportion_matched = mean(!is.na(invasiveness))
  )

print(progression_match_summary)


# Remaining South African serotypes without a progression-rate estimate
unmatched_serotypes <- south_africa_with_progression |>
  filter(is.na(invasiveness)) |>
  count(
    serotype_original,
    serotype_clean,
    serotype_progression,
    sort = TRUE,
    name = "n_observations"
  )

print(unmatched_serotypes, n = Inf)


# Check the manually reassigned serotypes
manual_match_check <- south_africa_with_progression |>
  filter(
    serotype_clean %in% names(serotype_progression_crosswalk)
  ) |>
  count(
    serotype_original,
    serotype_clean,
    serotype_progression,
    name = "n_observations"
  ) |>
  arrange(serotype_clean)

print(manual_match_check, n = Inf)


# Confirm that excluded serotypes are absent
south_africa_with_progression |>
  filter(serotype_clean %in% serotypes_to_exclude) |>
  count(serotype_clean)

# save the cleaned data so far
south_africa_progression_data <- south_africa_with_progression



# I download all SA samples from GPS between 2000 and 2014 and
# the carriage data in the Belman paper.

#match gps data
months_mono <- read.csv("data/monocle-metadata.csv")

south_africa_with_month <- south_africa_metadata_raw |>
  mutate(
    Sample_accession = stringr::str_trim(
      as.character(Sample_accession)
    )
  ) |>
  left_join(
    months_mono |>
      transmute(
        ERS = stringr::str_trim(as.character(ERS)),
        Month
      ),
    by = c("Sample_accession" = "ERS")
  )

# check how many have months 
south_africa_with_month |>
  summarise(
    n_total = n(),
    n_with_month = sum(!is.na(Month)),
    n_without_month = sum(is.na(Month)),
    proportion_with_month = mean(!is.na(Month))
  )

# check how many matched the accessions
match_summary <- tibble(
  n_ERS = nrow(months_mono),
  n_accessions = nrow(south_africa_metadata_raw),
  n_matched = nrow(
    inner_join(
      months_mono,
      south_africa_metadata_raw,
      by = c("ERS" = "Sample_accession")
    )
  ),
  n_unmatched = nrow(
    anti_join(
      months_mono,
      south_africa_metadata_raw,
      by = c("ERS" = "Sample_accession")
    )
  )
)

match_summary # hmm only 12 not matched

# now match the belman data based on "lane id"
# -- why this doesnt have accession i have no idea

belman <- read.csv("data/belman.csv")
head(belman)


south_africa_metadata <- south_africa_metadata_raw |>
  mutate(
    Lane = str_trim(as.character(Lane))
  ) |>
  left_join(
    belman |>
      transmute(
        Lane = str_trim(as.character(Lane)),
        Col_Month
      ),
    by = "Lane"
  )

south_africa_metadata |>
  summarise(
    n_total = n(),
    n_with_month = sum(!is.na(Col_Month)),
    n_without_month = sum(is.na(Col_Month)),
    proportion_matched = mean(!is.na(Col_Month))
  )

# perfect now we have months for everything. but there is one without a month

south_africa_metadata |>
  filter(Col_Month == "")
#lets just call it Jan of the year. its sero 9N so not affected by vax
south_africa_metadata$Col_Month[
  south_africa_metadata$Col_Month == ""
] <- "JAN"

# now merge the months with the progression data
south_africa_data <- south_africa_progression_data |>
  left_join(
    south_africa_metadata |>
      select(Sample_accession, Col_Month),
    by = "Sample_accession"
  )

# ---------------------------
# AGE at collection
# ---------------------------

south_africa_data_age <- south_africa_data |>
  mutate(
    Age_Year_num = readr::parse_number(as.character(Age_Year)),
    Age_Month_num = readr::parse_number(na_if(as.character(Age_Month), ".")),
    Age_Day_num = readr::parse_number(na_if(as.character(Age_Day), ".")),

    age_at_collection_years = case_when(
      !is.na(Age_Year_num)  ~ Age_Year_num,
      !is.na(Age_Month_num) ~ Age_Month_num / 12,
      !is.na(Age_Day_num)   ~ Age_Day_num / 365.25,
      TRUE                  ~ NA_real_
    )
  ) |>
  select(-Age_Year_num, -Age_Month_num, -Age_Day_num)


south_africa_data_age |>
  filter(is.na(age_at_collection_years)) |>
  select(
    Sample_accession,
    Age_Year,
    Age_Month,
    Age_Day,
    Col_Year,
    Col_Month,
    Source,
    Region
  )

# there are 272 observations with no age data at all.
# we can assume they are all kids

# now lets get an indicator for being under 6 when their
# sample was taken (the 272 with no age can be flagged as under 6)

south_africa_data_age <- south_africa_data_age |>
  mutate(
    child_under6 = is.na(age_at_collection_years) |
      age_at_collection_years < 6
  )

# ok now we are nice and squeaky clean :)
write_csv(south_africa_data_age, "data/meta_age_rates_clean.csv")


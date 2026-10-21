library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)
library(lubridate) # may not be needed

###########
# --- readme
###########

# this script refines the metadata and NO LONGER creates the group column for the South Africa data.
# it also groups rare things in a way I hope is sensible. (it wasn't, but for subtle reasons kind of) 

# outputs:
#   data/SouthAfrica/meta_clean.csv  one row per isolate; serotype and gpsc are the real
#                                    values, group is the (coarsened) unit used for x(t),
#                                    vaccine_type is PCV7 / PCV13 / NVT


###########
# ----- read the metadata csv file and harmonise serotype names
###########

meta <- read_csv("data/SouthAfrica/metadata.csv") 

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
glimpse(meta)


###########
# ------- ages, sample type, year; South Africa 
###########

child_cutoff <- 6      # 5 and under will be kids, consistent with previous work

## fixing the age year thing -- basically, age is recorded in number of years or number of days sometimes, those can be converted to year. there are still some missing

meta_clean <- meta %>%
    filter(Country == "South Africa") %>%
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


###########
# ------- vaccine type and grouping
###########
# Rule (raw isolate counts over the whole South Africa dataset, all years and sample types):
#   * a serotype-GPSC pair keeps its own group if it has >= min_gpsc_n isolates,
#     otherwise it becomes <serotype>_other
#   * non-vaccine serotypes with < min_serotype_n isolates in total are pooled into one
#     group, "rare_NVT". Vaccine serotypes are never pooled.
# serotype and gpsc keep their real values, so invasiveness is still applied per serotype
# in script 2; only `group` is coarsened.
# (previous rule: serotype-GPSC pairs with >= 10 isolates, others -> <serotype>_other,
#  and gpsc itself was overwritten with "other" but this is better - 80+, but not 140+, groups to model)
# with fewer groups that have 0 counts in some of the years - which could create problems for the modellign later 


# min_gpsc_n     <- 20   # isolates needed for a serotype-GPSC pair to be its own group
# min_serotype_n <- 40   # NVT serotypes with fewer isolates than this go to "rare_NVT"

pcv7_types  <- c("4", "6B", "9V", "14", "18C", "19F", "23F")
pcv13_extra <- c("1", "3", "5", "6A", "7F", "19A")   # the 6 serotypes PCV13 added to PCV7

meta_clean <- meta_clean %>%
    mutate(vaccine_type = case_when(serotype %in% pcv7_types  ~ "PCV7",
                                    serotype %in% pcv13_extra ~ "PCV13",
                                    TRUE                      ~ "NVT")) %>%
    add_count(serotype, gpsc, name = "n_sero_gpsc") %>%
    add_count(serotype,       name = "n_sero") %>%
    select(
        Sample_accession,
        year,
        serotype,
        gpsc,
#        group,
        vaccine_type,
        Age_Year,
        Age_Month,
        Age_Day,
        Age_Year_num,
        Age_Month_num,
        Age_Day_num,
        age_years,
        age_group,
        collection_type,
        Source
    )



# every group belongs to exactly one vaccine type (so VT / NVT can always be recovered)
# stopifnot(all(count(distinct(meta_clean, group, vaccine_type), group)$n == 1))

# message(n_distinct(meta_clean$group), " groups; ",
  #      sum(meta_clean$group == "rare_NVT"), " isolates in rare_NVT (",
  #      n_distinct(meta_clean$serotype[meta_clean$group == "rare_NVT"]), " serotypes); ",
  #     sum(grepl("_other$", meta_clean$group)), " isolates in <serotype>_other groups")


###########
# ------- checks
###########

# how many ages of each type. there are still NAs (272 with no age info at all)
meta_clean %>%
    summarise(
        year_only   = sum(!is.na(Age_Year_num) &  is.na(Age_Month_num) &  is.na(Age_Day_num)),
        month_only  = sum( is.na(Age_Year_num) & !is.na(Age_Month_num) &  is.na(Age_Day_num)),
        day_only    = sum( is.na(Age_Year_num) &  is.na(Age_Month_num) & !is.na(Age_Day_num)),
        multiple    = sum(rowSums(cbind(!is.na(Age_Year_num), !is.na(Age_Month_num), !is.na(Age_Day_num))) >  1),
        all_missing = sum(rowSums(cbind(!is.na(Age_Year_num), !is.na(Age_Month_num), !is.na(Age_Day_num))) == 0)
    )

# isolates by vaccine type and sample type
meta_clean %>% count(vaccine_type, collection_type) %>%
    pivot_wider(names_from = collection_type, values_from = n)


###########
# ------- save
###########

write_csv(meta_clean, "data/SouthAfrica/meta_clean.csv")



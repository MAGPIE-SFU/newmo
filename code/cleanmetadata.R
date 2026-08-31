library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)
library(lubridate) # may not be needed 

###########
# --- readme  
###########


# this script refines the metadata and creates the sc-st group column for the South Africa data 


###########
# ----- read the data do a bit of cleaning, checking, merge some serotypes etc 
###########

G <- as.matrix(read.csv("data/SouthAfrica/Ggrouped.csv", row.names=1))

meta <- read_csv("data/SouthAfrica/metadata.csv") # all of these are in the same order, though 

meta$serotype <- meta$In_Silico_Serotype
meta <- meta %>%
    mutate(serotype = recode(serotype,
                             "11B" = "11BC",
                             "11C" = "11BC",
                             "10X" = "10B",
                             "23B1" = "23",
                             "33A/33F" = "33F"  ))

glimpse(meta)

###########
# ------- create meta_clean, with group information, removing rare groups 
###########

child_cutoff <- 6      # 5 and under will be kids, consistent with previous work

# what are the sc-st combinations and how many isolates were there ever? 
group_totals <- meta %>%
    mutate(
        year = Col_Year,
        gpsc = new_GPSC
    ) %>%
    filter(Country == "South Africa") %>%
    count(serotype, gpsc, name = "n_total")
sum(group_totals$n_total<=10) # there are 253 groups with <= 10 observations ever 
# 
sum(group_totals$n_total[group_totals$n_total<=10]) # and they comprise about 1/10 of the data .
# it would be better to include them somehow .. merging some perhaps by their accessory genomes


common_groups <- group_totals %>%
    filter(n_total >= 10) %>%
    select(serotype, gpsc)

# FIX FIX ! PLEASE -- unknown Age_Year may be a known age given in months
#meta_clean <- meta %>%
#    mutate(
#        year = Col_Year,
#        gpsc = new_GPSC,
#        age_group = case_when(
#            is.na(Age_Year) ~ "unknown",
#            Age_Year < child_cutoff ~ "child",
#            TRUE ~ "adult"
#        ),
#        collection_type = if_else(Source == "NP", "carriage", "IPD")
#    ) %>%
#    filter(Country == "South Africa") %>%
#    left_join(
#        common_groups %>% mutate(keep_group = TRUE),
#        by = c("serotype", "gpsc")
#    ) %>%
#    mutate(
#        gpsc = if_else(is.na(keep_group), "other", as.character(gpsc)),
#        group = paste(serotype, gpsc, sep = "_")
#    ) %>%
#    select(
#        Sample_accession, year, serotype, gpsc, group,
#        Age_Year, age_group, collection_type, Source
#    )



## fixing the age year thing from above -- basically, age is recorded in number of years or number of days sometimes, those can be converted to year. there are still some missing 

meta_clean <- meta %>%
    mutate(
        # Convert "." and other non-numeric values to NA
        Age_Year_num = parse_number(
            na_if(as.character(Age_Year), ".")
        ),
        Age_Month_num = parse_number(
            na_if(as.character(Age_Month), ".")
        ),
        Age_Day_num = parse_number(
            na_if(as.character(Age_Day), ".")
        ),

        # Prioritize years, then months, then days
        age_years = coalesce(
            Age_Year_num,
            Age_Month_num / 12,
            Age_Day_num / 365.25
        ),

        year = Col_Year,
        gpsc = new_GPSC,

        age_group = case_when(
            is.na(age_years)       ~ "unknown",
            age_years < child_cutoff ~ "child",
            TRUE                   ~ "adult"
        ),

        collection_type = if_else(
            Source == "NP",
            "carriage",
            "IPD"
        )
    ) %>%
    filter(Country == "South Africa") %>%
    left_join(
        common_groups %>%
            mutate(keep_group = TRUE),
        by = c("serotype", "gpsc")
    ) %>%
    mutate(
        gpsc = if_else(
            is.na(keep_group),
            "other",
            as.character(gpsc)
        ),
        group = paste(serotype, gpsc, sep = "_")
    ) %>%
    select(
        Sample_accession,
        year,
        serotype,
        gpsc,
        group,
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

all(rownames(G) == meta_clean$Sample_accession) # TRUE good 

# check how many ages of each type. there are still NAs (272 with no age info at all)
meta_clean %>%
    summarise(
        year_only = sum(
            !is.na(Age_Year_num) &
            is.na(Age_Month_num) &
            is.na(Age_Day_num)
        ),

        month_only = sum(
            is.na(Age_Year_num) &
            !is.na(Age_Month_num) &
            is.na(Age_Day_num)
        ),

        day_only = sum(
            is.na(Age_Year_num) &
            is.na(Age_Month_num) &
            !is.na(Age_Day_num)
        ),

        multiple = sum(
            rowSums(
                cbind(
                    !is.na(Age_Year_num),
                    !is.na(Age_Month_num),
                    !is.na(Age_Day_num)
                )
            ) > 1
        ),

        all_missing = sum(
            rowSums(
                cbind(
                    !is.na(Age_Year_num),
                    !is.na(Age_Month_num),
                    !is.na(Age_Day_num)
                )
            ) == 0
        )
    )



 write_csv(meta_clean, "data/SouthAfrica/meta_clean.csv") 
#




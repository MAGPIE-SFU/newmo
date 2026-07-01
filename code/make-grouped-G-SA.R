library(dplyr)
library(ggplot2)
library(tidyr)
library(lubridate) # may not be needed 

###########
# --- readme  
###########


# this script creates a grouped version of the matrix G for the South Africa data 
# it uses the metadata and a "raw" matrix G. 
# In a "raw" matrix G, each row corresponds to an isolate, and the corresponding 
# metadata tells us information about that isolate.
# the elements of the raw G are 0 if isolate i has the gene at locus l and 0 otherwise. 
# Graw's entries all 0 and 1s 

# in a grouped matrix G, each row corresponds to a gpsc-st combination. 
# the elements of the grouped G, Gjl, are the fraction of isolates with the 
# gpsc and serotype combo that define strain j that have the gene at locus l 
# Gg entries are in [0,1]. 

###########
# ----- read the data, reorder G, do a bit of cleaning, checking, merge some serotypes etc 
###########


G <- as.matrix(read.csv("data/SouthAfrica/gene_matrix.csv",row.names=1) )
meta <- read_csv("data/SouthAfrica/metadata.csv") # all of these are in the same order, though 
ind = match(meta$Sample_accession, rownames(G)) # no NAs, good 

G = as.matrix(G[match(meta$Sample_accession, rownames(G)),]) # G, reordered . 
all(rownames(G) == meta$Sample_accession) # should be TRUE, good 

meta$serotype= meta$In_Silico_Serotype
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

meta_clean <- meta %>%
    mutate(
        year = Col_Year,
        gpsc = new_GPSC,
        age_group = case_when(
            is.na(Age_Year) ~ "unknown",
            Age_Year < child_cutoff ~ "child",
            TRUE ~ "adult"
        ),
        collection_type = if_else(Source == "NP", "carriage", "IPD")
    ) %>%
    filter(Country == "South Africa") %>%
    left_join(
        common_groups %>% mutate(keep_group = TRUE),
        by = c("serotype", "gpsc")
    ) %>%
    mutate(
        gpsc = if_else(is.na(keep_group), "other", as.character(gpsc)),
        group = paste(serotype, gpsc, sep = "_")
    ) %>%
    select(
        Sample_accession, year, serotype, gpsc, group,
        Age_Year, age_group, collection_type, Source
    )

all(rownames(G) == meta_clean$Sample_accession) # TRUE good 

# write_csv(meta_clean, "data/SouthAfrica/meta_clean.csv") 
# I left this line here because this script explains why meta_clean differs from the original metadata 



###########
# ----- now create the grouped version of G ------ 
###########


group_ids <- meta_clean$group

G_sum <- rowsum(G, group = group_ids)
groupcounts <- table(group_ids)[rownames(G_sum)]

# We want Gg[j,l] to be the mean prevalence of each locus; values in [0,1]
Gg <- G_sum / as.numeric(groupcounts)

# a spot check: 
hist(Gg[sample(1:nrow(Gg), 1), ])

# a visualization (but not a useful one) 
heatmap(Gg)

# write.csv(Gg, "data/SouthAfrica/Ggrouped.csv")




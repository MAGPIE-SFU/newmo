# --- readme ----
# this file creates the grouped G matrix (Ggrouped.csv), with columns sorted according to 
# the stability of the locus. The most stable columns will be first. 

# In a "raw" matrix G, each row corresponds to an isolate, and the corresponding 
# metadata tells us information about that isolate.
# the elements of the raw G are 0 if isolate i has the gene at locus l and 0 otherwise. 
# Graw's entries all 0 and 1s 

# in a grouped matrix G, each row corresponds to a gpsc-st combination. 
# the elements of the grouped G, Gjl, are the fraction of isolates with the 
# gpsc and serotype combo that define strain j that have the gene at locus l 
# Gg entries are in [0,1]. The loci are in order of increasing "delta statistic". 

library(readr)
library(dplyr)
library(tibble)

# G <- as.matrix(read.csv("data/SouthAfrica/gene_matrix.csv",row.names=1) )
G <- as.matrix(read.csv("data/SouthAfrica/pa_matrix_final.csv",row.names=1) )
meta_clean <- read_csv("data/SouthAfrica/meta_clean.csv", show_col_types = FALSE)

# reorder rows of G 
G <- as.matrix(G[match(meta_clean$Sample_accession, rownames(G)),]) # G, reordered .

all(rownames(G) == meta_clean$Sample_accession) # should be true 

###########
# ----- now create the grouped version of G ------ 
###########

group_ids <- meta_clean$group

G_sum <- rowsum(G, group = group_ids)
groupcounts <- table(group_ids)[rownames(G_sum)]

# We want Gg[j,l] to be the mean prevalence of each locus; values in [0,1]
Gg <- G_sum / as.numeric(groupcounts)

# a spot check: 
# hist(Gg[sample(1:nrow(Gg), 1), ])

# a visualization (but not a useful one) 
# heatmap(Gg)

# 

###########
# ----- now reorder the columns of G using the locus frequency data 
###########

locus_by_year <- read_csv("data/SouthAfrica/locus_freqs_deltas.csv", show_col_types = FALSE) |>
    rename(locus = `...1`) |>
    column_to_rownames("locus")

# want to order Gg's columns from lowest to highest delta 
# (Note that this order checking is redundant now, because i went back and sorted them 
# for the locus_freqs_deltas.csv file). 
locsinorder <- rownames(locus_by_year)[order(locus_by_year$delta)]


stopifnot(all(!is.na(match(locsinorder, colnames(Gg)))))

Greord <- Gg[, match(locsinorder, colnames(Gg))]

stopifnot(identical(colnames(Greord), locsinorder))


# sanity check: (same for tail; verifies) 
# locus_by_year$delta[ match(head(locsinorder), rownames(locus_by_year))][1:5]

Greord <- Gg[, match(locsinorder, colnames(Gg))]
all(colnames(Greord)==locsinorder) # true. 

write.csv(Greord, "data/SouthAfrica/G_cols_ordered.csv")

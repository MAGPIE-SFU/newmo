# --- readme ----
# this file creates the grouped G matrix (Ggrouped.csv), with columns sorted according to 
# the stability of the locus. The most stable columns will be first. 

# Inputs: the gene matrix (0 and 1, with accessions that match the metadata) 
#         the metadata (meta_clean.csv) 
#         the xt data (xt_grouped). this is the synthetic carriage population and it 
#         gives us the gene frequencies, which we use to order the columns of G, for the NFDS model
# Output: G for the NFDS model, with columns in order of stability (stablest first) and 
#         i think with only those genes that are at intermediate frequency. 


# what is G? 
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
xt <- read_csv("data/SouthAfrica/xt_grouped.csv", show_col_types = FALSE) 


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



###########################################################################
# locus frequencies : compute them , then re-order G's columns accordingly 
###########################################################################
years <- sort(unique(xt$year))

locus_by_year <- sapply(years, function(y) {
  p_df  <- xt %>% filter(year == y) %>% select(group, frac)
  p_vec <- setNames(p_df$frac, p_df$group)[rownames(Gg)]
  p_vec[is.na(p_vec)] <- 0
  as.numeric(p_vec %*% as.matrix(Gg))
})

rownames(locus_by_year) <- colnames(Gg)
colnames(locus_by_year) <- paste("t", years, sep = "")
locus_by_year = as.data.frame(locus_by_year)

# here is a locus frequency plot comparing two times 
ggplot(locus_by_year, aes(x=t2000, y=t2014)) + geom_point(alpha=0.3)
# interesting that a few points are very dark, why? for example the lowest one after 0.25 - 


# Corander et al defined Delta as (f - e)^2 / (1 - e(1-e)) 
# where e is the frequency before vaccination and f is the frequency after vaccination 
# (all pre- and during vaccination samples for e, all post-vaccination samples for f). 

# here, to do approximately the same, i would take the pre-2009 and post-2011 
# mean frequencies of these 

# NOTE ALSO - need to figure out what to do about filtering for only intermediate
# frequency genes. The lowest values of delta may be for very rare genes since 
# the top quantity is squared. That is unlikely to be what we want to choose 
# for the NFDS model. 

locus_by_year <- locus_by_year %>%
  mutate(
    pre2009  = rowMeans(across(num_range("t", 2000:2008))),
    post2011 = rowMeans(across(num_range("t", 2012:2014))),
    delta = (post2011 - pre2009)^2/(1-pre2009*(1-pre2009))
  )

hist(locus_by_year$delta,breaks = 50)
ggplot(filter(locus_by_year, (pre2009 < 0.95 & pre2009 > 0.05)), aes(x=log10(pre2009)))+geom_histogram()
ggplot(locus_by_year,aes( x=pre2009, y = delta) ) + geom_point(alpha=0.3) + scale_x_log10() +
  scale_y_log10()
# no clear correlation overall  
# between frequency pre-vaccination and delta. but the VERY low deltas are very low-freq genes 

# ---- 
# what about not using delta? 
#### 
# i am not convince that i like delta. it will be low where f and e are both low. 
# went down rabbit hole trying to find evidence for NFDS in this dataset. Failed. 

# what if we think about the variance over years,
# but we make it less dependent on the frequency? 
n_year   <- table(meta_clean$year)
n_pre    <- sum(n_year[as.character(2000:2008)])
n_post   <- sum(n_year[as.character(2012:2014)])
n_allyrs <- as.numeric(n_year[as.character(2000:2014)])
n_preyrs <- as.numeric(n_year[as.character(2000:2008)])

rowvar <- function(m) { m <- as.matrix(m); rowSums((m - rowMeans(m))^2) / (ncol(m) - 1) }
rowexp <- function(m, n) { m <- as.matrix(m); rowMeans(sweep(m * (1 - m), 2, n, "/")) }

locus_by_year <- locus_by_year %>%
  mutate(
    locus    = rownames(locus_by_year),          # keep gene names through mutate
    pre2009  = rowMeans(across(num_range("t", 2000:2008))),
    post2011 = rowMeans(across(num_range("t", 2012:2014))),
    delta    = (post2011 - pre2009)^2 / (1 - pre2009 * (1 - pre2009)),
    
    # standardised pre-to-post change: how many standard errors the locus moved
    z = (post2011 - pre2009) /
      sqrt(pre2009 * (1 - pre2009) / n_pre + post2011 * (1 - post2011) / n_post),
    
    # year-to-year variation against sampling noise, all 15 years
    obs_var  = rowvar(across(num_range("t", 2000:2014))),
    exp_var  = rowexp(across(num_range("t", 2000:2014)), n_allyrs),
    overdisp = obs_var / exp_var,
    
    # the same, pre-vaccine only -- the noise-model check, should sit near 1
    obs_pre      = rowvar(across(num_range("t", 2000:2008))),
    exp_pre      = rowexp(across(num_range("t", 2000:2008)), n_preyrs),
    overdisp_pre = obs_pre / exp_pre
  )

all(rownames(locus_by_year) == colnames(Gg)) # TRUE 

locus_by_year <- locus_by_year %>%
  mutate(in_model = pre2009 > 0.05 & pre2009 < 0.95 & !is.na(z)) %>%
  arrange(abs(z))                    # NaN z sorts to the bottom

Greord <- Gg[, locus_by_year$locus[locus_by_year$in_model]]

stopifnot(identical(colnames(Gg_ordered),
                    locus_by_year$locus[locus_by_year$in_model]))



write_csv(locus_by_year, "data/SouthAfrica/locus_freqs_deltas_zs.csv")

write.csv(Greord, "data/SouthAfrica/G_cols_ordered.csv")

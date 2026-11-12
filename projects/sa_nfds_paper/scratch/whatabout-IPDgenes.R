# what about IPD associated genes ?

# Mantel-Haenszel

# The within-group invasive-vs-carriage
# log odds ratio per locus (43 groups), with logOR_ngrp for how
# many lineages there were for that group. Positive = more common in invasive isolates.

# Within-group invasive-vs-carriage log-odds ratio per locus (Mantel-Haenszel),
# pooled over groups with at least min_n isolates of each type.
#   G        isolates x loci, 0/1, rows aligned to `group` and `is_ipd`
#   returns  data.frame(locus, logOR, ngrp)
#            logOR > 0 = gene more common in invasive isolates of the same lineage
#            ngrp      = lineages polymorphic enough to contribute (low = weak estimate)

# need to read in G (the raw one, not grouped, using the isolate level here) and the metadata
# and the locus by year (we don't strictly need this but it has all the other locus level stuff in it
# so it's a natural place to add this)

meta_clean <- read_csv("data/SouthAfrica/meta_clean.csv")

G <- as.matrix(G[match(meta_clean$Sample_accession, rownames(G)),]) # G, reordered .
all(rownames(G) == meta_clean$Sample_accession) # should be true

locus_by_year <- read.csv("data/SouthAfrica/locus_freqs_deltas.csv", check.names = FALSE)
rownames(locus_by_year) <- locus_by_year$locus      # only if something downstream wants them



ipd_logOR <- function(G, group, is_ipd, min_n = 10) {
  lev  <- sort(unique(group))
  ni   <- table(factor(group[is_ipd],  levels = lev))
  nc   <- table(factor(group[!is_ipd], levels = lev))
  keep <- lev[ni >= min_n & nc >= min_n]
  stopifnot(length(keep) > 0)

  prev <- function(rows) {                       # prevalence per kept lineage
    gr <- group[rows]
    s  <- rowsum(G[rows, , drop = FALSE], gr)
    s[keep, , drop = FALSE] / as.numeric(table(gr)[keep])
  }
  pI <- prev(is_ipd); pC <- prev(!is_ipd)
  nI <- as.numeric(ni[keep]); nC <- as.numeric(nc[keep])

  a <- pI * nI + 0.5; b <- (1 - pI) * nI + 0.5   # +0.5 continuity correction
  d <- pC * nC + 0.5; e <- (1 - pC) * nC + 0.5
  n <- a + b + d + e

  data.frame(locus = colnames(G),
             logOR = log(colSums(a * e / n) / colSums(b * d / n)),
             ngrp  = colSums((pI*nI + pC*nC) >= 1 & ((1-pI)*nI + (1-pC)*nC) >= 1),
             row.names = NULL)
}

lor <- ipd_logOR(G, meta_clean$group, meta_clean$collection_type == "IPD")

locus_by_year <- locus_by_year %>%
  mutate(
    locus    = rownames(locus_by_year),
    ... ,
    logitdiff = abs(logit(pre2009) - logit(post2011)),

    logOR_ipd  = lor$logOR[match(locus, lor$locus)],
    logOR_ngrp = lor$ngrp [match(locus, lor$locus)]
  )

# i did not do anything with this and full disclosure claude wrote the function above
# and i did not go into its details

# would doing this without also having the serotype groups be interesting?
# are there genes that are IPD associated in some STs but not others, and is that interesting?

# --- readme ----
# this creates the grouped G matrix, with columns sorted according to
# the stability of the locus. The most stable columns will be first.

# Inputs: the gene matrix (0 and 1, with accessions that match the metadata)
#         the metadata with groups (meta_clean_groups)
#         the xt data (xt_grouped). this is the synthetic carriage population and it
#         gives us the gene frequencies, which we use to order the columns of G, for the NFDS model
# Output: G for the NFDS model, with columns in order of stability (stablest first) and
#         only those genes that are at intermediate frequency.


# what is G?
# In a "raw" matrix G, each row corresponds to an isolate, and the corresponding
# metadata tells us information about that isolate.
# the elements of the raw G are 1 if isolate i has the gene at locus l and 0 otherwise.
# Graw's entries all 0 and 1s

# in a grouped matrix G, each row corresponds to a gpsc-st combination.
# the elements of the grouped G, Gjl, are the fraction of isolates with the
# gpsc and serotype combo that define strain j that have the gene at locus l
# Gg entries are in [0,1]. The loci are in order of increasing "delta statistic".

# why is the NFDS G different from some other G?
# (1) it only has intermediate frequency genes (min_freq, max_freq)
# (2) these frequencies are computed from the synthetic carriage sample xt
# (3) G's columns are in order of decreasing stability: the first K are the ones to use in NFDS
# (4) G's rows line up with the rows of xt. each row is a ST-SC combination and G's entries are
# the fraction of the group that has that column's locus
# a "general" G is a presence-absence matrix in general (rows could be isolates, groups, etc)

#' Build the grouped G matrix for the NFDS model
#'
#' Groups the isolate gene presence/absence matrix by model group, computes
#' locus frequencies per year from x(t), computes several locus stability
#' statistics, keeps intermediate-frequency loci and orders them by `logitdiff`.
#'
#' @param G Raw presence/absence matrix, one row per isolate (row names are
#'   sample accessions), one column per locus.
#' @param meta_clean Metadata with a `group` column, e.g.
#'   `reweight_carriage()$meta_clean_groups`.
#' @param xt Grouped x(t), e.g. from [reweight_carriage_uncertainty()].
#' @param pre_years Pre-vaccine years (the "pre2009" frequency).
#' @param post_years Post-vaccine years (the "post2011" frequency).
#' @param all_years All years used for the year-to-year variance.
#' @param min_freq,max_freq Loci with a pre-vaccine frequency outside this
#'   range are left out of the model.
#'
#' @return A list with `G` (the grouped, filtered and ordered G for the NFDS
#'   model) and `locus_by_year` (locus frequencies per year and the statistics,
#'   with an `in_model` column).
#'
#' @importFrom dplyr %>% filter select mutate across num_range arrange rowwise
#'   c_across all_of ungroup
#' @importFrom stats setNames var
#' @export
build_gnfds <- function(G, meta_clean, xt,
                        pre_years, post_years, all_years,
                        min_freq = 0.05, max_freq = 0.95) {

  # reorder rows of G
  G <- as.matrix(G[match(meta_clean$Sample_accession, rownames(G)),]) # G, reordered .

  stopifnot(all(rownames(G) == meta_clean$Sample_accession)) # should be true

  ###########
  # ----- now create the grouped version of G ------
  ###########

  group_ids <- meta_clean$group

  G_sum <- rowsum(G, group = group_ids)
  groupcounts <- table(group_ids)[rownames(G_sum)]

  # We want Gg[j,l] to be the mean prevalence of each locus; values in [0,1]
  Gg <- G_sum / as.numeric(groupcounts)

  grp_order <- xt$group[xt$year == min(xt$year)]     # group order as xt has it
  Gg <- Gg[grp_order, , drop = FALSE]
  stopifnot(identical(rownames(Gg), grp_order))

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

  # Frequencies are sums of group shares; floating-point rounding can push loci
  # that every group carries just above 1 (e.g. 1 + 2e-16). Clamp to [0, 1].
  locus_by_year <- pmin(pmax(locus_by_year, 0), 1)

  rownames(locus_by_year) <- colnames(Gg)
  colnames(locus_by_year) <- paste("t", years, sep = "")
  locus_by_year = as.data.frame(locus_by_year)

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
      pre2009  = rowMeans(across(num_range("t", pre_years))),
      post2011 = rowMeans(across(num_range("t", post_years))),
      delta = (post2011 - pre2009)^2/(1-pre2009*(1-pre2009))
    )

  # looked into correlations
  # between frequency pre-vaccination and delta and the other measures.
  # see locus-freq-explore.R script.
  # but the VERY low deltas are  low-freq genes, and there is one at frequency 1 oddly

  # ----
  # what about not using delta?
  # i am not convince that i like delta. it will be low where f and e are both low.
  # went down rabbit hole trying to find evidence for NFDS in this dataset. Failed for now.

  # what if we think about the variance over years,
  # but we make it less dependent on the frequency?
  n_year   <- table(meta_clean$year)
  n_pre    <- sum(n_year[as.character(pre_years)])
  n_post   <- sum(n_year[as.character(post_years)])
  n_allyrs <- as.numeric(n_year[as.character(all_years)])
  n_preyrs <- as.numeric(n_year[as.character(pre_years)])

  rowvar <- function(m) { m <- as.matrix(m); rowSums((m - rowMeans(m))^2) / (ncol(m) - 1) }
  rowexp <- function(m, n) { m <- as.matrix(m); rowMeans(sweep(m * (1 - m), 2, n, "/")) }

  logit <- function(x) log(x/(1-x))

  locus_by_year <- locus_by_year %>%
    mutate(
      locus    = rownames(locus_by_year),          # keep gene names through mutate
      pre2009  = rowMeans(across(num_range("t", pre_years))),
      post2011 = rowMeans(across(num_range("t", post_years))),
      delta    = (post2011 - pre2009)^2 / (1 - pre2009 * (1 - pre2009)),

      # standardised pre-to-post change: how many standard errors the locus moved
      z = (post2011 - pre2009) /
        sqrt(pre2009 * (1 - pre2009) / n_pre + post2011 * (1 - post2011) / n_post),

      # year-to-year variation compared to whta is expected due to sampling noise with n-ish samples, all years
      obs_var  = rowvar(across(num_range("t", all_years))),
      exp_var  = rowexp(across(num_range("t", all_years)), n_allyrs),
      overdisp = obs_var / exp_var,

      # the same, pre-vaccine only -- kind of a noise model: variance compared to expected variance
      obs_pre      = rowvar(across(num_range("t", pre_years))),
      exp_pre      = rowexp(across(num_range("t", pre_years)), n_preyrs),
      overdisp_pre = obs_pre / exp_pre,

      # logit comparison - new choice for how to rank the genes!
      logitdiff  = abs(logit(pre2009)-logit(post2011)) ,
      # log ratio - another choice, similar to logitdiff but not quite
      lograt = abs(log(pre2009/post2011))
    )

  stopifnot(all(rownames(locus_by_year) == colnames(Gg))) # TRUE

  locus_by_year <- locus_by_year %>%
    mutate(in_model = pre2009 > min_freq & pre2009 < max_freq & !is.na(z)) %>%
    arrange(logitdiff)                    # NaN z sorts to the bottom

  Greord <- Gg[, locus_by_year$locus[locus_by_year$in_model]]

  stopifnot(identical(colnames(Greord),
                      locus_by_year$locus[locus_by_year$in_model]))

  # --- adding some extra info to locus_by_year for pre-vaccine gene orderings
  # to see if they work in the model
  pre_cols <- paste0("t", pre_years)
  locus_by_year <- locus_by_year %>%
      rowwise() %>%
      mutate(varlogitpre = var(logit(c_across(all_of(pre_cols))))) %>%
      ungroup()

  list(
    G             = Greord,
    locus_by_year = locus_by_year
  )
}

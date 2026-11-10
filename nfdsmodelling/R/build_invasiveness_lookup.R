# ==============================================================================
# Build invasiveness lookup with intervals
# ==============================================================================

# ------------------------------------------------------------------------------
# Empirical Bayes shrinkage on the log scale
# DerSimonian-Laird estimate of between-serotype variance
# ------------------------------------------------------------------------------

#' Empirical Bayes shrinkage of invasiveness estimates
#'
#' Shrinks log invasiveness towards the overall mean, using a DerSimonian-Laird
#' estimate of the between-serotype variance.
#'
#' @param tab Table with columns Serotype, invasiveness, invasiveness_lower and
#'   invasiveness_upper.
#' @return `tab` with added columns logmean, logsd and v_shrunk, and attributes
#'   overall_mean, mu and tau.
#' @importFrom dplyr %>% mutate
#' @export
shrink_invasiveness <- function(tab) {

  required <- c(
    "Serotype",
    "invasiveness",
    "invasiveness_lower",
    "invasiveness_upper"
  )

  missing <- setdiff(required, names(tab))

  if (length(missing)) {
    stop(
      "Missing columns: ",
      paste(missing, collapse = ", ")
    )
  }

  tab <- tab %>%
    mutate(Serotype = trimws(as.character(Serotype)))

  if (anyDuplicated(tab$Serotype)) {
    stop("Duplicate serotypes in invasiveness table")
  }

  valid <- with(
    tab,
    is.finite(invasiveness) &
      invasiveness > 0 &
      is.finite(invasiveness_lower) &
      invasiveness_lower > 0 &
      is.finite(invasiveness_upper) &
      invasiveness_upper > invasiveness_lower
  )

  if (!all(valid)) {
    stop("Invalid invasiveness estimates or CI endpoints")
  }

  y <- log(tab$invasiveness)

  se <- (
    log(tab$invasiveness_upper) -
      log(tab$invasiveness_lower)
  ) / (2 * 1.96)

  w <- 1 / se^2

  mu0 <- sum(w * y) / sum(w)

  Q <- sum(w * (y - mu0)^2)

  denom <- sum(w) - sum(w^2) / sum(w)

  tau2 <- if (denom > 0) {
    max(0, (Q - (length(y) - 1)) / denom)
  } else {
    0
  }

  mu <- sum(y / (se^2 + tau2)) /
    sum(1 / (se^2 + tau2))

  B <- tau2 / (tau2 + se^2)

  message(
    sprintf(
      "Shrinkage: exp(mu) = %.3e, tau = %.3f",
      exp(mu), sqrt(tau2)
    )
  )

  result <- tab %>%
    mutate(
      logmean  = mu + B * (y - mu),   # posterior mean of log v
      logsd    = sqrt(B) * se,        # posterior sd of log v
      v_shrunk = exp(mu + B * (y - mu))
    )

  attr(result, "overall_mean") <- exp(mu)
  attr(result, "mu")           <- mu
  attr(result, "tau")          <- sqrt(tau2)
  result
}

# ------------------------------------------------------------------------------
# Extract serogroup
# Examples: 6A -> 6, 33D -> 33, 24F -> 24
# ------------------------------------------------------------------------------

serogroup_of <- function(s) {
  sub("^([0-9]+).*$", "\\1", as.character(s))
}

# ------------------------------------------------------------------------------
# Build lookup for one age category
# ------------------------------------------------------------------------------

#' @importFrom dplyr bind_rows tibble
#' @importFrom stats setNames
lookup_one_age <- function(tab, other_tab, serotypes,
                           serotype_substitutions, non_typeable) {

  t1 <- shrink_invasiveness(tab)
  t2 <- shrink_invasiveness(other_tab)

  # Current age table
  v1 <- setNames(t1$v_shrunk, t1$Serotype)
  raw1 <- setNames(t1$invasiveness, t1$Serotype)
  lower1 <- setNames(t1$invasiveness_lower, t1$Serotype)
  upper1 <- setNames(t1$invasiveness_upper, t1$Serotype)
  lm1 <- setNames(t1$logmean, t1$Serotype);  ls1 <- setNames(t1$logsd, t1$Serotype)
  # Other age table
  v2 <- setNames(t2$v_shrunk, t2$Serotype)
  raw2 <- setNames(t2$invasiveness, t2$Serotype)
  lower2 <- setNames(t2$invasiveness_lower, t2$Serotype)
  upper2 <- setNames(t2$invasiveness_upper, t2$Serotype)
  lm2 <- setNames(t2$logmean, t2$Serotype);  ls2 <- setNames(t2$logsd, t2$Serotype)

  # Shrinkage constants, used as the uncertainty for imputed serotypes
  mu_1  <- attr(t1, "mu")
  tau_1 <- attr(t1, "tau")
  tau_2 <- attr(t2, "tau")

  # Create an output row
  make_row <- function(s, donor, v, raw, lower, upper, source,
                       logmean = NA_real_, logsd = NA_real_) {

    tibble(
      serotype = s,
      v = as.numeric(v),
      v_raw = as.numeric(raw),
      v_lower = as.numeric(lower),
      v_upper = as.numeric(upper),
      v_logmean = as.numeric(logmean),
      v_logsd = as.numeric(logsd),
      v_source = source,
      v_donor = donor
    )
  }

  # Copy all values directly from a donor row
  copy_donor <- function(s, donor, age_table, source) {

    if (age_table == 1) {

      make_row(
        s = s,
        donor = donor,
        v = v1[[donor]],
        raw = raw1[[donor]],
        lower = lower1[[donor]],
        upper = upper1[[donor]],
        source = source,
        logmean = lm1[[donor]],
        logsd = ls1[[donor]]
      )

    } else {

      make_row(
        s = s,
        donor = donor,
        v = v2[[donor]],
        raw = raw2[[donor]],
        lower = lower2[[donor]],
        upper = upper2[[donor]],
        source = source,
        logmean = lm2[[donor]],
        logsd = ls2[[donor]]
      )
    }
  }

  # Geometric mean for ordinary serogroup fallback
  # Does not invent confidence intervals
  sg_mean <- function(v, sg) {

    x <- v[serogroup_of(names(v)) == sg]

    if (length(x)) {
      exp(mean(log(x)))
    } else {
      NA_real_
    }
  }

  bind_rows(lapply(serotypes, function(s) {

    s <- as.character(s)
    sg <- serogroup_of(s)

    # --------------------------------------------------------------------------
    # 1. Non-typeables
    # --------------------------------------------------------------------------

    if (s %in% non_typeable) {

      return(make_row(
        s, NA_character_,
        NA_real_, NA_real_, NA_real_, NA_real_,
        "non-typeable"
      ))
    }

    # --------------------------------------------------------------------------
    # 2. Hard-coded substitutions
    # Copy donor estimate and BOTH original CI endpoints.
    # --------------------------------------------------------------------------

    if (s %in% names(serotype_substitutions)) {

      donor <- unname(serotype_substitutions[[s]])

      if (donor %in% names(v1)) {

        return(copy_donor(
          s, donor, 1,
          paste("substituted", donor)
        ))
      }

      if (donor %in% names(v2)) {

        return(copy_donor(
          s, donor, 2,
          paste("substituted", donor, "other age table")
        ))
      }

      stop(
        "Hard-coded substitution donor not found: ",
        s, " -> ", donor
      )
    }

    # --------------------------------------------------------------------------
    # 3. Exact serotype match in current age table
    # --------------------------------------------------------------------------

    if (s %in% names(v1)) {

      return(copy_donor(
        s, s, 1, "exact"
      ))
    }

    # --------------------------------------------------------------------------
    # 6. Exact serotype in other age table
    # --------------------------------------------------------------------------

    if (s %in% names(v2)) {

      return(copy_donor(
        s, s, 2, "other age table"
      ))
    }


    # --------------------------------------------------------------------------
    # 4. Actual serogroup row in current age table
    # --------------------------------------------------------------------------

    if (sg %in% names(v1)) {

      return(copy_donor(
        s, sg, 1,
        paste("serogroup", sg)
      ))
    }

    # --------------------------------------------------------------------------
    # 5. Serogroup geometric mean in current age table
    # --------------------------------------------------------------------------

    m <- sg_mean(v1, sg)

    if (!is.na(m)) {

      return(make_row(
        s, paste0("serogroup ", sg),
        m, NA_real_, NA_real_, NA_real_,
        paste("serogroup mean", sg),
        logmean = log(m), logsd = tau_1
      ))
    }

    # --------------------------------------------------------------------------
    # 7. Actual serogroup row in other age table
    # --------------------------------------------------------------------------

    if (sg %in% names(v2)) {

      return(copy_donor(
        s, sg, 2,
        paste("other age table, serogroup", sg)
      ))
    }

    # --------------------------------------------------------------------------
    # 8. Serogroup geometric mean in other age table
    # --------------------------------------------------------------------------

    m <- sg_mean(v2, sg)

    if (!is.na(m)) {

      return(make_row(
        s, paste0("serogroup ", sg),
        m, NA_real_, NA_real_, NA_real_,
        paste("other age table, serogroup mean", sg),
        logmean = log(m), logsd = tau_2
      ))
    }

    # --------------------------------------------------------------------------
    # 9. Overall mean fallback
    # --------------------------------------------------------------------------

    make_row(
      s, NA_character_,
      attr(t1, "overall_mean"),
      NA_real_, NA_real_, NA_real_,
      "overall mean",
      logmean = mu_1, logsd = tau_1
    )
  }))
}

# ------------------------------------------------------------------------------
# Verify that donor intervals were copied from the intended source table
# ------------------------------------------------------------------------------

#' @importFrom dplyr %>% mutate filter
check_substitutions <- function(result, primary, secondary,
                                prefix, serotype_substitutions) {

  primary <- primary %>%
    mutate(Serotype = as.character(Serotype))

  secondary <- secondary %>%
    mutate(Serotype = as.character(Serotype))

  for (s in intersect(
    names(serotype_substitutions),
    result$serotype
  )) {

    donor <- unname(serotype_substitutions[[s]])

    source_row <- primary %>%
      filter(Serotype == donor)

    if (nrow(source_row) == 0) {
      source_row <- secondary %>%
        filter(Serotype == donor)
    }

    if (nrow(source_row) != 1) {
      stop("Cannot verify donor for ", s)
    }

    actual <- result %>%
      filter(serotype == s)

    stopifnot(
      isTRUE(all.equal(
        actual[[paste0(prefix, "_lower")]],
        as.numeric(source_row$invasiveness_lower)
      )),
      isTRUE(all.equal(
        actual[[paste0(prefix, "_upper")]],
        as.numeric(source_row$invasiveness_upper)
      ))
    )
  }
}

# ==============================================================================
# Main entry point: child and adult lookups, combined
# ==============================================================================

#' Build the invasiveness lookup with intervals
#'
#' For every serotype in `meta_clean`, finds child and adult invasiveness
#' estimates: from hard-coded substitutions, exact matches, serogroups or the
#' overall mean, in that order of priority (see the comments in the code).
#'
#' @param meta_clean Cleaned metadata, e.g. from [clean_metadata()].
#' @param inv_children,inv_adults Invasiveness tables for children and adults,
#'   with columns Serotype, invasiveness, invasiveness_lower, invasiveness_upper.
#' @param serotype_substitutions Named character vector: serotype = donor.
#'   These take priority over exact matches and serogroup averages.
#' @param non_typeable Serotype labels that get no invasiveness estimate.
#' @return One row per serotype, with v_child_* and v_adult_* columns.
#' @importFrom dplyr %>% rename left_join
#' @importFrom stats na.omit
#' @export
build_invasiveness_lookup <- function(meta_clean, inv_children, inv_adults,
                                      serotype_substitutions, non_typeable) {

  serotypes <- sort(unique(
    na.omit(as.character(meta_clean$serotype))
  ))

  message("Children:")
  child <- lookup_one_age(
    inv_children,
    inv_adults,
    serotypes,
    serotype_substitutions,
    non_typeable
  )

  message("Adults:")
  adult <- lookup_one_age(
    inv_adults,
    inv_children,
    serotypes,
    serotype_substitutions,
    non_typeable
  )

  invasiveness_lookup <- child %>%
    rename(
      v_child = v,
      v_child_raw = v_raw,
      v_child_lower = v_lower,
      v_child_upper = v_upper,
      v_child_logmean = v_logmean,
      v_child_logsd = v_logsd,
      v_child_source = v_source,
      v_child_donor = v_donor
    ) %>%
    left_join(
      adult %>%
        rename(
          v_adult = v,
          v_adult_raw = v_raw,
          v_adult_lower = v_lower,
          v_adult_upper = v_upper,
          v_adult_logmean = v_logmean,
          v_adult_logsd = v_logsd,
          v_adult_source = v_source,
          v_adult_donor = v_donor
        ),
      by = "serotype"
    )

  check_substitutions(
    invasiveness_lookup,
    inv_children, inv_adults,
    "v_child",
    serotype_substitutions
  )

  check_substitutions(
    invasiveness_lookup,
    inv_adults, inv_children,
    "v_adult",
    serotype_substitutions
  )

  invasiveness_lookup
}
# Build invasiveness lookup for South African pneumococcal data
# Inputs:
#   data/SouthAfrica/meta_clean.csv
#   data/Invasiveness/nick-inv-sts-children.xlsx
#   data/Invasiveness/nick-inv-sts-adults.xlsx
# Output:
#   data/Invasiveness/invasiveness_lookup.csv
#
# Each serotype receives child and adult invasiveness estimates. Exact matches
# retain the original published 95% CI endpoints (invasiveness_lower/upper).
# For serogroup/cross-age/overall-mean imputations, CI endpoints are NA: a
# published interval for the imputed estimate is unavailable.
# Note: v_child and v_adult are empirical-Bayes shrunken point estimates;
# the saved v_*_lower/upper are ORIGINAL (unshrunken) published intervals.

library(dplyr)
library(readr)
library(openxlsx)

meta_clean <- read_csv("data/SouthAfrica/meta_clean.csv", show_col_types = FALSE)
nickkids   <- read.xlsx("data/Invasiveness/nick-inv-sts-children.xlsx")
nickadults <- read.xlsx("data/Invasiveness/nick-inv-sts-adults.xlsx")

non_typeable <- c("untypable", "Swiss_NT")

# Empirical-Bayes shrinkage on the log scale, with DerSimonian-Laird tau^2.
shrink_invasiveness <- function(tab) {
  required <- c("Serotype", "invasiveness", "invasiveness_lower", "invasiveness_upper")
  missing <- setdiff(required, names(tab))
  if (length(missing)) stop("Missing columns in invasiveness table: ", paste(missing, collapse = ", "))

  valid <- with(tab,
    is.finite(invasiveness) & invasiveness > 0 &
    is.finite(invasiveness_lower) & invasiveness_lower > 0 &
    is.finite(invasiveness_upper) & invasiveness_upper > invasiveness_lower
  )
  if (!all(valid)) stop("Invalid or missing invasiveness estimates / CI limits in input table")
  if (anyDuplicated(tab$Serotype)) stop("Duplicate serotypes in invasiveness input table")

  y  <- log(tab$invasiveness)
  se <- (log(tab$invasiveness_upper) - log(tab$invasiveness_lower)) / (2 * 1.96)
  if (any(se <= 0)) stop("Non-positive log-scale standard errors")
  w <- 1 / se^2
  mu0 <- sum(w * y) / sum(w)
  Q <- sum(w * (y - mu0)^2)
  denom <- sum(w) - sum(w^2) / sum(w)
  tau2 <- if (denom > 0) max(0, (Q - (length(y) - 1)) / denom) else 0
  mu <- sum(y / (se^2 + tau2)) / sum(1 / (se^2 + tau2))
  B <- tau2 / (tau2 + se^2)
  message(sprintf("  shrinkage: exp(mu) = %.2e, tau = %.2f (log scale)", exp(mu), sqrt(tau2)))
  result <- tab %>% mutate(v_shrunk = exp(mu + B * (y - mu)))
  attr(result, "overall_mean") <- exp(mu)
  result
}

# Extract the initial numeric serogroup: '33D' -> '33'.
serogroup_of <- function(s) sub("^([0-9]+).*$", "\\1", s)

lookup_one_age <- function(tab, other_tab, serotypes) {
  t1 <- shrink_invasiveness(tab)
  t2 <- shrink_invasiveness(other_tab)

  v1 <- setNames(t1$v_shrunk, t1$Serotype)
  v2 <- setNames(t2$v_shrunk, t2$Serotype)
  raw1 <- setNames(t1$invasiveness, t1$Serotype)
  lower1 <- setNames(t1$invasiveness_lower, t1$Serotype)
  upper1 <- setNames(t1$invasiveness_upper, t1$Serotype)

  sg_mean <- function(v, sg) {
    x <- v[serogroup_of(names(v)) == sg]
    if (length(x)) exp(mean(log(x))) else NA_real_
  }

  bind_rows(lapply(serotypes, function(s) {
    sg <- serogroup_of(s)
    row <- function(v, raw, lower, upper, src) {
      tibble(serotype = s, v = v, v_raw = raw,
             v_lower = lower, v_upper = upper, v_source = src)
    }

    if (s %in% non_typeable) {
      return(row(NA_real_, NA_real_, NA_real_, NA_real_, "non-typeable"))
    }
    if (s %in% names(v1)) {
      return(row(v1[[s]], raw1[[s]], lower1[[s]], upper1[[s]], "exact"))
    }
    if (!is.na(m <- sg_mean(v1, sg))) {
      return(row(m, NA_real_, NA_real_, NA_real_, paste("serogroup", sg)))
    }
    if (s %in% names(v2)) {
      return(row(v2[[s]], NA_real_, NA_real_, NA_real_, "other age table"))
    }
    if (!is.na(m <- sg_mean(v2, sg))) {
      return(row(m, NA_real_, NA_real_, NA_real_, paste("other age table, serogroup", sg)))
    }
    row(attr(t1, "overall_mean"), NA_real_, NA_real_, NA_real_, "overall mean")
  }))
}

serotypes <- sort(unique(na.omit(meta_clean$serotype)))
message("children:")
child <- lookup_one_age(nickkids, nickadults, serotypes)
message("adults:")
adult <- lookup_one_age(nickadults, nickkids, serotypes)

invasiveness_lookup <- child %>%
  rename(
    v_child = v, v_child_raw = v_raw,
    v_child_lower = v_lower, v_child_upper = v_upper,
    v_child_source = v_source
  ) %>%
  left_join(
    adult %>% rename(
      v_adult = v, v_adult_raw = v_raw,
      v_adult_lower = v_lower, v_adult_upper = v_upper,
      v_adult_source = v_source
    ),
    by = "serotype"
  )

out_path <- "data/SouthAfrica/processed-inputs/invasiveness_lookup_interval.csv"
write_csv(invasiveness_lookup, out_path)
message("Wrote ", nrow(invasiveness_lookup), " serotypes to ", out_path)

# Inspect values that required imputation or have no intervals.
invasiveness_lookup %>%
  filter(v_child_source != "exact" | v_adult_source != "exact") %>%
  select(serotype, v_child_source, v_adult_source) %>%
  print(n = Inf)

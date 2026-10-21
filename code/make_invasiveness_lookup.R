

library(dplyr)
library(readr)
library(openxlsx)

###########
# --- readme  
#
# Input:  data/SouthAfrica/meta_clean.csv            
#         data/Invasiveness/nick-inv-sts-children.xlsx
#         data/Invasiveness/nick-inv-sts-adults.xlsx
# Output: data/Invasiveness/invasiveness_lookup.csv
#         one row per serotype in meta_clean, with child and adult invasiveness.
#
# Two things:
# (1) Shrinkage. Nick's estimates are pulled toward the overall mean on the log
#     scale, more strongly the wider their 95% CI (this is some kind of empirical Bayes /
#     DerSimonian-Laird random-effects thingy). Precisely estimated serotypes have their inv unchanged; serotypes with near-zero IPD in Nick's data (e.g. 19B, 35A, 33A)
#     move up to ~10-fold, so a single IPD case can no longer be a high portion of the carriage pop in a year.

# Description: "Serotype invasiveness estimates were stabilised by empirical Bayes
# shrinkage on the log scale, using a DerSimonian–Laird random-effects model 
# with standard errors derived from the reported 95% intervals.". Recommended by AI. 
# equation: log( inv_i_shrunk ) = mu + Bi (log (inv_i) - mu) 
# when Bi is 1 the "shrunk" estimate is the inv_i from Nick's table for type i 
# when Bi is near 0, it's the mean mu of all of them (log of course) 
# Bi = tau^2 / (tau^2 + si^2) 
# where si is the standard error of the estimate for type i from Nick's thing (derived from his 95% CIs)
# (log upper − log lower) / (2 × 1.96). 
# $\tau^2$ is the true between-serotype variance estimated with DerSimonsonLaird, tau~1.27 
# mu is the mean of the log inv estimates 

# (2) Coverage. Every serotype in meta_clean gets a value:
#       exact match in the age table
#       first: want the geometric mean of its serogroup (same age table)
#       or if not, the same serotype in the other age table
#       or if neither of those, its serogroup in the other age table
#       finally the overall mean.
#     Non-typeables get NA: they are left out of the IPD-based components and
#     enter the synthetic population through carriage only.
#
# Re-run this whenever we change meta_clean (ie when we re-run cleanmetadata.R) 
###########



meta_clean <- read_csv("data/SouthAfrica/meta_clean.csv", show_col_types = FALSE)
nickkids   <- read.xlsx("data/Invasiveness/nick-inv-sts-children.xlsx")
nickadults <- read.xlsx("data/Invasiveness/nick-inv-sts-adults.xlsx")

non_typeable <- c("untypable", "Swiss_NT")

# ---------------------------------------------------------------------------
# (1) Shrinkage on the log scale as described in my notes above 
#   log v_i^shrunk = mu + B_i (log v_i - mu),  B_i = tau^2 / (tau^2 + s_i^2)
#   s_i  : log-scale SE from the 95% CI
#   mu   : precision-weighted mean of the log estimates (mean of the logs)
#   tau^2: between-serotype variance of true log invasiveness (DerSimonian-Laird)
# ---------------------------------------------------------------------------
shrink_invasiveness <- function(tab) {
  y  <- log(tab$invasiveness)
  se <- (log(tab$invasiveness_upper) - log(tab$invasiveness_lower)) / (2 * 1.96)
  w  <- 1 / se^2
  mu0  <- sum(w * y) / sum(w)
  Q    <- sum(w * (y - mu0)^2)
  tau2 <- max(0, (Q - (length(y) - 1)) / (sum(w) - sum(w^2) / sum(w)))
  mu   <- sum(y / (se^2 + tau2)) / sum(1 / (se^2 + tau2))
  B    <- tau2 / (tau2 + se^2)
  message(sprintf("  shrinkage: exp(mu) = %.2e, tau = %.2f (log scale)", exp(mu), sqrt(tau2)))
  tab %>% mutate(v_shrunk = exp(mu + B * (y - mu))) %>%
    structure(overall_mean = exp(mu))
}

# ---------------------------------------------------------------------------
# (2) invasiveness vals for the serotypes
# ---------------------------------------------------------------------------
serogroup_of <- function(s) sub("^([0-9]+).*$", "\\1", s)   # "33D" -> "33" , ie  gets the 'serogroup'

lookup_one_age <- function(tab, other_tab, serotypes) {
  t1 <- shrink_invasiveness(tab); t2 <- shrink_invasiveness(other_tab)
  v1   <- setNames(t1$v_shrunk, t1$Serotype)
  v2   <- setNames(t2$v_shrunk, t2$Serotype)
  raw1 <- setNames(t1$invasiveness, t1$Serotype)
  sg_mean <- function(v, sg) {
    x <- v[serogroup_of(names(v)) == sg]
    if (length(x)) exp(mean(log(x))) else NA
  }
  bind_rows(lapply(serotypes, function(s) {
    sg <- serogroup_of(s)
    row <- function(v, raw, src) tibble(serotype = s, v = v, v_raw = raw, v_source = src)
    if (s %in% non_typeable)          return(row(NA, NA, "non-typeable"))
    if (s %in% names(v1))             return(row(v1[[s]], raw1[[s]], "exact"))
    if (!is.na(m <- sg_mean(v1, sg))) return(row(m, NA, paste("serogroup", sg)))
    if (s %in% names(v2))             return(row(v2[[s]], NA, "other age table"))
    if (!is.na(m <- sg_mean(v2, sg))) return(row(m, NA, paste("other age table, serogroup", sg)))
    row(attr(t1, "overall_mean"), NA, "overall mean")
  }))
}

serotypes <- sort(unique(meta_clean$serotype))

message("children:"); child <- lookup_one_age(nickkids,   nickadults, serotypes)
message("adults:");   adult <- lookup_one_age(nickadults, nickkids,   serotypes)

invasiveness_lookup <- child %>%
  rename(v_child = v, v_child_raw = v_raw, v_child_source = v_source) %>%
  left_join(adult %>% rename(v_adult = v, v_adult_raw = v_raw, v_adult_source = v_source),
            by = "serotype")

write_csv(invasiveness_lookup, "data/Invasiveness/invasiveness_lookup.csv")

# what was not an exact match (worth a look after each run)
invasiveness_lookup %>%
  filter(v_child_source != "exact" | v_adult_source != "exact") %>%
  select(serotype, v_child_source, v_adult_source) %>%
  print(n = Inf)

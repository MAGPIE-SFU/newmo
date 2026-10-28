library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(tibble)
library(ggplot2)

# ============================================================
# this script simulates an expected IPD population from a simulated carriage pop
# of 100k per year, includes uncertainty propagation from invasiveness estimates
# NOTE: need to set the correct simulation parameters
# the scrpit plot_projected_ipd.R will plot the outputs of this script
# ============================================================

# Parameters ---------------------------------------------------------------
N_REPS <- 1000L
N_CARRIAGE <- 100000
NYEARS <- 6L
FIRST_SIM_CALENDAR_YEAR <- 2011L # corresponds to year.1 in the simulation
set.seed(1)

INVASIVE_FILEPATH <- "data/SouthAfrica/processed-inputs/invasiveness_lookup_interval.csv"
GROUP_FILEPATH <- "data/SouthAfrica/processed-inputs/grouplookup.csv"
META_FILEPATH <- "data/SouthAfrica/processed-inputs/meta_clean.csv"
OUT_DIR <- "data/SouthAfrica/processed-inputs"

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
stopifnot(all(file.exists(c(INVASIVE_FILEPATH, GROUP_FILEPATH, META_FILEPATH))))

# IMPORTANT: `samodel`, `modelpars`, and `simulate()` must be defined/sourced
# before this script is run; use your chosen NFDS simulation parameters.
stopifnot(exists("samodel"), exists("modelpars"), exists("simulate", mode = "function"))

# 1. Simulated carriage composition ---------------------------------------
sim <- simulate(samodel = samodel, modelpars = modelpars, nyears = NYEARS)
stopifnot(nrow(sim) == length(samodel$st),
          !is.null(rownames(sim)), !anyDuplicated(rownames(sim)))

sim_long <- sim %>%
  as.data.frame() %>%
  rownames_to_column("group") %>%
  pivot_longer(cols = starts_with("year."),
               names_to = "year", values_to = "x") %>%
  mutate(group = as.character(group),
         year = as.integer(sub("year.", "", year))) %>%
  group_by(year) %>%
  mutate(frac = x / sum(x), carriage = N_CARRIAGE * frac) %>%
  ungroup()
stopifnot(all(is.finite(sim_long$frac)), all(sim_long$frac >= 0))

# 2. The SAME child progression-rate sampling scheme as the revised
#    carriage uncertainty script ------------------------------------------
inv <- read_csv(INVASIVE_FILEPATH, show_col_types = FALSE)
inv_dist <- inv %>%
  transmute(serotype = as.character(serotype),
            v_mean = v_child,
            logmean = v_child_logmean,
            logsd = v_child_logsd) %>%
  mutate(v_mean = if_else(serotype %in% c("Swiss_NT", "untypable"), 0, v_mean))
stopifnot(!anyDuplicated(inv_dist$serotype))

# Each iteration uses one child-rate draw per serotype, held fixed across
# all simulation years and all groups. The same draw is used for IPD totals,
# group composition and serotype composition.
draw_v <- function() {
  z <- inv_dist
  zero_type <- z$serotype %in% c("Swiss_NT", "untypable")
  valid <- is.finite(z$logmean) & is.finite(z$logsd) &
    z$logsd >= 0 & is.finite(z$v_mean) & z$v_mean > 0
  invalid <- !zero_type & !valid
  if (any(invalid)) {
    stop("Missing/invalid child log-normal progression-rate parameters for: ",
         paste(z$serotype[invalid], collapse = ", "))
  }
  z$v_sample <- 0
  z$v_sample[valid & !zero_type] <- rlnorm(
    sum(valid & !zero_type),
    meanlog = z$logmean[valid & !zero_type],
    sdlog = z$logsd[valid & !zero_type]
  )
  stopifnot(all(is.finite(z$v_sample)), all(z$v_sample >= 0))
  z %>% select(serotype, v_sample)
}

# 3. Map unchanged simulation groups to serotypes ------------------------
group_lookup <- read_csv(GROUP_FILEPATH, show_col_types = FALSE) %>%
  mutate(group = as.character(group), serotype = as.character(serotype))
meta_clean <- read_csv(META_FILEPATH, show_col_types = FALSE) %>%
  mutate(serotype = as.character(serotype))

stopifnot(!anyDuplicated(group_lookup[c("serotype", "gpsc")]))

# Allocate each simulated group among its constituent serotypes.
# This weighting is only necessary for mixed-serotype groups (e.g. rare_NVT).
# Use observed child carriage counts across all years, as in the old script.
serotype_weights <- group_lookup %>%
  distinct(group, serotype) %>%
  left_join(
    meta_clean %>%
      filter(age_group == "child", collection_type == "carriage") %>%
      count(serotype, gpsc, name = "n") %>%
      left_join(group_lookup %>% distinct(serotype, gpsc, group),
                by = c("serotype", "gpsc"), relationship = "many-to-one") %>%
      group_by(group, serotype) %>%
      summarise(n = sum(n), .groups = "drop"),
    by = c("group", "serotype"), relationship = "one-to-one"
  ) %>%
  mutate(n = replace_na(n, 0)) %>%
  group_by(group) %>%
  mutate(n_serotypes = n(),
         weight = case_when(
           n_serotypes == 1 ~ 1,
           sum(n) > 0 ~ n / sum(n),
           TRUE ~ NA_real_
         )) %>%
  ungroup()

if (anyNA(serotype_weights$weight)) {
  print(serotype_weights %>% filter(is.na(weight)))
  stop("No observed child-carriage weights for one or more mixed-serotype groups")
}
stopifnot(all(abs(serotype_weights %>% group_by(group) %>%
                    summarise(s = sum(weight), .groups = "drop") %>%
                    pull(s) - 1) < 1e-8))

missing_groups <- setdiff(unique(sim_long$group), unique(serotype_weights$group))
if (length(missing_groups))
  stop("Unmatched simulated groups: ", paste(missing_groups, collapse = ", "))

carriage_serotype <- sim_long %>%
  left_join(serotype_weights %>% select(group, serotype, weight),
            by = "group", relationship = "many-to-many") %>%
  mutate(carriage_serotype = carriage * weight)

missing_v <- setdiff(unique(carriage_serotype$serotype), inv_dist$serotype)
if (length(missing_v))
  stop("Missing serotypes in invasiveness lookup: ", paste(missing_v, collapse = ", "))

check_car <- carriage_serotype %>%
  group_by(year) %>%
  summarise(n = sum(carriage_serotype), .groups = "drop")
stopifnot(all(abs(check_car$n - N_CARRIAGE) < 1e-6))

# 4. Monte Carlo: a SINGLE collection of draws for ALL outputs -----------
# This propagates uncertainty in estimated rates, holding simulated
# carriage fixed. `ipd_expected` is a Poisson mean, not a realized case count.
ipd_serotype_component_draws <- map_dfr(seq_len(N_REPS), function(b) {
  v_draw <- draw_v()
  carriage_serotype %>%
    left_join(v_draw, by = "serotype", relationship = "many-to-one") %>%
    mutate(ipd_expected = carriage_serotype * v_sample,
           draw = b) %>%
    select(draw, year, group, serotype, frac, carriage,
           carriage_serotype, ipd_expected)
})

# 5. Group-level IPD predictions ------------------------------------------
ipd_draws <- ipd_serotype_component_draws %>%
  group_by(draw, year, group) %>%
  summarise(frac = first(frac), carriage = first(carriage),
            ipd_expected = sum(ipd_expected), .groups = "drop")

ipd_summary <- ipd_draws %>%
  group_by(year, group) %>%
  summarise(frac = mean(frac), carriage = mean(carriage),
            ipd_mean = mean(ipd_expected),
            ipd_expected_sd = sd(ipd_expected),
            ipd_sd = sqrt(mean(ipd_expected) + var(ipd_expected)),
            ipd_lower = quantile(ipd_expected, 0.025),
            ipd_upper = quantile(ipd_expected, 0.975),
            .groups = "drop")

# 6. Total IPD: aggregate INSIDE draws ------------------------------------
ipd_yearly_draws <- ipd_draws %>%
  group_by(draw, year) %>%
  summarise(total_ipd = sum(ipd_expected), .groups = "drop")

ipd_yearly <- ipd_yearly_draws %>%
  group_by(year) %>%
  summarise(ipd_mean = mean(total_ipd),
            ipd_expected_sd = sd(total_ipd),
            ipd_sd = sqrt(mean(total_ipd) + var(total_ipd)),
            ipd_lower = quantile(total_ipd, 0.025),
            ipd_upper = quantile(total_ipd, 0.975),
            .groups = "drop")

# 7. Serotype-level IPD and composition using the SAME draws --------------
ipd_serotype_draws <- ipd_serotype_component_draws %>%
  group_by(draw, year, serotype) %>%
  summarise(ipd_expected = sum(ipd_expected), .groups = "drop") %>%
  group_by(draw, year) %>%
  mutate(ipd_frac = ipd_expected / sum(ipd_expected)) %>%
  ungroup()

ipd_serotype_summary <- ipd_serotype_draws %>%
  group_by(year, serotype) %>%
  summarise(frac = mean(ipd_frac), frac_sd = sd(ipd_frac),
            frac_lower = quantile(ipd_frac, 0.025),
            frac_upper = quantile(ipd_frac, 0.975),
            ipd_mean = mean(ipd_expected), .groups = "drop") %>%
  arrange(year, serotype)

stopifnot(all(abs(ipd_serotype_summary %>%
                    group_by(year) %>% summarise(s = sum(frac), .groups = "drop") %>%
                    pull(s) - 1) < 1e-8))

# 8. Save outputs ---------------------------------------------------------
write_csv(ipd_summary, file.path(OUT_DIR, "simulated_ipd_uncertainty.csv"))
write_csv(ipd_yearly, file.path(OUT_DIR, "simulated_ipd_yearly_uncertainty.csv"))
write_csv(ipd_yearly_draws, file.path(OUT_DIR, "simulated_ipd_yearly_draws.csv"))
write_csv(ipd_serotype_summary, file.path(OUT_DIR, "simulated_ipd_serotype_uncertainty.csv"))

# 9. Plot to check results ------------------------------------------------

# Observed vs predicted serotype composition: calendar-year alignment
observed_ipd <- meta_clean %>%
  filter(collection_type == "IPD") %>%
  count(year, serotype, name = "ipd_count") %>%
  group_by(year) %>%
  mutate(frac = ipd_count / sum(ipd_count)) %>%
  ungroup()

predicted_ipd <- ipd_serotype_summary %>%
  mutate(year = FIRST_SIM_CALENDAR_YEAR + year - 1L, source = "Predicted")

p_comparison <- ggplot() +
  geom_ribbon(data = predicted_ipd,
              aes(x = year, ymin = frac_lower, ymax = frac_upper),
              fill = "#C44E52", alpha = 0.2) +
  geom_line(data = predicted_ipd, aes(year, frac),
            colour = "#C44E52", linewidth = 0.8) +
  geom_point(data = observed_ipd, aes(year, frac),
             colour = "black", size = 1.8) +
  facet_wrap(~serotype, scales = "free_y") +
  labs(x = "Calendar year", y = "Proportion of IPD") +
  theme_bw()
print(p_comparison)
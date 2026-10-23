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


# ============================================================
# Parameters
# ============================================================

N_REPS <- 1000
N_CARRIAGE <- 100000
NYEARS <- 6
#samodel <-  DEFINE THESE!!!!
#modelpars <- 

set.seed(1)

INVASIVE_FILEPATH <-
  "data/SouthAfrica/processed-inputs/invasiveness_lookup_interval.csv"

GROUP_FILEPATH <-
  "data/SouthAfrica/processed-inputs/grouplookup.csv"

META_FILEPATH <-
  "data/SouthAfrica/processed-inputs/meta_clean.csv"

OUT_DIR <- "data/SouthAfrica/processed-inputs"

# ============================================================
# 1. Run simulator
# ============================================================
# fill this in currently CCmagic 
sim <- simulate(
  samodel = samodel,
  modelpars = modelpars,
  nyears = NYEARS
)

stopifnot(
  nrow(sim) == length(samodel$st),
  !is.null(rownames(sim)),
  !anyDuplicated(rownames(sim))
)

# Convert simulated abundances into carriage composition
sim_long <- sim %>%
  rownames_to_column("group") %>%
  pivot_longer(
    cols = starts_with("year."),
    names_to = "year",
    values_to = "x"
  ) %>%
  mutate(
    year = as.integer(sub("year.", "", year))
  ) %>%
  group_by(year) %>%
  mutate(
    frac = x / sum(x),
    carriage = N_CARRIAGE * frac
  ) %>%
  ungroup()

# ============================================================
# 2. Prepare invasiveness uncertainty
# Same lognormal parameterization as previous script
# ============================================================

inv <- read_csv(INVASIVE_FILEPATH, show_col_types = FALSE)

inv_dist <- inv %>%
  transmute(
    serotype = as.character(serotype),
    v_mean = v_child,
    v_lower = v_child_lower,
    v_upper = v_child_upper
  ) %>%
  mutate(
    valid_ci =
      is.finite(v_mean) &
      is.finite(v_lower) &
      is.finite(v_upper) &
      v_mean > 0 &
      v_lower > 0 &
      v_upper > v_lower,

    SE = if_else(
      valid_ci,
      (v_upper - v_lower) / (2 * 1.96),
      0
    ),

    logmean = log(
      v_mean^2 / sqrt(v_mean^2 + SE^2)
    ),

    logsd = sqrt(
      log(1 + SE^2 / v_mean^2)
    )
  )

 inv_dist <- inv_dist %>%
  mutate(
    v_mean = if_else(
      serotype %in% c("Swiss_NT", "untypable"),
      0, v_mean
    ),
    logmean = if_else(
      serotype %in% c("Swiss_NT", "untypable"),
      NA_real_, logmean
    ),
    logsd = if_else(
      serotype %in% c("Swiss_NT", "untypable"),
      0, logsd
    )
  ) 

stopifnot(!anyDuplicated(inv_dist$serotype))

# ============================================================
# 3. Map simulated groups to serotypes
# ============================================================

group_lookup <- read_csv(
  GROUP_FILEPATH, show_col_types = FALSE
)

meta_clean <- read_csv(
  META_FILEPATH, show_col_types = FALSE
)

group_lookup <- group_lookup %>%
  mutate(
    group = as.character(group),
    serotype = as.character(serotype)
  )

# Estimate within-group serotype weights using
# observed childhood carriage across all years.
# For single-serotype groups the weight is exactly 1.

serotype_weights <- group_lookup %>%
  distinct(group, serotype) %>%
  left_join(
    meta_clean %>%
      filter(
        age_group == "child",
        collection_type == "carriage"
      ) %>%
      count(serotype, gpsc, name = "n") %>%
      left_join(
        group_lookup %>%
          distinct(serotype, gpsc, group),
        by = c("serotype", "gpsc")
      ) %>%
      group_by(group, serotype) %>%
      summarise(n = sum(n), .groups = "drop"),
    by = c("group", "serotype")
  ) %>%
  mutate(n = replace_na(n, 0)) %>%
  group_by(group) %>%
  mutate(
    n_serotypes = n(),
    weight = case_when(
      n_serotypes == 1 ~ 1,
      sum(n) > 0 ~ n / sum(n),
      TRUE ~ NA_real_
    )
  ) %>%
  ungroup()

# Stop if there is no evidence for weights
# in a group containing multiple serotypes
if (any(is.na(serotype_weights$weight))) {
  print(
    serotype_weights %>%
      filter(is.na(weight))
  )
  stop("Cannot estimate within-group serotype weights")
}

missing_groups <- setdiff(
  unique(sim_long$group),
  unique(serotype_weights$group)
)

if (length(missing_groups)) {
  stop(
    "Unmatched simulated groups: ",
    paste(missing_groups, collapse = ", ")
  )
}

# Expand group carriage into serotype contributions
carriage_serotype <- sim_long %>%
  left_join(
    serotype_weights %>%
      select(group, serotype, weight),
    by = "group",
    relationship = "many-to-many"
  ) %>%
  mutate(
    carriage_serotype = carriage * weight
  )

# check 
serotype_weights %>%
  group_by(group) %>%
  summarise(
    n_serotypes = n(),
    total_weight = sum(weight),
    .groups = "drop"
  ) %>%
  filter(abs(total_weight - 1) > 1e-8 | is.na(total_weight)) # should return nothing

carriage_serotype %>%
  group_by(year) %>%
  summarise(
    total_carriage = sum(carriage_serotype),
    .groups = "drop"
  )


# Validate available invasiveness estimates
missing_v <- carriage_serotype %>%
  distinct(serotype) %>%
  anti_join(
    inv_dist %>%
      filter(is.finite(v_mean), v_mean >= 0),
    by = "serotype"
  )

if (nrow(missing_v) > 0) {
  print(missing_v)
  stop("Missing or invalid invasiveness estimates")
}


# ============================================================
# 4. Monte Carlo uncertainty propagation
# ============================================================

ipd_draws <- lapply(seq_len(N_REPS), function(b) {

  # Draw invasiveness ONCE per serotype per iteration
  # note that Swiss_NT and untypable are always exactly zero, 
  # so they will not be drawn from the lognormal distribution
v_draw <- inv_dist %>%
  mutate(v_sample = v_mean)

# Only sample serotypes with positive invasiveness
valid <- is.finite(v_draw$logmean) &
         is.finite(v_draw$logsd) &
         v_draw$v_mean > 0

v_draw$v_sample[valid] <- rlnorm(
  sum(valid),
  meanlog = v_draw$logmean[valid],
  sdlog = v_draw$logsd[valid]
)

# Non-typeables have exactly zero invasiveness
v_draw$v_sample[
  v_draw$serotype %in% c("Swiss_NT", "untypable")
] <- 0

# Check for any remaining invalid estimates
stopifnot(
  all(is.finite(v_draw$v_sample)),
  all(v_draw$v_sample >= 0)
)

v_draw <- v_draw %>%
  select(serotype, v_sample)

  # Convert simulated carriage to expected IPD
  carriage_serotype %>%
    left_join(
      v_draw,
      by = "serotype"
    ) %>%
    mutate(
      ipd_expected = carriage_serotype * v_sample
    ) %>%
    group_by(year, group) %>%
    summarise(
      frac = first(frac),
      carriage = first(carriage),
      ipd_expected = sum(ipd_expected),
      .groups = "drop"
    ) %>%
    mutate(draw = b)
})

ipd_draws <- bind_rows(ipd_draws)

# ============================================================
# 5. Aggregate uncertainty using total variance
# ============================================================

# Conditional Poisson variance = expected IPD
# Variance due to invasiveness = variance across draws

ipd_summary <- ipd_draws %>%
  group_by(year, group) %>%
  summarise(
    frac = mean(frac),
    carriage = mean(carriage),

    ipd_mean = mean(ipd_expected),

    # Uncertainty in expected IPD only
    ipd_expected_sd = sd(ipd_expected),

    # Law of total variance: Poisson + invasiveness
    ipd_sd = sqrt(
      mean(ipd_expected) + var(ipd_expected)
    ),

    ipd_lower = quantile(ipd_expected, 0.025),
    ipd_upper = quantile(ipd_expected, 0.975),

    .groups = "drop"
  )

# ============================================================
# 6. Total IPD per year
# Important: sum within each draw FIRST
# ============================================================

ipd_yearly_draws <- ipd_draws %>%
  group_by(draw, year) %>%
  summarise(
    total_ipd = sum(ipd_expected),
    .groups = "drop"
  )

ipd_yearly <- ipd_yearly_draws %>%
  group_by(year) %>%
  summarise(
    ipd_mean = mean(total_ipd),

    ipd_expected_sd = sd(total_ipd),

    ipd_sd = sqrt(
      mean(total_ipd) + var(total_ipd)
    ),

    ipd_lower = quantile(total_ipd, 0.025),
    ipd_upper = quantile(total_ipd, 0.975),

    .groups = "drop"
  )

# ============================================================
# 7. Save outputs
# ============================================================

write_csv(
  ipd_summary,
  file.path(OUT_DIR, "simulated_ipd_uncertainty.csv")
)

write_csv(
  ipd_yearly,
  file.path(OUT_DIR, "simulated_ipd_yearly_uncertainty.csv")
)

write_csv(
  ipd_yearly_draws,
  file.path(OUT_DIR, "simulated_ipd_yearly_draws.csv")
)

# ============================================================
# 8. Checks
# ============================================================

# Carriage should total 100,000 per year
carriage_serotype %>%
  group_by(year) %>%
  summarise(
    total_carriage = sum(carriage_serotype),
    .groups = "drop"
  ) %>%
  print()

# Look at annual IPD estimates
print(ipd_yearly)

# ============================================================
# 9. Plot
# ============================================================

ggplot(ipd_yearly, aes(x = year, y = ipd_mean)) +
  geom_ribbon(
    aes(
      ymin = pmax(0, ipd_mean - 1.96 * ipd_sd),
      ymax = ipd_mean + 1.96 * ipd_sd
    ),
    alpha = 0.2
  ) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  labs(
    x = "Simulation year",
    y = "Predicted IPD cases",
    title = "Expected IPD from 100,000 annual carriage events"
  ) +
  theme_bw()



# ============================================================
# Predicted serotype composition of IPD with uncertainty
# ============================================================

# Recalculate serotype-level IPD in each Monte Carlo iteration
# using the same invasiveness sampling scheme

set.seed(1)

ipd_serotype_draws <- purrr::map_dfr(
  seq_len(N_REPS),
  function(b) {

    v_draw <- inv_dist %>%
      mutate(v_sample = v_mean)

    valid <- is.finite(v_draw$logmean) &
             is.finite(v_draw$logsd) &
             v_draw$v_mean > 0

    v_draw$v_sample[valid] <- rlnorm(
      sum(valid),
      meanlog = v_draw$logmean[valid],
      sdlog = v_draw$logsd[valid]
    )

    v_draw$v_sample[
      v_draw$serotype %in% c("Swiss_NT", "untypable")
    ] <- 0

    stopifnot(
      all(is.finite(v_draw$v_sample)),
      all(v_draw$v_sample >= 0)
    )

    carriage_serotype %>%
      left_join(
        v_draw %>% select(serotype, v_sample),
        by = "serotype"
      ) %>%
      mutate(
        ipd_expected = carriage_serotype * v_sample
      ) %>%
      group_by(year, serotype) %>%
      summarise(
        ipd_expected = sum(ipd_expected),
        .groups = "drop"
      ) %>%
      group_by(year) %>%
      mutate(
        ipd_frac = ipd_expected / sum(ipd_expected)
      ) %>%
      ungroup() %>%
      mutate(draw = b)
  }
)

# ============================================================
# Summarize composition and uncertainty
# ============================================================

ipd_serotype_summary <- ipd_serotype_draws %>%
  group_by(year, serotype) %>%
  summarise(
    frac = mean(ipd_frac),
    frac_sd = sd(ipd_frac),
    frac_lower = quantile(ipd_frac, 0.025),
    frac_upper = quantile(ipd_frac, 0.975),
    ipd_mean = mean(ipd_expected),
    .groups = "drop"
  ) %>%
  arrange(year, serotype)

# Check that estimated serotype compositions sum to 1
ipd_serotype_summary %>%
  group_by(year) %>%
  summarise(total = sum(frac))


# this is here to compare with the observed IPD composition in the next section
carriage_serotype_summary <- carriage_serotype %>%
  group_by(year, serotype) %>%
  summarise(
    carriage = sum(carriage_serotype),
    .groups = "drop"
  ) %>%
  group_by(year) %>%
  mutate(
    frac = carriage / sum(carriage)
  ) %>%
  ungroup()


# observed ipd composition
observed_ipd <- meta_clean %>%
  filter(collection_type == "IPD") %>%
  count(year, serotype, name = "ipd_count") %>%
  group_by(year) %>%
  mutate(
    frac = ipd_count / sum(ipd_count)
  ) %>%
  ungroup()

# align simulation years to post vax years 
predicted_ipd <- ipd_serotype_summary %>%
  mutate(
    year = 2010 + year,
    source = "Predicted"
  )

observed_ipd <- observed_ipd %>%
  mutate(source = "Observed")

# plot observed vs predicted IPD composition by serotype with uncertainty

ggplot() +
  geom_ribbon(
    data = predicted_ipd,
    aes(
      x = year,
      ymin = frac_lower,
      ymax = frac_upper
    ),
    alpha = 0.2
  ) +
  geom_line(
    data = predicted_ipd,
    aes(x = year, y = frac),
    colour = "#C44E52",
    linewidth = 0.8
  ) +
  geom_point(
    data = observed_ipd,
    aes(x = year, y = frac),
    colour = "black",
    size = 1.8
  ) +
  facet_wrap(~ serotype, scales = "free_y") +
  labs(
    x = "Year",
    y = "Proportion of IPD"
  ) +
  theme_bw()


# ============================================================
# predicted serotype comp in the final year 
# predicted proportion of invasive disease attributable to 
# each serotype in the final year, ordered from most to least common.
# ===========================================================

final_year_composition <- ipd_serotype_summary %>%
  filter(year == max(year)) %>%
  select(
    serotype,
    frac,
    frac_sd,
    frac_lower,
    frac_upper,
    ipd_mean
  ) %>%
  arrange(desc(frac))

print(final_year_composition, n = Inf)

# as percentage 
final_year_composition <- ipd_serotype_summary %>%
  filter(year == max(year)) %>%
  transmute(
    serotype,
    percent = 100 * frac,
    lower_95 = 100 * frac_lower,
    upper_95 = 100 * frac_upper
  ) %>%
  arrange(desc(percent))

print(final_year_composition, n = Inf)

## plot composition in final year 

ggplot(final_year_composition,
       aes(x = reorder(serotype, percent), y = percent)) +
  geom_col(fill = "#397C83") +
  geom_errorbar(
    aes(ymin = lower_95, ymax = upper_95),
    width = 0.2
  ) +
  coord_flip() +
  labs(
    x = "Serotype",
    y = "Predicted proportion of IPD (%)",
    title = "Final-year predicted IPD serotype composition"
  ) +
  theme_bw()


# total IPD counts 


# Mean predicted IPD count per serotype and year
ipd_plot <- ipd_serotype_draws %>%
  group_by(year, serotype) %>%
  summarise(
    ipd_mean = mean(ipd_expected),
    .groups = "drop"
  )

# Order serotypes by overall IPD contribution
sero_order <- ipd_plot %>%
  group_by(serotype) %>%
  summarise(total = sum(ipd_mean), .groups = "drop") %>%
  arrange(desc(total)) %>%
  pull(serotype)

ipd_plot <- ipd_plot %>%
  mutate(serotype = factor(serotype, levels = sero_order))

# Total predicted IPD for each year
total_ipd <- ipd_plot %>%
  group_by(year) %>%
  summarise(
    total = sum(ipd_mean),
    .groups = "drop"
  )

# Plot
ggplot() +
  geom_area(
    data = ipd_plot,
    aes(x = year, y = ipd_mean, fill = serotype),
    position = "stack",
    alpha = 0.9
  ) +
  geom_line(
    data = total_ipd,
    aes(x = year, y = total),
    colour = "black",
    linewidth = 1
  ) +
  scale_x_continuous(breaks = sort(unique(ipd_plot$year))) +
  labs(
    x = "Simulation year",
    y = "Expected IPD cases per 100,000 carriage events",
    fill = "Serotype"
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )


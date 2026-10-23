# compare_fits.R

library(readr)
library(dplyr)
library(tidyr)
library(gt)

# Update this path if needed
output_dir <- "data/SouthAfrica/fit_comparison"

USE_SWITCH <- T
suffix <- if (USE_SWITCH) "_switch" else ""

summary_file <- file.path(output_dir, paste0("fit_summary", suffix, ".csv"))
parameters_file <- file.path(output_dir, paste0("all_fitted_parameters", suffix, ".csv"))

fit_summary <- read_csv(summary_file, show_col_types = FALSE)
parameters <- read_csv(parameters_file, show_col_types = FALSE)

# ------------------------------------------------------------
# Combined fit summary and parameter table
# ------------------------------------------------------------

parameter_cols <- setdiff(names(parameters), c("loss", "method"))

# Load CC_MAGIC baseline parameters
baseline <- readRDS(file.path(output_dir, "CC_MAGIC.rds"))

# Ensure baseline has the same parameter columns
baseline <- baseline |>
  select(all_of(parameter_cols)) |>
  slice(1)

# Create baseline rows, one for each loss function
baseline_rows <- expand_grid(
  loss = unique(fit_summary$loss)
) |>
  mutate(
    method = "CC_MAGIC",
    initial = NA_real_,
    fitted = NA_real_,
    improvement = NA_real_,
    pct_change = NA_real_,
    convergence = NA_integer_,
    runtime_seconds = NA_real_
  )

for (p in parameter_cols) {
  baseline_rows[[p]] <- baseline[[p]][1]
}

# Combine fit statistics and fitted parameters
combined_table <- fit_summary |>
  select(
    loss, method, initial, fitted, improvement,
    pct_change, convergence, runtime_seconds
  ) |>
  left_join(
    parameters |> select(loss, method, all_of(parameter_cols)),
    by = c("loss", "method")
  )

# Put baseline parameter values alongside the fit results
combined_table <- bind_rows(
  combined_table,
  baseline_rows
) |>
  arrange(loss, method != "CC_MAGIC", method)

# Fill baseline loss values using the initial loss for each loss function
combined_table <- combined_table |>
  group_by(loss) |>
  mutate(
    initial = if_else(
      method == "CC_MAGIC",
      first(initial[method != "CC_MAGIC"]),
      initial
    ),
    fitted = if_else(method == "CC_MAGIC", initial, fitted),
    improvement = if_else(method == "CC_MAGIC", 0, improvement),
    pct_change = if_else(method == "CC_MAGIC", 0, pct_change)
  ) |>
  ungroup()

# Format and save the combined table
combined_gt <- combined_table |>
  gt(groupname_col = "loss") |>
  tab_header(
    title = "Optimization results and fitted parameters",
    subtitle = "CC_MAGIC gives the initial parameter values and baseline loss"
  ) |>
  fmt_number(
    columns = c(initial, fitted, improvement),
    decimals = 4
  ) |>
  fmt_number(
    columns = pct_change,
    decimals = 2
  ) |>
  fmt_number(
    columns = c(runtime_seconds, all_of(parameter_cols)),
    decimals = 4
  ) |>
  cols_label(
    method = "Method",
    initial = "Initial loss",
    fitted = "Fitted loss",
    improvement = "Improvement",
    pct_change = "Improvement (%)",
    convergence = "Convergence",
    runtime_seconds = "Runtime (s)",
    v1 = "v1",
    v2 = "v2",
    sigma = "Sigma",
    migr = "Migration",
    portion = "Portion",
    whigh = "w_high",
    wlow = "w_low",
    tswitch = "Switch time"
  ) |>
  tab_options(table.font.size = px(11))

gtsave(
  combined_gt,
  file.path(output_dir, paste0("combined_fit_comparison", suffix, ".html"))
)

write_csv(
  combined_table,
  file.path(output_dir, paste0("combined_fit_comparison", suffix, ".csv"))
)

print(combined_gt)


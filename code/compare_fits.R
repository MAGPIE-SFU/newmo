# compare_fits.R

library(readr)
library(dplyr)
library(tidyr)
library(gt)

# Update this path if needed
output_dir <- "data/SouthAfrica/fit_comparison"

summary_file <- file.path(output_dir, "fit_summary.csv")
parameters_file <- file.path(output_dir, "all_fitted_parameters.csv")

fit_summary <- read_csv(summary_file, show_col_types = FALSE)
parameters <- read_csv(parameters_file, show_col_types = FALSE)

# ------------------------------------------------------------
# 1. Fit summary: improvement over CC_MAGIC baseline
# ------------------------------------------------------------

fit_summary <- fit_summary |>
  arrange(loss, fitted) |>
  mutate(
    improvement = initial - fitted,
    pct_change = 100 * improvement / abs(initial)
  )

summary_table <- fit_summary |>
  select(
    loss, method, initial, fitted, improvement,
    pct_change, convergence, runtime_seconds
  ) |>
  gt(groupname_col = "loss") |>
  tab_header(
    title = "Optimization results",
    subtitle = "Improvement relative to the initial CC_MAGIC loss"
  ) |>
  fmt_number(
    columns = c(initial, fitted, improvement),
    decimals = 4
  ) |>
  fmt_number(columns = pct_change, decimals = 2) |>
  fmt_number(columns = runtime_seconds, decimals = 2) |>
  cols_label(
    method = "Method",
    initial = "Initial loss",
    fitted = "Fitted loss",
    improvement = "Absolute improvement",
    pct_change = "Improvement (%)",
    convergence = "Convergence code",
    runtime_seconds = "Runtime (s)"
  ) |>
  data_color(
    columns = improvement,
    palette = c("white", "lightgreen")
  ) |>
  tab_options(table.font.size = px(12))

gtsave(summary_table, file.path(output_dir, "fit_summary.html"))
print(summary_table)


# ------------------------------------------------------------
# 2. Parameter comparison: one row per fit, parameters as columns
# ------------------------------------------------------------

parameter_cols <- setdiff(names(parameters), c("loss", "method"))

baseline <- readRDS(file.path(output_dir, "CC_MAGIC.rds")) |>
  select(all_of(parameter_cols)) |>
  slice(1)

if (nrow(baseline) != 1) {
  stop("Expected exactly one CC_MAGIC baseline row.")
}

# Add the baseline as a clearly labelled row for each loss function
baseline_rows <- expand_grid(
  loss = unique(parameters$loss[parameters$loss != "CC_MAGIC"])
) |>
  mutate(method = "CC_MAGIC")

for (p in parameter_cols) {
  baseline_rows[[p]] <- baseline[[p]][1]
}

parameter_table_data <- bind_rows(
  parameters |> filter(loss != "CC_MAGIC"),
  baseline_rows
) |>
  arrange(loss, method)

# Save the table data
write_csv(
  parameter_table_data,
  file.path(output_dir, "fitted_parameters_wide.csv")
)

# Display as a formatted table, grouped by loss function
parameter_table <- parameter_table_data |>
  gt(groupname_col = "loss") |>
  tab_header(
    title = "Fitted parameter values",
    subtitle = "Each optimizer compared with the CC_MAGIC baseline"
  ) |>
  fmt_number(
    columns = all_of(parameter_cols),
    decimals = 5
  ) |>
  cols_label(
    method = "Method",
    v1 = "v1",
    v2 = "v2",
    sigma = "Sigma",
    migr = "Migration",
    portion = "Portion",
    whigh = "w_high",
    wlow = "w_low",
    tswitch = "Switch time"
  ) |>
  tab_options(table.font.size = px(12))

gtsave(
  parameter_table,
  file.path(output_dir, "parameter_comparison.html")
)

print(parameter_table)


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
  file.path(output_dir, "combined_fit_comparison.html")
)

write_csv(
  combined_table,
  file.path(output_dir, "combined_fit_comparison.csv")
)

print(combined_gt)


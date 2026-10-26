# ============================================================
#  fit_report_functions.R
#  Turn fit results into a combined comparison table (gt + csv).
#  Works on results held in memory (from save_fit_results()) or
#  re-read from disk with read_fit_results().
# ============================================================

library(readr)
library(dplyr)
library(tidyr)
library(gt)

REPORT_ID_COLS <- c("fit_id", "loss", "method", "start_id")

# Re-load saved summary and parameter tables from disk.
read_fit_results <- function(output_dir, use_switch) {
  suffix <- if (use_switch) "_switch" else ""
  list(
    summary = read_csv(
      file.path(output_dir, paste0("fit_summary", suffix, ".csv")),
      show_col_types = FALSE
    ),
    parameters = read_csv(
      file.path(output_dir, paste0("all_fitted_parameters", suffix, ".csv")),
      show_col_types = FALSE
    )
  )
}

# Fit statistics + fitted parameters in one table, optionally with a
# CC_MAGIC baseline row per loss function.
#   baseline  : one-row data frame of parameters (e.g. from format_parameters())
#               or NULL for no baseline rows
#   best_only : keep only the best start per loss x method
build_combined_table <- function(fit_summary, parameters,
                                 baseline = NULL, best_only = FALSE) {
  parameter_cols <- setdiff(names(parameters), REPORT_ID_COLS)

  # Baseline loss = loss at the default start, taken before any filtering
  baseline_loss <- fit_summary |>
    filter(start_id == "default") |>
    group_by(loss) |>
    summarise(initial = first(initial), .groups = "drop")

  if (best_only) {
    fit_summary <- fit_summary |>
      filter(is.finite(fitted)) |>
      group_by(loss, method) |>
      slice_min(fitted, n = 1, with_ties = FALSE) |>
      ungroup()
  }

  fits <- fit_summary |>
    select(
      fit_id, loss, method, start_id, initial, fitted, improvement,
      pct_change, convergence, n_feval, runtime_seconds
    ) |>
    left_join(
      parameters |> select(fit_id, all_of(parameter_cols)),
      by = "fit_id"
    )

  if (is.null(baseline)) {
    return(fits |> arrange(loss, fitted))
  }

  baseline <- baseline |>
    select(all_of(parameter_cols)) |>
    slice(1)

  baseline_rows <- baseline_loss |>
    mutate(
      fit_id          = paste(loss, "CC_MAGIC", sep = "__"),
      method          = "CC_MAGIC",
      start_id        = NA_character_,
      fitted          = initial,
      improvement     = 0,
      pct_change      = 0,
      convergence     = NA_integer_,
      n_feval         = NA_integer_,
      runtime_seconds = NA_real_
    ) |>
    bind_cols(baseline)

  bind_rows(fits, baseline_rows) |>
    mutate(is_baseline = method == "CC_MAGIC") |>
    arrange(loss, !is_baseline, fitted) |>   # baseline first, then best fit first
    select(-is_baseline)
}

make_combined_gt <- function(combined_table, parameter_cols) {
  labels <- list(
    method = "Method", start_id = "Start",
    initial = "Initial loss", fitted = "Fitted loss",
    improvement = "Improvement", pct_change = "Improvement (%)",
    convergence = "Convergence", n_feval = "Evaluations",
    runtime_seconds = "Runtime (s)",
    v1 = "v1", v2 = "v2", sigma = "Sigma", migr = "Migration",
    portion = "Portion", whigh = "w_high", wlow = "w_low",
    tswitch = "Switch time"
  )
  labels <- labels[names(labels) %in% names(combined_table)]

  combined_table |>
    gt(groupname_col = "loss") |>
    tab_header(
      title = "Optimization results and fitted parameters",
      subtitle = "CC_MAGIC gives the initial parameter values and baseline loss"
    ) |>
    cols_hide(fit_id) |>
    fmt_number(columns = c(initial, fitted, improvement), decimals = 4) |>
    fmt_number(columns = pct_change, decimals = 2) |>
    fmt_number(columns = runtime_seconds, decimals = 1) |>
    fmt_number(columns = n_feval, decimals = 0) |>
    fmt_number(columns = any_of(setdiff(parameter_cols, "migr")), decimals = 4) |>
    fmt_scientific(columns = any_of("migr"), decimals = 2) |>   # migr is ~1e-5
    sub_missing(missing_text = "-") |>
    cols_label(.list = labels) |>
    tab_options(table.font.size = px(11))
}

# Build, save (html + csv) and return the table and gt object.
write_fit_report <- function(fit_summary, parameters, output_dir, use_switch,
                             baseline = NULL, best_only = FALSE, name = NULL) {
  suffix <- if (use_switch) "_switch" else ""
  if (is.null(name)) {
    name <- if (best_only) "combined_fit_comparison_best" else "combined_fit_comparison"
  }
  parameter_cols <- setdiff(names(parameters), REPORT_ID_COLS)

  combined_table <- build_combined_table(fit_summary, parameters, baseline, best_only)
  combined_gt    <- make_combined_gt(combined_table, parameter_cols)

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  gtsave(combined_gt, file.path(output_dir, paste0(name, suffix, ".html")))
  write_csv(combined_table, file.path(output_dir, paste0(name, suffix, ".csv")))

  invisible(list(table = combined_table, gt = combined_gt))
}
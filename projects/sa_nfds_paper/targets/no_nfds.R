# ============================================================
#  Supplementary: model fits without NFDS (sigma = 0)
#
#  With sigma = 0 the NFDS term in the model is zero, so `portion` (and the
#  locus weights) have no effect: only v1 is fitted. One fit per loss
#  function, on the same data and weights as the main fits (xt_fit).
#  v1 is scanned on a log grid within the main fitting bounds, then refined
#  with a 1D optimizer around the best grid point.
# ============================================================

NO_NFDS_N_GRID <- 200   # log-spaced v1 values scanned before refining

fit_no_nfds <- function(samodel, xt, loss_name, tswitch, nyears, use_switch, n_grid) {
  fit_data <- load_fit_data(xt, nyears)

  loss_at <- function(log_v1) {
    sim <- tryCatch(
      run_model(list(v1 = exp(log_v1), sigma = 0, portion = 0),
                use_switch, samodel, nyears, tswitch),
      error = function(e) NULL)
    if (is.null(sim)) return(Inf)
    v <- loss_functions[[loss_name]](sim, fit_data$data, fit_data$sd)
    if (is.finite(v)) v else Inf
  }

  # same bounds for v1 as the main fits
  lo <- log(PAR_LOWER[["v1"]])
  hi <- log(PAR_UPPER[["v1"]])

  grid <- seq(lo, hi, length.out = n_grid)
  vals <- vapply(grid, loss_at, numeric(1))
  i <- which.min(vals)

  # refine between the neighbours of the best grid point
  opt <- optimize(loss_at, c(grid[max(i - 1, 1)], grid[min(i + 1, n_grid)]))
  best_log_v1 <- if (opt$objective <= vals[i]) opt$minimum else grid[i]

  tibble::tibble(loss = loss_name, v1 = exp(best_log_v1), sigma = 0, portion = 0,
                 fitted = min(opt$objective, vals[i]))
}

no_nfds_simulate <- function(fits, samodel, nyears, tswitch, use_switch) {
  setNames(lapply(seq_len(nrow(fits)), function(i)
    run_model(list(v1 = fits$v1[i], sigma = 0, portion = 0),
              use_switch, samodel, nyears, tswitch)),
    paste("No NFDS:", fits$loss))
}

targets_no_nfds <- list(
  # one fit per loss function
  tar_target(no_nfds_loss, names(loss_functions)),
  tar_target(
    no_nfds_fit,
    fit_no_nfds(samodel, xt_fit, no_nfds_loss, T_SWITCH, NYEARS, FIT_USE_SWITCH, NO_NFDS_N_GRID),
    pattern = map(no_nfds_loss)
  ),

  # how much worse each loss is without NFDS than the best NFDS fit for that loss
  tar_target(
    no_nfds_comparison,
    dplyr::left_join(
      no_nfds_fit |> dplyr::select(loss, v1_no_nfds = v1, fitted_no_nfds = fitted),
      fit_selection$all_best |>
        dplyr::group_by(loss) |>
        dplyr::slice_min(fitted, n = 1, with_ties = FALSE) |>
        dplyr::ungroup() |>
        dplyr::select(loss, nfds_key = key, fitted_nfds = fitted),
      by = "loss") |>
      dplyr::mutate(loss_ratio = fitted_no_nfds / fitted_nfds)
  ),
  tar_target(no_nfds_comparison_csv,
             write_csv_output(no_nfds_comparison,
                              file.path(results_dir, "computed_data", "fits", "no_nfds_comparison.csv")),
             format = "file"),

  # figure: the overall best NFDS fit against the no-NFDS fits
  tar_target(sims_no_nfds,
             no_nfds_simulate(no_nfds_fit, samodel, NYEARS, T_SWITCH, FIT_USE_SWITCH)),
  tar_target(fig_no_nfds_pdf,
             save_plots_pdf(
               fit_figure_set(c(setNames(sim_overall_best, "With NFDS (overall best)"), sims_no_nfds),
                              samodel, xt_grouped_uncertainty,
                              "Fits with and without NFDS (sigma = 0)", FIG_SETTINGS),
               file.path(results_dir, "plots", "fits", "fit_figures_no_nfds.pdf"),
               width = 14, height = 10),
             format = "file")
)
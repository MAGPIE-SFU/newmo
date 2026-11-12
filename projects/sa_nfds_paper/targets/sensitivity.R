sens_dir <- file.path(results_dir, "plots", "sensitivity")

targets_sensitivity <- list(
  # Optimal fit per loss: the best across seeds
  tar_target(sens_opt_pars,
             fit_selection$all_best |>
               dplyr::group_by(loss) |>
               dplyr::slice_min(fitted, n = 1, with_ties = FALSE) |>
               dplyr::ungroup()),
  tar_target(sens_top_groups, sensitivity_top_groups(samodel, SENS_N_GROUPS)),
  tar_target(sensitivity,
             run_sensitivity(sens_opt_pars, samodel, tswitch = T_SWITCH, nyears = NYEARS,
                             sd_pars = SENS_SD_PARS, use_switch = FIT_USE_SWITCH,
                             n_sims = SENS_N_SIMS, n_cores = SENS_N_CORES)),
  tar_target(sensitivity_failed_csv,
             write_csv_output(sensitivity$failed,
                              file.path(results_dir, "computed_data", "sensitivity_failed.csv")),
             format = "file"),

  # The plots the script only displayed, as one multi-page PDF
  tar_target(sensitivity_overview_pdf,
             save_plots_pdf(list(
               plot_sensitivity_optimal(sensitivity, sens_top_groups),
               plot_sensitivity_draws(sensitivity),
               plot_sensitivity_ribbons(sensitivity, sens_top_groups, SENS_START_YEAR, SENS_QTILE),
               plot_sensitivity_trajectories(sensitivity, sens_top_groups, SENS_START_YEAR)
             ), file.path(sens_dir, "sensitivity_overview.pdf"), width = 9, height = 7),
             format = "file"),

  # The plots the script saved
  tar_target(sensitivity_by_loss_pdfs,
             unname(vapply(unique(sens_opt_pars$loss), function(l)
               save_plot_output(
                 plot_sensitivity_trajectories(sensitivity, sens_top_groups, SENS_START_YEAR, losses = l),
                 file.path(sens_dir, paste0("sim_trajectories_", l, ".pdf")), width = 9, height = 7),
               character(1))),
             format = "file"),
  tar_target(sensitivity_vs_data_pdf,
             save_plot_output(
               plot_sensitivity_vs_data(sensitivity, xt_grouped_uncertainty, sens_top_groups,
                                        SENS_START_YEAR, SENS_QTILE),
               file.path(sens_dir, paste0("sensitivity_top", SENS_N_GROUPS, "_ribbon", SENS_QTILE, ".pdf")),
               width = 9, height = 7),
             format = "file"),
  tar_target(sensitivity_vs_data_best_pdf,
             save_plot_output(
               plot_sensitivity_vs_data(sensitivity, xt_grouped_uncertainty, sens_top_groups,
                                        SENS_START_YEAR, SENS_QTILE, best_loss = SENS_BEST_LOSS),
               file.path(sens_dir, paste0("sensitivity_top", SENS_N_GROUPS, "_ribbon", SENS_QTILE,
                                          "_", SENS_BEST_LOSS, ".pdf")),
               width = 9, height = 7),
             format = "file")
)

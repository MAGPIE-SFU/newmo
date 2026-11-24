diag_dir  <- file.path(results_dir, "computed_data", "landscape")
diag_plot_height <- 3.2 * length(names(loss_functions))

targets_fit_diagnostics <- list(
  # 1D profiles through the best fit of every seed x loss
  tar_target(
    fit_profiles,
    profile_best_fits(fit_selection$all_best, samodel, xt_fit,
                      tswitch = T_SWITCH, nyears = NYEARS, use_switch = FIT_USE_SWITCH,
                      n_cores = DIAG_N_CORES, profile_decades = PROFILE_DECADES,
                      n_grid = PROFILE_N_GRID)
  ),
  tar_target(fit_profiles_rds,
             write_rds_output(fit_profiles, file.path(diag_dir, "profiles_best_by_seed.rds")),
             format = "file"),
  tar_target(fit_profiles_pdf,
             save_plot_output(plot_profiles(fit_profiles, relative_to = PROFILE_RELATIVE_TO,
                                            ymax = PROFILE_YMAX),
                              file.path(results_dir, "plots", "fits", "profiles_best_by_seed.pdf"),
                              width = 16, height = diag_plot_height),
             format = "file"),

  # Loss landscape from random starts
  tar_target(
    landscape_screen,
    screen_landscape(samodel, xt_fit, default_pars = DEFAULT_PARS,
                     seed = LANDSCAPE_SEED, tswitch = T_SWITCH, nyears = NYEARS,
                     use_switch = FIT_USE_SWITCH, n_screen = LANDSCAPE_N_SCREEN,
                     n_cores = DIAG_N_CORES)
  ),
  tar_target(landscape_screen_rds,
             write_rds_output(landscape_screen, file.path(diag_dir, "screen.rds")),
             format = "file"),
  tar_target(landscape_pdf,
             save_plot_output(plot_landscape(landscape_screen),
                              file.path(results_dir, "plots", "fits", "landscape_all.pdf"),
                              width = 16, height = diag_plot_height),
             format = "file"),
  tar_target(landscape_best_starts_csv,
             write_csv_output(landscape_best_starts(landscape_screen),
                              file.path(diag_dir, "landscape_best_starts.csv")),
             format = "file"),
  # Simulations
  tar_target(sim_cc_magic,
             simulate_nfds(samodel, cc_magic_modelpars(DEFAULT_PARS, T_SWITCH, FIT_USE_SWITCH), NYEARS)),
  tar_target(sims_all_best,
             simulate_fits(fit_selection$all_best, samodel, NYEARS, T_SWITCH, FIT_USE_SWITCH)),
  tar_target(sim_overall_best,
             simulate_fits(fit_selection$overall_best, samodel, NYEARS, T_SWITCH, FIT_USE_SWITCH)),
  # Figures
  tar_target(fig_cc_magic_pdf,
             save_plots_pdf(
               fit_figure_set(list("CC_MAGIC" = sim_cc_magic), samodel, xt_grouped_uncertainty,
                              "CC_MAGIC default parameters", FIG_SETTINGS),
               file.path(results_dir, "plots", "fits", "fit_figures_cc_magic.pdf"), width = 14, height = 10),
             format = "file"),
  tar_target(fig_all_fits_pdf,
             save_plots_pdf(
               all_fits_figure_set(sims_all_best, fit_selection$all_best, samodel,
                                   xt_grouped_uncertainty, FIG_SETTINGS),
               file.path(results_dir, "plots", "fits", "fit_figures_all_fits.pdf"), width = 14, height = 10),
             format = "file"),
  tar_target(fig_overall_best_pdf,
             save_plots_pdf(
               fit_figure_set(setNames(sim_overall_best, "Overall best"), samodel,
                              xt_grouped_uncertainty,
                              with(fit_selection$overall_best,
                                   sprintf("Overall best: seed %s, %s (%s, %s)",
                                           seed, loss, method, start_id)),
                              FIG_SETTINGS),
               file.path(results_dir, "plots", "fits", "fit_figures_overall_best.pdf"), width = 14, height = 10),
             format = "file")
)

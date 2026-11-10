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
                              file.path(results_dir, "plots", "profiles_best_by_seed.pdf"),
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
                              file.path(results_dir, "plots", "landscape_all.pdf"),
                              width = 16, height = diag_plot_height),
             format = "file"),
  tar_target(landscape_best_starts_csv,
             write_csv_output(landscape_best_starts(landscape_screen),
                              file.path(diag_dir, "landscape_best_starts.csv")),
             format = "file")
)

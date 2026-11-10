fits_dir <- file.path(results_dir, "computed_data", "fits")

targets_fits <- list(
  tar_target(fit_seed, FIT_SEEDS),
  tar_target(
    seq_fit,
    run_sequential_fits(
      samodel, xt_grouped_uncertainty,
      seed = fit_seed, default_pars = DEFAULT_PARS,
      tswitch = T_SWITCH, nyears = NYEARS,
      use_switch = FIT_USE_SWITCH, use_log = FIT_USE_LOG,
      chain_methods = FIT_CHAIN_METHODS, n_rounds = FIT_N_ROUNDS,
      start_from = FIT_START_FROM, n_screen = FIT_N_SCREEN,
      maxfeval = FIT_MAXFEVAL, rel_tol = FIT_REL_TOL,
      parallel = FIT_PARALLEL, n_workers = FIT_N_WORKERS
    ),
    pattern = map(fit_seed),
    iteration = "list"
  ),
  tar_target(seq_fit_files,
             write_seed_fit_outputs(seq_fit, fits_dir, FIT_USE_SWITCH),
             pattern = map(seq_fit), format = "file"),
  tar_target(
    fit_selection,
    select_best_fits(seq_fit, samodel, xt_grouped_uncertainty,
                     tswitch = T_SWITCH, nyears = NYEARS, use_switch = FIT_USE_SWITCH)
  ),
  tar_target(all_fits_rds,
             write_rds_output(fit_selection$all_parameters, file.path(fits_dir, "all_fitted_parameters.rds")),
             format = "file"),
  tar_target(all_best_fits_rds,
             write_rds_output(fit_selection$all_best, file.path(fits_dir, "all_best_fitted_parameters.rds")),
             format = "file"),
  tar_target(overall_best_fit_rds,
             write_rds_output(fit_selection$overall_best, file.path(fits_dir, "overall_best_fitted_parameters.rds")),
             format = "file")
)

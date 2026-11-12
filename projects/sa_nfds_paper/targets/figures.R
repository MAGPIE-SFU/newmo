fig_dir <- file.path(results_dir, "plots")

targets_figures <- list(
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
               file.path(fig_dir, "fit_figures_cc_magic.pdf"), width = 14, height = 10),
             format = "file"),
  tar_target(fig_all_fits_pdf,
             save_plots_pdf(
               all_fits_figure_set(sims_all_best, fit_selection$all_best, samodel,
                                   xt_grouped_uncertainty, FIG_SETTINGS),
               file.path(fig_dir, "fit_figures_all_fits.pdf"), width = 14, height = 10),
             format = "file"),
  tar_target(fig_overall_best_pdf,
             save_plots_pdf(
               fit_figure_set(setNames(sim_overall_best, "Overall best"), samodel,
                              xt_grouped_uncertainty,
                              with(fit_selection$overall_best,
                                   sprintf("Overall best: seed %s, %s (%s, %s)",
                                           seed, loss, method, start_id)),
                              FIG_SETTINGS),
               file.path(fig_dir, "fit_figures_overall_best.pdf"), width = 14, height = 10),
             format = "file")
)

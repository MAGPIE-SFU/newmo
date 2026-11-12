ipd_data_dir <- file.path(results_dir, "computed_data", "projected_ipd")
ipd_plot_dir <- file.path(results_dir, "plots", "projected_ipd")

targets_projected_ipd <- list(
  tar_target(
    projected_ipd,
    simulate_projected_ipd(sim_overall_best[[1]], samodel, invasiveness_lookup,
                           group_lookup, meta_clean, non_typeable = NON_TYPEABLE,
                           n_reps = IPD_N_REPS, n_carriage = IPD_N_CARRIAGE)
  ),
  tar_target(projected_ipd_csvs, c(
    write_csv_output(projected_ipd$ipd_summary,          file.path(ipd_data_dir, "simulated_ipd_uncertainty.csv")),
    write_csv_output(projected_ipd$ipd_yearly,           file.path(ipd_data_dir, "simulated_ipd_yearly_uncertainty.csv")),
    write_csv_output(projected_ipd$ipd_yearly_draws,     file.path(ipd_data_dir, "simulated_ipd_yearly_draws.csv")),
    write_csv_output(projected_ipd$ipd_serotype_summary, file.path(ipd_data_dir, "simulated_ipd_serotype_uncertainty.csv"))
  ), format = "file"),

  # Plots
  tar_target(ipd_sero_plot_data,
             prepare_ipd_serotypes(projected_ipd$ipd_serotype_summary, NON_TYPEABLE)),
  tar_target(ipd_check_pdf,
             save_plot_output(
               plot_projected_ipd_check(projected_ipd$ipd_serotype_summary, meta_clean,
                                        FIRST_SIM_CALENDAR_YEAR),
               file.path(ipd_plot_dir, "ipd_observed_vs_predicted.pdf"),
               width = 16, height = 12),
             format = "file"),
  tar_target(ipd_total_pdf,
             save_plot_output(
               plot_ipd_total(ipd_sero_plot_data),
               file.path(ipd_plot_dir, "ipd_serotype_total.pdf"),
               width = 11,
               height = max(6, 2 + 0.18 * dplyr::n_distinct(ipd_sero_plot_data$serotype))),
             format = "file"),
  tar_target(ipd_dotplot_pdf,
             save_plot_output(
               plot_ipd_dotplot(ipd_sero_plot_data, top_n = IPD_DOT_TOP_N),
               file.path(ipd_plot_dir, "ipd_serotype_dotplot.pdf"),
               width = 11,
               height = max(6, 3 + 0.31 * min(IPD_DOT_TOP_N,
                                              dplyr::n_distinct(ipd_sero_plot_data$serotype)))),
             format = "file")
)

vd_data_dir <- file.path(results_dir, "computed_data", "vaccine_design")

targets_vaccine_design <- list(
  tar_target(vd_modelpars, fit_modelpars(fit_selection$overall_best, tswitch = 0)),
  tar_target(
    designed_vaccines,
    design_vaccines(vd_modelpars, samodel, xt_grouped_uncertainty, invasiveness_lookup,
                    group_lookup, base_vaccine = VACCINE_TYPES[["PCV7"]],
                    prevax_years = PREVAX_YEARS, target_size = VD_TARGET_SIZE,
                    n_candidates = VD_N_CANDIDATES, years = VD_SCORE_YEARS,
                    n_cores = VD_N_CORES)
  ),
  tar_target(designed_vaccines_rds,
             write_rds_output(designed_vaccines, file.path(vd_data_dir, "vaccines.rds")),
             format = "file"),
    # Score the licensed vaccines the same way as the invented ones
  tar_target(
    licensed_vaccine_ipd,
    sapply(HEATMAP_VACCINES, function(v)
      score_vaccine(vd_modelpars, VACCINE_TYPES[[v]], samodel, invasiveness_lookup,
                    group_lookup, years = VD_SCORE_YEARS))
  ),
  tar_target(
    vaccine_heatmap,
    plot_vaccine_heatmap(VACCINE_TYPES[HEATMAP_VACCINES], licensed_vaccine_ipd,
                         designed_vaccines, samodel$st,
                         n_exhaustive = HEATMAP_N_EXHAUSTIVE,
                         method_labels = HEATMAP_METHOD_LABELS)
  ),
  tar_target(vaccine_heatmap_pdf,
             save_plot_output(vaccine_heatmap$plot,
                              file.path(results_dir, "plots", "vaccine_design", "vaccine_heatmap.pdf"),
                              width = vaccine_heatmap$width, height = vaccine_heatmap$height),
             format = "file")
)

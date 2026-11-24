vd_eval_data_dir <- file.path(
  results_dir,
  "computed_data",
  "vaccine_proj_ipd"
)

targets_vax_proj_ipd <- list(

  # Select vaccine compositions
  tar_target(
    selected_vaccines,
    select_best_vaccines(
      readRDS(designed_vaccines_rds),
      target_size = VD_TARGET_SIZE,
      licensed_vaccines = VACCINE_TYPES["PCV13"]
    ),
    iteration = "list"
  ),

  # Evaluate each vaccine
  tar_target(
    selected_vaccine_results,
    evaluate_vaccine_ipd(
      vaccine = selected_vaccines,
      modelpars = vd_modelpars,
      samodel = samodel,
      invasiveness_lookup = invasiveness_lookup,
      group_lookup = group_lookup,
      meta_clean = meta_clean,
      non_typeable = NON_TYPEABLE,
      years = 5L,
      n_reps = 1000L,
      n_carriage = 100000,
      seed = 1L
    ),
    pattern = map(selected_vaccines),
    iteration = "list"
  ),

  # Combine cumulative IPD results
  tar_target(
    vaccine_ipd_composition,
    dplyr::bind_rows(
      setNames(
        purrr::map(selected_vaccine_results, "cumulative"),
        names(selected_vaccines)
      ),
      .id = "vaccine"
    )
  ),

  # Combine annual IPD results
  tar_target(
    vaccine_ipd_yearly,
    dplyr::bind_rows(
      setNames(
        purrr::map(selected_vaccine_results, "yearly"),
        names(selected_vaccines)
      ),
      .id = "vaccine"
    )
  ),

  # Save cumulative IPD results
  tar_target(
    vaccine_ipd_composition_rds,
    write_rds_output(
      vaccine_ipd_composition,
      file.path(
        vd_eval_data_dir,
        "best_vax_proj_ipd.rds"
      )
    ),
    format = "file"
  ),

  # Save annual IPD results
  tar_target(
    vaccine_ipd_yearly_rds,
    write_rds_output(
      vaccine_ipd_yearly,
      file.path(
        vd_eval_data_dir,
        "best_vax_proj_ipd_yearly.rds"
      )
    ),
    format = "file"
  ),

  # Plot projected cumulative IPD composition
  tar_target(
    vaccine_ipd_dotplot,
    plot_vaccine_ipd_dotplot(
      vaccine_ipd_composition,
      invasiveness_lookup,
      top_n = IPD_DOT_TOP_N
    )
  ),

  # Save dot plot
  tar_target(
    vaccine_ipd_dotplot_pdf,
    save_plot_output(
      vaccine_ipd_dotplot,
      file.path(
        results_dir,
        "plots",
        "vax_proj_ipd",
        "vaccine_ipd_dotplot.pdf"
      ),
      width = 12,
      height = 11
    ),
    format = "file"
  ),

  # Plot annual total IPD and serotype-specific contributions
  tar_target(
    vaccine_total_ipd_plot,
    plot_vaccine_total_ipd(
      vaccine_ipd_yearly,
      invasiveness_lookup,
      top_n = IPD_DOT_TOP_N
    )
  ),

  # Save stacked-area plot
  tar_target(
    vaccine_total_ipd_pdf,
    save_plot_output(
      vaccine_total_ipd_plot,
      file.path(
        results_dir,
        "plots",
        "vax_proj_ipd",
        "vaccine_total_ipd.pdf"
      ),
      width = 13,
      height = 6
    ),
    format = "file"
  )
)
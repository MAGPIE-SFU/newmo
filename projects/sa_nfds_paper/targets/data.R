targets_data <- list(
  tar_target(meta_raw_file, META_RAW_PATH, format = "file"),
  tar_target(
    meta_clean,
    clean_metadata(
      readr::read_csv(meta_raw_file, show_col_types = FALSE),
      serotype_renames = SEROTYPE_RENAMES
    )
  ),
  tar_target(
    meta_clean_csv,
    write_csv_output(meta_clean, file.path(results_dir, "computed_data", "meta_clean.csv")),
    format = "file"
  ),
  tar_target(inv_children_file, INV_CHILDREN_PATH, format = "file"),
  tar_target(inv_adults_file,   INV_ADULTS_PATH,   format = "file"),
  tar_target(
    invasiveness_lookup,
    build_invasiveness_lookup(
      meta_clean,
      openxlsx::read.xlsx(inv_children_file),
      openxlsx::read.xlsx(inv_adults_file),
      INVASIVENESS_DONORS,
      NON_TYPEABLE
    )
  ),
  tar_target(
    invasiveness_lookup_csv,
    write_csv_output(invasiveness_lookup,
                     file.path(results_dir, "computed_data", "invasiveness_lookup_interval.csv")),
    format = "file"
  ),
    tar_target(
    reweighted,
    reweight_carriage(
      meta_clean, invasiveness_lookup,
      tauA = REWEIGHT_TAU_A, tauU = REWEIGHT_TAU_U, k = REWEIGHT_K,
      prevax_years = PREVAX_YEARS,
      min_sero_frac = MIN_SERO_FRAC, min_sg_frac = MIN_SG_FRAC
    )
  ),
  tar_target(group_lookup, reweighted$group_lookup),
  tar_target(group_lookup_csv,
             write_csv_output(group_lookup, file.path(results_dir, "computed_data", "group_lookup.csv")),
             format = "file"),
  tar_target(meta_clean_groups_csv,
             write_csv_output(reweighted$meta_clean_groups, file.path(results_dir, "computed_data", "meta_clean_groups.csv")),
             format = "file"),
  tar_target(
    xt_grouped_uncertainty,
    reweight_carriage_uncertainty(
      meta_clean, invasiveness_lookup, group_lookup,
      tauA = REWEIGHT_TAU_A, tauU = REWEIGHT_TAU_U, k = REWEIGHT_K,
      lambda = XT_LAMBDA, n_reps = XT_N_REPS
    )
  ),
  tar_target(xt_grouped_uncertainty_csv,
             write_csv_output(xt_grouped_uncertainty,
                              file.path(results_dir, "computed_data", "xt_grouped_uncertainty.csv")),
             format = "file"),
  tar_target(xt_top_groups_pdf,
             save_plot_output(plot_xt_top_groups(xt_grouped_uncertainty),
                              file.path(results_dir, "plots", "xt_top_groups.pdf"),
                              width = 12, height = 9),
             format = "file")
)

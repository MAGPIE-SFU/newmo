targets_data <- list(
  tar_target(meta_raw_file, META_RAW_PATH, format = "file"),
  tar_target(meta_clean, clean_metadata(readr::read_csv(meta_raw_file, show_col_types = FALSE))),
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
      SEROTYPE_SUBSTITUTIONS,
      NON_TYPEABLE
    )
  ),
  tar_target(
    invasiveness_lookup_csv,
    write_csv_output(invasiveness_lookup,
                     file.path(results_dir, "computed_data", "invasiveness_lookup_interval.csv")),
    format = "file"
  )
)


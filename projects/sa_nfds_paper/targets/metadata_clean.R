targets_metadata_clean <- list(
  tar_target(meta_raw_file, META_RAW_PATH, format = "file"),
  tar_target(meta_clean, clean_metadata(readr::read_csv(meta_raw_file, show_col_types = FALSE))),
  tar_target(
    meta_clean_csv,
    write_csv_output(meta_clean, file.path(results_dir, "computed_data", "meta_clean.csv")),
    format = "file"
  )
)

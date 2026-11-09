targets_vaccines <- list(
  tar_target(vaccines, VACCINE_TYPES),
  tar_target(
    vaccine_table_csv,
    write_vaccine_table(vaccines, file.path(results_dir, "computed_data", "vaccine_serotypes.csv")),
    format = "file"
  ),
  tar_target(
    vaccine_coverage_pdf,
    plot_vaccine_coverage(vaccines, file.path(results_dir, "plots", "vaccine_coverage.pdf")),
    format = "file"
  )
)
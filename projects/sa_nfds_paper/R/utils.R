write_csv_output <- function(x, out_file) {
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(x, out_file)
  out_file
}

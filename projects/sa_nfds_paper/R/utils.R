write_csv_output <- function(x, out_file) {
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(x, out_file)
  out_file
}

save_plot_output <- function(p, out_file, width, height) {
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(out_file, p, width = width, height = height, bg = "white")
  out_file
}

write_matrix_output <- function(x, out_file) {
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, out_file)
  out_file
}

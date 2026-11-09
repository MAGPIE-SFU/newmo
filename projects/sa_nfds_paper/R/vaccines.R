vaccines_long <- function(vaccines) {
  data.frame(
    vaccine  = rep(names(vaccines), lengths(vaccines)),
    serotype = unlist(vaccines, use.names = FALSE)
  )
}

write_vaccine_table <- function(vaccines, out_file) {
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(vaccines_long(vaccines), out_file, row.names = FALSE)
  out_file
}

plot_vaccine_coverage <- function(vaccines, out_file) {
  tab <- vaccines_long(vaccines)
  tab$serotype <- factor(tab$serotype, levels = unique(tab$serotype))
  tab$vaccine  <- factor(tab$vaccine, levels = rev(names(vaccines)))

  p <- ggplot(tab, aes(serotype, vaccine)) +
    geom_tile(fill = "#1F6F8B", colour = "white") +
    labs(x = "Serotype", y = NULL) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  ggsave(out_file, p, width = 10, height = 3, bg = "white")
  out_file
}
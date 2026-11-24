# Common look for all paper figures. A partial theme: it modifies each plot's
# own theme instead of replacing it, so plot-specific settings stay.
theme_paper <- function(base_size = 11) {
  ggplot2::theme(
    text             = ggplot2::element_text(size = base_size),
    plot.title       = ggplot2::element_text(face = "bold", size = base_size + 2),
    panel.background = ggplot2::element_rect(fill = "white", colour = NA),
    panel.border     = ggplot2::element_rect(colour = "grey30", fill = NA),
    # panel.grid.minor = ggplot2::element_blank(),  # turns off minor grid
    strip.background = ggplot2::element_rect(fill = "grey90", colour = "grey30"),
    legend.position  = "bottom"
  )
}

# Plots are ggplots or patchworks (heatmap, vaccine candidates); patchwork needs &
apply_theme_paper <- function(p) {
  if (inherits(p, "patchwork")) p & theme_paper() else p + theme_paper()
}

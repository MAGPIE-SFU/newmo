# Common look for all paper figures. A partial theme: it modifies each plot's
# own theme instead of replacing it, so plot-specific settings stay.
theme_paper <- function(base_size = 11) {
  ggplot2::theme(
    text = ggplot2::element_text(size = base_size),
    plot.title = ggplot2::element_text(face = "bold", size = base_size + 2),
    panel.background = ggplot2::element_rect(fill = "white", colour = NA),
    panel.border = ggplot2::element_rect(colour = "grey30", fill = NA),
    # panel.grid.minor = ggplot2::element_blank(),  # turns off minor grid
    strip.background = ggplot2::element_rect(
      fill = "grey90",
      colour = "grey30"
    ),
    legend.position = "bottom"
  )
}

# Plots are ggplots or patchworks (heatmap, vaccine candidates); patchwork needs &
apply_theme_paper <- function(p) {
  if (inherits(p, "patchwork")) p & theme_paper() else p + theme_paper()
}

#' Common plotting theme for the NFDS paper
#'
#' Provides consistent formatting for publication figures, including
#' typography, panel borders, gridlines, facet labels and legends.
#'
#' @param base_size Base font size in points.
#' @param base_family Base font family.
#' @param legend_position Position of the legend.
#'
#' @return A ggplot2 theme object.
#' @export
theme_paper2 <- function(
    base_size = 12,
    base_family = "sans",
    legend_position = "bottom"
) {

  ggplot2::theme_bw(
    base_size = base_size,
    base_family = base_family
  ) +
    ggplot2::theme(

      # Panel appearance
      panel.background = ggplot2::element_rect(
        fill = "white",
        colour = NA
      ),
      panel.border = ggplot2::element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.5
      ),
      panel.grid.major = ggplot2::element_line(
        colour = "grey90",
        linewidth = 0.3
      ),
      panel.grid.minor = ggplot2::element_blank(),
      panel.spacing = grid::unit(0.8, "lines"),

      # Axes
      axis.title = ggplot2::element_text(
        size = ggplot2::rel(1),
        colour = "black"
      ),
      axis.text = ggplot2::element_text(
        size = ggplot2::rel(0.9),
        colour = "black"
      ),
      axis.ticks = ggplot2::element_line(
        colour = "black",
        linewidth = 0.4
      ),
      axis.ticks.length = grid::unit(2, "pt"),

      # Facet strips
      strip.background = ggplot2::element_rect(
        fill = "grey95",
        colour = "black",
        linewidth = 0.5
      ),
      strip.text = ggplot2::element_text(
        size = ggplot2::rel(1),
        face = "bold",
        colour = "black",
        margin = ggplot2::margin(5, 5, 5, 5)
      ),

      # Legend
      legend.position = legend_position,
      legend.background = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      legend.key.width = grid::unit(0.65, "cm"),
      legend.key.height = grid::unit(0.4, "cm"),
      legend.spacing.x = grid::unit(0.4, "cm"),
      legend.title = ggplot2::element_text(
        size = ggplot2::rel(0.95)
      ),
      legend.text = ggplot2::element_text(
        size = ggplot2::rel(0.85)
      ),

      # Titles and captions
      plot.title = ggplot2::element_text(
        size = ggplot2::rel(1.2),
        face = "bold",
        hjust = 0
      ),
      plot.subtitle = ggplot2::element_text(
        size = ggplot2::rel(0.95),
        colour = "grey30",
        margin = ggplot2::margin(b = 8)
      ),
      plot.caption = ggplot2::element_text(
        size = ggplot2::rel(0.8),
        colour = "grey40",
        hjust = 0
      ),
      plot.tag = ggplot2::element_text(
        size = ggplot2::rel(1.2),
        face = "bold"
      ),

      # Overall spacing
      plot.background = ggplot2::element_rect(
        fill = "white",
        colour = NA
      ),
      plot.margin = ggplot2::margin(8, 8, 8, 8)
    )
}

#' Generate discrete blue colours for the NFDS paper
#'
#' Generates alternating light and dark shades of blue to improve
#' differentiation between adjacent categories.
#'
#' @param n Number of colours.
#'
#' @return Character vector of hexadecimal colours.
#' @export
palette_paper_blue <- function(n) {

  if (n == 0L) {
    return(character(0))
  }

  pal <- grDevices::colorRampPalette(
    c(
      "#D6EAF8",
      "#85B9DF",
      "#3182BD",
      "#08519C",
      "#08306B"
    )
  )(n)

  idx <- as.vector(rbind(
    seq_len(ceiling(n / 2)),
    rev(seq_len(n))[seq_len(ceiling(n / 2))]
  ))

  pal[unique(idx)]
}
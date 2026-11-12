# ============================================================
# Predicted South African IPD:
# total IPD with serotype contributions, and uncertainty dot plot
# ============================================================

#' Prepare the predicted IPD serotype composition for plotting
#'
#' Checks the columns, drops non-typeables and missing values, and orders the
#' serotypes by overall predicted IPD burden.
#'
#' @param ipd_serotype_summary From [simulate_projected_ipd()].
#' @param non_typeable Serotype labels to leave out.
#' @return The filtered table, with serotype as a factor (largest burden last).
#' @importFrom dplyr %>% mutate filter group_by summarise arrange desc pull
#' @export
prepare_ipd_serotypes <- function(ipd_serotype_summary, non_typeable) {

  required <- c(
    "year", "serotype", "ipd_mean",
    "frac", "frac_lower", "frac_upper"
  )

  missing_cols <- setdiff(required, names(ipd_serotype_summary))

  if (length(missing_cols) > 0) {
    stop(
      "Missing columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  ipd_serotype_summary <- ipd_serotype_summary %>%
    mutate(
      serotype = as.character(serotype),
      year = as.integer(year)
    ) %>%
    filter(
      !is.na(serotype),
      !is.na(year),
      is.finite(ipd_mean),
      !serotype %in% non_typeable
    )

  # Order serotypes by overall predicted IPD burden
  sero_order <- ipd_serotype_summary %>%
    group_by(serotype) %>%
    summarise(
      total = sum(ipd_mean),
      .groups = "drop"
    ) %>%
    arrange(desc(total)) %>%
    pull(serotype)

  ipd_serotype_summary %>%
    mutate(
      serotype = factor(
        serotype,
        levels = rev(sero_order)
      )
    )
}

# ============================================================
# Panel A: Total predicted IPD with serotype contributions
# ============================================================

#' Total predicted IPD with serotype contributions
#'
#' @param ipd_sero From [prepare_ipd_serotypes()].
#' @return A ggplot.
#' @importFrom dplyr %>% group_by summarise pull mutate if_else arrange desc
#' @importFrom ggplot2 ggplot geom_area aes geom_line geom_point scale_fill_viridis_d
#'   scale_x_continuous scale_y_continuous expansion labs theme_bw theme
#'   element_blank element_text
#' @export
plot_ipd_total <- function(ipd_sero) {

  # Identify serotypes contributing the most IPD overall
  top_seros <- ipd_sero %>%
    group_by(serotype) %>%
    summarise(
      total_ipd = sum(ipd_mean),
      .groups = "drop"
    ) %>%
    pull(serotype)

  # Prepare stacked IPD data
  ipd_stacked <- ipd_sero %>%
    mutate(
      serotype_plot = if_else(
        serotype %in% top_seros,
        as.character(serotype),
        "Other"
      )
    ) %>%
    group_by(year, serotype_plot) %>%
    summarise(
      ipd_mean = sum(ipd_mean),
      .groups = "drop"
    )

  # Order serotypes by total contribution
  sero_order <- ipd_stacked %>%
    group_by(serotype_plot) %>%
    summarise(
      total = sum(ipd_mean),
      .groups = "drop"
    ) %>%
    arrange(desc(total)) %>%
    pull(serotype_plot)

  ipd_stacked <- ipd_stacked %>%
    mutate(
      serotype_plot = factor(
        serotype_plot,
        levels = sero_order
      )
    )

  # Calculate the total predicted IPD for each year
  total_ipd <- ipd_stacked %>%
    group_by(year) %>%
    summarise(
      total = sum(ipd_mean),
      .groups = "drop"
    )

  ggplot() +
    geom_area(
      data = ipd_stacked,
      aes(
        x = year,
        y = ipd_mean,
        fill = serotype_plot
      ),
      position = "stack",
      alpha = 0.85
    ) +
    geom_line(
      data = total_ipd,
      aes(x = year, y = total),
      colour = "black",
      linewidth = 1.1
    ) +
    geom_point(
      data = total_ipd,
      aes(x = year, y = total),
      colour = "black",
      size = 2
    ) +
    scale_fill_viridis_d(
      option = "mako",
      begin = 0.15,
      end = 0.9
    ) +
    scale_x_continuous(
      breaks = sort(unique(ipd_stacked$year))
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.05))
    ) +
    labs(
      x = "Simulation year",
      y = "Expected IPD cases per 100,000 carriage events",
      fill = "Serotype",
      title = "Predicted invasive pneumococcal disease",
      subtitle = "Total IPD and serotype-specific contributions"
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "right",
      plot.title = element_text(face = "bold")
    )
}

# ============================================================
# PANEL B: Dot plot with 95% uncertainty intervals
# ============================================================

#' Predicted serotype contribution to IPD, with uncertainty intervals
#'
#' Intervals represent uncertainty in invasiveness, not Poisson predictive
#' intervals. All years in one panel, coloured and dodged by year.
#'
#' @param ipd_sero From [prepare_ipd_serotypes()].
#' @param top_n Number of serotypes (by total expected IPD) to show.
#' @return A ggplot.
#' @importFrom dplyr %>% group_by summarise arrange desc pull filter mutate
#' @importFrom ggplot2 ggplot aes geom_errorbar geom_point position_dodge
#'   scale_colour_viridis_d scale_x_continuous labs theme_bw theme element_blank
#'   element_text
#' @export
plot_ipd_dotplot <- function(ipd_sero, top_n = 20) {

  # Identify top serotypes by total expected IPD
  sero_order <- ipd_sero %>%
    group_by(serotype) %>%
    summarise(
      total_ipd = sum(ipd_mean, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(total_ipd)) %>%
    pull(serotype)

  top_seros <- head(sero_order, top_n)

  dot_data <- ipd_sero %>%
    filter(serotype %in% top_seros) %>%
    mutate(
      serotype = factor(serotype, levels = rev(top_seros)),
      year = factor(year, levels = sort(unique(year)))
    )

  dodge <- position_dodge(
    width = 0.75,
    reverse = TRUE
  )

  ggplot(
    dot_data,
    aes(
      x = 100 * frac,
      y = serotype,
      colour = year,
      group = year
    )
  ) +
    geom_errorbar(
      aes(
        xmin = 100 * frac_lower,
        xmax = 100 * frac_upper
      ),
      orientation = "y",
      position = dodge,
      width = 0.4,
      linewidth = 0.6,
      alpha = 0.35
    ) +
    geom_point(
      position = dodge,
      size = 2.2,
      alpha = 0.9
    ) +
    scale_colour_viridis_d(
      option = "mako",
      begin = 0.8,
      end = 0.1,
      direction = 1
    ) +
    scale_x_continuous(
      labels = scales::label_number(suffix = "%"),
      limits = c(0, NA)
    ) +
    labs(
      x = "Proportion of total predicted IPD (%)",
      y = "Serotype",
      colour = "Simulation year",
      title = "Predicted serotype contribution to IPD",
      subtitle = paste0(
        "Top ", length(top_seros),
        " serotypes by total expected IPD"
      )
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "right",
      plot.title = element_text(face = "bold")
    )
}

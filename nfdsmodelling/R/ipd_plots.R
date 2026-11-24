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
    "year",
    "serotype",
    "ipd_mean",
    "frac",
    "frac_lower",
    "frac_upper"
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


#' Predicted serotype contribution to IPD, with uncertainty intervals
#' Compares projected serotype-specific IPD composition across vaccines.
#' Serotypes are selected by total expected IPD and ordered by child
#' invasiveness.
#'
#' @param dat Data frame containing cumulative IPD projections by vaccine
#'   and serotype.
#' @param invasiveness_lookup Data frame containing serotype-specific
#'   child invasiveness estimates (`v_child`).
#' @param top_n Number of serotypes to display, selected by total
#'   expected IPD. If NULL, all serotypes are displayed.
#'
#' @return A ggplot object.
#' @export
plot_vaccine_ipd_dotplot <- function(
    dat,
    invasiveness_lookup,
    top_n = 20L
) {

  # Select serotypes by total expected IPD
  serotype_ranks <- dat %>%
    dplyr::group_by(serotype) %>%
    dplyr::summarise(
      total_ipd = sum(ipd_mean),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(total_ipd))

  top_seros <- serotype_ranks$serotype

  if (!is.null(top_n)) {
    top_seros <- head(top_seros, top_n)
  }

  # Order selected serotypes by child invasiveness
  serotype_order <- invasiveness_lookup %>%
    dplyr::filter(serotype %in% top_seros) %>%
    dplyr::arrange(v_child) %>%
    dplyr::pull(serotype)

  # Prepare plotting data
  dot_data <- dat %>%
    dplyr::filter(serotype %in% top_seros) %>%
    dplyr::mutate(
      serotype = factor(
        serotype,
        levels = serotype_order
      ),
      vaccine = dplyr::recode(
        vaccine,
        "best_greedy" = "Greedy",
        "best_exhaustive" = "Exhaustive"
      ),
      vaccine = factor(vaccine)
    )

  dodge <- ggplot2::position_dodge(width = 0.7)

  ggplot2::ggplot(
    dot_data,
    ggplot2::aes(
      x = 100 * frac,
      y = serotype,
      colour = vaccine,
      group = vaccine
    )
  ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        xmin = 100 * frac_lower,
        xmax = 100 * frac_upper
      ),
      orientation = "y",
      position = dodge,
      width = 0.4,
      linewidth = 0.6,
      alpha = 0.35
    ) +
    ggplot2::geom_point(
      position = dodge,
      size = 2.2,
      alpha = 0.9
    ) +
    ggplot2::scale_colour_viridis_d(
      option = "mako",
      begin = 0.8,
      end = 0.1,
      direction = 1
    ) +
    ggplot2::scale_x_continuous(
      labels = scales::label_number(suffix = "%"),
      limits = c(0, NA)
    ) +
    ggplot2::labs(
      x = "Proportion of total predicted IPD (%)",
      y = "Serotype",
      colour = "Vaccine",
      title = "Predicted serotype contribution to IPD",
      subtitle = paste0(
        "Top ",
        length(top_seros),
        " serotypes by expected IPD, ordered by invasiveness"
      )
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "right",
      plot.title = ggplot2::element_text(face = "bold")
    )
}


#' Plot projected total IPD across vaccine formulations
#'
#' Produces stacked-area plots of annual expected IPD, with
#' serotype-specific contributions and a line indicating total IPD.
#' Each vaccine is displayed in a separate panel with common axis scales.
#'
#' @param dat Data frame containing annual IPD projections by vaccine,
#'   year and serotype, including `ipd_mean`.
#' @param invasiveness_lookup Data frame containing serotype-specific
#'   child invasiveness estimates (`v_child`).
#' @param top_n Number of serotypes displayed individually, ranked by
#'   total expected IPD across vaccines and years. Remaining serotypes
#'   are combined into an "Other" category.
#'
#' @return A ggplot object.
#' @export
plot_vaccine_total_ipd <- function(
    dat,
    invasiveness_lookup,
    top_n = 20L
) {

  # 1. Select the top serotypes ------------------------------------------

  top_seros <- dat %>%
    dplyr::group_by(serotype) %>%
    dplyr::summarise(
      total_ipd = sum(ipd_mean),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(total_ipd)) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::pull(serotype)

  # 2. Order serotypes by invasiveness -----------------------------------

  serotype_order <- invasiveness_lookup %>%
    dplyr::filter(serotype %in% top_seros) %>%
    dplyr::arrange(dplyr::desc(v_child)) %>%
    dplyr::pull(serotype)

  serotype_order <- unique(c(
    serotype_order,
    setdiff(top_seros, serotype_order),
    "Other"
  ))

  # 3. Vaccine labels ---------------------------------------------------

  dat <- dat %>%
    dplyr::mutate(
      vaccine = dplyr::recode(
        vaccine,
        "best_greedy" = "Greedy",
        "best_exhaustive" = "Exhaustive"
      )
    )

  # 4. Stacked serotype-specific IPD -------------------------------------

  ipd_stacked <- dat %>%
    dplyr::mutate(
      serotype_plot = dplyr::if_else(
        serotype %in% top_seros,
        as.character(serotype),
        "Other"
      ),
      serotype_plot = factor(
        serotype_plot,
        levels = serotype_order
      )
    ) %>%
    dplyr::group_by(vaccine, year, serotype_plot) %>%
    dplyr::summarise(
      ipd_mean = sum(ipd_mean),
      .groups = "drop"
    )

  # 5. Total expected IPD ------------------------------------------------

  total_ipd <- dat %>%
    dplyr::group_by(vaccine, year) %>%
    dplyr::summarise(
      total = sum(ipd_mean),
      .groups = "drop"
    )

  # 6. Plot -------------------------------------------------------------

  ggplot2::ggplot() +
    ggplot2::geom_area(
      data = ipd_stacked,
      ggplot2::aes(
        x = year,
        y = ipd_mean,
        fill = serotype_plot
      ),
      position = "stack",
      alpha = 0.85
    ) +
    ggplot2::geom_line(
      data = total_ipd,
      ggplot2::aes(
        x = year,
        y = total
      ),
      colour = "black",
      linewidth = 1.1
    ) +
    ggplot2::geom_point(
      data = total_ipd,
      ggplot2::aes(
        x = year,
        y = total
      ),
      colour = "black",
      size = 2
    ) +
    ggplot2::facet_wrap(
      ~ vaccine,
      nrow = 1,
      scales = "fixed"
    ) +
    ggplot2::scale_fill_viridis_d(
      option = "mako",
      begin = 0.15,
      end = 0.9
    ) +
    ggplot2::scale_x_continuous(
      breaks = sort(unique(ipd_stacked$year))
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(
        mult = c(0, 0.05)
      )
    ) +
    ggplot2::labs(
      x = "Simulation year",
      y = "Expected IPD cases per 100,000 carriage events",
      fill = "Serotype",
      title = "Predicted invasive pneumococcal disease",
      subtitle = "Total IPD and serotype-specific contributions by vaccine"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "right",
      plot.title = ggplot2::element_text(face = "bold")
    )
}
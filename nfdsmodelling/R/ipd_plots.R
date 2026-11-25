# ============================================================
# Predicted South African IPD:
# total IPD with serotype contributions, and uncertainty dot plot
# ============================================================

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
  vaccine_colours,
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

  # Set vaccine ordering
  vaccine_order <- c(
    "PCV13",
    "Greedy",
    "Exhaustive"
  )

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
      vaccine = factor(
        vaccine,
        levels = vaccine_order
      )
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
    ggplot2::scale_colour_manual(
      values = vaccine_colours,
      breaks = vaccine_order
    ) +
    ggplot2::scale_x_continuous(
      labels = scales::label_number(suffix = "%"),
      limits = c(0, NA)
    ) +
    ggplot2::labs(
      x = "Proportion of total predicted IPD (%)",
      y = "Serotype",
      colour = "Vaccine"
    ) +
    theme_paper2(legend_position = "bottom") +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        nrow = 1,
        byrow = TRUE
      )
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

  # 3. Vaccine labels and panel ordering ---------------------------------

  dat <- dat %>%
    dplyr::mutate(
      vaccine = dplyr::recode(
        vaccine,
        "best_greedy" = "Greedy",
        "best_exhaustive" = "Exhaustive"
      ),
      vaccine = factor(
        vaccine,
        levels = c("PCV13", "Greedy", "Exhaustive")
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

  # Discrete blue palette
  blue_palette <- grDevices::colorRampPalette(
    c("#D6EAF8", "#85B9DF", "#3182BD", "#08519C", "#08306B")
  )(length(serotype_order))

  # Alternate light and dark shades
  n <- length(blue_palette)

  blue_palette <- blue_palette[
    c(rbind(seq_len(ceiling(n / 2)), rev(seq_len(n))[seq_len(ceiling(n / 2))]))
  ]

  blue_palette <- unique(blue_palette)

  serotype_colours <- setNames(
    blue_palette,
    serotype_order
  )

  serotype_colours["Other"] <- "grey85"

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
      ~vaccine,
      nrow = 1,
      scales = "fixed"
    ) +
    ggplot2::scale_fill_manual(
      values = serotype_colours,
      drop = FALSE
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(
        nrow = 2,
        byrow = TRUE
      )
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
    ) +
    theme_paper2() +
    ggplot2::guides(
      fill = ggplot2::guide_legend(
        nrow = 2,
        byrow = TRUE
      )
    )
}

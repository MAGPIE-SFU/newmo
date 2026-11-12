# ============================================================
#  vaccine_heatmap.R
#  Heatmap of licensed and invented vaccines: which serotypes each
#  contains, and the IPD each leaves.
# ============================================================

#' Vaccine heatmap
#'
#' One row per vaccine (licensed first, then the invented candidates), one
#' column per serotype, plus a column with the IPD each vaccine leaves.
#' Serotype columns: one block per licensed vaccine (its serotypes not in any
#' earlier vaccine), then everything else.
#'
#' @param known_serotypes Named list of licensed vaccines (serotype vectors),
#'   in display order.
#' @param known_ipd Named numeric vector: IPD score per licensed vaccine.
#' @param vaccine_data Invented vaccines, e.g. from [design_vaccines()]:
#'   columns ipd, vaccine (string with "-" or "," separators), method.
#' @param all_serotypes All serotypes in the model (e.g. samodel$st), so
#'   serotypes in no vaccine get a column too.
#' @param n_exhaustive Number of exhaustive vaccines to show, lowest IPD first.
#'   All non-exhaustive vaccines are always shown.
#' @param method_labels Named character vector: display name per method.
#'   Methods not listed keep their own name.
#' @param ipd_scale IPD is multiplied by this for display.
#' @param ipd_label Header of the IPD column. Keep in step with `ipd_scale`.
#' @return A list with `plot` (the combined plot), `width` and `height`
#'   (suggested size in inches).
#' @importFrom dplyr %>% tibble filter slice_min bind_rows mutate coalesce arrange
#'   group_by n row_number ungroup select distinct transmute if_else case_when n_distinct
#' @importFrom tidyr unnest_longer complete nesting
#' @importFrom stringr str_split str_extract
#' @importFrom ggplot2 ggplot aes geom_tile geom_vline geom_hline scale_fill_manual
#'   scale_x_discrete scale_y_discrete facet_grid labs theme_minimal theme
#'   element_blank unit element_text margin geom_text scale_fill_gradient
#'   scale_colour_identity
#' @importFrom patchwork plot_layout
#' @export
plot_vaccine_heatmap <- function(known_serotypes, known_ipd, vaccine_data, all_serotypes,
                                 n_exhaustive = 6,
                                 method_labels = character(0),
                                 ipd_scale = 1e5,
                                 ipd_label = "IPD per\n100,000") {

  # Stops here if the scores are anything other than one number each
  stopifnot(is.numeric(known_ipd), length(known_ipd) == length(known_serotypes))

  all_serotypes <- unique(as.character(all_serotypes))
  all_serotypes <- all_serotypes[!is.na(all_serotypes)]

  # Build vaccine table ----------------------------------------------------------

  known_long <- tibble(
    vaccine_label = names(known_serotypes),
    group = "Licensed",
    ipd = unname(known_ipd[names(known_serotypes)]),
    serotypes = unname(known_serotypes)
  )

  exhaustive_best <- vaccine_data %>%
    filter(method == "exhaustive") %>%
    slice_min(ipd, n = n_exhaustive, with_ties = FALSE)

  candidate_long <- vaccine_data %>%
    filter(method != "exhaustive") %>%
    bind_rows(exhaustive_best) %>%
    mutate(
      vaccine = as.character(vaccine),
      # Vaccine strings use both "-" and "," as separators
      serotypes = str_split(vaccine, "[-,]"),
      valency = lengths(serotypes),
      group = "Candidate",
      method_name = coalesce(unname(method_labels[method]), method)
    ) %>%
    arrange(method != "greedy", method != "exhaustive", method, valency, ipd) %>%
    group_by(method, valency) %>%
    # Number the rows only where a method repeats a valency, best IPD first
    mutate(
      vaccine_label = if (n() > 1) {
        sprintf("%s %d-valent #%d", method_name, valency, row_number())
      } else {
        sprintf("%s %d-valent", method_name, valency)
      }
    ) %>%
    ungroup() %>%
    select(vaccine_label, group, method, ipd, serotypes)

  vaccine_table <- bind_rows(known_long, candidate_long) %>%
    mutate(vaccine_label = factor(vaccine_label, levels = rev(vaccine_label)))

  vaccine_long <- vaccine_table %>%
    select(-ipd) %>%
    unnest_longer(serotypes, values_to = "serotype") %>%
    mutate(serotype = as.character(serotype))

  # Column order -----------------------------------------------------------------
  sort_serotypes <- function(x) {
    x <- unique(as.character(x))
    x[order(as.numeric(str_extract(x, "^\\d+")), x, na.last = TRUE)]
  }

  # One block per licensed vaccine, holding its serotypes that are not in any
  # earlier vaccine, then everything else
  serotype_blocks <- list()
  seen_serotypes <- character(0)
  for (v in names(known_serotypes)) {
    serotype_blocks[[v]] <- sort_serotypes(setdiff(known_serotypes[[v]], seen_serotypes))
    seen_serotypes <- union(seen_serotypes, known_serotypes[[v]])
  }

  other_serotypes <- setdiff(
    union(all_serotypes, vaccine_long$serotype),
    seen_serotypes
  )

  serotype_levels <- c(
    unlist(serotype_blocks, use.names = FALSE),
    sort_serotypes(other_serotypes)
  )

  # Separator lines --------------------------------------------------------------

  line_colour <- "#1F2933"
  line_width <- 0.6

  # Vertical: after each licensed vaccine's block
  serotype_lines <- cumsum(lengths(serotype_blocks)) + 0.5

  # Horizontal: between candidate methods. Rows are counted from the bottom.
  candidate_methods <- vaccine_table$method[vaccine_table$group == "Candidate"]
  n_candidates <- length(candidate_methods)
  method_changes <- which(candidate_methods[-1] != candidate_methods[-n_candidates])

  method_lines <- data.frame(
    group = factor("Candidate", levels = c("Licensed", "Candidate")),
    y = n_candidates - method_changes + 0.5
  )

  # Serotype heatmap -------------------------------------------------------------

  plot_data <- vaccine_long %>%
    distinct(vaccine_label, group, serotype) %>%
    mutate(included = TRUE) %>%
    complete(
      nesting(vaccine_label, group),
      serotype = serotype_levels,
      fill = list(included = FALSE)
    ) %>%
    mutate(
      serotype = factor(serotype, levels = serotype_levels),
      group = factor(group, levels = c("Licensed", "Candidate"))
    )

  serotype_plot <- ggplot(plot_data, aes(serotype, vaccine_label, fill = included)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_vline(xintercept = serotype_lines, colour = line_colour, linewidth = line_width) +
    geom_hline(
      data = method_lines, aes(yintercept = y), inherit.aes = FALSE,
      colour = line_colour, linewidth = line_width
    ) +
    scale_fill_manual(
      values = c(`TRUE` = "#1F6F8B", `FALSE` = "#ECEFF1"),
      labels = c(`TRUE` = "Included", `FALSE` = "Not included"),
      breaks = c(TRUE, FALSE),
      name = NULL
    ) +
    scale_x_discrete(position = "top", expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    facet_grid(group ~ ., scales = "free_y", space = "free_y") +
    labs(x = "Serotype", y = NULL) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      panel.spacing.y = unit(0.6, "lines"),
      axis.text.x.top = element_text(angle = 90, vjust = 0.5, hjust = 0, margin = margin(b = 2)),
      axis.text.y = element_blank(),
      strip.text.y = element_blank(),
      legend.position = "bottom"
    )

  # IPD column -------------------------------------------------------------------
  ipd_data <- vaccine_table %>%
    transmute(
      vaccine_label,
      group = factor(group, levels = c("Licensed", "Candidate")),
      column = ipd_label,
      ipd_scaled = ipd * ipd_scale
    )

  ipd_midpoint <- mean(range(ipd_data$ipd_scaled, na.rm = TRUE))

  ipd_data <- ipd_data %>%
    mutate(
      value_label = if_else(is.na(ipd_scaled), "n/a", sprintf("%.1f", ipd_scaled)),
      text_colour = case_when(
        is.na(ipd_scaled) ~ "#6B7280",
        ipd_scaled > ipd_midpoint ~ "white",
        TRUE ~ "#1F2933"
      )
    )

  ipd_plot <- ggplot(ipd_data, aes(column, vaccine_label, fill = ipd_scaled)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_hline(
      data = method_lines, aes(yintercept = y), inherit.aes = FALSE,
      colour = line_colour, linewidth = line_width
    ) +
    geom_text(aes(label = value_label, colour = text_colour), size = 3.2) +
    # A different hue from the heatmap, as this shows a value rather than yes/no
    scale_fill_gradient(low = "#FBE9DD", high = "#B4541A", na.value = "#ECEFF1") +
    scale_colour_identity() +
    scale_x_discrete(position = "top", expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    facet_grid(group ~ ., scales = "free_y", space = "free_y", switch = "y") +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      panel.spacing.y = unit(0.6, "lines"),
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0, face = "bold"),
      legend.position = "none"
    )

  # Combine ----------------------------------------------------------------------

  heatmap_plot <- ipd_plot + serotype_plot +
    plot_layout(widths = c(3, length(serotype_levels)))

  list(
    plot   = heatmap_plot,
    width  = max(8, 3.2 + 0.22 * length(serotype_levels)),
    height = 2 + 0.3 * n_distinct(plot_data$vaccine_label)
  )
}

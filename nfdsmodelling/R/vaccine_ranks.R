#' Rank plot of invented vaccines, with licensed vaccines for comparison
#'
#' Exhaustive vaccines ranked by the IPD they leave (rank 1 = lowest IPD):
#' the best `n_top`, then `n_sample` ranks evenly spread over the rest.
#' Vaccines from `other_methods` are placed at the rank they would have among
#' the exhaustive vaccines. Licensed vaccines are dashed horizontal lines.
#'
#' @param vaccine_data Invented vaccines, e.g. from [design_vaccines()].
#' @param known_ipd Named numeric vector: IPD score per licensed vaccine.
#' @param n_top Number of best exhaustive vaccines shown individually.
#' @param n_sample Number of evenly spaced ranks shown from the remaining ones.
#' @param other_methods Methods placed among the exhaustive ranks.
#' @param method_labels Named character vector: display name per method.
#' @param ipd_scale IPD is multiplied by this for display.
#' @param ipd_label y-axis label. Keep in step with `ipd_scale`.
#' @return A ggplot.
#' @importFrom dplyr %>% filter arrange mutate row_number if_else coalesce
#'   bind_rows tibble
#' @importFrom ggplot2 ggplot aes geom_hline geom_text geom_point scale_x_log10
#'   scale_colour_brewer labs theme_bw
#' @export
plot_vaccine_ranks <- function(vaccine_data, known_ipd,
                               n_top = 100,
                               n_sample = 1000,
                               other_methods = "greedy",
                               method_labels = character(0),
                               ipd_scale = 1e5,
                               ipd_label = "IPD per 100,000") {

  exhaustive <- vaccine_data %>%
    filter(method == "exhaustive") %>%
    arrange(ipd) %>%
    mutate(rank = row_number())
  n_exh <- nrow(exhaustive)

  # the best n_top, then n_sample ranks evenly spread over the rest
  sample_ranks <- if (n_exh > n_top) {
    unique(round(seq(n_top + 1, n_exh, length.out = n_sample)))
  } else integer(0)

  shown_exh <- exhaustive %>%
    filter(rank <= n_top | rank %in% sample_ranks) %>%
    mutate(part = if_else(rank <= n_top,
                          paste("Exhaustive: best", n_top),
                          "Exhaustive: rest (evenly spaced)"))

  # other methods at the rank they would have among the exhaustive vaccines:
  # 1 + the number of exhaustive vaccines with a lower IPD
  others <- vaccine_data %>%
    filter(method %in% other_methods) %>%
    mutate(
      rank = findInterval(ipd, exhaustive$ipd) + 1,
      part = coalesce(unname(method_labels[method]), method),
      valency = lengths(str_split(as.character(vaccine), "[-,]"))
    )

  plot_df  <- bind_rows(shown_exh, others) %>% mutate(ipd_scaled = ipd * ipd_scale)

  min_gap <- 0.12   # minimum horizontal distance between labels, in decades (log10 units)

  licensed <- tibble(vaccine    = names(known_ipd),
                     ipd_scaled = unname(known_ipd) * ipd_scale,
                     rank       = findInterval(unname(known_ipd), exhaustive$ipd) + 1) %>%
    arrange(rank)

  # push each label right until it is at least min_gap away from the one before
  lx <- log10(licensed$rank)
  for (i in seq_along(lx)[-1]) lx[i] <- max(lx[i], lx[i - 1] + min_gap)
  licensed$label_x <- 10^lx

  ggplot(plot_df, aes(rank, ipd_scaled, colour = part)) +
    # geom_hline(data = licensed, aes(yintercept = ipd_scaled),
    #            linetype = "dashed", colour = "grey30") +
    # geom_vline(data = licensed, aes(xintercept = rank),
    #            linetype = "dashed", colour = "grey30") +
    # geom_text(data = licensed, aes(x = Inf, y = ipd_scaled, label = vaccine),
    #           hjust = 1.1, vjust = -0.4, size = 3.2, colour = "grey30",
    #           inherit.aes = FALSE) +
        # horizontal: from the left up to the vaccine's rank
    geom_segment(data = licensed,
                 aes(x = 1, xend = rank, y = ipd_scaled, yend = ipd_scaled),
                 linetype = "dashed", colour = "grey30", inherit.aes = FALSE) +
    # vertical: from the bottom up to the vaccine's IPD
    geom_segment(data = licensed,
                 aes(x = rank, xend = rank, y = -Inf, yend = ipd_scaled),
                 linetype = "dashed", colour = "grey30", inherit.aes = FALSE) +
    # vaccine names below the x-axis, at their rank
    geom_text(data = licensed,
              aes(x = label_x, y = -Inf, label = vaccine),
              angle = 90, hjust = 1.6, vjust = 0, size = 3.2, colour = "grey30",
              inherit.aes = FALSE) +

    geom_text(data = ~ filter(.x, method != "exhaustive"),
              aes(label = valency), vjust = -1, size = 3, show.legend = FALSE) +

    geom_point(data = ~ filter(.x, method == "exhaustive"), size = 1.2, shape = 21, alpha = 0.4) +
    geom_point(data = ~ filter(.x, method != "exhaustive"), size = 2.6, shape = 17) +
    scale_x_log10(labels = scales::label_comma()) +
    scale_colour_brewer(palette = "Dark2", name = NULL) +
    labs(x = "Rank among exhaustive vaccines (log scale)", y = ipd_label,
         title = "Invented vaccines ranked by predicted IPD",
         subtitle = paste0("Exhaustive: best ", n_top, " and ", length(sample_ranks),
                           " evenly spaced of the remaining ",
                           format(n_exh - n_top, big.mark = ","),
                           "; other methods at their equivalent rank; dashed: licensed")) +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme(plot.margin  = margin(5.5, 5.5, 40, 5.5),
            axis.title.x = element_text(margin = margin(t = 40)))
}

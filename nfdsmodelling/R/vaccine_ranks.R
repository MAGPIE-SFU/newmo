#' Rank plot of invented vaccines, with licensed vaccines for comparison
#'
#' Invented vaccines sorted by the IPD they leave (rank 1 = lowest IPD),
#' coloured by method. Licensed vaccines are shown as dashed horizontal lines.
#'
#' @param vaccine_data Invented vaccines, e.g. from [design_vaccines()].
#' @param known_ipd Named numeric vector: IPD score per licensed vaccine.
#' @param n_top Number of best-ranked vaccines to show (Inf = all).
#' @param valencies Only rank vaccines with these valencies (NULL = all).
#' @param method_labels Named character vector: display name per method.
#' @param ipd_scale IPD is multiplied by this for display.
#' @param ipd_label y-axis label. Keep in step with `ipd_scale`.
#' @return A ggplot.
#' @importFrom dplyr %>% mutate filter arrange row_number coalesce tibble
#' @importFrom stringr str_split
#' @importFrom ggplot2 ggplot aes geom_hline geom_text geom_point
#'   scale_colour_brewer labs theme_bw
#' @export
plot_vaccine_ranks <- function(vaccine_data, known_ipd,
                               n_top = 100,
                               valencies = NULL,
                               method_labels = character(0),
                               ipd_scale = 1e5,
                               ipd_label = "IPD per 100,000") {

  ranked <- vaccine_data %>%
    mutate(valency = lengths(str_split(as.character(vaccine), "[-,]")))
  if (!is.null(valencies)) ranked <- ranked %>% filter(valency %in% valencies)
  n_all <- nrow(ranked)

  ranked <- ranked %>%
    arrange(ipd) %>%
    mutate(rank        = row_number(),
           method_name = coalesce(unname(method_labels[method]), method),
           ipd_scaled  = ipd * ipd_scale) %>%
    filter(rank <= n_top)

  licensed <- tibble(vaccine    = names(known_ipd),
                     ipd_scaled = unname(known_ipd) * ipd_scale)

  ggplot(ranked, aes(rank, ipd_scaled)) +
    geom_hline(data = licensed, aes(yintercept = ipd_scaled),
               linetype = "dashed", colour = "grey30") +
    geom_text(data = licensed, aes(x = Inf, y = ipd_scaled, label = vaccine),
              hjust = 1.1, vjust = -0.4, size = 3.2, colour = "grey30",
              inherit.aes = FALSE) +
    geom_point(aes(colour = method_name), size = 1.8) +
    scale_colour_brewer(palette = "Dark2", name = "Method") +
    labs(x = "Rank (lowest IPD first)", y = ipd_label,
         title = "Invented vaccines ranked by predicted IPD",
         subtitle = paste0("Best ", nrow(ranked), " of ", n_all,
                           " vaccines; dashed lines: licensed vaccines")) +
    theme_bw()
}

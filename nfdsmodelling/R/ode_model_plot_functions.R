# Data = grey, simulations = coloured (one colour per simulation).
# top_group_plot / top_sero_plot drop VT serotypes (modelpars$vaccine) when drop_vt = TRUE (frequencies are NOT renormalised).
# vt_nvt_plot shows the VT and NVT totals in separate panels.

# helper: sum frequencies within a grouping variable (e.g. serotype)
#' @importFrom dplyr distinct mutate if_else left_join group_by across all_of summarise
agg <- function(df, by, xt, vaccine) {
  ginfo <- xt |>
    distinct(group, serotype) |>
    mutate(vt = if_else(serotype %in% vaccine, "VT", "NVT"))

  df |>
    left_join(ginfo, by = "group") |>
    group_by(across(all_of(c(by, intersect(c("year", "sim"), names(df)))))) |>
    summarise(frac = sum(frac), .groups = "drop")
}

# colour scale for the simulations, from a ready-made palette
# palette: a RColorBrewer name ("Set1", "Dark2", "Paired", "Set2", ...)
#          or a viridis option ("viridis", "magma", "plasma", "inferno", "cividis", "turbo")
#' @importFrom ggplot2 scale_colour_viridis_d scale_colour_brewer
sim_colour_scale <- function(palette = "Dark2") {
  viridis_opts <- c("viridis", "magma", "plasma", "inferno", "cividis", "turbo")
  if (palette %in% viridis_opts) {
    scale_colour_viridis_d(option = palette, end = 0.9)
  } else {
    scale_colour_brewer(palette = palette)
  }
}

# plot instructions, this goes in the ggplot call
#' @importFrom ggplot2 geom_line aes geom_point geom_segment geom_vline labs theme unit
compare_layers <- function(d, s, b, vax_start, prevax_years, palette = "Dark2") {
  list(
    # data: neutral grey
    geom_line(data = d, aes(year, frac), colour = "grey40", alpha = 0.3),
    geom_point(data = d, aes(year, frac), colour = "grey20", size = 1.5),

    # simulations: coloured by sim
    geom_line(data = s, aes(year, frac, colour = sim), linewidth = 1),

    # x0 (pre-vaccine mean): neutral dashed segment
    geom_segment(data = b, aes(y = frac, yend = frac),
                 x = min(prevax_years), xend = vax_start - 0.5,
                 colour = "black", linetype = "dashed"),
    geom_vline(xintercept = vax_start - 0.5, linetype = "dotted"),

    sim_colour_scale(palette),
    labs(x = NULL, y = "Frequency", colour = "Simulation",
         caption = "Points: data. Coloured lines: simulations. Dashed: x0 (pre-vaccine mean)."),
    theme(legend.key.width = unit(1.5, "cm"))
  )
}

#' @importFrom tibble rownames_to_column tibble
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr starts_with mutate select bind_rows
make_sim_long <- function(sim, samodel, vax_start) {
  # model: rows = groups, columns = year.1 ... year.k (make sure each year sums to 1)
  sim_long <- sim |>
    sweep(2, colSums(sim), "/") |>
    as.data.frame() |>
    rownames_to_column("group") |>
    pivot_longer(starts_with("year."), names_to = "k",
                 names_prefix = "year.", values_to = "frac") |>
    mutate(year = vax_start - 1 + as.integer(k)) |>
    select(group, year, frac)

  # x0 as frequencies (the pre-vaccine mean the model starts from)
  x0_df <- tibble(group = rownames(samodel$G),
                  frac  = as.numeric(samodel$x0) / sum(samodel$x0))
  bind_rows(x0_df |> mutate(year = vax_start - 0.5), sim_long)
}

# stack a list of simulations -- so we can compare multiple simulations on the same plot
# sims should be a list, like list(baseline = mysim1, lowerv = mysim2) and so on
#' @importFrom dplyr bind_rows mutate
stack_sims <- function(sims, samodel, vax_start) {
  if (is.data.frame(sims)) sims <- list(model = sims)   # allow a single simulation too
  bind_rows(lapply(sims, make_sim_long, samodel = samodel, vax_start = vax_start), .id = "sim") |>
    mutate(sim = factor(sim, levels = unique(sim)))
}


#' Plot VT and NVT totals, model and data together
#'
#' Data = the synthetic carriage sample x(t). VT and NVT go in separate panels
#' (colour is used for the simulations).
#'
#' @param vax_start First vaccine year.
#' @param prevax_years Pre-vaccine years (for the x0 segment).
#' @param samodel The model structure.
#' @param xt Grouped x(t).
#' @param sims A simulation, or a named list of simulations.
#' @param palette RColorBrewer name or viridis option for the simulations.
#' @param vaccine Vaccine serotypes (defines VT vs NVT).
#' @return A ggplot.
#' @importFrom dplyr select
#' @importFrom tibble tibble
#' @importFrom ggplot2 ggplot facet_wrap ggtitle
#' @export
vt_nvt_plot <- function(vax_start = 2010, prevax_years = 2000:2009,
                        samodel, xt, sims, palette = "Dark2", vaccine) {
  sim_long  <- stack_sims(sims, samodel, vax_start)
  data_long <- xt |> select(group, year, frac, frac_sd)
  x0_df     <- tibble(group = rownames(samodel$G),
                      frac  = as.numeric(samodel$x0) / sum(samodel$x0))

  ggplot() +
    compare_layers(agg(data_long, "vt", xt, vaccine), agg(sim_long, "vt", xt, vaccine), agg(x0_df, "vt", xt, vaccine),
                   vax_start = vax_start, prevax_years = prevax_years,
                   palette = palette) +
    facet_wrap(~ vt, scales = "free_y") +
    ggtitle("Vaccine-type vs non-vaccine-type frequency")
}


#' Plot the most common ST-SC groups, model and data together
#'
#' Strip colour: VT = red, NVT = blue.
#'
#' @inheritParams vt_nvt_plot
#' @param n_groups Number of top ST-SC groups to show.
#' @param drop_vt Drop vaccine-type groups (frequencies are NOT renormalised).
#' @return A ggplot.
#' @importFrom dplyr select distinct filter pull group_by summarise slice_max mutate if_else arrange
#' @importFrom tibble tibble
#' @importFrom ggplot2 ggplot geom_errorbar aes ggtitle
#' @export
top_group_plot <- function(vax_start = 2010, prevax_years = 2000:2009,
                           n_groups = 12, samodel, xt, sims, drop_vt = TRUE, palette = "Dark2", vaccine) {
  sim_long  <- stack_sims(sims, samodel, vax_start)
  data_long <- xt |> select(group, year, frac, frac_sd)
  x0_df     <- tibble(group = rownames(samodel$G),
                      frac  = as.numeric(samodel$x0) / sum(samodel$x0))

  # optionally drop vaccine-type groups
  if (drop_vt) {
    keep <- xt |> distinct(group, serotype) |>
      filter(!serotype %in% vaccine) |> pull(group)
    data_long <- data_long |> filter(group %in% keep)
  }

  # top groups (in the data)
  top_groups <- data_long |>
    group_by(group) |>
    summarise(m = mean(frac)) |>
    slice_max(m, n = n_groups) |>
    pull(group)

  # facet order + strip text colour: VT = red, NVT = blue
  grp_levels <- sort(top_groups)
  strip_cols <- xt |>
    distinct(group, serotype) |>
    filter(group %in% grp_levels) |>
    mutate(col = if_else(serotype %in% vaccine, "red", "blue")) |>
    arrange(match(group, grp_levels)) |>
    pull(col)

  d_g <- data_long |> filter(group %in% top_groups) |> mutate(group = factor(group, levels = grp_levels))
  s_g <- sim_long  |> filter(group %in% top_groups) |> mutate(group = factor(group, levels = grp_levels))
  b_g <- x0_df     |> filter(group %in% top_groups) |> mutate(group = factor(group, levels = grp_levels))

  ggplot() +
    geom_errorbar(data = d_g, aes(year, ymin = pmax(frac - frac_sd, 0),
                                  ymax = frac + frac_sd),
                  colour = "grey40", width = 0, alpha = 0.4) +
    compare_layers(d_g, s_g, b_g, vax_start = vax_start, prevax_years = prevax_years,
                   palette = palette) +
    ggh4x::facet_wrap2(~ group, scales = "free_y",
                       strip = ggh4x::strip_themed(
                         text_x = ggh4x::elem_list_text(colour = strip_cols))) +
    ggtitle(paste("The", n_groups, "most common ST-SC groups"))
}

#' Plot the most common serotypes, model and data together
#'
#' Strip colour: VT = red, NVT = blue.
#'
#' @inheritParams vt_nvt_plot
#' @param n_serotypes Number of top serotypes to show.
#' @param drop_vt Drop vaccine-type serotypes (frequencies are NOT renormalised).
#' @return A ggplot.
#' @importFrom dplyr select filter group_by summarise slice_max pull mutate
#' @importFrom tibble tibble
#' @importFrom ggplot2 ggplot ggtitle
#' @export
top_sero_plot <- function(vax_start = 2010, prevax_years = 2000:2009,
                          n_serotypes = 12, samodel, xt, sims, drop_vt = TRUE, palette = "Dark2", vaccine) {
  sim_long  <- stack_sims(sims, samodel, vax_start)
  data_long <- xt |> select(group, year, frac, frac_sd)
  x0_df     <- tibble(group = rownames(samodel$G),
                      frac  = as.numeric(samodel$x0) / sum(samodel$x0))

  d_st <- agg(data_long, "serotype", xt, vaccine)
  if (drop_vt) d_st <- d_st |> filter(!serotype %in% vaccine)
  top_st <- d_st |>
    group_by(serotype) |>
    summarise(m = mean(frac)) |>
    slice_max(m, n = n_serotypes) |>
    pull(serotype)

  # facet order + strip text colour: VT = red, NVT = blue
  sero_levels <- sort(top_st)
  strip_cols  <- ifelse(sero_levels %in% vaccine, "red", "blue")

  fac <- function(df) df |> filter(serotype %in% top_st) |>
    mutate(serotype = factor(serotype, levels = sero_levels))

  ggplot() +
    compare_layers(fac(d_st),
                   fac(agg(sim_long, "serotype", xt, vaccine)),
                   fac(agg(x0_df,    "serotype", xt, vaccine)),
                   vax_start = vax_start, prevax_years = prevax_years,
                   palette = palette) +
    ggh4x::facet_wrap2(~ serotype, scales = "free_y",
                       strip = ggh4x::strip_themed(
                         text_x = ggh4x::elem_list_text(colour = strip_cols))) +
    ggtitle(paste("The", n_serotypes, "most common serotypes"))
}

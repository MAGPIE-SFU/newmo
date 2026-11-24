# ============================================================
#  Assess sensitivity of the NFDS model to its input parameters by perturbing
#  inputs near the optimized parameter settings.
# ============================================================
# Display names for the loss functions in the plots
LOSS_LABELS <- c(
  SSE = "SSE",
  weighted_composition_distance = "Weighted composition",
  weighted_SSE = "Weighted SSE"
)
loss_labeller <- function(x) {
  ifelse(x %in% names(LOSS_LABELS), LOSS_LABELS[x], x)
}

# Common sensitivity figure styling. Colours can be overridden via config
# when calling add_sensitivity_ribbon(); loss colours are distinct blues.
SENS_LOSS_COLOURS <- c(
  SSE = "#084594",
  weighted_composition_distance = "#4292C6",
  weighted_SSE = "#9ECAE1"
)

sensitivity_group_labeller <- function(x, pcv13 = c(
    "1", "3", "4", "5", "6A", "6B", "7F", "9V", "14",
    "18C", "19A", "19F", "23F"),
    vt_colour = "#A94C22", nvt_colour = "#084594") {
  serotype <- sub("_.*$", "", x)
  colour <- ifelse(serotype %in% pcv13, vt_colour, nvt_colour)
  paste0("<span style='color:", colour, "'><b>", x, "</b></span>")
}

sensitivity_paper_theme <- function(coloured_strips = TRUE) {
  ptheme <- theme_paper2(legend_position = "bottom")
  if (coloured_strips) {
    ptheme + ggplot2::theme(
      strip.text = ggtext::element_markdown(face = "bold")
    )
  } else {
    ptheme + ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold")
    )
  }
}

#' Largest groups at the start of the simulation
#'
#' @param samodel The model structure.
#' @param n_groups Number of groups.
#' @return Data frame with group, IC and frac, largest first.
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr mutate arrange desc slice_head
#' @export
sensitivity_top_groups <- function(samodel, n_groups) {
  data.frame(IC = samodel$x0) |>
    rownames_to_column(var = "group") |>
    mutate(frac = IC / sum(IC)) |>
    arrange(desc(frac)) |>
    slice_head(n = n_groups)
}

#' Sensitivity of the model to its parameters around the fitted optima
#'
#' For each optimal fit (one row of `opt_pars` per loss function), draws
#' `n_sims` parameter sets from normal distributions around the fitted values
#' (log scale for v1 and sigma, truncated to the bounds lo and hi) and simulates each.
#'
#' @param opt_pars Optimal fits, one row per loss, with a `loss` column and PAR_NAMES.
#' @param samodel,tswitch,nyears,use_switch As in [run_model()].
#' @param sd_pars Data frame with columns par (logv1, logsigma, portion), SD,
#'   lo and hi (bounds on the sampling scale).
#' @param n_sims Number of parameter draws per loss function.
#' @param n_cores Cores for the simulations (> 1 uses mclapply, Linux/macOS only).
#' @return A list with `opt_traj` (trajectories of the optimal fits),
#'   `par_df` (the draws, sampling scale), `opt_settings` (the optimal values,
#'   sampling scale), `sim_out` (trajectories of all draws) and `failed`
#'   (draws whose simulation failed).
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr bind_rows bind_cols mutate left_join select between everything
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom purrr list_rbind keep discard
#' @importFrom stats runif pnorm qnorm
#' @export
run_sensitivity <- function(
  opt_pars,
  samodel,
  tswitch,
  nyears,
  sd_pars,
  use_switch = TRUE,
  n_sims = 100,
  n_cores = 1
) {
# pull out initial conditions
  IC <- data.frame(value = samodel$x0) |>
    rownames_to_column(var = "group") |>
    bind_cols(year = 0)
# Optimal trajectories
# ==============================================================================
# process results for each of the optimal parameter settings
  opt_traj <- lapply(seq_len(nrow(opt_pars)), function(i) {
    run_model(
      as.list(unlist(opt_pars[i, PAR_NAMES])),
      use_switch,
      samodel,
      nyears,
      tswitch
    ) |>
      rownames_to_column(var = "group") |>
      pivot_longer(
        -group,
        names_to = "year",
        names_prefix = "year.",
        names_transform = as.numeric
      ) |>
      bind_rows(IC) |>
      bind_cols(loss = opt_pars$loss[i])
  }) |>
    list_rbind()
# Perturbations
# ==============================================================================
# For each parameter simulate n_sims draws from a normal distribution with mean
#  equal to the optimized values and standard deviation loosely based on the
#  loss function profiles
  par_lookup <- lapply(seq_len(nrow(opt_pars)), function(i) {
    as.data.frame(opt_pars[i, ]) |>
      mutate(
        logv1 = log(v1),
        logsigma = log(sigma),
        portion = portion,
        .keep = "none"
      ) |>
      pivot_longer(everything(), names_to = "par") |>
      left_join(sd_pars, by = "par")
  })
# draw from a normal truncated to [lo, hi] via the inverse CDF
  sim_pars <- Vectorize(
    function(lookup, par_name) {
      row <- lookup[lookup$par == par_name, ]
      p <- runif(
        1,
        pnorm(row$lo, row$value, row$SD),
        pnorm(row$hi, row$value, row$SD)
      )
      qnorm(p, row$value, row$SD)
    },
    vectorize.args = "par_name"
  )
  par_df <- lapply(seq_len(nrow(opt_pars)), function(j) {
    lapply(1:n_sims, function(i) {
      data.frame(par = sd_pars$par) |>
        mutate(val = sim_pars(par_lookup[[j]], par)) |>
        bind_cols(rep = i)
    }) |>
      list_rbind() |>
      bind_cols(loss = opt_pars$loss[j])
  }) |>
    list_rbind()
# pull out parameters used in the optimal fits
  opt_settings <- lapply(seq_along(par_lookup), function(i) {
    par_lookup[[i]] |>
      select(par, value) |>
      bind_cols(loss = opt_pars$loss[i])
  }) |>
    list_rbind()
  sim_par_sets <- par_df |>
    pivot_wider(names_from = "par", values_from = "val") |>
    mutate(sigma = exp(logsigma), v1 = exp(logv1), .keep = "unused")
  stopifnot(all(between(sim_par_sets$portion, 0, 1)))
# run the model for each perturbed parameter set, keeping the rep/loss labels
#  alongside the trajectories so the results can be summarised and plotted
  run_one <- function(k) {
    pars_k <- sim_par_sets[k, ]
# a bad parameter draw shouldn't kill the whole batch; drop it and record why
    tryCatch(
      run_model(
        as.list(select(pars_k, -rep, -loss)),
        use_switch,
        samodel,
        nyears,
        tswitch
      ) |>
        rownames_to_column(var = "group") |>
        pivot_longer(
          -group,
          names_to = "year",
          names_prefix = "year.",
          names_transform = as.numeric
        ) |>
        bind_rows(IC) |>
        bind_cols(select(pars_k, rep, loss)),
      error = function(e) {
        data.frame(
          rep = pars_k$rep,
          loss = pars_k$loss,
          error = conditionMessage(e)
        )
      }
    )
  }
# ODE solve times vary a lot across parameter draws, so hand out one job at a
#  time (no prescheduling) rather than one fixed block per core
  sim_res <- if (n_cores > 1) {
    parallel::mclapply(
      seq_len(nrow(sim_par_sets)),
      run_one,
      mc.cores = n_cores,
      mc.preschedule = FALSE
    )
  } else {
    lapply(seq_len(nrow(sim_par_sets)), run_one)
  }
  failed <- keep(sim_res, \(x) "error" %in% names(x)) |> list_rbind()
  if (NROW(failed) > 0) {
    warning(nrow(failed), " parameter sets failed; see `failed`")
  } else {
    failed <- data.frame(
      rep = integer(),
      loss = character(),
      error = character()
    )
  }
  sim_out <- discard(sim_res, \(x) "error" %in% names(x)) |>
    list_rbind()
  list(
    opt_traj = opt_traj,
    par_df = par_df,
    opt_settings = opt_settings,
    sim_out = sim_out,
    failed = failed
  )
}

# fitted trajectories of the top groups, on the calendar scale
#' @importFrom dplyr filter mutate
sensitivity_opt_fit <- function(sens, top_groups, start_year) {
  sens$opt_traj |>
    filter(group %in% top_groups$group) |>
    mutate(year = start_year + year)
}

# simulated trajectory summaries
#' @importFrom dplyr filter group_by summarise mutate
#' @importFrom stats quantile
sensitivity_summary <- function(sens, top_groups, start_year, qtile) {
  sens$sim_out |>
    filter(group %in% top_groups$group) |>
    group_by(group, year, loss) |>
    summarise(
      L = quantile(value, probs = (1 - qtile) / 2),
      med = quantile(value, probs = 0.5),
      U = quantile(value, probs = 1 - (1 - qtile) / 2),
      .groups = "drop"
    ) |>
    mutate(year = start_year + year)
}

#' Sensitivity plots
#'
#' @param sens Result of [run_sensitivity()].
#' @param top_groups Result of [sensitivity_top_groups()].
#' @param start_year Year prior to start of model evaluation (i.e., fitted
#'   trajectories start at start_year + 1).
#' @param qtile Summary quantile for the ribbons.
#' @param losses Loss functions to show (NULL = all).
#' @param xt Grouped x(t) (the data).
#' @param best_loss Loss function of the chosen optimal fit, or NULL for all.
#' @return A ggplot.
#' @name sensitivity_plots
NULL
#' @describeIn sensitivity_plots Trajectories of the optimal fits.
#' @importFrom dplyr filter
#' @importFrom ggplot2 ggplot aes geom_point geom_line facet_wrap theme
#' @export
plot_sensitivity_optimal <- function(sens, top_groups) {
  sens$opt_traj |>
    filter(group %in% top_groups$group) |>
    ggplot(aes(x = year, y = value, colour = loss)) +
    geom_point() +
    geom_line() +
    facet_wrap(~group, ncol = 4, scales = "free_y",
               labeller = ggplot2::labeller(group = sensitivity_group_labeller)) +
    ggplot2::scale_colour_manual(values = SENS_LOSS_COLOURS, labels = loss_labeller) +
    sensitivity_paper_theme()
}

#' @describeIn sensitivity_plots Range of the simulated parameters (sampling scale).
#' @importFrom ggplot2 ggplot aes geom_histogram facet_grid geom_vline
#' @export
plot_sensitivity_draws <- function(sens) {
  sens$par_df |>
    ggplot(aes(x = val)) +
    geom_histogram(aes(fill = loss), colour = "black") +
    facet_grid(loss ~ par, scales = "free") +
    geom_vline(
      data = sens$opt_settings,
      aes(xintercept = value),
      linewidth = 1,
      linetype = "dashed"
    ) +
    ggplot2::scale_fill_manual(values = SENS_LOSS_COLOURS, labels = loss_labeller) +
    sensitivity_paper_theme(coloured_strips = FALSE)
}

#' @describeIn sensitivity_plots Ribbons with the optimal trajectories (dashed).
#' @importFrom ggplot2 ggplot aes geom_ribbon geom_line facet_wrap labs
#'   theme
#' @export
plot_sensitivity_ribbons <- function(sens, top_groups, start_year, qtile) {
  opt_fit <- sensitivity_opt_fit(sens, top_groups, start_year)
  sensitivity_summary(sens, top_groups, start_year, qtile) |>
    ggplot(aes(x = year, colour = loss)) +
    geom_ribbon(
      aes(ymin = L, ymax = U, fill = loss),
      alpha = 0.25,
      colour = NA
    ) +
    geom_line(aes(y = med), linewidth = 1) +
    geom_line(
      data = opt_fit,
      aes(y = value),
      linetype = "dashed",
      linewidth = 1
    ) +
    facet_wrap(~group, ncol = 4, scales = "free_y",
               labeller = ggplot2::labeller(group = sensitivity_group_labeller)) +
    labs(
      x = "Year",
      y = "ST-SC group frequency",
      colour = "Loss function",
      fill = "Loss function"
    ) +
    ggplot2::scale_colour_manual(values = SENS_LOSS_COLOURS, labels = loss_labeller) +
    ggplot2::scale_fill_manual(values = SENS_LOSS_COLOURS, labels = loss_labeller) +
    sensitivity_paper_theme()
}

#' @describeIn sensitivity_plots All simulated trajectories, for all or some losses.
#' @importFrom dplyr filter mutate
#' @importFrom ggplot2 ggplot aes geom_line facet_wrap labs
#'   theme
#' @export
plot_sensitivity_trajectories <- function(
  sens,
  top_groups,
  start_year,
  losses = NULL
) {
  opt_fit <- sensitivity_opt_fit(sens, top_groups, start_year)
  d <- sens$sim_out |> filter(group %in% top_groups$group)
  if (!is.null(losses)) {
    d <- d |> filter(loss %in% losses)
    opt_fit <- opt_fit |> filter(loss %in% losses)
  }
  d |>
    mutate(year = start_year + year) |>
    ggplot(aes(x = year, colour = loss)) +
    geom_line(
      aes(y = value, colour = loss, group = interaction(rep, loss)),
      alpha = 0.1
    ) +
    geom_line(
      data = opt_fit,
      aes(y = value),
      linetype = "dashed",
      linewidth = 1
    ) +
    facet_wrap(~group, ncol = 4, scales = "free_y",
               labeller = ggplot2::labeller(group = sensitivity_group_labeller)) +
    labs(
      x = "Year",
      y = "ST-SC group frequency",
      colour = "Loss function",
      fill = "Loss function"
    ) +
    ggplot2::scale_colour_manual(values = SENS_LOSS_COLOURS, labels = loss_labeller) +
    ggplot2::scale_fill_manual(values = SENS_LOSS_COLOURS, labels = loss_labeller) +
    sensitivity_paper_theme()
}

#' @describeIn sensitivity_plots The data with the optimal fits and the
#'   sensitivity ribbons, for all losses or only `best_loss`.
#' @importFrom dplyr filter
#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_segment geom_ribbon
#'   labs facet_wrap theme
#' @export
plot_sensitivity_vs_data <- function(
    sens,
    xt,
    top_groups,
    start_year,
    qtile,
    best_loss = NULL
) {

  opt_fit <- sensitivity_opt_fit(
    sens, top_groups, start_year
  )

  sim_summ <- sensitivity_summary(
    sens, top_groups, start_year, qtile
  )

  obs <- xt |>
    filter(group %in% top_groups$group)

  # Initialise plot
  p <- ggplot(obs, aes(x = year))

  if (is.null(best_loss)) {

    # Sensitivity ribbons for all loss functions
    p <- p +
      geom_ribbon(
        data = sim_summ,
        aes(ymin = L, ymax = U, fill = loss),
        colour = NA,
        alpha = 0.2
      ) +
      ggplot2::scale_colour_manual(
        values = SENS_LOSS_COLOURS,
        labels = loss_labeller
      ) +
      ggplot2::scale_fill_manual(
        values = SENS_LOSS_COLOURS,
        labels = loss_labeller
      )

  } else {

    # Sensitivity ribbon for the selected fit
    p <- p +
      geom_ribbon(
        data = filter(sim_summ, loss == best_loss),
        aes(ymin = L, ymax = U),
        fill = "#6BAED6",
        colour = NA,
        alpha = 0.25
      )
  }

  # Observed uncertainty intervals and trajectories
  p <- p +
    geom_segment(
      aes(
        y = frac - frac_sd * 1.96,
        yend = frac + frac_sd * 1.96,
        xend = year
      ),
      colour = "grey70",
      linewidth = 0.5
    ) +
    geom_line(
      aes(y = frac),
      linewidth = 0.5,
      linetype = "dashed",
      colour = "grey65"
    )

  # Add fitted trajectories above the ribbons
  if (is.null(best_loss)) {

    p <- p +
      geom_line(
        data = opt_fit,
        aes(y = value, colour = loss),
        linewidth = 1
      )

  } else {

    p <- p +
      geom_line(
        data = filter(opt_fit, loss == best_loss),
        aes(y = value),
        colour = "#2171B5",
        linewidth = 1
      )
  }

  # Add observed points last so they remain fully visible
  p +
    geom_point(
      aes(y = frac),
      size = 1.5,
      colour = "black",
      alpha = 1
    ) +
    labs(
      x = "Year",
      y = "ST-SC group frequency",
      colour = "Loss function",
      fill = "Loss function"
    ) +
    facet_wrap(
      ~group,
      scales = "free_y",
      labeller = ggplot2::labeller(
        group = sensitivity_group_labeller
      )
    ) +
    sensitivity_paper_theme()
}

#' Add a sensitivity band to a top_group_plot()
#'
#' Adds the quantile band of the sensitivity draws for one fit underneath a
#' plot made with [top_group_plot()], on the same frequency scale and time axis.
#'
#' @param p A plot from [top_group_plot()].
#' @param sens Result of [run_sensitivity()].
#' @param key Id (`loss` column in `sens`) of the fit to show the band for.
#' @param vax_start The `vax_start` used in `p`.
#' @param qtile Quantile range of the band.
#' @param fill,alpha Fill colour and transparency of the band.
#' @param pcv13 Character vector of PCV13 serotypes.
#' @param vt_colour,nvt_colour Colours for vaccine and non-vaccine facet labels.
#' @param line_colour Colour for the single selected simulation trajectory.
#' @return The plot with the band underneath all other layers.
#' @importFrom dplyr %>% filter group_by mutate ungroup summarise
#' @importFrom ggplot2 geom_ribbon aes
#' @importFrom stats quantile
#' @export
add_sensitivity_ribbon <- function(
  p,
  sens,
  key,
  vax_start,
  qtile = 0.95,
  fill = "#6BAED6",
  alpha = 0.25,
  pcv13 = c("1", "3", "4", "5", "6A", "6B", "7F", "9V", "14",
            "18C", "19A", "19F", "23F"),
  vt_colour = "#A94C22",
  nvt_colour = "#084594",
  line_colour = "#2171B5"
) {
# the groups shown in the plot (the first layer of top_group_plot() is the data error bars)
  plot_groups <- levels(p$layers[[1]]$data$group)
# top_group_plot() shows simulations as frequencies, with model year k at vax_start - 1 + k
  ribbon <- sens$sim_out %>%
    filter(loss == key, year > 0) %>% # year 0 = x0, shown separately
    group_by(rep, year) %>%
    mutate(frac = value / sum(value)) %>% # normalise within each draw and year
    ungroup() %>%
    filter(group %in% plot_groups) %>%
    group_by(group, year) %>%
    summarise(
      L = quantile(frac, (1 - qtile) / 2),
      U = quantile(frac, 1 - (1 - qtile) / 2),
      .groups = "drop"
    ) %>%
    mutate(
      year = vax_start - 1 + year,
      group = factor(group, levels = plot_groups)
    )
  p$layers <- c(
    list(geom_ribbon(
      data = ribbon,
      aes(x = year, ymin = L, ymax = U),
      fill = fill,
      alpha = alpha,
      inherit.aes = FALSE
    )),
    p$layers
  )
  # Preserve the existing facet layout, scales and ordering.
  p$facet <- ggplot2::ggproto(NULL, p$facet)
  p$facet$params$labeller <- ggplot2::labeller(
    group = function(x) sensitivity_group_labeller(
      x, pcv13 = pcv13, vt_colour = vt_colour, nvt_colour = nvt_colour
    )
  )

  # top_group_plot() supplies the trajectory colour mapping. Replace its
  # scale for this single-fit figure, preserving the simulation legend.
  p +
    ggplot2::scale_colour_manual(values = line_colour) +
    sensitivity_paper_theme() +
    ggplot2::labs(title = NULL, subtitle = NULL, caption = NULL)
}

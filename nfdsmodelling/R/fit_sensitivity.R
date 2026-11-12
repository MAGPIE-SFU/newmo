# ============================================================
#  Assess sensitivity of the NFDS model to its input parameters by perturbing
#  inputs near the optimized parameter settings.
# ============================================================

# Display names for the loss functions in the plots
LOSS_LABELS <- c(
  SSE                           = "SSE",
  weighted_composition_distance = "Weighted composition",
  weighted_SSE                  = "Weighted SSE"
)
loss_labeller <- function(x) ifelse(x %in% names(LOSS_LABELS), LOSS_LABELS[x], x)

# default theme for the plots
#' @importFrom ggplot2 theme element_rect element_line
clean_theme <- function() {
  theme(panel.background = element_rect(fill="white", colour="black"),
        panel.grid.major = element_line("grey85"),
        panel.grid.minor = element_line('grey95'))
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
    rownames_to_column(var="group") |>
    mutate(frac = IC / sum(IC)) |>
    arrange(desc(frac)) |>
    slice_head(n=n_groups)
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
run_sensitivity <- function(opt_pars, samodel, tswitch, nyears, sd_pars,
                            use_switch = TRUE, n_sims = 100, n_cores = 1) {

  # pull out initial conditions
  IC <- data.frame(value = samodel$x0) |>
    rownames_to_column(var="group") |>
    bind_cols(year=0)

  # Optimal trajectories
  # ==============================================================================

  # process results for each of the optimal parameter settings
  opt_traj <- lapply(seq_len(nrow(opt_pars)), function(i) {
    run_model(as.list(unlist(opt_pars[i, PAR_NAMES])), use_switch, samodel, nyears, tswitch) |>
      rownames_to_column(var="group") |>
      pivot_longer(-group, names_to="year", names_prefix="year.", names_transform = as.numeric) |>
      bind_rows(IC) |>
      bind_cols(loss = opt_pars$loss[i])
  }) |> list_rbind()

  # Perturbations
  # ==============================================================================

  # For each parameter simulate n_sims draws from a normal distribution with mean
  #  equal to the optimized values and standard deviation loosely based on the
  #  loss function profiles

  par_lookup <- lapply(seq_len(nrow(opt_pars)), function(i) {as.data.frame(opt_pars[i,]) |>
      mutate(logv1 = log(v1),
             logsigma = log(sigma),
             portion = portion,
             .keep = "none") |>
      pivot_longer(everything(), names_to="par") |>
      left_join(sd_pars, by="par")
    }
  )

  # draw from a normal truncated to [lo, hi] via the inverse CDF
  sim_pars <- Vectorize(function(lookup, par_name){
    row <- lookup[lookup$par == par_name, ]
    p <- runif(1, pnorm(row$lo, row$value, row$SD), pnorm(row$hi, row$value, row$SD))
    qnorm(p, row$value, row$SD)},
    vectorize.args="par_name")

  par_df <- lapply(seq_len(nrow(opt_pars)), function(j) {
    lapply(1:n_sims, function(i){
      data.frame(par = sd_pars$par) |>
        mutate(val = sim_pars(par_lookup[[j]], par)) |>
        bind_cols(rep=i)
    }) |>
      list_rbind() |>
      bind_cols(loss=opt_pars$loss[j])
  }) |>
    list_rbind()

  # pull out parameters used in the optimal fits
  opt_settings <- lapply(seq_along(par_lookup), function(i) {
    par_lookup[[i]] |>
      select(par, value) |>
      bind_cols(loss = opt_pars$loss[i])
  }) |> list_rbind()

  sim_par_sets <- par_df |>
    pivot_wider(names_from="par", values_from = "val") |>
    mutate(sigma = exp(logsigma),
           v1 = exp(logv1),
           .keep = "unused")

  stopifnot(all(between(sim_par_sets$portion, 0, 1)))

  # run the model for each perturbed parameter set, keeping the rep/loss labels
  #  alongside the trajectories so the results can be summarised and plotted
  run_one <- function(k) {
    pars_k <- sim_par_sets[k, ]
    # a bad parameter draw shouldn't kill the whole batch; drop it and record why
    tryCatch(
      run_model(as.list(select(pars_k, -rep, -loss)), use_switch, samodel, nyears, tswitch) |>
        rownames_to_column(var="group") |>
        pivot_longer(-group, names_to="year", names_prefix="year.", names_transform = as.numeric) |>
        bind_rows(IC) |>
        bind_cols(select(pars_k, rep, loss)),
      error = function(e) {
        data.frame(rep = pars_k$rep, loss = pars_k$loss, error = conditionMessage(e))
      }
    )
  }

  # ODE solve times vary a lot across parameter draws, so hand out one job at a
  #  time (no prescheduling) rather than one fixed block per core
  sim_res <- if (n_cores > 1) {
    parallel::mclapply(seq_len(nrow(sim_par_sets)), run_one,
                       mc.cores = n_cores, mc.preschedule = FALSE)
  } else {
    lapply(seq_len(nrow(sim_par_sets)), run_one)
  }

  failed <- keep(sim_res, \(x) "error" %in% names(x)) |> list_rbind()
  if (NROW(failed) > 0) {
    warning(nrow(failed), " parameter sets failed; see `failed`")
  } else {
    failed <- data.frame(rep = integer(), loss = character(), error = character())
  }

  sim_out <- discard(sim_res, \(x) "error" %in% names(x)) |>
    list_rbind()

  list(
    opt_traj     = opt_traj,
    par_df       = par_df,
    opt_settings = opt_settings,
    sim_out      = sim_out,
    failed       = failed
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
    summarise(L = quantile(value, probs=(1-qtile)/2),
              med = quantile(value, probs=0.5),
              U = quantile(value, probs=1-(1-qtile)/2),
              .groups = "drop") |>
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
#' @importFrom ggplot2 ggplot aes geom_point geom_line facet_wrap theme scale_colour_brewer
#' @export
plot_sensitivity_optimal <- function(sens, top_groups) {
  sens$opt_traj |>
    filter(group %in% top_groups$group) |>
    ggplot(aes(x=year, y=value, colour=loss)) +
    geom_point() +
    geom_line() +
    facet_wrap(~group, ncol=4, scales="free_y") +
    scale_colour_brewer(palette="Dark2", labels = loss_labeller) +
    clean_theme() +
    theme(legend.position="bottom")
}

#' @describeIn sensitivity_plots Range of the simulated parameters (sampling scale).
#' @importFrom ggplot2 ggplot aes geom_histogram facet_grid geom_vline
#' @export
plot_sensitivity_draws <- function(sens) {
  sens$par_df |>
    ggplot(aes(x=val)) +
    geom_histogram(aes(fill=loss), colour="black") +
    facet_grid(loss~par, scales="free") +
    geom_vline(data=sens$opt_settings, aes(xintercept=value), linewidth=1, linetype="dashed") +
    clean_theme()
}

#' @describeIn sensitivity_plots Ribbons with the optimal trajectories (dashed).
#' @importFrom ggplot2 ggplot aes geom_ribbon geom_line facet_wrap labs
#'   scale_colour_brewer scale_fill_brewer theme
#' @export
plot_sensitivity_ribbons <- function(sens, top_groups, start_year, qtile) {
  opt_fit <- sensitivity_opt_fit(sens, top_groups, start_year)
  sensitivity_summary(sens, top_groups, start_year, qtile) |>
    ggplot(aes(x=year, colour=loss)) +
    geom_ribbon(aes(ymin=L, ymax=U, fill=loss), alpha=0.25, colour=NA) +
    geom_line(aes(y=med), linewidth=1) +
    geom_line(data=opt_fit, aes(y=value), linetype="dashed", linewidth=1) +
    facet_wrap(~group, ncol=4, scales="free_y") +
    labs(x="Year", y="ST-SC group frequency", colour="Loss function", fill="Loss function") +
    scale_colour_brewer(palette="Dark2", labels = loss_labeller) +
    scale_fill_brewer(palette="Dark2", labels = loss_labeller) +
    clean_theme() +
    theme(legend.position="bottom")
}

#' @describeIn sensitivity_plots All simulated trajectories, for all or some losses.
#' @importFrom dplyr filter mutate
#' @importFrom ggplot2 ggplot aes geom_line facet_wrap labs scale_colour_brewer
#'   scale_fill_brewer theme
#' @export
plot_sensitivity_trajectories <- function(sens, top_groups, start_year, losses = NULL) {
  opt_fit <- sensitivity_opt_fit(sens, top_groups, start_year)
  d <- sens$sim_out |> filter(group %in% top_groups$group)
  if (!is.null(losses)) {
    d <- d |> filter(loss %in% losses)
    opt_fit <- opt_fit |> filter(loss %in% losses)
  }
  d |>
    mutate(year = start_year + year) |>
    ggplot(aes(x=year, colour=loss)) +
    geom_line(aes(y=value, colour=loss, group=interaction(rep,loss)), alpha=0.1) +
    geom_line(data=opt_fit, aes(y=value), linetype="dashed", linewidth=1) +
    facet_wrap(~group, ncol=4, scales="free_y") +
    labs(x="Year", y="ST-SC group frequency", colour="Loss function", fill="Loss function") +
    scale_colour_brewer(palette="Dark2", labels = loss_labeller) +
    scale_fill_brewer(palette="Dark2", labels = loss_labeller) +
    clean_theme() +
    theme(legend.position="bottom")
}

#' @describeIn sensitivity_plots The data with the optimal fits and the
#'   sensitivity ribbons, for all losses or only `best_loss`.
#' @importFrom dplyr filter
#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_segment geom_ribbon
#'   labs scale_colour_brewer scale_fill_brewer facet_wrap theme
#' @export
plot_sensitivity_vs_data <- function(sens, xt, top_groups, start_year, qtile, best_loss = NULL) {
  opt_fit  <- sensitivity_opt_fit(sens, top_groups, start_year)
  sim_summ <- sensitivity_summary(sens, top_groups, start_year, qtile)
  n_groups <- nrow(top_groups)

  base <- xt |>
    filter(group %in% top_groups$group) |>
    ggplot(aes(x=year)) +
    geom_point(aes(y=frac), size=1) +
    geom_line(aes(y=frac), linewidth=0.5, linetype="dashed") +
    geom_segment(aes(y=frac-frac_sd*1.96, yend=frac+frac_sd*1.96))

  if (is.null(best_loss)) {
    # summarise sensitivity for all loss functions
    base +
      geom_line(data=opt_fit, aes(y=value, colour=loss), linewidth=1) +
      geom_ribbon(data=sim_summ, aes(ymin=L, ymax=U, fill=loss), alpha=0.2) +
      labs(x="Year", y="ST-SC group frequency", colour="Loss function", fill="Loss function",
           title=paste0("Top ",n_groups, " ST-GPSC groups. Ribbons show ", qtile*100, "% trajectory band.")) +
      scale_colour_brewer(palette="Dark2", labels = loss_labeller) +
      scale_fill_brewer(palette="Dark2", labels = loss_labeller) +
      facet_wrap(~group, scales="free_y") +
      clean_theme() +
      theme(legend.position="bottom")
  } else {
    # summarise sensitivity for just the chosen loss function
    base +
      geom_line(data=filter(opt_fit, loss==best_loss), aes(y=value), linewidth=1) +
      geom_ribbon(data=filter(sim_summ, loss==best_loss), aes(ymin=L, ymax=U), alpha=0.2) +
      labs(x="Year", y="ST-SC group frequency",
           title=paste0("Top ",n_groups, " ST-GPSC groups. Ribbons show ", qtile*100, "% range."),
           subtitle=paste("Based on the optimal fit from the", best_loss, "function.")) +
      facet_wrap(~group, scales="free_y") +
      clean_theme() +
      theme(legend.position="bottom")
  }
}

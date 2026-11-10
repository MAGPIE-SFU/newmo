# ============================================================
#  fit_diagnostics.R
#  Diagnostic plots for the fits:
#  (1) 1D loss profiles through the best fit of every seed x loss
#      (was profile_best_fits.R)
#  (2) the loss landscape from many random starts
#      (was visualize_landscape.R)
# ============================================================

# ------------------------------------------------------------
# Shared plotting settings
# ------------------------------------------------------------

# Plot these parameters on a log10 axis (match the log-scale starts in make_starts)
LOG_AXIS <- c(v1 = TRUE, sigma = TRUE, portion = FALSE)

axis_labels <- function() {
  setNames(
    ifelse(LOG_AXIS, paste0("log(", names(LOG_AXIS), ")"), names(LOG_AXIS)),
    names(LOG_AXIS)
  )
}

#' @importFrom ggplot2 theme_bw theme element_rect element_blank unit element_text
theme_panels <- function(legend_position) {
  theme_bw() +
    theme(
      panel.border     = element_rect(colour = "grey30", fill = NA, linewidth = 0.8),
      panel.grid.minor = element_blank(),
      panel.spacing    = unit(0.9, "lines"),
      strip.background = element_rect(fill = "grey90", colour = "grey30"),
      strip.text       = element_text(face = "bold"),
      legend.position  = legend_position
    )
}

# n_cores > 1 uses parallel::mclapply (Linux/macOS only)
par_lapply <- function(X, FUN, n_cores) {
  if (n_cores > 1) parallel::mclapply(X, FUN, mc.cores = n_cores) else lapply(X, FUN)
}

# ============================================================
# (1) 1D loss profiles through the best fit of every seed x loss.
#     Vary one parameter at a time, hold the others at the fitted values.
# ============================================================

#' 1D loss profiles through the best fits
#'
#' For the best fit of every seed x loss, varies one parameter at a time over
#' a grid (holding the others at the fitted values) and evaluates the loss.
#'
#' @param best Best fit per seed x loss, e.g. `select_best_fits()$all_best`.
#' @param samodel,xt,tswitch,nyears,use_switch As in [run_sequential_fits()].
#' @param losses Loss functions to profile.
#' @param n_cores Cores for the profiles (> 1 uses mclapply, Linux/macOS only).
#' @param profile_decades Log-scale parameters: +/- this many decades around the fit.
#' @param n_grid Grid points per parameter.
#' @return One row per (seed, loss, parameter, grid value).
#' @importFrom dplyr filter select all_of bind_rows n_distinct
#' @export
profile_best_fits <- function(best, samodel, xt, tswitch, nyears,
                              use_switch      = TRUE,
                              losses          = names(loss_functions),
                              n_cores         = 1,
                              profile_decades = 1.5,
                              n_grid          = 50) {

  fit_data <- load_fit_data(xt, nyears)

  # ------------------------------------------------------------
  # Best fit per seed x loss
  # ------------------------------------------------------------
  best <- best |>
    filter(loss %in% losses) |>
    select(seed, loss, fit_id, start_id, method, fitted, all_of(PAR_NAMES))
  stopifnot(nrow(best) > 0, !anyNA(best[PAR_NAMES]))
  message("Profiling ", nrow(best), " fits (", n_distinct(best$seed), " seeds x ",
          n_distinct(best$loss), " losses)")

  # ------------------------------------------------------------
  # Evaluation helpers
  # ------------------------------------------------------------
  eval_par <- function(p, loss_name) {
    sim <- tryCatch(run_model(as.list(p), use_switch, samodel, nyears, tswitch),
                    error = function(e) e)
    if (inherits(sim, "error")) return(NA_real_)
    v <- tryCatch(loss_functions[[loss_name]](sim, fit_data$data, fit_data$sd),
                  error = function(e) NA_real_)
    if (is.finite(v)) v else NA_real_
  }

  profile_grid <- function(par_name, center_value) {
    lo <- PAR_LOWER[[par_name]] + 1e-6     # stay strictly inside the bounds
    hi <- PAR_UPPER[[par_name]] - 1e-6
    cv <- min(max(center_value, lo), hi)
    if (LOG_AXIS[[par_name]]) {
      lo <- max(lo, cv / 10^profile_decades)
      hi <- min(hi, cv * 10^profile_decades)
      g  <- 10^seq(log10(lo), log10(hi), length.out = n_grid)
    } else {
      g  <- seq(lo, hi, length.out = n_grid)
    }
    sort(unique(c(g, cv)))
  }

  centers <- lapply(seq_len(nrow(best)), function(i) unlist(best[i, PAR_NAMES]))

  # Fail loudly if the model can't even run at the fitted points
  best$center_loss <- unlist(par_lapply(seq_len(nrow(best)), function(i)
    eval_par(centers[[i]], best$loss[i]), n_cores))
  if (anyNA(best$center_loss)) {
    bad <- which(is.na(best$center_loss))
    print(as.data.frame(best[bad, c("seed", "loss", "fit_id", PAR_NAMES)]))
    i <- bad[1]
    run_model(as.list(centers[[i]]), use_switch, samodel, nyears, tswitch)   # rerun unprotected to show the error
    stop("Loss is NA at the fitted parameters for the fits listed above")
  }

  # ------------------------------------------------------------
  # Profiles: one job per (fit, parameter)
  # ------------------------------------------------------------
  jobs <- expand.grid(i = seq_len(nrow(best)), par = PAR_NAMES, stringsAsFactors = FALSE)

  t_prof <- system.time(
    prof <- bind_rows(par_lapply(seq_len(nrow(jobs)), function(j) {
      i  <- jobs$i[j]; nm <- jobs$par[j]
      center <- centers[[i]]
      g <- profile_grid(nm, center[[nm]])
      vals <- vapply(g, function(x) { p <- center; p[[nm]] <- x; eval_par(p, best$loss[i]) },
                     numeric(1))
      data.frame(seed = best$seed[i], loss = best$loss[i], parameter = nm,
                 par_value = g, loss_value = vals,
                 center_value = center[[nm]], center_loss = best$center_loss[i])
    }, n_cores))
  )[["elapsed"]]
  message("Profiles done in ", round(t_prof), " s; ",
          sum(is.na(prof$loss_value)), " of ", nrow(prof), " evaluations failed")

  prof
}

#' Plot the 1D loss profiles
#'
#' Rows = loss, columns = parameter, one line per seed.
#'
#' @param prof Profiles from [profile_best_fits()].
#' @param losses Loss functions, in plotting order.
#' @param relative_to "best_seed": ratio to the lowest fitted loss across seeds
#'   (shows which seed fit best). "own": ratio to each seed's own fitted loss
#'   (every curve touches 1).
#' @param ymax Cap on the loss ratio, so one bad region doesn't flatten the rest.
#' @return A ggplot.
#' @importFrom dplyr group_by mutate ungroup distinct
#' @importFrom ggplot2 ggplot aes geom_hline geom_vline geom_line geom_point
#'   scale_colour_brewer facet_grid labeller scale_y_log10 labs
#' @export
plot_profiles <- function(prof, losses = names(loss_functions),
                          relative_to = "best_seed", ymax = 10) {

  prof <- prof |>
    group_by(loss) |>
    mutate(ref_loss = if (relative_to == "own") center_loss else min(center_loss)) |>
    ungroup() |>
    mutate(
      parameter    = factor(parameter, levels = PAR_NAMES),
      loss         = factor(loss, levels = losses),
      seed         = factor(seed),
      is_log       = LOG_AXIS[as.character(parameter)],
      x            = ifelse(is_log, log10(par_value), par_value),
      x_center     = ifelse(is_log, log10(center_value), center_value),
      ratio        = pmin(loss_value / ref_loss, ymax),
      center_ratio = pmin(center_loss / ref_loss, ymax)
    )

  centers_df <- distinct(prof, seed, loss, parameter, x_center, center_ratio)

  ggplot(prof, aes(x, ratio, colour = seed)) +
    geom_hline(yintercept = 1, linetype = "dotted", colour = "grey50") +
    geom_vline(data = centers_df, aes(xintercept = x_center, colour = seed),
               linetype = "dashed", alpha = 0.4) +
    geom_line(alpha = 0.8, na.rm = TRUE) +
    geom_point(data = centers_df, aes(x_center, center_ratio), size = 2.2) +
    scale_colour_brewer(palette = "Dark2", name = "Seed") +
    facet_grid(loss ~ parameter, scales = "free_x",
               labeller = labeller(parameter = axis_labels())) +
    scale_y_log10() +
    labs(
      title    = "1D loss profiles through the best fit of each seed",
      subtitle = paste0(
        "dots / dashed lines = fitted values; other parameters held at that seed's fit; ",
        if (relative_to == "own") "loss relative to each seed's own fit"
        else "loss relative to the best fit across seeds",
        "; capped at ", ymax, "x"),
      x = NULL,
      y = if (relative_to == "own") "loss / own fitted loss (log)" else "loss / best fitted loss (log)"
    ) +
    theme_panels("bottom")
}

# ============================================================
# (2) Loss landscape: screen many random starts, evaluate every loss
#     function on each, and plot loss against each parameter.
# ============================================================

#' Screen random starts for the loss landscape
#'
#' @inheritParams run_sequential_fits
#' @param n_screen Number of random starts (the default start is added too).
#' @param losses Loss functions to evaluate.
#' @param n_cores Cores (> 1 uses mclapply, Linux/macOS only).
#' @return One row per start: start_id, the parameters, and one column per loss.
#' @export
screen_landscape <- function(samodel, xt, default_pars, seed, tswitch, nyears,
                             use_switch = TRUE,
                             n_screen   = 1000,
                             losses     = names(loss_functions),
                             n_cores    = 1) {

  fit_data <- load_fit_data(xt, nyears)

  # ------------------------------------------------------------
  # Screen: simulate once per start, apply every loss to that simulation
  # ------------------------------------------------------------
  set.seed(seed)

  starts <- make_starts(n_random = n_screen, default_pars = default_pars)      # includes "default"

  eval_start <- function(s) {
    sim <- tryCatch(run_model(as.list(s), use_switch, samodel, nyears, tswitch),
                    error = function(e) NULL)
    vapply(loss_functions[losses], function(f) {
      if (is.null(sim)) return(Inf)
      v <- tryCatch(f(sim, fit_data$data, fit_data$sd), error = function(e) Inf)
      if (is.finite(v)) v else Inf
    }, numeric(1))
  }

  elapsed <- system.time(
    results <- if (n_cores > 1) {
      parallel::mclapply(starts, eval_start, mc.cores = n_cores)
    } else {
      lapply(starts, eval_start)
    }
  )[["elapsed"]]
  message("Screened ", length(starts), " starts in ", round(elapsed), " s")

  screen <- data.frame(
    start_id = names(starts),
    do.call(rbind, starts),
    do.call(rbind, results),
    row.names = NULL,
    check.names = FALSE
  )

  # Failures (errors or non-finite losses)
  for (l in losses) {
    message(l, ": ", sum(!is.finite(screen[[l]])), " of ", nrow(screen), " starts failed")
  }

  screen
}

#' Plot the loss landscape
#'
#' Loss against each parameter, one row per loss function.
#'
#' @param screen Screen from [screen_landscape()].
#' @param losses Loss functions, in plotting order.
#' @return A ggplot.
#' @importFrom dplyr filter mutate group_by ungroup all_of
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_point scale_colour_brewer facet_grid
#'   labeller scale_y_log10 labs
#' @export
plot_landscape <- function(screen, losses = names(loss_functions)) {

  n_screen <- sum(screen$start_id != "default")

  # ------------------------------------------------------------
  # Long format for plotting
  # ------------------------------------------------------------
  long <- screen |>
    pivot_longer(all_of(PAR_NAMES), names_to = "parameter", values_to = "par_value") |>
    pivot_longer(all_of(losses),    names_to = "loss",      values_to = "loss_value") |>
    filter(is.finite(loss_value), loss_value > 0) |>
    mutate(
      parameter = factor(parameter, levels = PAR_NAMES),
      loss      = factor(loss, levels = losses),
      x         = ifelse(LOG_AXIS[as.character(parameter)], log10(par_value), par_value)
    )

  # ------------------------------------------------------------
  # Plot: loss against each parameter, one row per loss function
  # ------------------------------------------------------------
  best <- long |>
    group_by(loss) |>
    filter(loss_value == min(loss_value)) |>
    ungroup()

  ggplot(long, aes(x, loss_value)) +
    geom_point(aes(colour = loss), alpha = 0.3, size = 0.7) +
    scale_colour_brewer(palette = "Dark2") +
    geom_point(data = filter(long, start_id == "default"), colour = "red", size = 2.5) +
    geom_point(data = best, colour = "dodgerblue3", size = 2.5) +
    facet_grid(loss ~ parameter, scales = "free",
               labeller = labeller(parameter = axis_labels())) +
    scale_y_log10() +
    labs(
      title = paste0("Loss vs parameter: ", n_screen, " random starts"),
      subtitle = "red = default start, blue = best start for that loss",
      x = NULL, y = "loss (log scale)"
    ) +
    theme_panels("none")
}

#' Best screened start per loss
#'
#' @param screen Screen from [screen_landscape()].
#' @param losses Loss functions.
#' @return One row per loss: the best start, its loss, the default start's loss
#'   and the parameters.
#' @export
landscape_best_starts <- function(screen, losses = names(loss_functions)) {
  do.call(rbind, lapply(losses, function(l) {
    i <- which.min(screen[[l]])
    data.frame(loss = l, start_id = screen$start_id[i], best_loss = screen[[l]][i],
               default_loss = screen[[l]][screen$start_id == "default"],
               screen[i, PAR_NAMES])
  }))
}
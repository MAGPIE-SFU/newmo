# ============================================================
#  fitting.R
#  Everything needed to fit the ODE model: parameter bookkeeping,
#  data, model wrapper, loss functions, starting points, optimizers,
#  summaries, the sequential fits for one seed, and the selection of
#  the best fit across seeds.
# ============================================================

# ------------------------------------------------------------
# 1. Parameter bookkeeping
# ------------------------------------------------------------

#' Names of the fitted parameters
#' @export
PAR_NAMES <- c("v1", "sigma", "portion")

FIXED_MIGR <- 0.00001
FIXED_WLOW <- 0.0
FIXED_WHIGH <- 1

# PAR_LOWER <- c(v1 = 0.001, sigma = 0,    migr = 0, portion = 0, whigh = 0.002, wfrac = 0)
# PAR_UPPER <- c(v1 = 0.99,  sigma = 1000, migr = 1, portion = 1, whigh = 1,     wfrac = 1)
PAR_LOWER <- c(v1 = 0.001, sigma = 0, portion = 0)
PAR_UPPER <- c(v1 = 1, sigma = 10000, portion = 1)

# Keep a start strictly inside the bounds (the bounded optimizers dislike starting exactly on a boundary).
clamp_to_bounds <- function(
  par,
  lower = PAR_LOWER,
  upper = PAR_UPPER,
  eps = 1e-6
) {
  pmin(pmax(par, lower + eps), upper - eps)
}

# ------------------------------------------------------------
# 2. Data
# ------------------------------------------------------------

#' @importFrom dplyr %>% select all_of
#' @importFrom tidyr pivot_wider
#' @importFrom tibble column_to_rownames
pivot_to_matrix <- function(xt, value_col) {
  xt %>%
    select(group, year, all_of(value_col)) %>%
    pivot_wider(
      names_from = year,
      values_from = all_of(value_col),
      names_glue = "year.{year}"
    ) %>%
    column_to_rownames("group") %>%
    as.matrix()
}

#' Fit data from x(t)
#'
#' @param xt Grouped x(t) with columns group, year, frac, frac_sd.
#' @param nyears Number of (final) years to fit to.
#' @param drop_groups Groups to leave out, or NULL.
#' @return list(data = <group x year matrix>, sd = <matching SD matrix>).
#' @export
load_fit_data <- function(
  xt,
  nyears,
  drop_groups = NULL
) {
  data <- pivot_to_matrix(xt, "frac")
  sd <- pivot_to_matrix(xt, "frac_sd")

  data <- data[, (ncol(data) - nyears + 1):ncol(data), drop = FALSE]
  sd <- sd[, (ncol(sd) - nyears + 1):ncol(sd), drop = FALSE]

  stopifnot(identical(colnames(data), colnames(sd)))

  # Zero SDs would give infinite weights; use FIT_SD_FLOOR (or lambda > 0) to avoid them
  stopifnot("Zero SDs in the fit data" = all(sd > 0))

  if (!is.null(drop_groups)) {
    data <- data[!rownames(data) %in% drop_groups, , drop = FALSE]
    sd <- sd[!rownames(sd) %in% drop_groups, , drop = FALSE]
  }

  list(data = data, sd = sd)
}

# ------------------------------------------------------------
# 3. Model
# ------------------------------------------------------------

#' Simulate the model for a set of fitted parameters
#'
#' Simulation only: no data, no loss. Adds the fixed parameters and the
#' vaccines (PCV7, then PCV13) to `model_pars`.
#'
#' @param model_pars Named list or vector with the fitted parameters (PAR_NAMES).
#' @param use_switch Use the vaccine switch at `tswitch`.
#' @param samodel The model structure, e.g. from [build_samodel()].
#' @param nyears Number of years to simulate.
#' @param tswitch Time of the vaccine switch (used if `use_switch`).
#' @export
run_model <- function(
  model_pars = list(),
  use_switch,
  samodel,
  nyears,
  tswitch
) {
  pars <- as.list(model_pars)

  # Fixed Parameters are added:
  pars$whigh <- FIXED_WHIGH
  pars$wlow <- FIXED_WLOW
  pars$v2 <- pars$v1
  pars$migr <- FIXED_MIGR
  pars$tswitch <- if (use_switch) tswitch else 0

  pars <- c(pars,
            list(vaccine1 = VACCINE_TYPES[["PCV7"]],
                 vaccine2 = VACCINE_TYPES[["PCV13"]]),
            list(K = sum(samodel$x0)))

  simulate_nfds(samodel, pars, nyears)
}

# ------------------------------------------------------------
# 4. Loss functions: all have signature (sim, data, sd)
# ------------------------------------------------------------

#' Loss functions for fitting
#'
#' A named list; each has signature (sim, data, sd).
#' @export
loss_functions <- list(
  SSE = function(sim, data, sd) {
    sum((sim - data)^2)
  },

  weighted_SSE = function(sim, data, sd) {
    sum((1 / sd^2) * (sim - data)^2)
  },

  weighted_composition_distance = function(sim, data, sd, gamma = 1) {
    weights <- 1 / sd^2
    total <- 0

    for (j in seq_len(ncol(sim))) {
      x1 <- data[, j]^gamma
      x2 <- sim[, j]^gamma

      y1 <- x1 / sum(x1)
      y2 <- x2 / sum(x2)

      d_sq <- (y1 - y2)^2
      c_vec <- y1 + y2

      idx <- d_sq != 0
      d_sq[idx] <- d_sq[idx] / c_vec[idx]

      total <- total + sum(weights[, j] * d_sq)
    }
    (1 / gamma) * sqrt(2 * nrow(sim)) * sqrt(total)
  }
)

#' Turn a loss function into an objective of the parameter vector alone
#' @inheritParams run_model
#' @param data,sd Data and SD matrices, from [load_fit_data()].
#' @param loss_fn One of [loss_functions].
#' @export
make_loss <- function(
  data,
  sd,
  loss_fn,
  use_switch,
  samodel,
  nyears,
  tswitch
) {
  force(data)
  force(sd)
  force(loss_fn)
  force(use_switch)
  force(samodel)
  force(nyears)
  force(tswitch)

  function(par) {
    names(par) <- PAR_NAMES
    sim <- run_model(as.list(par), use_switch, samodel, nyears, tswitch)
    loss_fn(sim, data, sd)
  }
}

# Wrap a loss so the optimizer minimizes log(loss). Non-finite or negative
# values get a large finite penalty (log-scale losses are small, so 1e10 is
# plenty).
log_loss <- function(fn, penalty = 1e10) {
  force(fn)
  function(par) {
    v <- fn(par)
    if (!is.finite(v) || v < 0) {
      return(penalty)
    }
    log(max(v, .Machine$double.xmin))
  }
}

# ------------------------------------------------------------
# 5. Starting points
# ------------------------------------------------------------

# Ranges and sampling scale for RANDOM STARTS only (the fitting bounds
# PAR_LOWER/PAR_UPPER are unchanged). Log-scale needs a positive lower end,
# so sigma and migr can't start at exactly 0.
# These ranges are my placeholders: please set them to what's plausible.
# START_LOWER <- c(v1 = 0.001, sigma = 0.001, migr = 1e-8, portion = 0.01, whigh = 0.002, wfrac = 1e-3)
# START_UPPER <- c(v1 = 0.99,  sigma = 10,    migr = 1e-2, portion = 0.99, whigh = 1,     wfrac = 1)

START_LOWER <- c(v1 = 0.001, sigma = 0.001, portion = 0.01)
START_UPPER <- c(v1 = 1, sigma = 10000, portion = 0.99)

START_LOG <- c(v1 = TRUE, sigma = TRUE, portion = FALSE
)

# Latin hypercube on [0,1]^k: one point in each of n strata per dimension
#' @importFrom stats runif
lhs_unit <- function(n, k) {
  matrix(
    unlist(lapply(seq_len(k), function(j) (sample.int(n) - runif(n)) / n)),
    nrow = n
  )
}

# default_pars: the CC_MAGIC parameters, used as the "default" start
make_starts <- function(
  n_random = 0,
  default_pars,
  lower = PAR_LOWER,
  upper = PAR_UPPER,
  start_lower = START_LOWER,
  start_upper = START_UPPER,
  start_log = START_LOG
) {
  stopifnot(
    identical(names(lower), PAR_NAMES),
    identical(names(upper), PAR_NAMES),
    identical(names(start_lower), PAR_NAMES),
    identical(names(start_upper), PAR_NAMES),
    identical(names(start_log), PAR_NAMES),
    all(start_lower >= lower),
    all(start_upper <= upper),
    all(start_lower > 0)
  )

  base <- clamp_to_bounds(unlist(default_pars)[PAR_NAMES], lower, upper)
  stopifnot(!anyNA(base))
  starts <- list(default = base)

  if (n_random > 0) {
    u <- lhs_unit(n_random, length(PAR_NAMES))
    for (i in seq_len(n_random)) {
      s <- ifelse(
        start_log,
        exp(log(start_lower) + u[i, ] * (log(start_upper) - log(start_lower))),
        start_lower + u[i, ] * (start_upper - start_lower)
      )
      names(s) <- PAR_NAMES
      starts[[sprintf("start_%02d", i)]] <- clamp_to_bounds(s, lower, upper)
    }
  }
  starts
}

# ------------------------------------------------------------
# 6. Optimizers
# ------------------------------------------------------------

# L-BFGS-B with finite-difference gradients. `portion` is held fixed: it is
# rounded to a whole number of loci inside simulate(), so the loss is a step
# function of it and has no usable gradient.
GRADIENT_FIXED <- "portion"

#' @importFrom stats optim
lbfgsb <- function(par, fn, lower, upper, control = list(maxfeval = 20000)) {
  free <- !(names(par) %in% GRADIENT_FIXED)
  full <- function(p) {
    x <- par
    x[free] <- p
    x
  }
  n_free <- sum(free)

  r <- optim(
    par[free],
    function(p) fn(full(p)),
    method = "L-BFGS-B",
    lower = lower[free],
    upper = upper[free],
    control = list(
      # iterations allowed within the evaluation budget (1 + 2 * n_free per iteration)
      maxit = control$maxfeval %/% (1 + 2 * n_free),
      parscale = pmax(abs(par[free]), 1e-6), # relative steps despite very different scales
      ndeps = rep(1e-3, n_free)
    )
  )

  list(
    par = full(r$par),
    value = r$value,
    feval = r$counts[["function"]] * (1 + 2 * n_free), # counts include gradient evals
    convergence = r$convergence,
    message = r$message
  )
}

# Wrapped in functions so dfoptim is looked up when called, not when the package is built
OPTIMIZERS <- list(
  nmkb   = function(...) dfoptim::nmkb(...),
  hjkb   = function(...) dfoptim::hjkb(...),
  lbfgsb = function(...) lbfgsb(...)
)

# Run one optimizer and standardize its output. `value` is always the
# loss on the ORIGINAL scale (re-evaluated at the optimum), so fits are
# comparable whether or not use_log was used.
run_optimizer <- function(
  method,
  fn,
  initial,
  lower,
  upper,
  use_log = TRUE,
  maxfeval = 10000
) {
  optimizer <- OPTIMIZERS[[method]]
  if (is.null(optimizer)) {
    stop("Unknown method: ", method)
  }

  fn_used <- if (use_log) log_loss(fn) else fn

  r <- optimizer(
    par = initial,
    fn = fn_used,
    lower = lower,
    upper = upper,
    control = list(maxfeval = maxfeval)
  )

  list(
    par = r$par,
    value = fn(r$par),
    convergence = r$convergence,
    n_feval = if (is.null(r$feval)) NA_integer_ else r$feval,
    message = if (is.null(r$message)) NA_character_ else r$message
  )
}

# One fit = one (loss, method, start). Never throws: errors are recorded.
fit_one <- function(
  fn,
  loss_name,
  method,
  start_id,
  initial,
  lower = PAR_LOWER,
  upper = PAR_UPPER,
  use_log = TRUE,
  maxfeval = 20000
) {
  initial_loss <- tryCatch(fn(initial), error = function(e) NA_real_)

  elapsed <- system.time({
    result <- tryCatch(
      run_optimizer(method, fn, initial, lower, upper, use_log, maxfeval),
      error = function(e) {
        list(
          par = rep(NA_real_, length(PAR_NAMES)),
          value = NA_real_,
          convergence = NA_integer_,
          n_feval = NA_integer_,
          message = conditionMessage(e)
        )
      }
    )
  })[["elapsed"]]

  c(
    result,
    list(
      fit_id = paste(loss_name, method, start_id, sep = "__"),
      loss = loss_name,
      method = method,
      start_id = start_id,
      use_log = use_log,
      initial_par = initial,
      initial_loss = initial_loss,
      runtime_seconds = elapsed
    )
  )
}

# ------------------------------------------------------------
# 7. Summaries
# ------------------------------------------------------------

# Fitted parameters in the format the model uses downstream.
#' @importFrom stats setNames
format_parameters <- function(par, use_switch, tswitch) {
  p <- as.data.frame(as.list(setNames(par, PAR_NAMES)))

  # Should not be there, but keeping it for now.
  p$wfrac <- NULL

  p$whigh <- FIXED_WHIGH
  p$wlow <- FIXED_WLOW
  p$migr <- FIXED_MIGR

  p$tswitch <- if (use_switch) tswitch else 0
  p$v2 <- p$v1
  p
}

summarise_fits <- function(fits) {
  rows <- lapply(fits, function(f) {
    data.frame(
      fit_id = f$fit_id,
      loss = f$loss,
      method = f$method,
      start_id = f$start_id,
      use_log = f$use_log,
      initial = f$initial_loss,
      fitted = f$value,
      improvement = f$initial_loss - f$value,
      pct_change = 100 * (f$initial_loss - f$value) / abs(f$initial_loss),
      convergence = f$convergence,
      n_feval = f$n_feval,
      runtime_seconds = f$runtime_seconds,
      message = f$message,
      row.names = NULL
    )
  })
  tbl <- do.call(rbind, rows)
  tbl <- tbl[order(tbl$loss, tbl$fitted), ]
  row.names(tbl) <- NULL
  tbl
}

parameter_table <- function(fits, use_switch, tswitch) {
  rows <- lapply(fits, function(f) {
    if (!all(is.finite(f$par))) {
      return(NULL)
    }
    data.frame(
      fit_id = f$fit_id,
      loss = f$loss,
      method = f$method,
      start_id = f$start_id,
      format_parameters(f$par, use_switch, tswitch)
    )
  })
  tbl <- do.call(rbind, rows)
  row.names(tbl) <- NULL
  tbl
}

# Best (lowest fitted loss) fit for each loss x method.
best_per_group <- function(summary_table) {
  ok <- summary_table[is.finite(summary_table$fitted), ]
  ok <- ok[order(ok$loss, ok$method, ok$fitted), ]
  ok[!duplicated(ok[c("loss", "method")]), ]
}

# Previously save_fit_results(): now returns the tables instead of writing them.
collect_fit_results <- function(fits, use_switch, tswitch) {
  summary_tbl <- summarise_fits(fits)
  params_tbl <- parameter_table(fits, use_switch, tswitch)
  best <- best_per_group(summary_tbl)

  list(summary = summary_tbl, parameters = params_tbl, best = best)
}

# ------------------------------------------------------------
# 8. Sequential fits for one seed (was sequential_fits.R)
#    Chained refits: each round starts from the previous round's optimum.
#    Checks whether the fit is stable or keeps drifting / improving.
# ------------------------------------------------------------

#' Sequential (chained) fits for one seed
#'
#' Screens `n_screen` random starts (plus the default start), then for each
#' loss function runs `n_rounds` fits, each starting from the previous
#' round's optimum, cycling through `chain_methods`.
#'
#' @param samodel The model structure, e.g. from [build_samodel()].
#' @param xt Grouped x(t) to fit to.
#' @param seed Random seed for this run.
#' @param default_pars The CC_MAGIC parameters, used as the default start.
#' @param tswitch Time of the vaccine switch.
#' @param nyears Number of (final) years to fit.
#' @param use_switch,use_log Use the vaccine switch / minimise log(loss).
#' @param losses Names of the loss functions to fit.
#' @param chain_methods Optimizers, cycled through the rounds.
#' @param n_rounds Rounds per chain.
#' @param start_from "screen", "random", "default", or a numeric start vector.
#' @param n_screen Number of random starts to screen.
#' @param maxfeval Evaluation budget per fit.
#' @param rel_tol Stop a chain when the relative improvement falls below
#'   this, or NULL to always run all rounds.
#' @param parallel,n_workers Parallelise the start screen.
#'
#' @return A list with `seed`, `summary` and `parameters` (both with a
#'   "round_00" start row per loss and a seed column), `screen` (losses x
#'   starts, or NULL) and `fits` (all fit objects).
#'
#' @importFrom dplyr filter transmute bind_rows mutate arrange
#' @importFrom tibble as_tibble_row
#' @importFrom stats setNames
#' @export
run_sequential_fits <- function(samodel, xt, seed, default_pars, tswitch, nyears,
                                use_switch    = TRUE,
                                use_log       = TRUE,
                                losses        = names(loss_functions),
                                chain_methods = c("nmkb", "lbfgsb", "hjkb"),
                                n_rounds      = 6,
                                start_from    = "screen",
                                n_screen      = 20000,
                                maxfeval      = 20000,
                                rel_tol       = NULL,
                                parallel      = TRUE,
                                n_workers     = max(1, parallelly::availableCores() - 1)) {

  # ------------------------------------------------------------
  # Setup
  # ------------------------------------------------------------
  seed_value <- seed
  set.seed(seed)
  fit_data <- load_fit_data(xt, nyears)
  screen <- NULL

  if (identical(start_from, "screen")) {

    screen_starts <- make_starts(n_random = n_screen, default_pars = default_pars)   # includes "default"

    # One start -> vector of losses (one per loss function)
    # simulate once per start, apply every loss to that simulation
    screen_one <- function(s) {
      sim <- tryCatch(run_model(as.list(s), use_switch, samodel, nyears, tswitch),
                      error = function(e) NULL)
      vapply(loss_functions, function(f) {
        if (is.null(sim)) return(Inf)
        v <- tryCatch(f(sim, fit_data$data, fit_data$sd), error = function(e) Inf)
        if (is.finite(v)) v else Inf
      }, numeric(1))
    }

    # losses x starts
    if (parallel && n_workers > 1) {
      # forking where possible (Linux/macOS terminal, Docker); fresh sessions otherwise
      use_fork <- parallelly::supportsMulticore()
      if (use_fork) {
        future::plan(future::multicore, workers = n_workers)
      } else {
        future::plan(future::multisession, workers = n_workers)
      }

      t_screen <- system.time(
        screen <- future.apply::future_sapply(
          screen_starts, screen_one,
          future.seed       = TRUE,
          future.chunk.size = ceiling(length(screen_starts) / (n_workers * 4))
        )
      )[["elapsed"]]

      future::plan(future::sequential)
      message("Screened ", length(screen_starts), " starts on ", n_workers, " workers (",
              if (use_fork) "multicore" else "multisession", ") in ", round(t_screen), " s")
    } else {
      t_screen <- system.time(
        screen <- sapply(screen_starts, screen_one)
      )[["elapsed"]]
      message("Screened ", length(screen_starts), " starts sequentially in ",
              round(t_screen), " s")
    }

    best_idx <- apply(screen, 1, which.min)    # best start per loss function

  } else {
    start <- if (is.numeric(start_from)) {
      clamp_to_bounds(start_from[PAR_NAMES])
    } else {
      starts <- make_starts(n_random = as.integer(start_from == "random"),
                            default_pars = default_pars)
      starts[[length(starts)]]
    }
  }

  # ------------------------------------------------------------
  # Chains: round r starts where round r-1 ended
  # ------------------------------------------------------------
  fits         <- list()
  chain_starts <- list()   # starting parameters for each loss (for the report)

  for (loss_name in losses) {
    fn <- make_loss(fit_data$data, fit_data$sd, loss_functions[[loss_name]],
                    use_switch, samodel, nyears, tswitch)
    current <- if (identical(start_from, "screen")) {
        message("  best screened start: ", names(screen_starts)[best_idx[[loss_name]]],
                " (loss ", signif(screen[loss_name, best_idx[[loss_name]]], 4),
                "; default start: ", signif(screen[loss_name, "default"], 4), ")")
        screen_starts[[best_idx[[loss_name]]]]
      } else start

    # remember where this chain starts
    chain_starts[[loss_name]] <- clamp_to_bounds(current)

    for (r in seq_len(n_rounds)) {
      method <- chain_methods[(r - 1) %% length(chain_methods) + 1]
      message("Fitting: ", loss_name, " / round ", r, " / ", method)

      fit <- fit_one(
        fn, loss_name, method, sprintf("round_%02d", r),
        initial  = current,
        use_log  = use_log,
        maxfeval = maxfeval
      )
      fits[[fit$fit_id]] <- fit

      if (!all(is.finite(fit$par))) {
        message("  fit failed, stopping this chain: ", fit$message)
        break
      }

      current <- setNames(fit$par, PAR_NAMES)
      current <- clamp_to_bounds(current)

      if (!is.null(rel_tol) && is.finite(fit$initial_loss)) {
        rel_gain <- (fit$initial_loss - fit$value) / abs(fit$initial_loss)
        if (rel_gain < rel_tol) {
          message("  improvement below rel_tol, stopping this chain")
          break
        }
      }
    }
  }

  # ------------------------------------------------------------
  # Collect
  # ------------------------------------------------------------
  stopifnot("No fits were produced - check the chain loop above" = length(fits) > 0)

  results <- collect_fit_results(fits, use_switch, tswitch)

  # ------------------------------------------------------------
  # add each chain's starting point as a "round_00" row
  # ------------------------------------------------------------
  # Loss at the start = initial loss of round 1 (same point, same loss function)
  start_summary <- results$summary |>
    filter(start_id == "round_01") |>
    transmute(
      fit_id      = paste0(loss, "__start"),
      loss,
      start_id    = "round_00",
      method      = "start",
      initial,
      fitted      = initial,
      convergence = NA,
      n_feval     = 0L
    )

  start_pars <- bind_rows(lapply(names(chain_starts), function(l) {
    as_tibble_row(as.list(chain_starts[[l]][PAR_NAMES])) |>
      mutate(fit_id = paste0(l, "__start"), loss = l,
             method = "start", start_id = "round_00", .before = 1)
  }))

  summary_with_start    <- bind_rows(start_summary, results$summary)    |> arrange(loss, start_id)
  parameters_with_start <- bind_rows(start_pars,    results$parameters) |> arrange(loss, start_id)

  # Tag every row with the seed, so runs can be stacked and compared later
  summary_with_start    <- summary_with_start    |> mutate(seed = seed_value, .before = 1)
  parameters_with_start <- parameters_with_start |> mutate(seed = seed_value, .before = 1)

  list(
    seed       = seed_value,
    summary    = summary_with_start,
    parameters = parameters_with_start,
    screen     = screen,
    fits       = fits
  )
}

# ------------------------------------------------------------
# 9. Best fit across seeds (was the seed-comparison script)
# ------------------------------------------------------------

#' Pick the best fits across seeds
#'
#' Takes the best fit per seed x loss (excluding the round_00 starting
#' points), simulates every candidate, evaluates it on every loss function,
#' and scores it as the sum over loss functions of (loss / best loss across
#' all candidates). Lower is better; the minimum possible score is the number
#' of loss functions.
#'
#' @param seq_fits List of results from [run_sequential_fits()], one per seed.
#' @inheritParams run_sequential_fits
#'
#' @return A list with `all_parameters` and `all_summary` (all fits of all
#'   seeds), `all_best` (best fit per seed x loss, with key and score) and
#'   `overall_best` (the one row with the lowest score).
#'
#' @importFrom dplyr bind_rows filter group_by slice_min ungroup select left_join mutate
#' @importFrom stats setNames
#' @export
select_best_fits <- function(seq_fits, samodel, xt, tswitch, nyears, use_switch = TRUE) {

  summary_all <- bind_rows(lapply(seq_fits, `[[`, "summary"))
  params_all  <- bind_rows(lapply(seq_fits, `[[`, "parameters"))

  # ---- Best fit per seed x loss (excluding the round_00 starting points) --------
  best <- summary_all |>
    filter(method != "start", is.finite(fitted)) |>
    group_by(seed, loss) |>
    slice_min(fitted, n = 1, with_ties = FALSE) |>
    ungroup() |>
    select(seed, fit_id, loss, method, start_id, fitted) |>
    left_join(params_all |> select(-loss, -method, -start_id),
              by = c("seed", "fit_id")) |>
    mutate(key = paste0("seed", seed, "_", fit_id))   # fit_id repeats across seeds

  # ---- Simulate every candidate --------------------------------------------------
  # Same parameters as parameter_generator() built before: the fitted PAR_NAMES
  # plus the fixed values, tswitch and PCV7 -> PCV13.
  simulations <- setNames(lapply(seq_len(nrow(best)), function(i) {
    run_model(as.list(unlist(best[i, PAR_NAMES])), use_switch, samodel, nyears, tswitch)
  }), best$key)

  # ---- Cross-evaluate every candidate on every loss -----------------------------
  data_sd <- load_fit_data(xt, nyears)
  errors <- t(sapply(simulations, function(sim)
    sapply(loss_functions, function(f) f(sim, data_sd$data, data_sd$sd))
  ))
  rel <- sweep(errors, 2, apply(errors, 2, min), "/")   # relative to best across ALL seeds
  best$score <- rowSums(rel)[best$key]

  best_params <- best |> slice_min(score, n = 1, with_ties = FALSE)

  list(
    all_parameters = params_all,
    all_summary    = summary_all,
    all_best       = best,
    overall_best   = best_params
  )
}

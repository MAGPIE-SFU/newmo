# ============================================================
#  parameter_fitting_functions.R
#  Everything needed to fit the ODE model. Nothing is run here;
#  see run_fits.R for the fitting process.
# ============================================================

library(deSolve)
library(readr)
library(dplyr)
library(tibble)
library(tidyr)
library(dfoptim)

common <- new.env()
sys.source(here::here("code", "common.R"), envir = common)
source(here::here("code", "simulate_ode.R"))

# ------------------------------------------------------------
# 1. Parameter bookkeeping
# ------------------------------------------------------------

PAR_NAMES <- c("v1", "sigma", "migr", "portion", "whigh", "wfrac")

DEFAULT_PARS <- list(
  v1      = 0.065,
  sigma   = 0.13,
  migr    = 0.00001,
  portion = 0.3,
  whigh   = 0.1363,
  wfrac   = 0.0023 / 0.1363
)

FIXED_PARS <- list(
  vaccine1 = common$VACCINE_TYPES[["PCV7"]],
  vaccine2 = common$VACCINE_TYPES[["PCV13"]]
)

PAR_LOWER <- c(v1 = 0.001, sigma = 0,    migr = 0, portion = 0, whigh = 0.002, wfrac = 0)
PAR_UPPER <- c(v1 = 0.99,  sigma = 1000, migr = 1, portion = 1, whigh = 1,     wfrac = 1)

NYEARS <- common$NYEARS

# Keep a start strictly inside the bounds (the bounded optimizers dislike
# starting exactly on a boundary).
clamp_to_bounds <- function(par, lower = PAR_LOWER, upper = PAR_UPPER, eps = 1e-6) {
  pmin(pmax(par, lower + eps), upper - eps)
}

# ------------------------------------------------------------
# 2. Data
# ------------------------------------------------------------

pivot_to_matrix <- function(xt, value_col) {
  xt %>%
    select(group, year, all_of(value_col)) %>%
    pivot_wider(
      names_from  = year,
      values_from = all_of(value_col),
      names_glue  = "year.{year}"
    ) %>%
    column_to_rownames("group") %>%
    as.matrix()
}

# Returns list(data = <group x year matrix>, sd = <matching SD matrix>).
load_fit_data <- function(xt_path = common$XT_PATH,
                          nyears = NYEARS,
                          drop_groups = NULL) {
  xt <- readr::read_csv(xt_path, show_col_types = FALSE)

  data <- pivot_to_matrix(xt, "frac")
  sd   <- pivot_to_matrix(xt, "frac_sd")

  data <- data[, (ncol(data) - nyears + 1):ncol(data), drop = FALSE]
  sd   <- sd[, (ncol(data) - nyears + 1):ncol(data), drop = FALSE]

  # Avoid zero SDs (they would give infinite weights)
  min_nonzero <- min(sd[sd > 0])
  sd[sd == 0] <- min_nonzero / 2

  if (!is.null(drop_groups)) {
    data <- data[!rownames(data) %in% drop_groups, , drop = FALSE]
    sd   <- sd[!rownames(sd) %in% drop_groups, , drop = FALSE]
  }

  list(data = data, sd = sd)
}

# ------------------------------------------------------------
# 3. Model
# ------------------------------------------------------------

# Cached per path, so readRDS() runs once instead of once per evaluation.
get_samodel <- local({
  cache <- list()
  function(path = common$SA_MODEL_PATH) {
    if (is.null(cache[[path]])) cache[[path]] <<- readRDS(path)
    cache[[path]]
  }
})

# Simulation only: no data, no loss.
run_model <- function(model_pars = list(),
                      use_switch,
                      sa_model_path = common$SA_MODEL_PATH,
                      nyears = NYEARS) {
  samodel <- get_samodel(sa_model_path)

  pars <- DEFAULT_PARS
  pars[names(model_pars)] <- model_pars

  # wfrac = wlow / whigh, so wlow <= whigh always holds
  if ("wfrac" %in% names(pars)) {
    pars$wlow  <- pars$whigh * pars$wfrac
    pars$wfrac <- NULL
  }
  pars$v2      <- pars$v1
  pars$tswitch <- if (use_switch) common$T_SWITCH else 0

  pars <- c(pars, FIXED_PARS, list(K = sum(samodel$x0)))

  simulate(samodel, pars, nyears)
}

# ------------------------------------------------------------
# 4. Loss functions: all have signature (sim, data, sd)
# ------------------------------------------------------------

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

      d_sq  <- (y1 - y2)^2
      c_vec <- y1 + y2

      idx <- d_sq != 0
      d_sq[idx] <- d_sq[idx] / c_vec[idx]

      total <- total + sum(weights[, j] * d_sq)
    }
    (1 / gamma) * sqrt(2 * nrow(sim)) * sqrt(total)
  }
)

# Turns a loss function into an objective of the parameter vector alone.
make_loss <- function(data, sd, loss_fn, use_switch,
                      sa_model_path = common$SA_MODEL_PATH,
                      nyears = NYEARS) {
  force(data); force(sd); force(loss_fn); force(use_switch)

  function(par) {
    names(par) <- PAR_NAMES
    sim <- run_model(as.list(par), use_switch, sa_model_path, nyears)
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
    if (!is.finite(v) || v < 0) return(penalty)
    log(max(v, .Machine$double.xmin))
  }
}

# ------------------------------------------------------------
# 5. Starting points
# ------------------------------------------------------------

# Returns a named list of starting vectors: "default" plus `n_random` extras.
#   type = "jitter":  default * exp(N(0, jitter_sd)), i.e. multiplicative noise
#   type = "uniform": uniform within start_lower / start_upper (tighten these
#                     if the bounds are very wide, e.g. sigma up to 1000)
make_starts <- function(n_random = 0,
                        type = c("jitter", "uniform"),
                        jitter_sd = 0.5,
                        lower = PAR_LOWER,
                        upper = PAR_UPPER,
                        start_lower = lower,
                        start_upper = upper) {
  type <- match.arg(type)
  stopifnot(
    identical(names(lower), PAR_NAMES),
    identical(names(upper), PAR_NAMES)
  )

  base <- clamp_to_bounds(unlist(DEFAULT_PARS)[PAR_NAMES], lower, upper)
  stopifnot(!anyNA(base))

  starts <- list(default = base)
  for (i in seq_len(n_random)) {
    s <- switch(type,
      jitter  = base * exp(rnorm(length(base), 0, jitter_sd)),
      uniform = runif(length(base), start_lower, start_upper)
    )
    names(s) <- PAR_NAMES
    starts[[sprintf("start_%02d", i)]] <- clamp_to_bounds(s, lower, upper)
  }
  starts
}

# ------------------------------------------------------------
# 6. Optimizers
# ------------------------------------------------------------

# Add new methods here; run_optimizer() picks them up by name.
OPTIMIZERS <- list(
  nmkb = dfoptim::nmkb,
  hjkb = dfoptim::hjkb
)

# Run one optimizer and standardize its output. `value` is always the
# loss on the ORIGINAL scale (re-evaluated at the optimum), so fits are
# comparable whether or not use_log was used.
run_optimizer <- function(method, fn, initial, lower, upper,
                          use_log = TRUE, maxfeval = 20000) {
  optimizer <- OPTIMIZERS[[method]]
  if (is.null(optimizer)) stop("Unknown method: ", method)

  fn_used <- if (use_log) log_loss(fn) else fn

  r <- optimizer(
    par = initial, fn = fn_used, lower = lower, upper = upper,
    control = list(maxfeval = maxfeval)
  )

  list(
    par         = r$par,
    value       = fn(r$par),
    convergence = r$convergence,
    n_feval     = if (is.null(r$feval)) NA_integer_ else r$feval,
    message     = if (is.null(r$message)) NA_character_ else r$message
  )
}

# One fit = one (loss, method, start). Never throws: errors are recorded.
fit_one <- function(fn, loss_name, method, start_id, initial,
                    lower = PAR_LOWER, upper = PAR_UPPER,
                    use_log = TRUE, maxfeval = 20000) {
  initial_loss <- tryCatch(fn(initial), error = function(e) NA_real_)

  elapsed <- system.time({
    result <- tryCatch(
      run_optimizer(method, fn, initial, lower, upper, use_log, maxfeval),
      error = function(e) list(
        par         = rep(NA_real_, length(PAR_NAMES)),
        value       = NA_real_,
        convergence = NA_integer_,
        n_feval     = NA_integer_,
        message     = conditionMessage(e)
      )
    )
  })[["elapsed"]]

  c(result, list(
    fit_id          = paste(loss_name, method, start_id, sep = "__"),
    loss            = loss_name,
    method          = method,
    start_id        = start_id,
    use_log         = use_log,
    initial_par     = initial,
    initial_loss    = initial_loss,
    runtime_seconds = elapsed
  ))
}

# ------------------------------------------------------------
# 7. Summaries and saving
# ------------------------------------------------------------

# Fitted parameters in the format the model uses downstream.
format_parameters <- function(par, use_switch) {
  p <- as.data.frame(as.list(setNames(par, PAR_NAMES)))

  p$wlow  <- p$whigh * p$wfrac
  p$wfrac <- NULL

  p$tswitch <- if (use_switch) common$T_SWITCH else 0
  p$v2      <- p$v1
  p
}

summarise_fits <- function(fits) {
  rows <- lapply(fits, function(f) {
    data.frame(
      fit_id          = f$fit_id,
      loss            = f$loss,
      method          = f$method,
      start_id        = f$start_id,
      use_log         = f$use_log,
      initial         = f$initial_loss,
      fitted          = f$value,
      improvement     = f$initial_loss - f$value,
      pct_change      = 100 * (f$initial_loss - f$value) / abs(f$initial_loss),
      convergence     = f$convergence,
      n_feval         = f$n_feval,
      runtime_seconds = f$runtime_seconds,
      message         = f$message,
      row.names       = NULL
    )
  })
  tbl <- do.call(rbind, rows)
  tbl <- tbl[order(tbl$loss, tbl$fitted), ]
  row.names(tbl) <- NULL
  tbl
}

parameter_table <- function(fits, use_switch) {
  rows <- lapply(fits, function(f) {
    if (!all(is.finite(f$par))) return(NULL)
    data.frame(
      fit_id   = f$fit_id,
      loss     = f$loss,
      method   = f$method,
      start_id = f$start_id,
      format_parameters(f$par, use_switch)
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

save_fit_results <- function(fits, use_switch, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  suffix <- if (use_switch) "_switch" else ""

  summary_tbl <- summarise_fits(fits)
  params_tbl  <- parameter_table(fits, use_switch)

  saveRDS(summary_tbl, file.path(output_dir, paste0("fit_summary", suffix, ".rds")))
  write.csv(summary_tbl, file.path(output_dir, paste0("fit_summary", suffix, ".csv")),
            row.names = FALSE)
  saveRDS(params_tbl, file.path(output_dir, paste0("all_fitted_parameters", suffix, ".rds")))
  write.csv(params_tbl, file.path(output_dir, paste0("all_fitted_parameters", suffix, ".csv")),
            row.names = FALSE)
  saveRDS(fits, file.path(output_dir, paste0("all_fit_objects", suffix, ".rds")))

  # One parameter file per loss x method (the best start), as before
  best <- best_per_group(summary_tbl)
  for (i in seq_len(nrow(best))) {
    pars <- params_tbl[params_tbl$fit_id == best$fit_id[i], ]
    saveRDS(pars, file.path(
      output_dir,
      paste0(make.names(best$loss[i]), "_", best$method[i], suffix, ".rds")
    ))
  }

  invisible(list(summary = summary_tbl, parameters = params_tbl, best = best))
}

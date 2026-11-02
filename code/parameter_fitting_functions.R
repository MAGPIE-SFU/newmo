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

PAR_NAMES <- c("v1", "sigma", "portion")

FIXED_MIGR <- 0.00001
FIXED_WLOW <- 0.0
FIXED_WHIGH <- 1

# These are the CC-magic parameters, hard coded here, so no need to save them to disk
DEFAULT_PARS <- common$DEFAULT_PARS

FIXED_PARS <- list(
  vaccine1 = common$VACCINE_TYPES[["PCV7"]],
  vaccine2 = common$VACCINE_TYPES[["PCV13"]]
)

# PAR_LOWER <- c(v1 = 0.001, sigma = 0,    migr = 0, portion = 0, whigh = 0.002, wfrac = 0)
# PAR_UPPER <- c(v1 = 0.99,  sigma = 1000, migr = 1, portion = 1, whigh = 1,     wfrac = 1)
PAR_LOWER <- c(v1 = 0.001, sigma = 0, portion = 0)
PAR_UPPER <- c(v1 = 1, sigma = 10000, portion = 1)

NYEARS <- common$NYEARS

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

# Returns list(data = <group x year matrix>, sd = <matching SD matrix>).
load_fit_data <- function(
  xt_path = common$XT_PATH,
  nyears = NYEARS,
  drop_groups = NULL
) {
  xt <- readr::read_csv(xt_path, show_col_types = FALSE)

  data <- pivot_to_matrix(xt, "frac")
  sd <- pivot_to_matrix(xt, "frac_sd")

  data <- data[, (ncol(data) - nyears + 1):ncol(data), drop = FALSE]
  sd <- sd[, (ncol(data) - nyears + 1):ncol(data), drop = FALSE]

  # Avoid zero SDs (they would give infinite weights)
  min_nonzero <- min(sd[sd > 0])
  sd[sd == 0] <- min_nonzero / 2

  if (!is.null(drop_groups)) {
    data <- data[!rownames(data) %in% drop_groups, , drop = FALSE]
    sd <- sd[!rownames(sd) %in% drop_groups, , drop = FALSE]
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
    if (is.null(cache[[path]])) {
      cache[[path]] <<- readRDS(path)
    }
    cache[[path]]
  }
})

# Simulation only: no data, no loss.
run_model <- function(
  model_pars = list(),
  use_switch,
  sa_model_path = common$SA_MODEL_PATH,
  nyears = NYEARS
) {
  samodel <- get_samodel(sa_model_path)

  pars <- as.list(model_pars)

  # Fixed Parameters are added:
  pars$whigh <- FIXED_WHIGH
  pars$wlow <- FIXED_WLOW
  pars$v2 <- pars$v1
  pars$migr <- FIXED_MIGR
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

      d_sq <- (y1 - y2)^2
      c_vec <- y1 + y2

      idx <- d_sq != 0
      d_sq[idx] <- d_sq[idx] / c_vec[idx]

      total <- total + sum(weights[, j] * d_sq)
    }
    (1 / gamma) * sqrt(2 * nrow(sim)) * sqrt(total)
  }
)

# Turns a loss function into an objective of the parameter vector alone.
make_loss <- function(
  data,
  sd,
  loss_fn,
  use_switch,
  sa_model_path = common$SA_MODEL_PATH,
  nyears = NYEARS
) {
  force(data)
  force(sd)
  force(loss_fn)
  force(use_switch)

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
lhs_unit <- function(n, k) {
  matrix(
    unlist(lapply(seq_len(k), function(j) (sample.int(n) - runif(n)) / n)),
    nrow = n
  )
}

make_starts <- function(
  n_random = 0,
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

  base <- clamp_to_bounds(unlist(DEFAULT_PARS)[PAR_NAMES], lower, upper)
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

OPTIMIZERS <- list(nmkb = dfoptim::nmkb, hjkb = dfoptim::hjkb, lbfgsb = lbfgsb)
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
# 7. Summaries and saving
# ------------------------------------------------------------

# Fitted parameters in the format the model uses downstream.
format_parameters <- function(par, use_switch) {
  p <- as.data.frame(as.list(setNames(par, PAR_NAMES)))

  # Should not be there, but keeping it for now.
  p$wfrac <- NULL

  p$whigh <- FIXED_WHIGH
  p$wlow <- FIXED_WLOW
  p$migr <- FIXED_MIGR

  p$tswitch <- if (use_switch) common$T_SWITCH else 0
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

parameter_table <- function(fits, use_switch) {
  rows <- lapply(fits, function(f) {
    if (!all(is.finite(f$par))) {
      return(NULL)
    }
    data.frame(
      fit_id = f$fit_id,
      loss = f$loss,
      method = f$method,
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
  params_tbl <- parameter_table(fits, use_switch)

  saveRDS(
    summary_tbl,
    file.path(output_dir, paste0("fit_summary", suffix, ".rds"))
  )
  write.csv(
    summary_tbl,
    file.path(output_dir, paste0("fit_summary", suffix, ".csv")),
    row.names = FALSE
  )
  saveRDS(
    params_tbl,
    file.path(output_dir, paste0("all_fitted_parameters", suffix, ".rds"))
  )
  write.csv(
    params_tbl,
    file.path(output_dir, paste0("all_fitted_parameters", suffix, ".csv")),
    row.names = FALSE
  )
  saveRDS(
    fits,
    file.path(output_dir, paste0("all_fit_objects", suffix, ".rds"))
  )

  # One parameter file per loss x method (the best start), as before
  best <- best_per_group(summary_tbl)
  for (i in seq_len(nrow(best))) {
    pars <- params_tbl[params_tbl$fit_id == best$fit_id[i], ]
    saveRDS(
      pars,
      file.path(
        output_dir,
        paste0(make.names(best$loss[i]), "_", best$method[i], suffix, ".rds")
      )
    )
  }

  invisible(list(summary = summary_tbl, parameters = params_tbl, best = best))
}

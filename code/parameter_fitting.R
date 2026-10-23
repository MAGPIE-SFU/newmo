# ============================================================
#  Fitting script
# ============================================================

library(deSolve)
library(readr)
library(dplyr)
library(tibble)
library(tidyr)
library(psych)
library(dfoptim)

common <- new.env()
sys.source(here::here("code", "common.R"), envir = common)

source(here::here("code", "simulate_ode.R"))

# Parameter bookkeeping
PAR_NAMES <- c("v1", "sigma", "migr", "portion", "whigh", "wfrac")

DEFAULT_PARS <- list(
    v1       = 0.065,
    # v2       = 0.065,
    sigma   = 0.13,
    migr    = 0.00001,
    portion = 0.3,
    whigh   = 0.1363,
    wfrac   = 0.0023 / 0.1363
    # tswitch = common$T_SWITCH
)

FIXED_PARS <- list(
    vaccine = common$VACCINE_TYPES[["PCV13"]]
)

NYEARS <- common$NYEARS

# Input data
xt <- readr::read_csv(common$XT_PATH, show_col_types = FALSE)

# Make a nice matrix from xt values:
XT_MATRIX <- xt %>%
    select(group, year, frac) %>%
    pivot_wider(
        names_from = year,
        values_from = frac,
        names_glue = "year.{year}"
    ) %>%
    column_to_rownames("group") %>%
    as.matrix()
XT_MATRIX <- XT_MATRIX[, (ncol(XT_MATRIX) - NYEARS + 1):ncol(XT_MATRIX), drop = FALSE]

# ------------------------------------
# DEBUG START
# ------------------------------------
# library(ggplot2)
# library(dplyr)
# library(tidyr)
# library(tibble)

# # Convert XT_MATRIX back to long format
# xt_matrix_long <- as.data.frame(XT_MATRIX) %>%
#     rownames_to_column("group") %>%
#     pivot_longer(
#         cols = -group,
#         names_to = "year",
#         values_to = "frac"
#     ) %>%
#     mutate(year = as.numeric(sub("year\\.", "", year)))

# # Plot original xt values and overlay XT_MATRIX values
# ggplot(xt, aes(x = year, y = frac)) +
#     geom_line(linewidth = 0.7, alpha = 0.7) +
#     geom_point(
#         data = xt_matrix_long,
#         aes(x = year, y = frac),
#         size = 2
#     ) +
#     facet_wrap(~ group, scales = "free_y") +
#     labs(
#         x = "Year",
#         y = "XT value"
#     ) +
#     theme_minimal()
# ------------------------------------
#  DEBUG END
# ------------------------------------


# Building the matrix of weights for the weighted SSE
XT_SD_MATRIX <- xt %>%
    select(group, year, frac_sd) %>%
    pivot_wider(
        names_from = year,
        values_from = frac_sd,
        names_glue = "year.{year}"
    ) %>%
    column_to_rownames("group") %>%
    as.matrix()
XT_SD_MATRIX <- XT_SD_MATRIX[, seq_len(NYEARS), drop=FALSE]
# Doing crazy adding to avoid 0 in SD matrix
min_nonzero <- min(XT_SD_MATRIX[XT_SD_MATRIX > 0])
XT_SD_MATRIX[XT_SD_MATRIX == 0] <- min_nonzero / 2

# ------------------------------------------------------------
# 2. Model object cache
# ------------------------------------------------------------

# Without this, readRDS() runs once per objective evaluation —
# thousands of deserialisations over a full fit.
get_samodel <- local({
    cache <- NULL
    function(path = common$SA_MODEL_PATH) {
        if (is.null(cache)) cache <<- readRDS(path)
        cache
    }
})

# ------------------------------------------------------------
# 3. Simulation only — no data, no loss
# ------------------------------------------------------------

run_model <- function(model_pars = list(),
                      sa_model_path = common$SA_MODEL_PATH,
                      nyears = NYEARS, use_switch) {

  samodel <- get_samodel(sa_model_path)

  pars <- DEFAULT_PARS
  pars[names(model_pars)] <- model_pars

  # wfrac = wlow / whigh, so wlow <= whigh always holds (no clamp, no dead zone).
  if ("wfrac" %in% names(pars)) {
    pars$wlow  <- pars$whigh * pars$wfrac
    pars$wfrac <- NULL
  }
  pars$v2 <- pars$v1
  pars$tswitch <- if (use_switch) common$T_SWITCH else 0

  pars <- c(pars, FIXED_PARS, list(K = sum(samodel$x0)))

  simulate(samodel, pars, nyears)
}

# ------------------------------------------------------------
# 4. Objective for optim()
# ------------------------------------------------------------

# ------------------------------------------------------------
# Loss functions
# ------------------------------------------------------------

loss_functions <- list(
    SSE = function(sim, data) {
        sum((sim - data)^2)
    },
    weighted_SSE = function(sim, data) {
        sum((1 / XT_SD_MATRIX^2) * (sim - data)^2)
    },
    weighted_composition_distance = function(sim, data) {
    gamma <- 1
    weights <- 1 / XT_SD_MATRIX^2

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

compare <- function(sim, data, loss_fn) {
    loss_fn(sim, data)
}

# ------------------------------------------------------------
# Loss factory
# ------------------------------------------------------------

make_loss <- function(data, loss_fn, use_switch, sa_model_path = common$SA_MODEL_PATH,
                      nyears = NYEARS) {
    force(data)

    function(par) {
        names(par) <- PAR_NAMES
        sim <- run_model(as.list(par), sa_model_path, nyears, use_switch)
        compare(sim, data, loss_fn)
    }
}

# ------------------------------------------------------------
#  DEBUG START
# ------------------------------------------------------------
# sim <- run_model(as.list(unlist(DEFAULT_PARS)[PAR_NAMES]))
# data <- XT_MATRIX

# for (name in names(loss_functions)) {
#     loss <- loss_functions[[name]](sim, XT_MATRIX)
#     cat(name, ":", loss, "\n")
# }
# ------------------------------------------------------------
#  DEBUG END
# ------------------------------------------------------------

# ------------------------------------------------------------
# 6. Fit
# ------------------------------------------------------------
# Compare all optimization methods, excluding multistart
USE_SWITCH <- F
set.seed(1337)

methods <- c("nlminb", "nmkb", "hjkb")

initial_raw <- unlist(DEFAULT_PARS)[PAR_NAMES]
lower <- c(v1 = 0.001, sigma = 0, migr = 0, portion = 0, whigh = 0.002, wfrac = 0)
upper <- c(v1 = 0.1, sigma = 100, migr = 1, portion = 1, whigh = 1, wfrac = 1)
initial <- pmin(pmax(initial_raw, lower + 1e-6), upper - 1e-6)

stopifnot(
  !anyNA(initial),
  identical(names(initial), names(lower)),
  identical(names(lower), names(upper))
)

# Output directory for comparison results
output_dir <- file.path(
  dirname(common$FITTED_PARAMETERS_PATH),
  "fit_comparison"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Run one optimizer and standardize its output
run_optimizer <- function(method, fn, initial, lower, upper) {
  switch(method,
    "nlminb" = {
      r <- nlminb(
        start = initial, objective = fn,
        lower = lower, upper = upper
      )
      list(par = r$par, value = r$objective,
           convergence = r$convergence, message = r$message)
    },
    "nmkb" = {
      r <- dfoptim::nmkb(
        par = initial, fn = fn, lower = lower, upper = upper,
        control = list(maxfeval = 20000)
      )
      list(par = r$par, value = r$value,
           convergence = r$convergence, message = NA_character_)
    },
    "hjkb" = {
      r <- dfoptim::hjkb(
        par = initial, fn = fn, lower = lower, upper = upper,
        control = list(maxfeval = 20000)
      )
      list(par = r$par, value = r$value,
           convergence = r$convergence, message = NA_character_)
    })
}

# Transform fitted parameters into the format used by the model
format_parameters <- function(par) {
  p <- as.data.frame(as.list(setNames(par, PAR_NAMES)))

  if ("wfrac" %in% names(p)) {
    p$wlow <- p$whigh * p$wfrac
    p$wfrac <- NULL
  }

  p$tswitch <- if (USE_SWITCH) common$T_SWITCH else 0
  p$v2 <- p$v1
  p
}

fits <- list()
summary_rows <- list()
parameter_rows <- list()

for (loss_name in names(loss_functions)) {
  fn <- make_loss(XT_MATRIX, loss_functions[[loss_name]], USE_SWITCH)
  initial_loss <- fn(initial_raw)

  for (method_name in methods) {
    fit_id <- paste(loss_name, method_name, sep = "__")
    message("Fitting: ", fit_id)

    elapsed <- system.time({
      result <- tryCatch(
        run_optimizer(method_name, fn, initial, lower, upper),
        error = function(e) list(
            par = rep(NA_real_, length(PAR_NAMES)),
            value = NA_real_,
            convergence = NA_integer_,
            message = conditionMessage(e)
        )
        )
    })[["elapsed"]]

    fits[[fit_id]] <- result
    fits[[fit_id]]$method <- method_name
    fits[[fit_id]]$loss <- loss_name
    fits[[fit_id]]$initial <- initial_loss
    fits[[fit_id]]$runtime_seconds <- elapsed

    summary_rows[[fit_id]] <- data.frame(
      loss = loss_name,
      method = method_name,
      initial = initial_loss,
      fitted = result$value,
      improvement = initial_loss - result$value,
      pct_change = 100 * (initial_loss - result$value) /
        abs(initial_loss),
      convergence = result$convergence,
      runtime_seconds = elapsed,
      message = if (is.null(result$message)) NA_character_ else
        result$message,
      row.names = NULL
    )

    # Save this fit's parameters separately
    if (all(is.finite(result$par))) {
      pars <- format_parameters(result$par)
      pars$loss <- loss_name
      pars$method <- method_name
      parameter_rows[[fit_id]] <- pars

      saveRDS(
        pars,
        file.path(
          output_dir,
          paste0(make.names(loss_name), "__", method_name, ".rds")
        )
      )
    }
  }
}

# Combined summary, sorted by loss then best fitted objective
summary_table <- do.call(rbind, summary_rows)
summary_table <- summary_table[
  order(summary_table$loss, summary_table$fitted),
]
row.names(summary_table) <- NULL

print(summary_table, row.names = FALSE)

# Combined parameter table
parameter_table <- do.call(rbind, parameter_rows)
row.names(parameter_table) <- NULL
parameter_table <- parameter_table[
  , c("loss", "method", setdiff(names(parameter_table), c("loss", "method")))
]

print(parameter_table, row.names = FALSE)

# Save combined outputs and the full optimizer results
saveRDS(summary_table, file.path(output_dir, "fit_summary.rds"))
write.csv(
  summary_table,
  file.path(output_dir, "fit_summary.csv"),
  row.names = FALSE
)
saveRDS(parameter_table, file.path(output_dir, "all_fitted_parameters.rds"))
write.csv(
  parameter_table,
  file.path(output_dir, "all_fitted_parameters.csv"),
  row.names = FALSE
)
saveRDS(fits, file.path(output_dir, "all_fit_objects.rds"))

# Also save the original CC_MAGIC parameters for comparison
cc_magic <- format_parameters(initial_raw)
cc_magic$loss <- "CC_MAGIC"
cc_magic$method <- "initial"
saveRDS(cc_magic, file.path(output_dir, "CC_MAGIC.rds"))

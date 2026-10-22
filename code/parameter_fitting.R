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
PAR_NAMES <- c("v", "sigma", "migr", "portion", "whigh", "wfrac")

DEFAULT_PARS <- list(
    v       = 0.065,
    sigma   = 0.13,
    migr    = 0.00001,
    portion = 0.3,
    whigh   = 0.1363,
    wfrac   = 0.0023 / 0.1363
)

FIXED_PARS <- list(
    vaccine = common$VACCINE_TYPES[["PCV13"]]
)

NYEARS <- 6

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
XT_SD_MATRIX <- XT_SD_MATRIX + min_nonzero / 2

# ------------------------------------------------------------
# 2. Model object cache
# ------------------------------------------------------------

# Without this, readRDS() runs once per objective evaluation —
# thousands of deserialisations over a full fit.
get_samodel <- local({
    cache <- NULL
    function(path = SA_MODEL_PATH) {
        if (is.null(cache)) cache <<- readRDS(path)
        cache
    }
})

# ------------------------------------------------------------
# 3. Simulation only — no data, no loss
# ------------------------------------------------------------

run_model <- function(model_pars = list(),
                      sa_model_path = SA_MODEL_PATH,
                      nyears = NYEARS) {

  samodel <- get_samodel(sa_model_path)

  pars <- DEFAULT_PARS
  pars[names(model_pars)] <- model_pars

  # wfrac = wlow / whigh, so wlow <= whigh always holds (no clamp, no dead zone).
  if ("wfrac" %in% names(pars)) {
    pars$wlow  <- pars$whigh * pars$wfrac
    pars$wfrac <- NULL
  }

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

make_loss <- function(data, loss_fn, sa_model_path = SA_MODEL_PATH,
                      nyears = NYEARS) {
    force(data)

    function(par) {
        names(par) <- PAR_NAMES
        sim <- run_model(as.list(par), sa_model_path, nyears)
        compare(sim, data, loss_fn)
    }
}

# ------------------------------------------------------------
#  DEBUG START
# ------------------------------------------------------------
sim <- run_model(as.list(unlist(DEFAULT_PARS)[PAR_NAMES]))
data <- XT_MATRIX

for (name in names(loss_functions)) {
    loss <- loss_functions[[name]](sim, XT_MATRIX)
    cat(name, ":", loss, "\n")
}
# ------------------------------------------------------------
#  DEBUG END
# ------------------------------------------------------------

# ------------------------------------------------------------
# 6. Fit
# ------------------------------------------------------------

initial <- unlist(DEFAULT_PARS)[PAR_NAMES]
lower <- c(v = 0.001, sigma = 0, migr = 0, portion = 0, whigh = 0.002, wfrac = 0)
upper <- c(v = 1,     sigma = 100, migr = 200, portion = 1, whigh = 1,     wfrac = 1)
stopifnot(!anyNA(initial), identical(names(initial), names(lower)))

initial_raw <- unlist(DEFAULT_PARS)[PAR_NAMES]  # values achieved my CC magic

initial <- pmin(pmax(initial, lower + 1e-6), upper - 1e-6)

set.seed(1337)

# "L-BFGS-B", "nlminb", "nmkb", "hjkb", "multistart"
method <- "nlminb"

fits <- lapply(loss_functions, function(loss_fn) {
  fn <- make_loss(XT_MATRIX, loss_fn)

  switch(method,
    "L-BFGS-B" = optim(initial, fn, method = "L-BFGS-B", lower = lower, upper = upper),

    "nlminb" = {
      r <- nlminb(initial, fn, lower = lower, upper = upper)
      list(par = r$par, value = r$objective, convergence = r$convergence, message = r$message)
    },

    "nmkb" = {
      r <- dfoptim::nmkb(initial, fn, lower = lower, upper = upper,
                         control = list(maxfeval = 5000))
      list(par = r$par, value = r$value, convergence = r$convergence)
    },

    "hjkb" = {
      r <- dfoptim::hjkb(initial, fn, lower = lower, upper = upper,
                         control = list(maxfeval = 5000))
      list(par = r$par, value = r$value, convergence = r$convergence)
    },

    # global-ish: nlminb from 20 random starts inside the bounds, keep the best
    "multistart" = {
      starts <- c(list(initial),
                  lapply(1:19, function(i) {
                    s <- lower + runif(length(lower)) * (upper - lower)
                    pmin(pmax(s, lower + 1e-6), upper - 1e-6)
                  }))
      runs <- lapply(starts, function(s) nlminb(s, fn, lower = lower, upper = upper))
      r <- runs[[which.min(sapply(runs, `[[`, "objective"))]]
      list(par = r$par, value = r$objective, convergence = r$convergence)
    }
  )
})

# computing loss for CC magic values
initial_values <- sapply(loss_functions, function(loss_fn) {
  make_loss(XT_MATRIX, loss_fn)(initial_raw)
})

# Compare the fitted parameters
results <- data.frame(
  loss        = names(fits),
  initial     = initial_values[names(fits)],
  fitted      = sapply(fits, `[[`, "value"),
  convergence = sapply(fits, `[[`, "convergence")
)
results$improvement <- results$initial - results$fitted
results$pct_change  <- 100 * results$improvement / abs(results$initial)

print(results, row.names = FALSE)

parameters <- as.data.frame(
    do.call(rbind, lapply(fits, `[[`, "par"))
)

parameters$loss <- names(fits)
parameters <- rbind(parameters,
                    data.frame(as.list(initial_raw), loss = "CC_MAGIC"))

# Convert wfrac to wlow for the saved parameter sets
if ("wfrac" %in% names(parameters)) {
    parameters$wlow <- parameters$whigh * parameters$wfrac
    parameters$wfrac <- NULL
}

print(parameters, row.names=FALSE)

# This should be looked at but currently not necessary to save it...
# saveRDS(fits, "data/SouthAfrica/parameter_fitting_results.rds")

saveRDS(parameters, common$FITTED_PARAMETERS_PATH)

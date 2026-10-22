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

source("code/simulate_ode.R")   # provides simulate()

# ------------------------------------------------------------
# 1. Parameter bookkeeping
# ------------------------------------------------------------

# Order is load-bearing: optim() passes par as an unnamed vector,
# so this is what puts the names back.
PAR_NAMES <- c("v", "sigma", "migr", "portion", "whigh", "wfrac")

DEFAULT_PARS <- list(
    v       = 0.065,
    sigma   = 0.13,
    migr    = 0.00001,
    portion = 0.3,
    whigh   = 0.1363,
    wfrac = 0.0023 / 0.1363
    # wlow    = 0.0023
    # delta   = 0.134 # whigh - wlow
)

FIXED_PARS <- list(
    vaccine = c("1","3","4","5","6A","6B","7F","9V","14","18C","19A","19F","23F")
)

SA_MODEL_PATH <- "data/SouthAfrica/samodel.rds"
NYEARS        <- 6

# Imputed input data to fit against
xt <- read_csv("data/SouthAfrica/processed-inputs/xt_grouped.csv", show_col_types = FALSE)
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
XT_MATRIX <- XT_MATRIX[, seq_len(NYEARS), drop = FALSE]

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
    # This does not work becuase division by 0...
    # kl_distance = function (sim, data) {
    #   return(sum((sim - data) * log(sim/data)))
    # },
    # This does not work because we have 0s as the data_gm and run into problems [division by 0]
    # log_ratio_transform = function(sim, data) {
    #     sim_gm  <- apply(sim, 2, psych::geometric.mean)
    #     data_gm <- apply(data, 2, psych::geometric.mean)

    #     sqrt(sum(
    #         (log(sim / sim_gm) - log(data / data_gm))^2
    #     ))
    # }
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
stopifnot(!anyNA(initial), identical(names(initial), names(lower)))
lower <- c(v = 0.001, sigma = 0, migr = 0, portion = 0, whigh = 0.002, wfrac = 0)
upper <- c(v = 1,     sigma = 100, migr = 200, portion = 1, whigh = 1,     wfrac = 1)

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

print(parameters, row.names=FALSE)

saveRDS(fits, "data/SouthAfrica/parameter_fitting_results.rds")

# ------------------------------------------------
# 3. Plot the 12 most abundant serotypes
# ------------------------------------------------

X <- as.matrix(XT_MATRIX)
big <- order(rowMeans(X), decreasing = TRUE)[1:min(12, nrow(X))]

op <- par(mfrow = c(5, 3), mar = c(3, 3, 2, 1))

for (i in big) {
  yr <- seq_len(ncol(X))
  offset <- 0.08

  ylim <- range(X[i, ], sim[i, ], sim_c[i, ], finite = TRUE)

  plot(yr, X[i, ], type = "n",
       main = rownames(X)[i],
       xlab = "Time index", ylab = "Abundance",
       ylim = ylim,
       xlim = c(0.7, ncol(X) + 0.3),
       cex.main = 0.8)

  # Observed data
  lines(yr - offset, X[i, ], col = "black", lwd = 2)
  points(yr - offset, X[i, ], col = "black", pch = 16)

  # Simulated data
  lines(yr, sim[i, ], col = "blue", lwd = 2, lty = 2)
  points(yr, sim[i, ], col = "blue", pch = 17, cex = 1.1)

  # Alternative simulation
  lines(yr + offset, sim_c[i, ], col = "red", lwd = 2, lty = 3)
  points(yr + offset, sim_c[i, ], col = "red", pch = 15, cex = 1.1)
}

# Empty panel for the shared legend
# Empty panel for a cleaner shared legend
plot.new()
legend("center",
       legend = c("Observed data", "sim", "sim_c"),
       col = c("black", "blue", "red"),
       lty = c(1, 2, 3),
       pch = c(16, 17, 15),
       lwd = 2,
       horiz = FALSE,
       bty = "n",
       cex = 1.1,
       x.intersp = 1.2,
       seg.len = 2)
par(op)

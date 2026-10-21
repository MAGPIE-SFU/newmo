# ============================================================
#  Fitting script
# ============================================================

library(deSolve)
library(readr)
library(dplyr)
library(tibble)
library(tidyr)
library(psych)

source("code/simulate_ode.R")   # provides simulate()

# ------------------------------------------------------------
# 1. Parameter bookkeeping
# ------------------------------------------------------------

# Order is load-bearing: optim() passes par as an unnamed vector,
# so this is what puts the names back.
PAR_NAMES <- c("v", "sigma", "migr", "portion", "whigh", "delta")

DEFAULT_PARS <- list(
    v       = 0.065,
    sigma   = 0.13,
    migr    = 0.00001,
    portion = 0.3,
    whigh   = 0.1363,
    # wlow    = 0.0023
    delta   = 0.134 # whigh - wlow
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

    if ("delta" %in% names(pars) && !"wlow" %in% names(pars)) {
      pars$wlow <- max(0, pars$whigh - pars$delta)
    }

    if ("delta" %in% names(pars)) {
      pars$wlow <- max(0, pars$whigh - pars$delta)
      pars$delta <- NULL
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
lower   <- c(v = 0.001, sigma = 0, migr = 0, portion = 0, whigh = 0.002, delta = 0.001)
upper   <- c(v = 1, sigma = 1, migr = 2, portion = 1, whigh = 1, delta = 0.998)

fits <- lapply(loss_functions, function(loss_fn) {
    optim(
        par    = initial,
        fn     = make_loss(XT_MATRIX, loss_fn),
        method = "L-BFGS-B",
        lower  = lower,
        upper  = upper
    )
})


# Compare the fitted parameters
results <- data.frame(
    loss       = names(fits),
    value      = sapply(fits, `[[`, "value"),
    convergence = sapply(fits, `[[`, "convergence")
)

print(results)

parameters <- as.data.frame(
    do.call(rbind, lapply(fits, `[[`, "par"))
)

parameters$loss <- names(fits)

print(parameters)

saveRDS(fits, "data/SouthAfrica/parameter_fitting_results.rds")

library(tidyverse)

source("simulate_ode.R")   # defines simulate(); also loads `samodel` via readRDS()

# -----------------------------
# Base parameters
# -----------------------------

modelpars <- list()

modelpars$v <- 0.061
modelpars$sigma <- 0.18
modelpars$migr <- 0.001
modelpars$portion <- 0.3

modelpars$vaccine <- c(
    "1", "3", "4", "5", "6A",
    "6B", "7F", "9V",
    "14", "18C", "19F", "23F"
)

modelpars$whigh <- 0.1363
# modelpars$wlow <- 0.0023
modelpars$wlow <- 0
modelpars$K <- sum(samodel$x0)

# match simulate_ode.R defaults: 6 years, yearly time steps
tspan <- seq(0, 72, by = 12)   # 0, 12, 24, 36, 48, 60, 72

# -----------------------------
# Parameter settings to test
# -----------------------------

parameter_grid <- expand_grid(
    v = c(0.03, 0.061, 0.1),
    sigma = c(0.1, 0.18, 0.3),
    migr = c(0.0001, 0.001, 0.01)
)

# -----------------------------
# Run simulations
# -----------------------------

all_simulations <- list()

for (i in seq_len(nrow(parameter_grid))) {

    cat("Running simulation", i, "\n")

    pars <- modelpars
    pars$v     <- parameter_grid$v[i]
    pars$sigma <- parameter_grid$sigma[i]
    pars$migr  <- parameter_grid$migr[i]

    sim <- simulate(samodel, pars, tspan)   # rows = strains, cols = "month.<t>"

    sim_df <- sim |>
        rownames_to_column("serotype") |>
        pivot_longer(
            cols = starts_with("month"),
            names_to = "month",
            values_to = "value"
        ) |>
        mutate(month = as.numeric(sub("^month\\.", "", month)))

    sim_df$simulation_id <- i
    sim_df$v     <- pars$v
    sim_df$sigma <- pars$sigma
    sim_df$migr  <- pars$migr

    all_simulations[[i]] <- sim_df
}

final_df <- bind_rows(all_simulations)

write.csv(
    final_df,
    "../data/SouthAfrica/parameter_exploration.csv",
    row.names = FALSE
)
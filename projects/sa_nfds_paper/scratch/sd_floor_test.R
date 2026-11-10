
Sys.setenv(TAR_PROJECT = "sa_nfds_paper")
pkgload::load_all("nfdsmodelling")
source("projects/sa_nfds_paper/config.R")
library(dplyr)

xt      <- targets::tar_read(xt_grouped_uncertainty)
samodel <- targets::tar_read(samodel)

# Look at the sd values first, to choose a floor
quantile(load_fit_data(xt, NYEARS)$sd, c(0, 0.01, 0.05, 0.1, 0.25, 0.5))

sd_floor <- 0.004   # pick from the quantiles above
xt_floor <- xt |> mutate(frac_sd = pmax(frac_sd, sd_floor))

fit <- run_sequential_fits(
  samodel, xt_floor, seed = 42, default_pars = DEFAULT_PARS,
  tswitch = T_SWITCH, nyears = NYEARS,
  losses = "weighted_SSE", n_screen = 2000, n_rounds = 3
)
fit$parameters      # does sigma stay reasonable now?

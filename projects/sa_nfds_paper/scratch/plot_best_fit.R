# TEMPORARY: plots the overall best fit from the pipeline with the old plot functions.
# Run from the repo root after tar_make().

source("code/plot_functions.R")

library(tidyr)
library(tibble)

Sys.setenv(TAR_PROJECT = "sa_nfds_paper")
pkgload::load_all("nfdsmodelling", quiet = TRUE)
source("projects/sa_nfds_paper/config.R")

# Load information from the pipeline
samodel  <- targets::tar_read(samodel)
xt       <- targets::tar_read(xt_grouped_uncertainty)
best_fit <- targets::tar_read(fit_selection)$overall_best

nyears <- NYEARS

# Same parameters as the fits use: fitted v1, sigma, portion + fixed values, tswitch, PCV7 -> PCV13
simulations <- list(
  "Best Parameter Fit" = run_model(
    as.list(unlist(best_fit[PAR_NAMES])),
    use_switch = FIT_USE_SWITCH, samodel = samodel,
    nyears = nyears, tswitch = T_SWITCH
  )
)

# plot the top groups
top_group_plot(
  vax_start = 2010, prevax_years = 2000:2009,
  n_groups = 20,
  samodel,
  xt,
  simulations,
  drop_vt = T,
  vaccine = VACCINE_TYPES[["PCV13"]]
)

# plot the top serotypes
top_sero_plot(vax_start = 2010, prevax_years = 2000:2009,
              n_serotypes = 18, samodel, xt, simulations,
              drop_vt = F,
              vaccine = VACCINE_TYPES[["PCV13"]])

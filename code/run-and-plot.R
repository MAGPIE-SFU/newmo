source("code/plot_functions.R")
source("code/simulate_ode.R")
source("code/parameter_generator.R")

library(tidyr)
library(tibble)

common <- new.env()
sys.source(here::here("code", "common.R"), envir = common)

# Load information
samodel <- readRDS(common$SA_MODEL_PATH)
xt <- read_csv(
    common$XT_PATH,
    show_col_types = FALSE
)

nyears <- common$NYEARS

# seq_seed <- 42
# rds_params_fit_path <- file.path(here::here("data", "SouthAfrica", "fit_sequential", paste0("seed_", seq_seed)), "all_fitted_parameters_switch.rds")

# rds_params_fit_path <- common$ALL_BEST_FITS_PATH
# rds_params_fit_path <- common$OVERALL_BEST_FITS_PATH
# all_fits <- readRDS(rds_params_fit_path)

# Optional filter, since 3 losses x 3 starts x 2 optimizers is a lot of curves.
# Examples: one optimizer, or one loss function:
# all_fits <- all_fits[all_fits$method == "hjkb", ]
# all_fits <- all_fits[all_fits$loss == "weighted_composition_distance", ]
# all_fits <- all_fits[all_fits$start_id == "round_06", ]

# modelpars_list <- setNames(
#   Map(function(id, loss, method) {
#     parameter_generator(
#       how_numerical = loss, optimizer = method, fit_id = id,
#       path = rds_params_fit_path,
#       vaccine_type1 = "PCV7", vaccine_type2 = "PCV13",
#       how_gene_order = "default",
#       use_switch = TRUE
#     )
#   }, all_fits$fit_id, all_fits$loss, all_fits$method),
#   all_fits$fit_id
# )

best_fit <- readRDS(common$OVERALL_BEST_FITS_PATH)
modelpars_list <- setNames(list(load_best_fit_pars()), "Best Paramter Fit")


# keep the baseline for comparison
# modelpars_list[["CC_MAGIC"]] <- parameter_generator(
#   how_numerical = "CC_MAGIC", use_switch = TRUE
# )

# Inspect parameters for each method
# str(modelpars_list)
# modelpars_list[["weighted_SSE"]]  # example: inspect one method

# Run simulations using the stored parameters
simulations <- lapply(modelpars_list, function(modelpars) {
  simulate(samodel, modelpars, nyears)
})

# plot the top groups
top_group_plot(
  vax_start=2010, prevax_years = 2000:2009,
  n_groups=20,
  samodel,
  xt,
  simulations,
  drop_vt = F,
  vaccine=common$VACCINE_TYPES[["PCV13"]]
)

# plot the top serotypes
top_sero_plot(vax_start=2010, prevax_years = 2000:2009,
              n_serotypes=18, samodel, xt, simulations,
              drop_vt = F,
              vaccine=common$VACCINE_TYPES[["PCV13"]])

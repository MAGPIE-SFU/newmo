library(here)

# Paths
SA_MODEL_PATH <- here::here("data", "SouthAfrica", "samodel.rds")

XT_PATH <- here::here(
    "data", "SouthAfrica", "processed-inputs", "xt_grouped_uncertainty.csv"
)

FIT_DIR            <- here::here("data", "SouthAfrica", "fit_comparison")
FIT_SEQUENTIAL_DIR <- here::here("data", "SouthAfrica", "fit_sequential")

# one row per loss x method (written by save_fit_results, see my earlier snippet)
FITTED_PARAMETERS_PATH        <- file.path(FIT_DIR, "all_fitted_parameters.rds")
FITTED_PARAMETERS_PATH_SWITCH <- file.path(FIT_DIR, "all_fitted_parameters_switch.rds")
FITTED_PARAMETERS_PATH_SEQ <- file.path(FIT_SEQUENTIAL_DIR, "all_fitted_parameters_switch.rds")

ALL_BEST_FITS_PATH <- file.path(here::here("data", "SouthAfrica"), "all_best_fitted_parameters.rds")
OVERALL_BEST_FITS_PATH <- file.path(here::here("data", "SouthAfrica"), "overall_best_fitted_parameters.rds")

# Default (CC_MAGIC) parameters: fit starting values and the baseline
DEFAULT_PARS <- list(
#   v1      = 0.065,
    v1      = -log(1 - 0.065),
#   sigma   = 0.13,
    sigma   = log(1 + 0.13),
  migr    = 0.00001,
  portion = 0.3,
  whigh   = 0.1363,
  wfrac   = 0.0023 / 0.1363
)

LOCUS_FILE_PATH <- here::here(
    "data", "SouthAfrica", "processed-inputs", "locus_freqs_deltas.csv"
)

FITTING_METHODS <- c(
    "CC_MAGIC",
    "weighted_composition_distance",
    "SSE",
    "weighted_SSE"
)

# Time point for the model with two vaccines
# 2 years between vaccine, PCV7 in 2009 and PCV13 in 2011
NYEARS         <- 6
UNITS_PER_YEAR <- 12

# Hard switch:
SWITCH_YEARS_AFTER_START <- 2
T_SWITCH <- SWITCH_YEARS_AFTER_START * UNITS_PER_YEAR

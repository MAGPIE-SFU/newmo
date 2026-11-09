library(here)

# Paths
SA_MODEL_PATH <- here::here("data", "SouthAfrica", "samodel.rds")

XT_PATH <- here::here(
    "data", "SouthAfrica", "processed-inputs", "xt_grouped_uncertainty.csv"
)

# Vaccine definitions
VACCINE_TYPES <- list(
  PCV7 = c(
    "4", "6B", "9V", "14", "18C", "19F", "23F"
  ),

  PCV13 = c(
    "1", "3", "4", "5", "6A", "6B", "7F", "9V", "14", "18C",
    "19A", "19F", "23F"
  ),

  PCV15 = c(
    "1", "3", "4", "5", "6A", "6B", "7F", "9V", "14", "18C",
    "19A", "19F", "22F", "23F", "33F"
  ),

  PCV21 = c(
    "3", "6A", "7F", "8", "9N", "10A", "11A", "12F", "15A", "15B",
    "16F", "17F", "19A", "20A", "22F", "23A", "23B", "24F", "31", "33F",
    "35B"
  ),

  PCV23 = c(
    "1", "2", "3", "4", "5", "6B", "7F", "8", "9N", "9V",
    "10A", "11A", "12F", "14", "15B", "17F", "18C", "19A", "19F", "20",
    "22F", "23F", "33F"
  )
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

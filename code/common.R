library(here)

# Paths
SA_MODEL_PATH <- here::here("data", "SouthAfrica", "samodel.rds")

XT_PATH <- here::here(
    "data", "SouthAfrica", "processed-inputs", "xt_grouped.csv"
)

# Vaccine definitions
VACCINE_TYPES <- list(
    PCV13 = c(
        "1", "3", "4", "5", "6A", "6B", "7F",
        "9V", "14", "18C", "19A", "19F", "23F"
    ),
    PCV7 = c(
        "4", "6B", "9V", "14", "18C", "19F", "23F"
    )
)

FITTED_PARAMETERS_PATH <- here::here(
    "data", "SouthAfrica", "fitted_parameters.rds"
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

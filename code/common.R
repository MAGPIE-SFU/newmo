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

# Cached model loader
get_samodel <- local({
    cache <- NULL

    function(path = SA_MODEL_PATH) {
        if (is.null(cache)) cache <<- readRDS(path)
        cache
    }
})

FITTED_PARAMETERS_PATH <- here::here(
    "data", "SouthAfrica", "fitted_parameters.rds"
)



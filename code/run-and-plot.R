library(deSolve)
library(readr)
library(tibble)
library(dplyr)
library(tidyr)
library(ggplot2)

source("code/plot_functions.R")
source("code/simulate_ode.R")
source("code/parameter_generator.R")

common <- new.env()
sys.source(here::here("code", "common.R"), envir = common)

# Load information
samodel <- readRDS(common$SA_MODEL_PATH)
xt <- read_csv(
    common$XT_PATH,
    show_col_types = FALSE
)

nyears <- common$NYEARS

simulations <- setNames(
    lapply(common$FITTING_METHODS, function(method) {
        modelpars <- parameter_generator(
            how_numerical = method,
            vaccine_type1 = "PCV7",
            vaccine_type2 = "PCV13",
            how_gene_order = "default",
            use_switch = FALSE
        )

        simulate(samodel, modelpars, nyears)
    }),
    common$FITTING_METHODS
)


# plot the top groups
top_group_plot(vax_start=2010, prevax_years = 2000:2009,
               n_groups=20, samodel, xt, simulations, drop_vt = FALSE,
              vaccine=common$VACCINE_TYPES[["PCV13"]])

# plot the top serotypes
top_sero_plot(vax_start=2010, prevax_years = 2000:2009,
              n_serotypes=18, samodel, xt, simulations, vaccine=common$VACCINE_TYPES[["PCV13"]])

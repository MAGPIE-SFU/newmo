library(deSolve)
library(readr)
library(tibble)
library(dplyr)
library(tidyr)
library(ggplot2)

source("code/plot_functions.R")
source("code/simulate_ode.R")
source("code/parameter_generator.R")

PCV13_VACCINE <- c("1", "3", "4", "5", "6A", "6B", "7F", "9V", "14", "18C", "19A", "19F", "23F")

# Load information
samodel <- readRDS("data/SouthAfrica/samodel.rds")
xt <- read_csv(
    "data/SouthAfrica/processed-inputs/xt_grouped.csv",
    show_col_types = FALSE
)

numerical_methods <- c(
    "CC_MAGIC",
    "weighted_composition_distance",
    "SSE",
    "weighted_SSE"
)

nyears <- 6

simulations <- setNames(
    lapply(numerical_methods, function(method) {
        modelpars <- parameter_generator(
            how_numerical = method,
            vaccine_type = "PCV13",
            how_gene_order = "sdmove"
        )

        simulate(samodel, modelpars, nyears)
    }),
    numerical_methods
)


# plot the top groups
top_group_plot(vax_start=2010, prevax_years = 2000:2009,
               n_groups=20, samodel, xt, simulations,drop_vt = FALSE,
              vaccine=PCV13_VACCINE)

# plot the top serotypes
top_sero_plot(vax_start=2010, prevax_years = 2000:2009,
              n_serotypes=18, samodel, xt, simulations, vaccine=PCV13_VACCINE)

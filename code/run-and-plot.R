library(deSolve)
library(readr)
library(tibble)
library(dplyr)
library(tidyr)
library(ggplot2)

source("code/plot_functions.R")
source("code/simulate_ode.R")
source("code/parameter_generator.R")

# Load information
samodel <- readRDS("data/SouthAfrica/samodel.rds")
xt <- read_csv(
    "data/SouthAfrica/processed-inputs/xt_grouped.csv",
    show_col_types = FALSE
)

# Set model parameters
modelpars <- parameter_generator(
    how_numerical = "CC_MAGIC",
    vaccine_type = "PCV13",
    how_gene_order = "default"
)

nyears <- 6

# run the model
mysim = simulate(samodel, modelpars, nyears)

# plot the top groups
top_group_plot(vax_start=2010, prevax_years = 2000:2009,
               n_groups=18, samodel, xt, mysim)

# plot the top serotypes
top_sero_plot(vax_start=2010, prevax_years = 2000:2009,
              n_serotypes=18, samodel, xt, mysim)


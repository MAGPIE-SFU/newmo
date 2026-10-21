library(deSolve)
library(readr)
library(tibble)
library(dplyr)

source("code/plot_functions.R") 
source("code/simulate_ode.R")

# load information 
samodel = readRDS("data/SouthAfrica/samodel.rds")
xt <- read_csv("data/SouthAfrica/xt_grouped.csv", show_col_types = FALSE)

# set model parameters
modelpars <-list()
modelpars$v <- 0.065 # there is a reference to 0.1 in my notes, too..
modelpars$sigma <- 0.13 # from previous NFDS models
modelpars$migr <- 0.00001 # migration parameter
modelpars$portion <- 0.3 # fraction under NFDS
modelpars$vaccine = c("1", "3", "4", "5", "6A","6B","7F","9V","14","18C","19A", "19F","23F") # PCV13 i think
modelpars$whigh = 0.1363 # from MA model and previous paper
modelpars$wlow = 0.0023 # from MA model too; could be 0 (Leonie)
modelpars$K=sum(samodel$x0)
nyears = 6 

# run the model 
mysim = simulate(samodel, modelpars, nyears)

# plot the top groups 
top_group_plot(vax_start=2010, prevax_years = 2000:2009,
               n_groups=18, samodel, xt, mysim)

# plot the top serotypes 
top_sero_plot(vax_start=2010, prevax_years = 2000:2009,
              n_serotypes=18, samodel, xt, mysim)




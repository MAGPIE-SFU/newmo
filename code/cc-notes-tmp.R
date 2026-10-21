# cc simulation notes 
library(deSolve)
library(readr)
library(tibble)
library(dplyr)

source("code/plot_functions.R") 
source("code/simulate_ode.R")
samodel = readRDS("data/SouthAfrica/samodel.rds")

xt <- read_csv("data/SouthAfrica/xt_grouped.csv", show_col_types = FALSE)

# model parameters 
 modelpars <-list()
 modelpars$v <- 0.065 # there is a reference to 0.1 in my notes, too..
 modelpars$sigma <- 0.13 # from previous NFDS models
 modelpars$migr <- 0.00001 # migration parameter
 modelpars$portion <- 0.3 # fraction under NFDS
 modelpars$vaccine = c("1", "3", "4", "5", "6A","6B","7F","9V","14","18C","19A", "19F","23F") # PCV13 i think
 modelpars$whigh = 0.1363 # from MA model and previous paper
 modelpars$wlow = 0.0023 # from MA model too; could be 0 (Leonie)
#
 modelpars$K=sum(samodel$x0) # ok but that's 1.. 
nyears = 6 
mp2 = modelpars; mp2$sigma = 0; 
# note that all 3 of sigma, portion and whigh (and wlow) give a strength of NFDS - should probably simplify 
# more wlow (0.07 only change) : 19A and 35B are worse, but rareNVT a bit better , 15A is worse (than baseline above) 
# more wlow (0.04) and lower sigma (0.1) 

 mysim = simulate(samodel, modelpars, nyears)
sim2 = simulate(samodel, mp2, nyears)

sims = list(baseline = mysim, noNFDS = sim2) 
 # -- check out these plots! see plot_functions.R 
 vt_nvt_plot(vax_start=2010, prevax_years = 2000:2009,
              samodel, xt, sims) # ha this works 
 # changing the fraction under NFDS doesn't change the vt totals; goo. 
 # well, it does it you put it all the way up to 0.5 with a 0.2 whigh, and the fit is not as good 
 
 top_group_plot(vax_start=2010, prevax_years = 2000:2009,
 n_groups=18, samodel, xt, sims)
 
 top_sero_plot(vax_start=2010, prevax_years = 2000:2009,
                n_serotypes=18, samodel, xt, sims)
 
 # at baseline we have too much 19A and rare NVT, and not enough 16F, 23A , 35B. 35B is very common in the data, so it's a bit problem 
 # huge nfds does not fix these things and worsens the VT to NVT match eventually though you need a lot of it 
 
 
 # -- some other checks 
 
 
 # is x0 actually an equilib of the model with v = 0 ? 
 mp0 <- modelpars; mp0$v <- 0; mp0$migr = 0.00001
 sim0 <- simulate(samodel, mp0, nyears = 6)
 sim0 <- sweep(sim0, 2, colSums(sim0), "/")
 tapply(sim0[, 1], samodel$st %in% modelpars$vaccine, sum)   # VT share, year 1, no vaccine
 tapply(sim0[, 6], samodel$st %in% modelpars$vaccine, sum)   # ... and year 6
 # yes unless migration is too high; therefore have set it very low . 
 
 fx0 <- as.numeric(t(samodel$G) %*% samodel$x0) / sum(samodel$x0)
 summary(fx0 - samodel$el)          # should be ~0 everywhere if x0 is an equilibrium
 
 
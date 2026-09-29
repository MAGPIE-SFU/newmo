# cc simulation notes 
library(deSolve)
library(readr)
# model parameters 
 modelpars <-list()
 modelpars$v <- 0.065 # there is a reference to 0.1 in my notes, too..
 modelpars$sigma <- 0.13 # from previous NFDS models
 modelpars$migr <- 0.00001 # migration parameter
 modelpars$portion <- 0.3 # fraction under NFDS
 modelpars$vaccine = c("1", "3", "4", "5", "6A","6B","7F","9V","14","18C","19F","23F") # PCV13 i think
 modelpars$whigh = 0.1363 # from MA model and previous paper
 modelpars$wlow = 0.0023 # from MA model too; could be 0 (Leonie)
#
 modelpars$K=sum(samodel$x0) # ok but that's 1.. 

 
 
 mysim = simulate(samodel, modelpars, nyears)

 
 # is x0 actually an equilib of the model with v = 0 ? 
 mp0 <- modelpars; mp0$v <- 0; mp0$migr = 0.00001
 sim0 <- simulate(samodel, mp0, nyears = 6)
 sim0 <- sweep(sim0, 2, colSums(sim0), "/")
 tapply(sim0[, 1], samodel$st %in% modelpars$vaccine, sum)   # VT share, year 1, no vaccine
 tapply(sim0[, 6], samodel$st %in% modelpars$vaccine, sum)   # ... and year 6
 # yes unless migration is too high; therefore have set it very low . 
 
 fx0 <- as.numeric(t(samodel$G) %*% samodel$x0) / sum(samodel$x0)
 summary(fx0 - samodel$el)          # should be ~0 everywhere if x0 is an equilibrium
 
# functions for the NFDS-based vaccine design 


getIPD <- function(simulation, inv) {} 

score_vaccine <- function(modelpars, vseros) {
  modelpars$vaccine1=vseros
  mysim = simulate(samodel, modelpars, nyears=6)
  ipd  <- get_ipd(mysim)                    # per group per year
  sum(ipd)                                # objective: cumulative IPD
}
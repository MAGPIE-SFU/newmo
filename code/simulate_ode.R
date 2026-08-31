library(deSolve)
library(readr)

# simulate ode version of the dynamics for SA
samodel <- readRDS("data/SouthAfrica/samodel.rds") # this is the samodel; has G, st, x0, el, delta

# model parameters 
# modelpars <-list()
# modelpars$v <- 0.061 # there is a reference to 0.1 in my notes, too..
# modelpars$sigma <- 0.18 # from previous NFDS models
# modelpars$migr <- 0.001 # migration parameter
# modelpars$portion <- 0.3 # fraction under NFDS
# modelpars$vaccine = c("1", "3", "4", "5", "6A","6B","7F","9V","14","18C","19F","23F") # the serotypes in the PCV7 vaccines
# modelpars$whigh = 0.1363 # from MA model and previous paper
# modelpars$wlow = 0.0023 # from MA model too; could be 0 (Leonie)
#
# modelpars$K=sum(samodel$x0)

# pcv 13 (1, 3, 4, 5, 6A, 6B, 7F, 9V, 14, 18C, 19A, 19F, and 23F) 
# note pcv7 is  c("4","6A","6B","9V","14","18C","19F","23F") 
tmax   <- 6 # years 
tspan  <- seq(0, tmax * 12, by = 12) # time steps are monthly 

simulate <- function(samodel, modelpars, tspan ) {
    # preliminary to running a simulation, but depends on model parameters: 
    
    # create logged versions to translate to continuous time
    rho <- log(1 + modelpars$sigma)
    
    # create vt: 0 or parameter v, one value for each strain (row of G) 
    vt <-  modelpars$v * (samodel$st %in% modelpars$vaccine)
    rr  <- -log(1 - vt)
    
    # create weight vector w (one value for each locus) 
    w <- rep(1, length(samodel$el))
    tmpind <- round(modelpars$portion*length(samodel$el))
    w[1:tmpind] <- modelpars$whigh
    w[(tmpind+1):length(w)] <- modelpars$wlow

    # create matrix multiples to simplify pi calculation 
    modelpars$Gwe <- samodel$G %*% (w * samodel$el) # this is sum_l Gil wl el , the first part of pi
                                        # called "weighteddev" in earlier models
    modifiedGT  <- t(samodel$G) *w # this is Gw
    modelpars$GtimesModG  <- as.matrix(samodel$G) %*% modifiedGT

    # define the differential equation
    diffxt <- function(t, x, modelpars) {
        with(modelpars, { 
            totalx <- sum(x)
        dx <- x * (
            log(K / totalx) -
                rr +
                rho * (Gwe - GtimesModG %*% (x / totalx))
        ) + migr
        list(as.numeric(dx))
        })
    }
    
    # now run the simulation 
    out <- ode(
        y = as.numeric(samodel$x0),
        times = tspan,
        func = diffxt,
        parms = modelpars,
        method = "lsoda",
        atol = 1e-8,
        rtol = 1e-5
    ) 
    # out's first column is time. Each row gives the strain numbers at that time. Weird. 
    
    # collect the results in a convenient shape
    outdf <- as.data.frame(t(out[, -1, drop = FALSE]))
    colnames(outdf) <- paste("month", out[, "time"] , sep=".")       # columns = time points
    rownames(outdf) <- rownames(samodel$G)          # rows = strain names
    return(outdf)
}

# mysim = simulate(samodel, modelpars, tspan)
#
# # modified dataframe saved to csv, path sensitive...
# library(tibble)
# mod_df <- rownames_to_column(mysim, var="id")
# write.csv(mod_df, "../data/SouthAfrica/simulation_data.csv", row.names=F)

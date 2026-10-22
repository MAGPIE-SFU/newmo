library(deSolve)
library(readr)

simulate <- function(samodel, modelpars, nyears, stepspy = 24 ) {
    # preliminary to running a simulation, but depends on model parameters:

    K <- sum(samodel$x0)

    # create logged versions to translate to continuous time
    rho <- log(1 + modelpars$sigma)

    # create weight vector w (one value for each locus)
    if (is.null(modelpars[["gene_order"]])) modelpars$gene_order = seq_len(length(samodel$el))

    w <- rep(modelpars$wlow, length(samodel$el))
    tmpind <- round(modelpars$portion*length(samodel$el))
    w[modelpars$gene_order[1:tmpind]] <- modelpars$whigh

    # create matrix multiples to simplify pi calculation
    modelpars$Gwe <- samodel$G %*% (w * samodel$el) # this is sum_l Gil wl el , the first part of pi
                                        # called "weighteddev" in earlier models
    modifiedGT  <- t(samodel$G) *w # this is Gw
    modelpars$GtimesModG  <- as.matrix(samodel$G) %*% modifiedGT

    # define the differential equation
    diffxt <- function(t, x, modelpars) {
        # create vt: 0 or parameter v, one value for each strain (row of G)

        # time dependent vaccines, where we switch from vaccine1 to vaccine2 at point tswitch
        vt <- if (t < modelpars$tswitch) {
            modelpars$v1 * (samodel$st %in% modelpars$vaccine1)
        } else {
            modelpars$v2 * (samodel$st %in% modelpars$vaccine2)
        }
        rr  <- -log(1 - vt)

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
    tfine <- seq(0, nyears * 12, length.out = nyears * stepspy + 1)
    # now run the simulation
    out <- ode(
        y = as.numeric(samodel$x0),
        times = tfine,
        func = diffxt,
        parms = modelpars,
        method = "lsoda",
        atol = 1e-8,
        rtol = 1e-5
    )
    # out's first column is time. Each row gives the strain numbers at that time. Weird.
    X <- out[, -1, drop = FALSE]

    # average over each year
    avg <- sapply(seq_len(nyears), function(k) {
      idx <- ((k - 1) * stepspy + 1):(k * stepspy + 1) # rows for year k
      XX <- X[idx, , drop = FALSE] # all strains, but only year k
      n <- nrow(XX)
      colSums(XX[-1, , drop = FALSE] + XX[-n, , drop = FALSE]) / (2 * (n - 1)) # trapezoid rule
      # which is (x_1/2 + x_2 + … + x_{n−1} + x_n/2) / (n−1) where the endpoints are shared with the other year
    })

    # collect the results in a convenient shape (sort of. though we reshape it again for plotting so..??)

    outdf <- as.data.frame(matrix(avg, ncol = nyears))
    rownames(outdf) <- rownames(samodel$G)             # rows = strains
    colnames(outdf) <- paste("year", seq_len(nyears), sep = ".")   # year k = months 12(k-1) to 12k
    return(outdf)
}

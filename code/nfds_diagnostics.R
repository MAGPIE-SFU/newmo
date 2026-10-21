# Diagnose which NVT groups the NFDS model favours after vaccination, and why.
# Needs: samodel (G, st, x0, el, delta), modelpars (sigma, portion, whigh, wlow, vaccine)
library(dplyr)
library(tibble)

# 0. the same weights the simulation uses
w <- rep(modelpars$wlow, length(samodel$el))

w[1:round(modelpars$portion * length(w))] <- modelpars$whigh
rho <- log(1 + modelpars$sigma)

# 1. locus "deficits" just after the vaccine removes VT
G   <- samodel$G
x0  <- as.numeric(samodel$x0)
vt  <- samodel$st %in% modelpars$vaccine
xpost <- ifelse(vt, 0, x0)
f_post <- as.numeric(t(G) %*% xpost) / sum(xpost)
d <- samodel$el - f_post      # > 0: locus now under-represented, so nfds pushes up

# 2. NFDS fitness of each NVT group in that post-vaccine population
fit <- tibble(group  = rownames(G),
              st     = samodel$st,
              vt     = vt,
              x0freq = x0 / sum(x0),
              fitness = rho * as.numeric(G %*% (w * d)),
              fit_high = rho * as.numeric(G %*% (w * d * (w == modelpars$whigh))),   # from high-weight loci
              fit_low  = rho * as.numeric(G %*% (w * d * (w == modelpars$wlow)))) |>
  filter(!vt) |>
  mutate(postfreq = x0freq / sum(x0freq),                          # frequency among NVT once VT are gone
         rel_fitness = fitness - sum(postfreq * fitness),           # frequency grows at this rate (per month)
         months_to_5pct = ifelse(rel_fitness > 0 & postfreq < 0.05,
                                 log(0.05 / postfreq) / rel_fitness, NA)) |>
  arrange(desc(rel_fitness))
print(fit, n = 15)

# 3. which loci make one group fitter than another
compare_groups <- function(a, b, top = 15) {
  tibble(locus = colnames(G),
         weight = w, deficit = d,
         in_a = G[a, ], in_b = G[b, ],
         contrib_diff = rho * (G[a, ] - G[b, ]) * w * d) |>
    filter(in_a != in_b) |>
    arrange(desc(abs(contrib_diff))) |>
    head(top)
}
print(compare_groups("19A_17", "35B_other"))

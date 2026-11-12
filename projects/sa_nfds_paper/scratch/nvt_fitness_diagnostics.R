# Diagnose which NVT groups the NFDS model favours right after vaccination, and why.
# Approximation: the NFDS selection term at the moment all VT groups are removed
# (ignores vaccination dynamics, the log(K/total) term and migration).
# Run from the repo root after tar_make().

library(dplyr)
library(tibble)

Sys.setenv(TAR_PROJECT = "sa_nfds_paper")
pkgload::load_all("nfdsmodelling", quiet = TRUE)

samodel  <- targets::tar_read(samodel)
best_fit <- targets::tar_read(fit_selection)$overall_best

# Same parameters as run_model() uses for the overall best fit
sigma   <- best_fit$sigma
portion <- best_fit$portion
whigh   <- FIXED_WHIGH
wlow    <- FIXED_WLOW
vaccine <- VACCINE_TYPES[["PCV13"]]   # population after the switch to PCV13

# 0. the same weights and rho the simulation uses
w <- rep(wlow, length(samodel$el))
w[seq_len(round(portion * length(w)))] <- whigh   # default gene order
rho <- sigma                                       # simulate_nfds() uses sigma directly

# 1. locus "deficits" just after the vaccine removes VT
G     <- samodel$G
x0    <- as.numeric(samodel$x0)
vt    <- samodel$st %in% vaccine
xpost <- ifelse(vt, 0, x0)
f_post <- as.numeric(t(G) %*% xpost) / sum(xpost)
d <- samodel$el - f_post      # > 0: locus now under-represented, so NFDS pushes up

# 2. NFDS fitness of each NVT group in that post-vaccine population
fit <- tibble(group    = rownames(G),
              st       = samodel$st,
              vt       = vt,
              x0freq   = x0 / sum(x0),
              fitness  = rho * as.numeric(G %*% (w * d)),
              fit_high = rho * as.numeric(G %*% (w * d * (w == whigh))),   # from high-weight loci
              fit_low  = rho * as.numeric(G %*% (w * d * (w == wlow)))) |> # from low-weight loci
  filter(!vt) |>
  mutate(postfreq       = x0freq / sum(x0freq),                # frequency among NVT once VT are gone
         rel_fitness    = fitness - sum(postfreq * fitness),   # frequency grows at this rate (per month)
         months_to_5pct = ifelse(rel_fitness > 0 & postfreq < 0.05,
                                 log(0.05 / postfreq) / rel_fitness, NA)) |>
  arrange(desc(rel_fitness))
print(fit, n = 15)

# 3. which loci make one group fitter than another
compare_groups <- function(a, b, top = 15) {
  missing <- setdiff(c(a, b), rownames(G))
  if (length(missing)) stop("Unknown group(s): ", paste(missing, collapse = ", "),
                            ". See rownames(samodel$G).")
  tibble(locus = colnames(G),
         weight = w, deficit = d,
         in_a = G[a, ], in_b = G[b, ],
         contrib_diff = rho * (G[a, ] - G[b, ]) * w * d) |>
    filter(in_a != in_b) |>
    arrange(desc(abs(contrib_diff))) |>
    head(top)
}

# e.g. the fittest NVT group against the least fit one
print(compare_groups(fit$group[1], fit$group[nrow(fit)]))

# ============================================================
#  sequential_fits.R
#  Chained refits: each round starts from the previous round's optimum.
#  Checks whether the fit is stable or keeps drifting / improving.
#  Uses the same functions as run_fits.R.
# ============================================================

source(here::here("code", "parameter_fitting_functions.R"))
source(here::here("code", "fit_report_functions.R"))

# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------
USE_SWITCH    <- TRUE
USE_LOG       <- TRUE
LOSSES        <- names(loss_functions)
CHAIN_METHODS <- c("nmkb", "lbfgsb", "hjkb")
N_ROUNDS      <- 6
START_FROM    <- "screen"
N_SCREEN      <- 1000
MAXFEVAL      <- 20000
SEED          <- 1337
REL_TOL       <- NULL                      # e.g. 1e-6 to stop a chain early when a round
                                           # improves the loss by less than this (relative)

# Kept separate so run_fits.R output is not overwritten
OUTPUT_DIR <- common$FIT_SEQUENTIAL_DIR

# ------------------------------------------------------------
# Setup
# ------------------------------------------------------------
set.seed(SEED)
fit_data <- load_fit_data()

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

if (identical(START_FROM, "screen")) {

  screen_starts <- make_starts(n_random = N_SCREEN)   # includes "default"

  # losses x starts: simulate once per start, apply every loss to that simulation
  t_screen <- system.time(
    screen <- sapply(screen_starts, function(s) {
      sim <- tryCatch(run_model(as.list(s), USE_SWITCH), error = function(e) NULL)
      vapply(loss_functions, function(f) {
        if (is.null(sim)) return(Inf)
        v <- tryCatch(f(sim, fit_data$data, fit_data$sd), error = function(e) Inf)
        if (is.finite(v)) v else Inf
      }, numeric(1))
    })
  )[["elapsed"]]
  message("Screened ", length(screen_starts), " starts in ", round(t_screen), " s")

  saveRDS(screen, file.path(OUTPUT_DIR, "screen_losses.rds"))
  best_idx <- apply(screen, 1, which.min)    # best start per loss function

} else {
  start <- if (is.numeric(START_FROM)) {
    clamp_to_bounds(START_FROM[PAR_NAMES])
  } else {
    starts <- make_starts(n_random = as.integer(START_FROM == "random"))
    starts[[length(starts)]]
  }
}

# ------------------------------------------------------------
# Chains: round r starts where round r-1 ended
# ------------------------------------------------------------
fits <- list()

for (loss_name in LOSSES) {
  fn <- make_loss(fit_data$data, fit_data$sd, loss_functions[[loss_name]], USE_SWITCH)
  current <- if (identical(START_FROM, "screen")) {
      message("  best screened start: ", names(screen_starts)[best_idx[[loss_name]]],
              " (loss ", signif(screen[loss_name, best_idx[[loss_name]]], 4),
              "; default start: ", signif(screen[loss_name, "default"], 4), ")")
      screen_starts[[best_idx[[loss_name]]]]
    } else start

  for (r in seq_len(N_ROUNDS)) {
    method <- CHAIN_METHODS[(r - 1) %% length(CHAIN_METHODS) + 1]
    message("Fitting: ", loss_name, " / round ", r, " / ", method)

    fit <- fit_one(
      fn, loss_name, method, sprintf("round_%02d", r),
      initial  = current,
      use_log  = USE_LOG,
      maxfeval = MAXFEVAL
    )
    fits[[fit$fit_id]] <- fit

    if (!all(is.finite(fit$par))) {
      message("  fit failed, stopping this chain: ", fit$message)
      break
    }

    current <- setNames(fit$par, PAR_NAMES)
    current <- clamp_to_bounds(current)

    if (!is.null(REL_TOL) && is.finite(fit$initial_loss)) {
      rel_gain <- (fit$initial_loss - fit$value) / abs(fit$initial_loss)
      if (rel_gain < REL_TOL) {
        message("  improvement below REL_TOL, stopping this chain")
        break
      }
    }
  }
}

# ------------------------------------------------------------
# Save and report
# ------------------------------------------------------------
results <- save_fit_results(fits, USE_SWITCH, OUTPUT_DIR)

# Round-by-round trajectory: does the loss keep falling? do parameters settle?
trajectory <- results$summary |>
  select(fit_id, loss, start_id, method, initial, fitted, convergence, n_feval) |>
  left_join(
    results$parameters |> select(-loss, -method, -start_id),
    by = "fit_id"
  ) |>
  arrange(loss, start_id) |>
  select(-fit_id)

print(as.data.frame(trajectory), row.names = FALSE)

report <- write_fit_report(
  results$summary, results$parameters,
  output_dir = OUTPUT_DIR,
  use_switch = USE_SWITCH,
  baseline   = NULL,
  name       = "sequential_fit_comparison"
)
print(report$gt)

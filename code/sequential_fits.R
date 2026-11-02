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
N_SCREEN      <- 20000
MAXFEVAL      <- 20000
SEED          <- 42
REL_TOL       <- NULL # 1e-6
PARALLEL      <- TRUE                      # parallelise the start screen
N_WORKERS     <- max(1, parallelly::availableCores() - 1)

# Kept separate so run_fits.R output is not overwritten.
# One subfolder per seed, so runs with different seeds don't overwrite each other.
SEED_TAG   <- sprintf("seed_%s", SEED)
OUTPUT_DIR <- file.path(common$FIT_SEQUENTIAL_DIR, SEED_TAG)

# ------------------------------------------------------------
# Setup
# ------------------------------------------------------------
set.seed(SEED)
fit_data <- load_fit_data()

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

if (identical(START_FROM, "screen")) {

  screen_starts <- make_starts(n_random = N_SCREEN)   # includes "default"

  # One start -> vector of losses (one per loss function)
  # simulate once per start, apply every loss to that simulation
  screen_one <- function(s) {
    sim <- tryCatch(run_model(as.list(s), USE_SWITCH), error = function(e) NULL)
    vapply(loss_functions, function(f) {
      if (is.null(sim)) return(Inf)
      v <- tryCatch(f(sim, fit_data$data, fit_data$sd), error = function(e) Inf)
      if (is.finite(v)) v else Inf
    }, numeric(1))
  }

  # losses x starts
  if (PARALLEL && N_WORKERS > 1) {
    # forking where possible (Linux/macOS terminal, Docker); fresh sessions otherwise
    use_fork <- parallelly::supportsMulticore()
    if (use_fork) {
      future::plan(future::multicore, workers = N_WORKERS)
    } else {
      future::plan(future::multisession, workers = N_WORKERS)
    }

    t_screen <- system.time(
      screen <- future.apply::future_sapply(
        screen_starts, screen_one,
        future.seed       = TRUE,
        future.chunk.size = ceiling(length(screen_starts) / (N_WORKERS * 4))
      )
    )[["elapsed"]]

    future::plan(future::sequential)
    message("Screened ", length(screen_starts), " starts on ", N_WORKERS, " workers (",
            if (use_fork) "multicore" else "multisession", ") in ", round(t_screen), " s")
  } else {
    t_screen <- system.time(
      screen <- sapply(screen_starts, screen_one)
    )[["elapsed"]]
    message("Screened ", length(screen_starts), " starts sequentially in ",
            round(t_screen), " s")
  }

  saveRDS(screen, file.path(OUTPUT_DIR, paste0("screen_losses_", SEED_TAG, ".rds")))
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
fits         <- list()
chain_starts <- list()   # starting parameters for each loss (for the report)

for (loss_name in LOSSES) {
  fn <- make_loss(fit_data$data, fit_data$sd, loss_functions[[loss_name]], USE_SWITCH)
  current <- if (identical(START_FROM, "screen")) {
      message("  best screened start: ", names(screen_starts)[best_idx[[loss_name]]],
              " (loss ", signif(screen[loss_name, best_idx[[loss_name]]], 4),
              "; default start: ", signif(screen[loss_name, "default"], 4), ")")
      screen_starts[[best_idx[[loss_name]]]]
    } else start

  # NEW: remember where this chain starts (only this one line goes in the loop)
  chain_starts[[loss_name]] <- clamp_to_bounds(current)

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
# Save
# ------------------------------------------------------------
stopifnot("No fits were produced - check the chain loop above" = length(fits) > 0)

results <- save_fit_results(fits, USE_SWITCH, OUTPUT_DIR)

# ------------------------------------------------------------
# NEW: add each chain's starting point as a "round_00" row
# (runs after the loop, once `results` exists)
# ------------------------------------------------------------
# Loss at the start = initial loss of round 1 (same point, same loss function)
start_summary <- results$summary |>
  filter(start_id == "round_01") |>
  transmute(
    fit_id      = paste0(loss, "__start"),
    loss,
    start_id    = "round_00",
    method      = "start",
    initial,
    fitted      = initial,
    convergence = NA,
    n_feval     = 0L
  )

start_pars <- bind_rows(lapply(names(chain_starts), function(l) {
  tibble::as_tibble_row(as.list(chain_starts[[l]][PAR_NAMES])) |>
    mutate(fit_id = paste0(l, "__start"), loss = l,
           method = "start", start_id = "round_00", .before = 1)
}))

summary_with_start    <- bind_rows(start_summary, results$summary)    |> arrange(loss, start_id)
parameters_with_start <- bind_rows(start_pars,    results$parameters) |> arrange(loss, start_id)

# Tag every row with the seed and save, so runs can be stacked and compared later
summary_with_start    <- summary_with_start    |> mutate(seed = SEED, .before = 1)
parameters_with_start <- parameters_with_start |> mutate(seed = SEED, .before = 1)

saveRDS(summary_with_start,    file.path(OUTPUT_DIR, paste0("summary_with_start_",    SEED_TAG, ".rds")))
saveRDS(parameters_with_start, file.path(OUTPUT_DIR, paste0("parameters_with_start_", SEED_TAG, ".rds")))

# ------------------------------------------------------------
# Report
# ------------------------------------------------------------
# Round-by-round trajectory: does the loss keep falling? do parameters settle?
trajectory <- summary_with_start |>
  select(fit_id, loss, start_id, method, initial, fitted, convergence, n_feval) |>
  left_join(
    parameters_with_start |> select(-seed, -loss, -method, -start_id),
    by = "fit_id"
  ) |>
  arrange(loss, start_id) |>
  select(-fit_id)

print(as.data.frame(trajectory), row.names = FALSE)

report <- write_fit_report(
  summary_with_start, parameters_with_start,
  output_dir = OUTPUT_DIR,
  use_switch = USE_SWITCH,
  baseline   = NULL,
  name       = paste0("sequential_fit_comparison_", SEED_TAG)
)
print(report$gt)

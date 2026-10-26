# ============================================================
#  run_fits.R
#  The fitting process. Functions live in parameter_fitting_functions.R
# ============================================================

source(here::here("code", "parameter_fitting_functions.R"))
source(here::here("code", "fit_report_functions.R"))

# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------
USE_SWITCH      <- TRUE
USE_LOG         <- TRUE                    # optimize log(loss)
METHODS         <- names(OPTIMIZERS)       # or e.g. c("nmkb")
LOSSES          <- names(loss_functions)   # or e.g. c("weighted_SSE")
N_RANDOM_STARTS <- 2                       # extra starts besides the default
MAXFEVAL        <- 20000
SEED            <- 1337

OUTPUT_DIR <- common$FIT_DIR

# ------------------------------------------------------------
# Setup
# ------------------------------------------------------------
set.seed(SEED)

fit_data <- load_fit_data()
starts   <- make_starts(
  n_random  = N_RANDOM_STARTS
)

# ------------------------------------------------------------
# Fit: loss x start x method
# ------------------------------------------------------------
fits <- list()

for (loss_name in LOSSES) {
  fn <- make_loss(fit_data$data, fit_data$sd, loss_functions[[loss_name]], USE_SWITCH)

  for (start_id in names(starts)) {
    for (method in METHODS) {
      message("Fitting: ", loss_name, " / ", method, " / ", start_id)

      fit <- fit_one(
        fn, loss_name, method, start_id,
        initial = starts[[start_id]],
        use_log = USE_LOG,
        maxfeval = MAXFEVAL
      )
      fits[[fit$fit_id]] <- fit
    }
  }
}

# ------------------------------------------------------------
# Save and report
# ------------------------------------------------------------
results <- save_fit_results(fits, USE_SWITCH, OUTPUT_DIR)

# ------------------------------------------------------------
# Report: combined table of fit statistics and fitted parameters
# ------------------------------------------------------------
# CC_MAGIC baseline = the default parameters, built in memory
# (no separate CC_MAGIC.rds needed any more)
baseline <- format_parameters(unlist(DEFAULT_PARS)[PAR_NAMES], USE_SWITCH)

report <- write_fit_report(
  results$summary, results$parameters,
  output_dir = OUTPUT_DIR,
  use_switch = USE_SWITCH,
  baseline   = baseline,
  best_only  = FALSE      # TRUE: one row per loss x method
)

print(report$gt)

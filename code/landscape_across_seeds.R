# ============================================================
#  profile_best_fits.R
#  1D loss profiles through the best fit of every seed x loss.
#  Vary one parameter at a time, hold the others at the fitted values.
# ============================================================

source(here::here("code", "parameter_fitting_functions.R"))
library(dplyr)
library(ggplot2)

# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------
USE_SWITCH      <- TRUE
LOSSES          <- names(loss_functions)
N_CORES         <- 1        # > 1 uses parallel::mclapply (Linux/macOS only)
PROFILE_DECADES <- 1.5      # log-scale parameters: +/- this many decades around the fit
N_GRID          <- 50       # grid points per parameter
PROFILE_YMAX    <- 10       # cap on the loss ratio so one bad region doesn't flatten the rest
# "best_seed": ratio to the lowest fitted loss across seeds (shows which seed fit best)
# "own":       ratio to each seed's own fitted loss (every curve touches 1)
RELATIVE_TO     <- "best_seed"

LOG_AXIS <- c(v1 = TRUE, sigma = TRUE, portion = FALSE)
X_LABEL  <- setNames(
  ifelse(LOG_AXIS, paste0("log(", names(LOG_AXIS), ")"), names(LOG_AXIS)),
  names(LOG_AXIS)
)

theme_panels <- theme_bw() +
  theme(
    panel.border     = element_rect(colour = "grey30", fill = NA, linewidth = 0.8),
    panel.grid.minor = element_blank(),
    panel.spacing    = unit(0.9, "lines"),
    strip.background = element_rect(fill = "grey90", colour = "grey30"),
    strip.text       = element_text(face = "bold"),
    legend.position  = "bottom"
  )

par_lapply <- function(X, FUN) {
  if (N_CORES > 1) parallel::mclapply(X, FUN, mc.cores = N_CORES) else lapply(X, FUN)
}

OUTPUT_DIR <- file.path(dirname(common$SA_MODEL_PATH), "landscape")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

fit_data <- load_fit_data()

# ------------------------------------------------------------
# Best fit per seed x loss, from every seed_* folder
# ------------------------------------------------------------
seed_dirs <- list.dirs(common$FIT_SEQUENTIAL_DIR, recursive = FALSE)
seed_dirs <- seed_dirs[grepl("^seed_", basename(seed_dirs))]
stopifnot("No seed_* folders found" = length(seed_dirs) > 0)

read_seed <- function(dir, prefix) {
  f <- file.path(dir, paste0(prefix, basename(dir), ".rds"))
  if (file.exists(f)) readRDS(f) else { message("Skipping ", basename(dir)); NULL }
}
summary_all <- bind_rows(lapply(seed_dirs, read_seed, prefix = "summary_with_start_"))
params_all  <- bind_rows(lapply(seed_dirs, read_seed, prefix = "parameters_with_start_"))

best <- summary_all |>
  filter(loss %in% LOSSES, method != "start", is.finite(fitted)) |>
  group_by(seed, loss) |>
  slice_min(fitted, n = 1, with_ties = FALSE) |>
  ungroup() |>
  select(seed, loss, fit_id, start_id, method, fitted) |>
  left_join(params_all |> select(seed, fit_id, all_of(PAR_NAMES)),
            by = c("seed", "fit_id"))
stopifnot(nrow(best) > 0, !anyNA(best[PAR_NAMES]))
message("Profiling ", nrow(best), " fits (", n_distinct(best$seed), " seeds x ",
        n_distinct(best$loss), " losses)")

# ------------------------------------------------------------
# Evaluation helpers
# ------------------------------------------------------------
eval_par <- function(p, loss_name) {
  sim <- tryCatch(run_model(as.list(p), USE_SWITCH), error = function(e) e)
  if (inherits(sim, "error")) return(NA_real_)
  v <- tryCatch(loss_functions[[loss_name]](sim, fit_data$data, fit_data$sd),
                error = function(e) NA_real_)
  if (is.finite(v)) v else NA_real_
}

profile_grid <- function(par_name, center_value) {
  lo <- PAR_LOWER[[par_name]] + 1e-6     # stay strictly inside the bounds
  hi <- PAR_UPPER[[par_name]] - 1e-6
  cv <- min(max(center_value, lo), hi)
  if (LOG_AXIS[[par_name]]) {
    lo <- max(lo, cv / 10^PROFILE_DECADES)
    hi <- min(hi, cv * 10^PROFILE_DECADES)
    g  <- 10^seq(log10(lo), log10(hi), length.out = N_GRID)
  } else {
    g  <- seq(lo, hi, length.out = N_GRID)
  }
  sort(unique(c(g, cv)))
}

centers <- lapply(seq_len(nrow(best)), function(i) unlist(best[i, PAR_NAMES]))

# Fail loudly if the model can't even run at the fitted points
best$center_loss <- unlist(par_lapply(seq_len(nrow(best)), function(i)
  eval_par(centers[[i]], best$loss[i])))
if (anyNA(best$center_loss)) {
  bad <- which(is.na(best$center_loss))
  print(as.data.frame(best[bad, c("seed", "loss", "fit_id", PAR_NAMES)]))
  i <- bad[1]
  run_model(as.list(centers[[i]]), USE_SWITCH)   # rerun unprotected to show the error
  stop("Loss is NA at the fitted parameters for the fits listed above")
}

# ------------------------------------------------------------
# Profiles: one job per (fit, parameter)
# ------------------------------------------------------------
jobs <- expand.grid(i = seq_len(nrow(best)), par = PAR_NAMES, stringsAsFactors = FALSE)

t_prof <- system.time(
  prof <- bind_rows(par_lapply(seq_len(nrow(jobs)), function(j) {
    i  <- jobs$i[j]; nm <- jobs$par[j]
    center <- centers[[i]]
    g <- profile_grid(nm, center[[nm]])
    vals <- vapply(g, function(x) { p <- center; p[[nm]] <- x; eval_par(p, best$loss[i]) },
                   numeric(1))
    data.frame(seed = best$seed[i], loss = best$loss[i], parameter = nm,
               par_value = g, loss_value = vals,
               center_value = center[[nm]], center_loss = best$center_loss[i])
  }))
)[["elapsed"]]
message("Profiles done in ", round(t_prof), " s; ",
        sum(is.na(prof$loss_value)), " of ", nrow(prof), " evaluations failed")
saveRDS(prof, file.path(OUTPUT_DIR, "profiles_best_by_seed.rds"))

# ------------------------------------------------------------
# Plot: rows = loss, columns = parameter, one line per seed
# ------------------------------------------------------------
prof <- prof |>
  group_by(loss) |>
  mutate(ref_loss = if (RELATIVE_TO == "own") center_loss else min(center_loss)) |>
  ungroup() |>
  mutate(
    parameter    = factor(parameter, levels = PAR_NAMES),
    loss         = factor(loss, levels = LOSSES),
    seed         = factor(seed),
    is_log       = LOG_AXIS[as.character(parameter)],
    x            = ifelse(is_log, log10(par_value), par_value),
    x_center     = ifelse(is_log, log10(center_value), center_value),
    ratio        = pmin(loss_value / ref_loss, PROFILE_YMAX),
    center_ratio = pmin(center_loss / ref_loss, PROFILE_YMAX)
  )

centers_df <- distinct(prof, seed, loss, parameter, x_center, center_ratio)

p_prof <- ggplot(prof, aes(x, ratio, colour = seed)) +
  geom_hline(yintercept = 1, linetype = "dotted", colour = "grey50") +
  geom_vline(data = centers_df, aes(xintercept = x_center, colour = seed),
             linetype = "dashed", alpha = 0.4) +
  geom_line(alpha = 0.8, na.rm = TRUE) +
  geom_point(data = centers_df, aes(x_center, center_ratio), size = 2.2) +
  scale_colour_brewer(palette = "Dark2", name = "Seed") +
  facet_grid(loss ~ parameter, scales = "free_x",
             labeller = labeller(parameter = X_LABEL)) +
  scale_y_log10() +
  labs(
    title    = "1D loss profiles through the best fit of each seed",
    subtitle = paste0(
      "dots / dashed lines = fitted values; other parameters held at that seed's fit; ",
      if (RELATIVE_TO == "own") "loss relative to each seed's own fit"
      else "loss relative to the best fit across seeds",
      "; capped at ", PROFILE_YMAX, "x"),
    x = NULL,
    y = if (RELATIVE_TO == "own") "loss / own fitted loss (log)" else "loss / best fitted loss (log)"
  ) +
  theme_panels

ggsave(file.path(OUTPUT_DIR, "profiles_best_by_seed.pdf"), p_prof,
       width = 16, height = 3.2 * n_distinct(prof$loss), dpi = 150)
print(p_prof)

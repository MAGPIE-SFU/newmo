# ============================================================
#  visualize_landscape.R
#  Screen many random starts, evaluate every loss function on each,
#  and plot loss against each parameter. Standalone: it only needs
#  parameter_fitting_functions.R for the model, data and losses.
# ============================================================

source(here::here("code", "parameter_fitting_functions.R"))
library(ggplot2)

# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------
USE_SWITCH <- TRUE
N_SCREEN   <- 1000                    # random starts (the default start is added too)
LOSSES     <- names(loss_functions)
SEED       <- 1337
N_CORES    <- 1                       # > 1 uses parallel::mclapply (Linux/macOS only)
RECOMPUTE  <- T
TOP_FRAC   <- 0.1                     # second plot zooms on the best fraction per loss

# Plot these parameters on a log10 axis (match the log-scale starts in make_starts)
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
    legend.position  = "none"
  )

# 1D profiles through the fitted optimum
# FITS_PATH       <- common$ALL_BEST_FITS_PATH  # fits table; the last round per loss is used
# seq_seed <- 1337
# FITS_PATH <- file.path(here::here("data", "SouthAfrica", "fit_sequential", paste0("seed_", seq_seed)), "all_fitted_parameters_switch.rds")
# PROFILE_DECADES <- 1.5     # log-scale parameters: +/- this many decades around the fit
# N_GRID          <- 50      # grid points per parameter
# PROFILE_YMAX    <- 10      # cap on loss / fitted loss, so one bad region doesn't flatten the rest

par_lapply <- function(X, FUN) {
  if (N_CORES > 1) parallel::mclapply(X, FUN, mc.cores = N_CORES) else lapply(X, FUN)
}

OUTPUT_DIR <- file.path(dirname(common$SA_MODEL_PATH), "landscape")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
SCREEN_FILE <- file.path(OUTPUT_DIR, "screen.rds")

fit_data <- load_fit_data()

# ------------------------------------------------------------
# Screen: simulate once per start, apply every loss to that simulation
# ------------------------------------------------------------
set.seed(SEED)

if (!RECOMPUTE && file.exists(SCREEN_FILE)) {
  screen <- readRDS(SCREEN_FILE)
  message("Loaded saved screen: ", nrow(screen), " starts")
} else {
  starts <- make_starts(n_random = N_SCREEN)      # includes "default"

  eval_start <- function(s) {
    sim <- tryCatch(run_model(as.list(s), USE_SWITCH), error = function(e) NULL)
    vapply(loss_functions[LOSSES], function(f) {
      if (is.null(sim)) return(Inf)
      v <- tryCatch(f(sim, fit_data$data, fit_data$sd), error = function(e) Inf)
      if (is.finite(v)) v else Inf
    }, numeric(1))
  }

  invisible(get_samodel())    # load the model once, before any forking

  elapsed <- system.time(
    results <- if (N_CORES > 1) {
      parallel::mclapply(starts, eval_start, mc.cores = N_CORES)
    } else {
      lapply(starts, eval_start)
    }
  )[["elapsed"]]
  message("Screened ", length(starts), " starts in ", round(elapsed), " s")

  screen <- data.frame(
    start_id = names(starts),
    do.call(rbind, starts),
    do.call(rbind, results),
    row.names = NULL,
    check.names = FALSE
  )
  saveRDS(screen, SCREEN_FILE)
}

# Failures (errors or non-finite losses)
for (l in LOSSES) {
  message(l, ": ", sum(!is.finite(screen[[l]])), " of ", nrow(screen), " starts failed")
}

# ------------------------------------------------------------
# Long format for plotting
# ------------------------------------------------------------
long <- screen |>
  pivot_longer(all_of(PAR_NAMES), names_to = "parameter", values_to = "par_value") |>
  pivot_longer(all_of(LOSSES),    names_to = "loss",      values_to = "loss_value") |>
  filter(is.finite(loss_value), loss_value > 0) |>
  mutate(
    parameter = factor(parameter, levels = PAR_NAMES),
    loss      = factor(loss, levels = LOSSES),
    x         = ifelse(LOG_AXIS[as.character(parameter)], log10(par_value), par_value)
  )

# ------------------------------------------------------------
# Plot: loss against each parameter, one row per loss function
# ------------------------------------------------------------
plot_landscape <- function(d, title) {
  best <- d |>
    group_by(loss) |>
    filter(loss_value == min(loss_value)) |>
    ungroup()

  ggplot(d, aes(x, loss_value)) +
    # geom_point(alpha = 0.25, size = 0.7) +
    geom_point(aes(colour = loss), alpha = 0.3, size = 0.7) +
    scale_colour_brewer(palette = "Dark2") +
    geom_point(data = filter(d, start_id == "default"), colour = "red", size = 2.5) +
    geom_point(data = best, colour = "dodgerblue3", size = 2.5) +
    facet_grid(loss ~ parameter, scales = "free",
               labeller = labeller(parameter = X_LABEL)) +
    scale_y_log10() +
    labs(
      title = title,
      subtitle = "red = default start, blue = best start for that loss",
      x = NULL, y = "loss (log scale)"
    ) +
    theme_panels
    # theme_minimal() +
    # theme(strip.text = element_text(face = "bold"))
}

p_all <- plot_landscape(long, paste0("Loss vs parameter: ", N_SCREEN, " random starts"))

top <- long |>
  group_by(loss) |>
  filter(loss_value <= quantile(loss_value, TOP_FRAC)) |>
  ungroup()
# p_top <- plot_landscape(top, paste0("Best ", TOP_FRAC * 100, "% of starts per loss"))

ggsave(file.path(OUTPUT_DIR, "landscape_all.pdf"), p_all,
       width = 16, height = 3.2 * length(LOSSES), dpi = 150)
# ggsave(file.path(OUTPUT_DIR, "landscape_top.pdf"), p_top,
#        width = 16, height = 3.2 * length(LOSSES), dpi = 150)

print(p_all)
# print(p_top)

# Best screened start per loss (feed these into the fits)
best_starts <- do.call(rbind, lapply(LOSSES, function(l) {
  i <- which.min(screen[[l]])
  data.frame(loss = l, start_id = screen$start_id[i], best_loss = screen[[l]][i],
             default_loss = screen[[l]][screen$start_id == "default"],
             screen[i, PAR_NAMES])
}))
print(best_starts, row.names = FALSE)

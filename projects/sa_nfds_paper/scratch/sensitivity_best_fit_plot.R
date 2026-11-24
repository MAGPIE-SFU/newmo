# ============================================================
#  Scratch: overall best fit plot (top_group_plot style) with the
#  sensitivity ribbon added
#
#  Reads the precomputed results from the pipeline (nothing is refitted
#  or resimulated), builds the same plot as fit_figures_overall_best.pdf
#  with top_group_plot(), and adds the sensitivity band underneath.
#  Run from the repo root after tar_make().
# ============================================================

library(dplyr)
library(ggplot2)

Sys.setenv(TAR_PROJECT = "sa_nfds_paper")
pkgload::load_all("nfdsmodelling", quiet = TRUE)
source("projects/sa_nfds_paper/config.R")

# ---- Pipeline results ---------------------------------------------------------
samodel     <- targets::tar_read(samodel)
xt          <- targets::tar_read(xt_grouped_uncertainty)
sim_best    <- targets::tar_read(sim_overall_best)          # named list, one simulation
sensitivity <- targets::tar_read(sensitivity)
best_key    <- targets::tar_read(fit_selection)$overall_best$key

# ---- Settings (defaults from config.R; change freely here) -------------------
s      <- FIG_SETTINGS     # vax_start, prevax_years, n_groups, vaccine
qtile  <- SENS_QTILE
sim_label   <- "Overall best"
ribbon_fill <- "#1B9E77"   # first Dark2 colour = colour of the first simulation line

# ---- The overall best fit plot, as in the pipeline ---------------------------
p <- top_group_plot(
  vax_start = s$vax_start, prevax_years = s$prevax_years,
  n_groups = s$n_groups,
  samodel, xt,
  setNames(sim_best, sim_label),
  drop_vt = FALSE,
  vaccine = s$vaccine
)

# the groups shown in the plot (the first layer is the data error bars)
plot_groups <- levels(p$layers[[1]]$data$group)

# ---- Sensitivity ribbon, on the same scale and time axis as the plot ---------
# top_group_plot() shows each simulation as frequencies (normalised per year)
# and puts model year k at vax_start - 1 + k. Do the same for every draw.
ribbon <- sensitivity$sim_out |>
  filter(loss == best_key, year > 0) |>              # year 0 = x0, shown separately
  group_by(rep, year) |>
  mutate(frac = value / sum(value)) |>               # normalise within each draw and year
  ungroup() |>
  filter(group %in% plot_groups) |>
  group_by(group, year) |>
  summarise(L = quantile(frac, (1 - qtile) / 2),
            U = quantile(frac, 1 - (1 - qtile) / 2),
            .groups = "drop") |>
  mutate(year  = s$vax_start - 1 + year,
         group = factor(group, levels = plot_groups))

ribbon_layer <- geom_ribbon(data = ribbon, aes(x = year, ymin = L, ymax = U),
                            fill = ribbon_fill, alpha = 0.25, inherit.aes = FALSE)

# put the ribbon underneath everything else
p$layers <- c(list(ribbon_layer), p$layers)

p <- p + labs(subtitle = paste0("Overall best fit (", best_key, ") with ",
                                qtile * 100, "% sensitivity band"))
print(p)

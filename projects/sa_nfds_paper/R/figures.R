# CC_MAGIC baseline parameters from DEFAULT_PARS (config.R).
# Same as parameter_generator(how_numerical = "CC_MAGIC"): wlow = whigh * wfrac,
# v2 = v1, plus tswitch and the vaccines (PCV7 -> PCV13).
cc_magic_modelpars <- function(default_pars, tswitch, use_switch = TRUE) {
  modelpars <- default_pars
  modelpars$wlow  <- modelpars$whigh * modelpars$wfrac
  modelpars$wfrac <- NULL
  modelpars$v2    <- modelpars$v1
  modelpars$tswitch  <- if (use_switch) tswitch else 0
  modelpars$vaccine1 <- VACCINE_TYPES[["PCV7"]]
  modelpars$vaccine2 <- VACCINE_TYPES[["PCV13"]]
  modelpars
}

# The three plots for one set of simulations
fit_figure_set <- function(sims, samodel, xt, title, s) {
  list(
    top_group_plot(s$vax_start, s$prevax_years, s$n_groups, samodel, xt, sims,
                   drop_vt = FALSE, vaccine = s$vaccine) + ggplot2::labs(subtitle = title),
    top_sero_plot(s$vax_start, s$prevax_years, s$n_serotypes, samodel, xt, sims,
                  drop_vt = FALSE, vaccine = s$vaccine) + ggplot2::labs(subtitle = title),
    vt_nvt_plot(s$vax_start, s$prevax_years, samodel, xt, sims,
                vaccine = s$vaccine) + ggplot2::labs(subtitle = title)
  )
}

# Best fit per seed x loss: one set of pages per loss, one line per seed
all_fits_figure_set <- function(sims, all_best, samodel, xt, s) {
  unlist(lapply(unique(all_best$loss), function(l) {
    rows <- all_best[all_best$loss == l, ]
    fit_figure_set(setNames(sims[rows$key], paste("seed", rows$seed)),
                   samodel, xt, paste("Best fit per seed:", l), s)
  }), recursive = FALSE)
}

# Simulate a table of fits (one row per fit, with PAR_NAMES and key)
simulate_fits <- function(fits, samodel, nyears, tswitch, use_switch) {
  setNames(lapply(seq_len(nrow(fits)), function(i) {
    run_model(as.list(unlist(fits[i, PAR_NAMES])), use_switch, samodel, nyears, tswitch)
  }), fits$key)
}

# Several plots into one multi-page PDF
save_plots_pdf <- function(plots, out_file, width, height) {
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  grDevices::pdf(out_file, width = width, height = height)
  on.exit(grDevices::dev.off())
  for (p in plots) print(p)
  out_file
}

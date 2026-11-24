sens_dir <- file.path(results_dir, "plots", "sensitivity")

targets_sensitivity <- list(
  # Optimal fit per loss: the best across seeds
  # One sensitivity run per fit (best per seed x loss); `loss` is the fit's id
  tar_target(
    sens_opt_pars,
    fit_selection$all_best |> dplyr::mutate(loss = key)
  ),
  tar_target(sens_top_groups, sensitivity_top_groups(samodel, SENS_N_GROUPS)),
  tar_target(
    sensitivity,
    run_sensitivity(
      sens_opt_pars,
      samodel,
      tswitch = T_SWITCH,
      nyears = NYEARS,
      sd_pars = SENS_SD_PARS,
      use_switch = FIT_USE_SWITCH,
      n_sims = SENS_N_SIMS,
      n_cores = SENS_N_CORES
    )
  ),
  tar_target(
    sensitivity_failed_csv,
    write_csv_output(
      sensitivity$failed,
      file.path(results_dir, "computed_data", "sensitivity_failed.csv")
    ),
    format = "file"
  ),

  # Overall best fit
  tar_target(
    sensitivity_best_pdf,
    save_plot_output(
      plot_sensitivity_vs_data(
        sensitivity,
        xt_grouped_uncertainty,
        sens_top_groups,
        SENS_START_YEAR,
        SENS_QTILE,
        best_loss = fit_selection$overall_best$key
      ),
      file.path(sens_dir, "sensitivity_overall_best.pdf"),
      width = 9,
      height = 7
    ),
    format = "file"
  ),

  # Every fit, one page each
  tar_target(
    sensitivity_per_fit_pdf,
    save_plots_pdf(
      lapply(sens_opt_pars$key, function(k) {
        plot_sensitivity_vs_data(
          sensitivity,
          xt_grouped_uncertainty,
          sens_top_groups,
          SENS_START_YEAR,
          SENS_QTILE,
          best_loss = k
        )
      }),
      file.path(sens_dir, "sensitivity_per_fit.pdf"),
      width = 9,
      height = 7
    ),
    format = "file"
  ),
  tar_target(
    sensitivity_by_loss_pdfs,
    unname(vapply(
      unique(sens_opt_pars$loss),
      function(l) {
        save_plot_output(
          plot_sensitivity_trajectories(
            sensitivity,
            sens_top_groups,
            SENS_START_YEAR,
            losses = l
          ),
          file.path(sens_dir, paste0("sim_trajectories_", l, ".pdf")),
          width = 9,
          height = 7
        )
      },
      character(1)
    )),
    format = "file"
  ),
  tar_target(
    sensitivity_best_groups_pdf,
    save_plot_output(
      add_sensitivity_ribbon(
        top_group_plot(
          FIG_SETTINGS$vax_start,
          FIG_SETTINGS$prevax_years,
          FIG_SETTINGS$n_groups,
          samodel,
          xt_grouped_uncertainty,
          setNames(sim_overall_best, "Overall best"),
          drop_vt = FALSE,
          vaccine = FIG_SETTINGS$vaccine
        ),
        sensitivity,
        fit_selection$overall_best$key,
        vax_start = FIG_SETTINGS$vax_start,
        qtile = SENS_QTILE
      ) +
        ggplot2::labs(
          subtitle = paste0(
            "Overall best fit (",
            fit_selection$overall_best$key,
            ") with ",
            SENS_QTILE * 100,
            "% sensitivity band"
          )
        ),
      file.path(sens_dir, "overall_best_with_sensitivity.pdf"),
      width = 14,
      height = 10
    ),
    format = "file"
  ),
  tar_target(
    sensitivity_per_fit_groups_pdf,
    save_plots_pdf(
      lapply(fit_selection$all_best$key, function(k) {
        add_sensitivity_ribbon(
          top_group_plot(
            FIG_SETTINGS$vax_start,
            FIG_SETTINGS$prevax_years,
            FIG_SETTINGS$n_groups,
            samodel,
            xt_grouped_uncertainty,
            sims_all_best[k],
            drop_vt = FALSE,
            vaccine = FIG_SETTINGS$vaccine
          ),
          sensitivity,
          k,
          vax_start = FIG_SETTINGS$vax_start,
          qtile = SENS_QTILE
        ) +
          ggplot2::labs(
            subtitle = paste0(
              "Fit ",
              k,
              " with ",
              SENS_QTILE * 100,
              "% sensitivity band"
            )
          )
      }),
      file.path(sens_dir, "per_fit_with_sensitivity.pdf"),
      width = 14,
      height = 10
    ),
    format = "file"
  )
)

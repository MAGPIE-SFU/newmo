# Package index

## Data preparation

- [`clean_metadata()`](clean_metadata.md) : Clean isolate metadata
- [`build_invasiveness_lookup()`](build_invasiveness_lookup.md) : Build
  the invasiveness lookup with intervals
- [`shrink_invasiveness()`](shrink_invasiveness.md) : Empirical Bayes
  shrinkage of invasiveness estimates
- [`reweight_carriage()`](reweight_carriage.md) : Reweight carriage
  composition and define model groups
- [`reweight_carriage_uncertainty()`](reweight_carriage_uncertainty.md)
  : Grouped x(t) with invasiveness uncertainty
- [`plot_xt_top_groups()`](plot_xt_top_groups.md) : Plot the most common
  groups in x(t) with uncertainty
- [`build_gnfds()`](build_gnfds.md) : Build the grouped G matrix for the
  NFDS model
- [`build_samodel()`](build_samodel.md) : Build the NFDS model structure

## Model

- [`simulate_nfds()`](simulate_nfds.md) : Simulate the NFDS model
- [`run_model()`](run_model.md) : Simulate the model for a set of fitted
  parameters
- [`VACCINE_TYPES`](VACCINE_TYPES.md) : Serotypes covered by each
  pneumococcal vaccine
- [`vaccine_serotypes()`](vaccine_serotypes.md) : Get the serotypes of
  one vaccine

## Fitting

- [`run_sequential_fits()`](run_sequential_fits.md) : Sequential
  (chained) fits for one seed
- [`select_best_fits()`](select_best_fits.md) : Pick the best fits
  across seeds
- [`load_fit_data()`](load_fit_data.md) : Fit data from x(t)
- [`make_loss()`](make_loss.md) : Turn a loss function into an objective
  of the parameter vector alone
- [`loss_functions`](loss_functions.md) : Loss functions for fitting
- [`PAR_NAMES`](PAR_NAMES.md) : Names of the fitted parameters

## Fit diagnostics and sensitivity

- [`profile_best_fits()`](profile_best_fits.md) : 1D loss profiles
  through the best fits
- [`plot_profiles()`](plot_profiles.md) : Plot the 1D loss profiles
- [`screen_landscape()`](screen_landscape.md) : Screen random starts for
  the loss landscape
- [`plot_landscape()`](plot_landscape.md) : Plot the loss landscape
- [`landscape_best_starts()`](landscape_best_starts.md) : Best screened
  start per loss
- [`run_sensitivity()`](run_sensitivity.md) : Sensitivity of the model
  to its parameters around the fitted optima
- [`sensitivity_top_groups()`](sensitivity_top_groups.md) : Largest
  groups at the start of the simulation
- [`plot_sensitivity_optimal()`](sensitivity_plots.md)
  [`plot_sensitivity_draws()`](sensitivity_plots.md)
  [`plot_sensitivity_ribbons()`](sensitivity_plots.md)
  [`plot_sensitivity_trajectories()`](sensitivity_plots.md)
  [`plot_sensitivity_vs_data()`](sensitivity_plots.md) : Sensitivity
  plots

## Projected IPD and vaccine design

- [`simulate_projected_ipd()`](simulate_projected_ipd.md) : Simulate
  projected IPD with uncertainty from invasiveness estimates
- [`plot_projected_ipd_check()`](plot_projected_ipd_check.md) : Observed
  vs predicted IPD serotype composition
- [`prepare_ipd_serotypes()`](prepare_ipd_serotypes.md) : Prepare the
  predicted IPD serotype composition for plotting
- [`plot_ipd_total()`](plot_ipd_total.md) : Total predicted IPD with
  serotype contributions
- [`plot_ipd_dotplot()`](plot_ipd_dotplot.md) : Predicted serotype
  contribution to IPD, with uncertainty intervals
- [`getIPD()`](getIPD.md) : Predicted IPD of a simulation
- [`score_vaccine()`](score_vaccine.md) : Score a vaccine by the IPD it
  leaves
- [`compress_vstring()`](compress_vstring.md) : Vaccine as a single
  string
- [`design_vaccines()`](design_vaccines.md) : Invent vaccines
- [`plot_vaccine_heatmap()`](plot_vaccine_heatmap.md) : Vaccine heatmap

## Model vs data plots

- [`top_group_plot()`](top_group_plot.md) : Plot the most common ST-SC
  groups, model and data together
- [`top_sero_plot()`](top_sero_plot.md) : Plot the most common
  serotypes, model and data together
- [`vt_nvt_plot()`](vt_nvt_plot.md) : Plot VT and NVT totals, model and
  data together

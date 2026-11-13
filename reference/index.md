# Package index

## Data preparation

- [`clean_metadata()`](https://magpie-sfu.github.io/newmo/reference/clean_metadata.md)
  : Clean isolate metadata
- [`build_invasiveness_lookup()`](https://magpie-sfu.github.io/newmo/reference/build_invasiveness_lookup.md)
  : Build the invasiveness lookup with intervals
- [`shrink_invasiveness()`](https://magpie-sfu.github.io/newmo/reference/shrink_invasiveness.md)
  : Empirical Bayes shrinkage of invasiveness estimates
- [`reweight_carriage()`](https://magpie-sfu.github.io/newmo/reference/reweight_carriage.md)
  : Reweight carriage composition and define model groups
- [`reweight_carriage_uncertainty()`](https://magpie-sfu.github.io/newmo/reference/reweight_carriage_uncertainty.md)
  : Grouped x(t) with invasiveness uncertainty
- [`plot_xt_top_groups()`](https://magpie-sfu.github.io/newmo/reference/plot_xt_top_groups.md)
  : Plot the most common groups in x(t) with uncertainty
- [`build_gnfds()`](https://magpie-sfu.github.io/newmo/reference/build_gnfds.md)
  : Build the grouped G matrix for the NFDS model
- [`build_samodel()`](https://magpie-sfu.github.io/newmo/reference/build_samodel.md)
  : Build the NFDS model structure

## Model

- [`simulate_nfds()`](https://magpie-sfu.github.io/newmo/reference/simulate_nfds.md)
  : Simulate the NFDS model
- [`run_model()`](https://magpie-sfu.github.io/newmo/reference/run_model.md)
  : Simulate the model for a set of fitted parameters
- [`VACCINE_TYPES`](https://magpie-sfu.github.io/newmo/reference/VACCINE_TYPES.md)
  : Serotypes covered by each pneumococcal vaccine
- [`vaccine_serotypes()`](https://magpie-sfu.github.io/newmo/reference/vaccine_serotypes.md)
  : Get the serotypes of one vaccine

## Fitting

- [`run_sequential_fits()`](https://magpie-sfu.github.io/newmo/reference/run_sequential_fits.md)
  : Sequential (chained) fits for one seed
- [`select_best_fits()`](https://magpie-sfu.github.io/newmo/reference/select_best_fits.md)
  : Pick the best fits across seeds
- [`load_fit_data()`](https://magpie-sfu.github.io/newmo/reference/load_fit_data.md)
  : Fit data from x(t)
- [`make_loss()`](https://magpie-sfu.github.io/newmo/reference/make_loss.md)
  : Turn a loss function into an objective of the parameter vector alone
- [`loss_functions`](https://magpie-sfu.github.io/newmo/reference/loss_functions.md)
  : Loss functions for fitting
- [`PAR_NAMES`](https://magpie-sfu.github.io/newmo/reference/PAR_NAMES.md)
  : Names of the fitted parameters

## Fit diagnostics and sensitivity

- [`profile_best_fits()`](https://magpie-sfu.github.io/newmo/reference/profile_best_fits.md)
  : 1D loss profiles through the best fits
- [`plot_profiles()`](https://magpie-sfu.github.io/newmo/reference/plot_profiles.md)
  : Plot the 1D loss profiles
- [`screen_landscape()`](https://magpie-sfu.github.io/newmo/reference/screen_landscape.md)
  : Screen random starts for the loss landscape
- [`plot_landscape()`](https://magpie-sfu.github.io/newmo/reference/plot_landscape.md)
  : Plot the loss landscape
- [`landscape_best_starts()`](https://magpie-sfu.github.io/newmo/reference/landscape_best_starts.md)
  : Best screened start per loss
- [`run_sensitivity()`](https://magpie-sfu.github.io/newmo/reference/run_sensitivity.md)
  : Sensitivity of the model to its parameters around the fitted optima
- [`sensitivity_top_groups()`](https://magpie-sfu.github.io/newmo/reference/sensitivity_top_groups.md)
  : Largest groups at the start of the simulation
- [`plot_sensitivity_optimal()`](https://magpie-sfu.github.io/newmo/reference/sensitivity_plots.md)
  [`plot_sensitivity_draws()`](https://magpie-sfu.github.io/newmo/reference/sensitivity_plots.md)
  [`plot_sensitivity_ribbons()`](https://magpie-sfu.github.io/newmo/reference/sensitivity_plots.md)
  [`plot_sensitivity_trajectories()`](https://magpie-sfu.github.io/newmo/reference/sensitivity_plots.md)
  [`plot_sensitivity_vs_data()`](https://magpie-sfu.github.io/newmo/reference/sensitivity_plots.md)
  : Sensitivity plots

## Projected IPD and vaccine design

- [`simulate_projected_ipd()`](https://magpie-sfu.github.io/newmo/reference/simulate_projected_ipd.md)
  : Simulate projected IPD with uncertainty from invasiveness estimates
- [`plot_projected_ipd_check()`](https://magpie-sfu.github.io/newmo/reference/plot_projected_ipd_check.md)
  : Observed vs predicted IPD serotype composition
- [`prepare_ipd_serotypes()`](https://magpie-sfu.github.io/newmo/reference/prepare_ipd_serotypes.md)
  : Prepare the predicted IPD serotype composition for plotting
- [`plot_ipd_total()`](https://magpie-sfu.github.io/newmo/reference/plot_ipd_total.md)
  : Total predicted IPD with serotype contributions
- [`plot_ipd_dotplot()`](https://magpie-sfu.github.io/newmo/reference/plot_ipd_dotplot.md)
  : Predicted serotype contribution to IPD, with uncertainty intervals
- [`getIPD()`](https://magpie-sfu.github.io/newmo/reference/getIPD.md) :
  Predicted IPD of a simulation
- [`score_vaccine()`](https://magpie-sfu.github.io/newmo/reference/score_vaccine.md)
  : Score a vaccine by the IPD it leaves
- [`compress_vstring()`](https://magpie-sfu.github.io/newmo/reference/compress_vstring.md)
  : Vaccine as a single string
- [`design_vaccines()`](https://magpie-sfu.github.io/newmo/reference/design_vaccines.md)
  : Invent vaccines
- [`plot_vaccine_heatmap()`](https://magpie-sfu.github.io/newmo/reference/plot_vaccine_heatmap.md)
  : Vaccine heatmap

## Model vs data plots

- [`top_group_plot()`](https://magpie-sfu.github.io/newmo/reference/top_group_plot.md)
  : Plot the most common ST-SC groups, model and data together
- [`top_sero_plot()`](https://magpie-sfu.github.io/newmo/reference/top_sero_plot.md)
  : Plot the most common serotypes, model and data together
- [`vt_nvt_plot()`](https://magpie-sfu.github.io/newmo/reference/vt_nvt_plot.md)
  : Plot VT and NVT totals, model and data together

# Sensitivity plots

Sensitivity plots

## Usage

``` r
plot_sensitivity_optimal(sens, top_groups)

plot_sensitivity_draws(sens)

plot_sensitivity_ribbons(sens, top_groups, start_year, qtile)

plot_sensitivity_trajectories(sens, top_groups, start_year, losses = NULL)

plot_sensitivity_vs_data(
  sens,
  xt,
  top_groups,
  start_year,
  qtile,
  best_loss = NULL
)
```

## Arguments

- sens:

  Result of
  [`run_sensitivity()`](https://magpie-sfu.github.io/newmo/reference/run_sensitivity.md).

- top_groups:

  Result of
  [`sensitivity_top_groups()`](https://magpie-sfu.github.io/newmo/reference/sensitivity_top_groups.md).

- start_year:

  Year prior to start of model evaluation (i.e., fitted trajectories
  start at start_year + 1).

- qtile:

  Summary quantile for the ribbons.

- losses:

  Loss functions to show (NULL = all).

- xt:

  Grouped x(t) (the data).

- best_loss:

  Loss function of the chosen optimal fit, or NULL for all.

## Value

A ggplot.

## Functions

- `plot_sensitivity_optimal()`: Trajectories of the optimal fits.

- `plot_sensitivity_draws()`: Range of the simulated parameters
  (sampling scale).

- `plot_sensitivity_ribbons()`: Ribbons with the optimal trajectories
  (dashed).

- `plot_sensitivity_trajectories()`: All simulated trajectories, for all
  or some losses.

- `plot_sensitivity_vs_data()`: The data with the optimal fits and the
  sensitivity ribbons, for all losses or only `best_loss`.

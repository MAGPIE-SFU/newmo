# Sensitivity of the model to its parameters around the fitted optima

For each optimal fit (one row of `opt_pars` per loss function), draws
`n_sims` parameter sets from normal distributions around the fitted
values (log scale for v1 and sigma, truncated to the bounds lo and hi)
and simulates each.

## Usage

``` r
run_sensitivity(
  opt_pars,
  samodel,
  tswitch,
  nyears,
  sd_pars,
  use_switch = TRUE,
  n_sims = 100,
  n_cores = 1
)
```

## Arguments

- opt_pars:

  Optimal fits, one row per loss, with a `loss` column and PAR_NAMES.

- samodel, tswitch, nyears, use_switch:

  As in [`run_model()`](run_model.md).

- sd_pars:

  Data frame with columns par (logv1, logsigma, portion), SD, lo and hi
  (bounds on the sampling scale).

- n_sims:

  Number of parameter draws per loss function.

- n_cores:

  Cores for the simulations (\> 1 uses mclapply, Linux/macOS only).

## Value

A list with `opt_traj` (trajectories of the optimal fits), `par_df` (the
draws, sampling scale), `opt_settings` (the optimal values, sampling
scale), `sim_out` (trajectories of all draws) and `failed` (draws whose
simulation failed).

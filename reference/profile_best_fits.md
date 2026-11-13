# 1D loss profiles through the best fits

For the best fit of every seed x loss, varies one parameter at a time
over a grid (holding the others at the fitted values) and evaluates the
loss.

## Usage

``` r
profile_best_fits(
  best,
  samodel,
  xt,
  tswitch,
  nyears,
  use_switch = TRUE,
  losses = names(loss_functions),
  n_cores = 1,
  profile_decades = 1.5,
  n_grid = 50
)
```

## Arguments

- best:

  Best fit per seed x loss, e.g. `select_best_fits()$all_best`.

- samodel, xt, tswitch, nyears, use_switch:

  As in
  [`run_sequential_fits()`](https://magpie-sfu.github.io/newmo/reference/run_sequential_fits.md).

- losses:

  Loss functions to profile.

- n_cores:

  Cores for the profiles (\> 1 uses mclapply, Linux/macOS only).

- profile_decades:

  Log-scale parameters: +/- this many decades around the fit.

- n_grid:

  Grid points per parameter.

## Value

One row per (seed, loss, parameter, grid value).

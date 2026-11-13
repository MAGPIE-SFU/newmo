# Sequential (chained) fits for one seed

Screens `n_screen` random starts (plus the default start), then for each
loss function runs `n_rounds` fits, each starting from the previous
round's optimum, cycling through `chain_methods`.

## Usage

``` r
run_sequential_fits(
  samodel,
  xt,
  seed,
  default_pars,
  tswitch,
  nyears,
  use_switch = TRUE,
  use_log = TRUE,
  losses = names(loss_functions),
  chain_methods = c("nmkb", "lbfgsb", "hjkb"),
  n_rounds = 6,
  start_from = "screen",
  n_screen = 20000,
  maxfeval = 20000,
  rel_tol = NULL,
  parallel = TRUE,
  n_workers = max(1, parallelly::availableCores() - 1)
)
```

## Arguments

- samodel:

  The model structure, e.g. from [`build_samodel()`](build_samodel.md).

- xt:

  Grouped x(t) to fit to.

- seed:

  Random seed for this run.

- default_pars:

  The CC_MAGIC parameters, used as the default start.

- tswitch:

  Time of the vaccine switch.

- nyears:

  Number of (final) years to fit.

- use_switch:

  Use the vaccine switch.

- use_log:

  Minimise log(loss) instead of the loss.

- losses:

  Names of the loss functions to fit.

- chain_methods:

  Optimizers, cycled through the rounds.

- n_rounds:

  Rounds per chain.

- start_from:

  "screen", "random", "default", or a numeric start vector.

- n_screen:

  Number of random starts to screen.

- maxfeval:

  Evaluation budget per fit.

- rel_tol:

  Stop a chain when the relative improvement falls below this, or NULL
  to always run all rounds.

- parallel, n_workers:

  Parallelise the start screen.

## Value

A list with `seed`, `summary` and `parameters` (both with a "round_00"
start row per loss and a seed column), `screen` (losses x starts, or
NULL) and `fits` (all fit objects).

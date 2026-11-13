# Screen random starts for the loss landscape

Screen random starts for the loss landscape

## Usage

``` r
screen_landscape(
  samodel,
  xt,
  default_pars,
  seed,
  tswitch,
  nyears,
  use_switch = TRUE,
  n_screen = 1000,
  losses = names(loss_functions),
  n_cores = 1
)
```

## Arguments

- samodel:

  The model structure, e.g. from [`build_samodel()`](build_samodel.md).

- xt:

  Grouped x(t) to fit to.

- default_pars:

  The CC_MAGIC parameters, used as the default start.

- seed:

  Random seed for this run.

- tswitch:

  Time of the vaccine switch.

- nyears:

  Number of (final) years to fit.

- use_switch:

  Use the vaccine switch.

- n_screen:

  Number of random starts (the default start is added too).

- losses:

  Loss functions to evaluate.

- n_cores:

  Cores (\> 1 uses mclapply, Linux/macOS only).

## Value

One row per start: start_id, the parameters, and one column per loss.

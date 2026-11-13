# Simulate the model for a set of fitted parameters

Simulation only: no data, no loss. Adds the fixed parameters and the
vaccines (PCV7, then PCV13) to `model_pars`.

## Usage

``` r
run_model(model_pars = list(), use_switch, samodel, nyears, tswitch)
```

## Arguments

- model_pars:

  Named list or vector with the fitted parameters (PAR_NAMES).

- use_switch:

  Use the vaccine switch at `tswitch`.

- samodel:

  The model structure, e.g. from [`build_samodel()`](build_samodel.md).

- nyears:

  Number of years to simulate.

- tswitch:

  Time of the vaccine switch (used if `use_switch`).

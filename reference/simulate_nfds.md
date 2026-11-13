# Simulate the NFDS model

Simulate the NFDS model

## Usage

``` r
simulate_nfds(samodel, modelpars, nyears, stepspy = 24)
```

## Arguments

- samodel:

  Model object (needs x0, G, el, st).

- modelpars:

  Parameter list, e.g. from parameter_generator().

- nyears:

  Number of years to simulate.

- stepspy:

  ODE output steps per year.

## Value

Data frame of yearly average strain sizes: one row per strain, columns
year.1 ... year..

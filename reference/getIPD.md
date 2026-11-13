# Predicted IPD of a simulation

Predicted IPD of a simulation

## Usage

``` r
getIPD(simulation, inv, group_lookup, years = 6)
```

## Arguments

- simulation:

  A simulation (rows = groups, columns year.1 ...).

- inv:

  Invasiveness lookup (serotype, v_child).

- group_lookup:

  Group definitions (group, serotype).

- years:

  Simulation years to average over, e.g. 6:10.

## Value

A list with `total` (summed IPD) and `seroIPD` (IPD per serotype).

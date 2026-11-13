# Build the NFDS model structure

Build the NFDS model structure

## Usage

``` r
build_samodel(G, xt, locus_by_year, x0_years)
```

## Arguments

- G:

  Grouped G for the NFDS model, e.g. `build_gnfds()$G`.

- xt:

  Grouped x(t), e.g. from
  [`reweight_carriage_uncertainty()`](reweight_carriage_uncertainty.md).

- locus_by_year:

  Locus statistics, e.g. `build_gnfds()$locus_by_year`.

- x0_years:

  Years averaged for the starting composition x0.

## Value

A list with G, st (serotype of each row of G), x0 and el (equilibrium
locus frequencies).

# Prepare the predicted IPD serotype composition for plotting

Checks the columns, drops non-typeables and missing values, and orders
the serotypes by overall predicted IPD burden.

## Usage

``` r
prepare_ipd_serotypes(ipd_serotype_summary, non_typeable)
```

## Arguments

- ipd_serotype_summary:

  From [`simulate_projected_ipd()`](simulate_projected_ipd.md).

- non_typeable:

  Serotype labels to leave out.

## Value

The filtered table, with serotype as a factor (largest burden last).

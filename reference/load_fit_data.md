# Fit data from x(t)

Fit data from x(t)

## Usage

``` r
load_fit_data(xt, nyears, drop_groups = NULL)
```

## Arguments

- xt:

  Grouped x(t) with columns group, year, frac, frac_sd.

- nyears:

  Number of (final) years to fit to.

- drop_groups:

  Groups to leave out, or NULL.

## Value

list(data = , sd = ).

# Predicted serotype contribution to IPD, with uncertainty intervals

Intervals represent uncertainty in invasiveness, not Poisson predictive
intervals. All years in one panel, coloured and dodged by year.

## Usage

``` r
plot_ipd_dotplot(ipd_sero, top_n = 20)
```

## Arguments

- ipd_sero:

  From [`prepare_ipd_serotypes()`](prepare_ipd_serotypes.md).

- top_n:

  Number of serotypes (by total expected IPD) to show.

## Value

A ggplot.

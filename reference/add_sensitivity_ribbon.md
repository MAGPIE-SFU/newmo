# Add a sensitivity band to a top_group_plot()

Adds the quantile band of the sensitivity draws for one fit underneath a
plot made with
[`top_group_plot()`](https://magpie-sfu.github.io/newmo/reference/top_group_plot.md),
on the same frequency scale and time axis.

## Usage

``` r
add_sensitivity_ribbon(
  p,
  sens,
  key,
  vax_start,
  qtile = 0.95,
  fill = "#1B9E77",
  alpha = 0.25
)
```

## Arguments

- p:

  A plot from
  [`top_group_plot()`](https://magpie-sfu.github.io/newmo/reference/top_group_plot.md).

- sens:

  Result of
  [`run_sensitivity()`](https://magpie-sfu.github.io/newmo/reference/run_sensitivity.md).

- key:

  Id (`loss` column in `sens`) of the fit to show the band for.

- vax_start:

  The `vax_start` used in `p`.

- qtile:

  Quantile range of the band.

- fill, alpha:

  Fill colour and transparency of the band.

## Value

The plot with the band underneath all other layers.

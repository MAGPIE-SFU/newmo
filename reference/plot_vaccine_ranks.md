# Rank plot of invented vaccines, with licensed vaccines for comparison

Exhaustive vaccines ranked by the IPD they leave (rank 1 = lowest IPD):
the best `n_top`, then `n_sample` ranks evenly spread over the rest.
Vaccines from `other_methods` are placed at the rank they would have
among the exhaustive vaccines. Licensed vaccines are dashed horizontal
lines.

## Usage

``` r
plot_vaccine_ranks(
  vaccine_data,
  known_ipd,
  n_top = 100,
  n_sample = 1000,
  other_methods = "greedy",
  method_labels = character(0),
  ipd_scale = 1e+05,
  ipd_label = "IPD per 100,000"
)
```

## Arguments

- vaccine_data:

  Invented vaccines, e.g. from
  [`design_vaccines()`](https://magpie-sfu.github.io/newmo/reference/design_vaccines.md).

- known_ipd:

  Named numeric vector: IPD score per licensed vaccine.

- n_top:

  Number of best exhaustive vaccines shown individually.

- n_sample:

  Number of evenly spaced ranks shown from the remaining ones.

- other_methods:

  Methods placed among the exhaustive ranks.

- method_labels:

  Named character vector: display name per method.

- ipd_scale:

  IPD is multiplied by this for display.

- ipd_label:

  y-axis label. Keep in step with `ipd_scale`.

## Value

A ggplot.

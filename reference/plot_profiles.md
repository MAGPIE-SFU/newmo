# Plot the 1D loss profiles

Rows = loss, columns = parameter, one line per seed.

## Usage

``` r
plot_profiles(
  prof,
  losses = names(loss_functions),
  relative_to = "best_seed",
  ymax = 10
)
```

## Arguments

- prof:

  Profiles from [`profile_best_fits()`](profile_best_fits.md).

- losses:

  Loss functions, in plotting order.

- relative_to:

  "best_seed": ratio to the lowest fitted loss across seeds (shows which
  seed fit best). "own": ratio to each seed's own fitted loss (every
  curve touches 1).

- ymax:

  Cap on the loss ratio, so one bad region doesn't flatten the rest.

## Value

A ggplot.

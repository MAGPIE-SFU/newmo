# Plot the loss landscape

Loss against each parameter, one row per loss function.

## Usage

``` r
plot_landscape(screen, losses = names(loss_functions))
```

## Arguments

- screen:

  Screen from [`screen_landscape()`](screen_landscape.md).

- losses:

  Loss functions, in plotting order.

## Value

A ggplot.

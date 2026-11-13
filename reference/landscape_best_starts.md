# Best screened start per loss

Best screened start per loss

## Usage

``` r
landscape_best_starts(screen, losses = names(loss_functions))
```

## Arguments

- screen:

  Screen from [`screen_landscape()`](screen_landscape.md).

- losses:

  Loss functions.

## Value

One row per loss: the best start, its loss, the default start's loss and
the parameters.

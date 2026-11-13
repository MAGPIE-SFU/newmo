# Pick the best fits across seeds

Takes the best fit per seed x loss (excluding the round_00 starting
points), simulates every candidate, evaluates it on every loss function,
and scores it as the sum over loss functions of (loss / best loss across
all candidates). Lower is better; the minimum possible score is the
number of loss functions.

## Usage

``` r
select_best_fits(seq_fits, samodel, xt, tswitch, nyears, use_switch = TRUE)
```

## Arguments

- seq_fits:

  List of results from
  [`run_sequential_fits()`](https://magpie-sfu.github.io/newmo/reference/run_sequential_fits.md),
  one per seed.

- samodel:

  The model structure, e.g. from
  [`build_samodel()`](https://magpie-sfu.github.io/newmo/reference/build_samodel.md).

- xt:

  Grouped x(t) to fit to.

- tswitch:

  Time of the vaccine switch.

- nyears:

  Number of (final) years to fit.

- use_switch:

  Use the vaccine switch.

## Value

A list with `all_parameters` and `all_summary` (all fits of all seeds),
`all_best` (best fit per seed x loss, with key and score) and
`overall_best` (the one row with the lowest score).

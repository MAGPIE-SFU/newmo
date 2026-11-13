# Empirical Bayes shrinkage of invasiveness estimates

Shrinks log invasiveness towards the overall mean, using a
DerSimonian-Laird estimate of the between-serotype variance.

## Usage

``` r
shrink_invasiveness(tab)
```

## Arguments

- tab:

  Table with columns Serotype, invasiveness, invasiveness_lower and
  invasiveness_upper.

## Value

`tab` with added columns logmean, logsd and v_shrunk, and attributes
overall_mean, mu and tau.

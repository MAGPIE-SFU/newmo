# Grouped x(t) with invasiveness uncertainty

Repeats the carriage reweighting `n_reps` times, each time with one
invasiveness draw per serotype (from the log-normal in the invasiveness
lookup), builds x(t) for the given groups, and combines the repetitions.

## Usage

``` r
reweight_carriage_uncertainty(
  meta_clean,
  inv,
  group_lookup,
  tauA,
  tauU,
  k,
  lambda = 0.01,
  n_reps = 1000,
  seed = 1
)
```

## Arguments

- meta_clean:

  Cleaned metadata, e.g. from
  [`clean_metadata()`](https://magpie-sfu.github.io/newmo/reference/clean_metadata.md).

- inv:

  Invasiveness lookup, e.g. from
  [`build_invasiveness_lookup()`](https://magpie-sfu.github.io/newmo/reference/build_invasiveness_lookup.md).

- group_lookup:

  Group definitions, e.g. `reweight_carriage()$group_lookup`.

- tauA:

  Weight for adult IPD, should be very small.

- tauU:

  Weight for unknown-age IPD.

- k:

  Dirichlet prior strength, in equivalent child carriage swabs.

- lambda:

  Fraction of the prior mean that is uniform.

- n_reps:

  Number of repetitions for the invasiveness uncertainty.

- seed:

  Random seed, set at the start as in the original script.

## Value

x(t) with columns year, group, serotype, vaccine_type, frac, frac_sd.

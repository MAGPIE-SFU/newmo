# Simulate projected IPD with uncertainty from invasiveness estimates

Scales a simulated carriage composition to `n_carriage` per year, splits
each group over its serotypes, and multiplies by child invasiveness
drawn `n_reps` times from its log-normal distribution. Simulated
carriage is held fixed. `ipd_expected` is a Poisson mean, not a realised
case count.

## Usage

``` r
simulate_projected_ipd(
  sim,
  samodel,
  inv,
  group_lookup,
  meta_clean,
  non_typeable,
  n_reps = 1000L,
  n_carriage = 1e+05,
  seed = 1
)
```

## Arguments

- sim:

  A simulation from
  [`run_model()`](https://magpie-sfu.github.io/newmo/reference/run_model.md)
  /
  [`simulate_nfds()`](https://magpie-sfu.github.io/newmo/reference/simulate_nfds.md):
  rows = groups, columns year.1 ... year.n.

- samodel:

  The model structure the simulation was run with.

- inv:

  Invasiveness lookup, e.g. from
  [`build_invasiveness_lookup()`](https://magpie-sfu.github.io/newmo/reference/build_invasiveness_lookup.md).

- group_lookup:

  Group definitions, e.g. `reweight_carriage()$group_lookup`.

- meta_clean:

  Cleaned metadata, e.g. from
  [`clean_metadata()`](https://magpie-sfu.github.io/newmo/reference/clean_metadata.md).

- non_typeable:

  Serotype labels with zero invasiveness.

- n_reps:

  Number of Monte Carlo draws.

- n_carriage:

  Simulated carriage population per year.

- seed:

  Random seed, set at the start as in the original script.

## Value

A list with `ipd_summary`, `ipd_yearly`, `ipd_yearly_draws` and
`ipd_serotype_summary`.

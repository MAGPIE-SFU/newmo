# Reweight carriage composition and define model groups

Combines child carriage counts with IPD-derived carriage estimates
(child, adult and unknown-age IPD, divided by invasiveness) into a
synthetic child carriage population per year, defines the model groups
from the pre-vaccine composition, and builds x(t).

## Usage

``` r
reweight_carriage(
  meta_clean,
  inv,
  tauA,
  tauU,
  k,
  prevax_years,
  min_sero_frac,
  min_sg_frac
)
```

## Arguments

- meta_clean:

  Cleaned metadata, e.g. from
  [`clean_metadata()`](https://magpie-sfu.github.io/newmo/reference/clean_metadata.md).

- inv:

  Invasiveness lookup, e.g. from
  [`build_invasiveness_lookup()`](https://magpie-sfu.github.io/newmo/reference/build_invasiveness_lookup.md).

- tauA:

  Weight for adult IPD, should be very small.

- tauU:

  Weight for unknown-age IPD (proportion of the dataset that is of
  unknown age).

- k:

  Dirichlet prior strength: how many child carriage swabs the IPD-based
  estimate is worth.

- prevax_years:

  Pre-vaccine years used to define the groups.

- min_sero_frac:

  NVT serotypes below this share of pre-vaccine synthetic carriage go to
  "rare_NVT".

- min_sg_frac:

  Serotype-GPSC pairs below this share go to "\_other".

## Value

A list with `group_lookup`, `meta_clean_groups` (meta_clean with a group
column) and `final_yearly` (the synthetic carriage composition per
serotype-GPSC pair, for checks).

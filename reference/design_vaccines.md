# Invent vaccines

Method 1 ("ipdfrac", "mostinvasive"): the top `target_size` serotypes by
their contribution to IPD in the inferred pre-vaccine population, or by
child invasiveness. Method 2 ("greedy"): starting from `base_vaccine`,
add the candidate that lowers IPD most, one at a time, up to
`target_size`. Method 3 ("exhaustive"): every combination of candidates
added to `base_vaccine` up to `target_size`.

## Usage

``` r
design_vaccines(
  modelpars,
  samodel,
  xt,
  inv,
  group_lookup,
  base_vaccine,
  prevax_years = 2000:2009,
  target_size = 13,
  n_candidates = 18,
  years = 6:10,
  n_cores = 1
)
```

## Arguments

- modelpars:

  Model parameters for [`simulate_nfds()`](simulate_nfds.md) (vaccines
  are set per score).

- samodel:

  The model structure.

- xt:

  Grouped x(t), with a serotype column.

- inv:

  Invasiveness lookup.

- group_lookup:

  Group definitions.

- base_vaccine:

  Serotypes every greedy / exhaustive vaccine starts from (PCV7).

- prevax_years:

  Years for the inferred pre-vaccine composition.

- target_size:

  Valency of the invented vaccines.

- n_candidates:

  Number of most invasive non-base serotypes to choose from.

- years:

  Simulation years to average the IPD over.

- n_cores:

  Cores for the exhaustive search (\> 1 uses mclapply, Linux/macOS
  only).

## Value

One row per vaccine: ipd, vaccine (string), method.

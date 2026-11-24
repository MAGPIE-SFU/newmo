# Invent vaccines

Method 1 ("ipdfrac", "mostinvasive"): the top `direct_size` serotypes by
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
  candidate_vaccines,
  prevax_years = 2000:2009,
  direct_size = 13,
  target_size = 15,
  n_extra = 5,
  years = 6:10,
  n_cores = 1
)
```

## Arguments

- modelpars:

  Model parameters for
  [`simulate_nfds()`](https://magpie-sfu.github.io/newmo/reference/simulate_nfds.md)
  (vaccines are set per score).

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

- candidate_vaccines:

  Named list of vaccines whose serotypes are candidates, in priority
  order (see
  [`vaccine_candidate_data()`](https://magpie-sfu.github.io/newmo/reference/vaccine_candidate_data.md)).

- prevax_years:

  Years for the inferred pre-vaccine composition.

- direct_size:

  Valency of the method 1 vaccines.

- target_size:

  Valency of the greedy and exhaustive vaccines.

- n_extra:

  Number of top NVT serotypes (by IPD contribution, and by invasiveness)
  added to the candidates.

- years:

  Simulation years to average the IPD over.

- n_cores:

  Cores for the exhaustive search (\> 1 uses mclapply, Linux/macOS
  only).

## Value

One row per vaccine: ipd, vaccine (string), method.

## Details

Candidates for methods 2 and 3: every serotype in one of
`candidate_vaccines`, plus the NVT serotypes in the top `n_extra` by IPD
contribution or by child invasiveness, minus `base_vaccine`.

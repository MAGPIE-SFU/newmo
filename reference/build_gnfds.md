# Build the grouped G matrix for the NFDS model

Groups the isolate gene presence/absence matrix by model group, computes
locus frequencies per year from x(t), computes several locus stability
statistics, keeps intermediate-frequency loci and orders them by
`logitdiff`.

## Usage

``` r
build_gnfds(
  G,
  meta_clean,
  xt,
  pre_years,
  post_years,
  all_years,
  min_freq = 0.05,
  max_freq = 0.95
)
```

## Arguments

- G:

  Raw presence/absence matrix, one row per isolate (row names are sample
  accessions), one column per locus.

- meta_clean:

  Metadata with a `group` column, e.g.
  `reweight_carriage()$meta_clean_groups`.

- xt:

  Grouped x(t), e.g. from
  [`reweight_carriage_uncertainty()`](https://magpie-sfu.github.io/newmo/reference/reweight_carriage_uncertainty.md).

- pre_years:

  Pre-vaccine years (the "pre2009" frequency).

- post_years:

  Post-vaccine years (the "post2011" frequency).

- all_years:

  All years used for the year-to-year variance.

- min_freq, max_freq:

  Loci with a pre-vaccine frequency outside this range are left out of
  the model.

## Value

A list with `G` (the grouped, filtered and ordered G for the NFDS model)
and `locus_by_year` (locus frequencies per year and the statistics, with
an `in_model` column).

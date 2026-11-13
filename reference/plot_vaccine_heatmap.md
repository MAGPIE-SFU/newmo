# Vaccine heatmap

One row per vaccine (licensed first, then the invented candidates), one
column per serotype, plus a column with the IPD each vaccine leaves.
Serotype columns: one block per licensed vaccine (its serotypes not in
any earlier vaccine), then everything else.

## Usage

``` r
plot_vaccine_heatmap(
  known_serotypes,
  known_ipd,
  vaccine_data,
  all_serotypes,
  n_exhaustive = 6,
  method_labels = character(0),
  ipd_scale = 1e+05,
  ipd_label = "IPD per\n100,000"
)
```

## Arguments

- known_serotypes:

  Named list of licensed vaccines (serotype vectors), in display order.

- known_ipd:

  Named numeric vector: IPD score per licensed vaccine.

- vaccine_data:

  Invented vaccines, e.g. from
  [`design_vaccines()`](https://magpie-sfu.github.io/newmo/reference/design_vaccines.md):
  columns ipd, vaccine (string with "-" or "," separators), method.

- all_serotypes:

  All serotypes in the model (e.g. samodel\$st), so serotypes in no
  vaccine get a column too.

- n_exhaustive:

  Number of exhaustive vaccines to show, lowest IPD first. All
  non-exhaustive vaccines are always shown.

- method_labels:

  Named character vector: display name per method. Methods not listed
  keep their own name.

- ipd_scale:

  IPD is multiplied by this for display.

- ipd_label:

  Header of the IPD column. Keep in step with `ipd_scale`.

## Value

A list with `plot` (the combined plot), `width` and `height` (suggested
size in inches).

# Invasiveness and IPD contribution of each serotype, by vaccine

Illustrates which serotypes are in which vaccines, and their (1) child
invasiveness and (2) contribution to pre-vaccine IPD. Shows why the
candidates for the vaccine design are chosen.

## Usage

``` r
plot_vaccine_candidates(
  compodat,
  vax_cols = c(PCV7 = "#440154", PCV13 = "#414487", PCV15 = "#2a788e", PCV20 = "#22a884",
    PCV23 = "#7ad151", NVT = "grey70")
)
```

## Arguments

- compodat:

  From
  [`vaccine_candidate_data()`](https://magpie-sfu.github.io/newmo/reference/vaccine_candidate_data.md).

- vax_cols:

  Named colours per vaccine (and "NVT").

## Value

A patchwork of two plots.

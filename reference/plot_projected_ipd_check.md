# Observed vs predicted IPD serotype composition

Observed vs predicted IPD serotype composition

## Usage

``` r
plot_projected_ipd_check(
  ipd_serotype_summary,
  meta_clean,
  first_sim_calendar_year
)
```

## Arguments

- ipd_serotype_summary:

  From [`simulate_projected_ipd()`](simulate_projected_ipd.md).

- meta_clean:

  Cleaned metadata (for the observed IPD).

- first_sim_calendar_year:

  Calendar year that corresponds to year.1 in the simulation.

## Value

A ggplot.

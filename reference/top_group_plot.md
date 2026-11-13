# Plot the most common ST-SC groups, model and data together

Strip colour: VT = red, NVT = blue.

## Usage

``` r
top_group_plot(
  vax_start = 2010,
  prevax_years = 2000:2009,
  n_groups = 12,
  samodel,
  xt,
  sims,
  drop_vt = TRUE,
  palette = "Dark2",
  vaccine
)
```

## Arguments

- vax_start:

  First vaccine year.

- prevax_years:

  Pre-vaccine years (for the x0 segment).

- n_groups:

  Number of top ST-SC groups to show.

- samodel:

  The model structure.

- xt:

  Grouped x(t).

- sims:

  A simulation, or a named list of simulations.

- drop_vt:

  Drop vaccine-type groups (frequencies are NOT renormalised).

- palette:

  RColorBrewer name or viridis option for the simulations.

- vaccine:

  Vaccine serotypes (defines VT vs NVT).

## Value

A ggplot.

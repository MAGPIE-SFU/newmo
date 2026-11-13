# Plot the most common serotypes, model and data together

Strip colour: VT = red, NVT = blue.

## Usage

``` r
top_sero_plot(
  vax_start = 2010,
  prevax_years = 2000:2009,
  n_serotypes = 12,
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

- n_serotypes:

  Number of top serotypes to show.

- samodel:

  The model structure.

- xt:

  Grouped x(t).

- sims:

  A simulation, or a named list of simulations.

- drop_vt:

  Drop vaccine-type serotypes (frequencies are NOT renormalised).

- palette:

  RColorBrewer name or viridis option for the simulations.

- vaccine:

  Vaccine serotypes (defines VT vs NVT).

## Value

A ggplot.

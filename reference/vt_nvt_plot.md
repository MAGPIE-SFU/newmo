# Plot VT and NVT totals, model and data together

Data = the synthetic carriage sample x(t). VT and NVT go in separate
panels (colour is used for the simulations).

## Usage

``` r
vt_nvt_plot(
  vax_start = 2010,
  prevax_years = 2000:2009,
  samodel,
  xt,
  sims,
  palette = "Dark2",
  vaccine
)
```

## Arguments

- vax_start:

  First vaccine year.

- prevax_years:

  Pre-vaccine years (for the x0 segment).

- samodel:

  The model structure.

- xt:

  Grouped x(t).

- sims:

  A simulation, or a named list of simulations.

- palette:

  RColorBrewer name or viridis option for the simulations.

- vaccine:

  Vaccine serotypes (defines VT vs NVT).

## Value

A ggplot.

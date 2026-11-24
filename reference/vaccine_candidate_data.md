# Pre-vaccine serotype composition, invasiveness and vaccine membership

Serotype-level carriage fraction in the inferred pre-vaccine population,
child invasiveness, their product (the IPD contribution), and the first
vaccine in `vaccines` that contains each serotype ("NVT" if none).

## Usage

``` r
vaccine_candidate_data(xt, inv, vaccines, prevax_years = 2000:2009)
```

## Arguments

- xt:

  Grouped x(t), with a serotype column.

- inv:

  Invasiveness lookup (serotype, v_child).

- vaccines:

  Named list of vaccines, in priority order (e.g. PCV7, PCV13, PCV15,
  PCV20, PCV23): a serotype gets the first vaccine it is in.

- prevax_years:

  Years for the pre-vaccine composition.

## Value

One row per serotype: serotype, frac, v_child, ipdfrac, ipdshare, vax.

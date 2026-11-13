# Turn a loss function into an objective of the parameter vector alone

Turn a loss function into an objective of the parameter vector alone

## Usage

``` r
make_loss(data, sd, loss_fn, use_switch, samodel, nyears, tswitch)
```

## Arguments

- data, sd:

  Data and SD matrices, from [`load_fit_data()`](load_fit_data.md).

- loss_fn:

  One of [loss_functions](loss_functions.md).

- use_switch:

  Use the vaccine switch at `tswitch`.

- samodel:

  The model structure, e.g. from [`build_samodel()`](build_samodel.md).

- nyears:

  Number of years to simulate.

- tswitch:

  Time of the vaccine switch (used if `use_switch`).

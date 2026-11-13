# Score a vaccine by the IPD it leaves

Simulates with `vseros` as the only vaccine (vaccine1 = vaccine2) and
returns the total IPD, averaged over `years`.

## Usage

``` r
score_vaccine(modelpars, vseros, samodel, inv, group_lookup, years = 6)
```

## Arguments

- modelpars:

  Model parameters for [`simulate_nfds()`](simulate_nfds.md).

- vseros:

  Serotypes in the vaccine.

- samodel:

  The model structure.

- inv:

  Invasiveness lookup (serotype, v_child).

- group_lookup:

  Group definitions (group, serotype).

- years:

  Simulation years to average over, e.g. 6:10.

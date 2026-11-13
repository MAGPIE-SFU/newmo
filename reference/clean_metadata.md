# Clean isolate metadata

Harmonises serotype names, keeps one country, and derives age in years,
age group, collection type and vaccine type (PCV7 / PCV13 / NVT).

## Usage

``` r
clean_metadata(
  meta,
  serotype_renames,
  keep_country = "South Africa",
  child_cutoff = 6
)
```

## Arguments

- meta:

  Raw metadata, as read from metadata.csv.

- serotype_renames:

  Named character vector, old label = new label. Renamed isolates are
  treated as the new serotype from here on.

- keep_country:

  Country to keep.

- child_cutoff:

  Isolates with an age below this (in years) are "child".

## Value

One row per isolate.

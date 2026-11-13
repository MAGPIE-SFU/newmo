# Build the invasiveness lookup with intervals

For every serotype in `meta_clean`, finds child and adult invasiveness
estimates: from hard-coded substitutions, exact matches, serogroups or
the overall mean, in that order of priority (see the comments in the
code).

## Usage

``` r
build_invasiveness_lookup(
  meta_clean,
  inv_children,
  inv_adults,
  serotype_substitutions,
  non_typeable
)
```

## Arguments

- meta_clean:

  Cleaned metadata, e.g. from [`clean_metadata()`](clean_metadata.md).

- inv_children, inv_adults:

  Invasiveness tables for children and adults, with columns Serotype,
  invasiveness, invasiveness_lower, invasiveness_upper.

- serotype_substitutions:

  Named character vector: serotype = donor. These take priority over
  exact matches and serogroup averages.

- non_typeable:

  Serotype labels that get no invasiveness estimate.

## Value

One row per serotype, with v_child\_\* and v_adult\_\* columns.

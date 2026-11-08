# Within-between (Mundlak) decomposition of a numeric variable

Splits a numeric vector into a between-group part (the group mean,
broadcast back to each row) and a within-group part (the deviation from
that mean). Used to separate the within-unit component of an exposure –
which is robust to any unit-level confounder, observed or not – from the
between-unit component.

## Usage

``` r
pollin_within_between(x, group)
```

## Arguments

- x:

  Numeric vector.

- group:

  Grouping factor of the same length as `x`.

## Value

A list with elements `within` and `between`, each the length of `x`.

## Examples

``` r
wb <- pollin_within_between(c(1, 3, 10, 12), c("a", "a", "b", "b"))
wb$between
#> [1]  2  2 11 11
wb$within
#> [1] -1  1 -1  1
```

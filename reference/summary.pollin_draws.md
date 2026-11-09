# Summarise draws of the estimand

The median and percentile interval of the estimand draws, per column.
[`as.matrix()`](https://rdrr.io/r/base/matrix.html) returns the matrix
of estimand draws.

## Usage

``` r
# S3 method for class 'pollin_draws'
summary(object, conf_level = 0.95, ...)

# S3 method for class 'pollin_draws'
as.matrix(x, ...)
```

## Arguments

- object, x:

  A
  [`pollin_draws()`](https://max578.github.io/pollinatoR/reference/pollin_draws.md)
  object.

- conf_level:

  Level of the percentile interval.

- ...:

  Not used.

## Value

A data frame with one row per column of the draws and the columns
`response`, `level` (the column name: a `by` level, `"overall"` or
`"population"`), `median`, `conf.low` and `conf.high` (the percentile
interval), `family`, `scale`, `basis`, `method` and `n`.
[`as.matrix()`](https://rdrr.io/r/base/matrix.html) returns the
`estimate` matrix of the draws.

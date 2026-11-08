# Validate the structure and internal consistency of a study's data

Runs aggregate structural checks against a declared table: missingness
per declared column, expected factor levels, and any user-supplied
internal-consistency identities (for example that a total equals a sum
of parts, or that a ratio equals its definition). Returns counts and
pass/warn flags only – no raw rows – so it is safe on sensitive data.
[`pollin_data()`](https://max578.github.io/pollinatoR/reference/pollin_data.md)
runs these checks without levels or identities when it builds its
object.

## Usage

``` r
pollin_validate(x, expect_levels = NULL, identities = NULL)
```

## Arguments

- x:

  A
  [`pollin_data()`](https://max578.github.io/pollinatoR/reference/pollin_data.md)
  object, or a
  [`pollin_spec()`](https://max578.github.io/pollinatoR/reference/pollin_spec.md)
  object.

- expect_levels:

  Optional named list mapping column names to the expected set of factor
  levels.

- identities:

  Optional named list of functions; each takes the data frame and
  returns a logical vector that is `TRUE` where the identity holds. The
  number of violations is reported.

## Value

An object of class `pollin_validate` carrying a tidy `table`.

## Examples

``` r
df <- data.frame(y = c(1, 2, NA), trt = c("Cage", "Open", "Open"),
                 Site = c("a", "a", "b"))
pd <- pollin_data(df, response = "y", treatment = "trt", unit = "Site")
pollin_validate(pd, expect_levels = list(trt = c("Cage", "Open")))
#> <pollin_validate>
#>       check status    detail
#>        rows   INFO         3
#>        cols   INFO         3
#>        NA:y   WARN         1
#>  levels:trt   PASS Cage,Open
```

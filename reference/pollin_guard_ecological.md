# Guard against reading a between-unit slope as a within-unit one

An exposure slope estimated across sites (between units) is an
ecological association: it need not hold for plants within a site, where
unit-level confounders are held fixed. This guard splits the exposure
into its between-unit and within-unit parts with
[`pollin_within_between()`](https://max578.github.io/pollinatoR/reference/pollin_within_between.md)
and fits `response ~ within + between` by least squares, reporting the
pooled, within-unit and between-unit slopes. It warns when

## Usage

``` r
pollin_guard_ecological(
  x,
  exposure,
  response = NULL,
  unit = NULL,
  alpha = 0.05
)
```

## Arguments

- x:

  A data frame, or a
  [`pollin_data()`](https://max578.github.io/pollinatoR/reference/pollin_data.md)
  object.

- exposure:

  Name of the numeric exposure column (for example insect abundance).

- response:

  Name of the numeric response column. With a
  [`pollin_data()`](https://max578.github.io/pollinatoR/reference/pollin_data.md)
  object it may be left `NULL` when one response is declared.

- unit:

  The unit column(s) the exposure is split by. Required for a data frame
  (several columns are combined); with a
  [`pollin_data()`](https://max578.github.io/pollinatoR/reference/pollin_data.md)
  object, one of the unit columns, taken within its outer units, and
  `NULL` (the default) uses the outermost unit.

- alpha:

  Significance level for the third check.

## Value

An object of class `pollin_ecological`: a list holding the slope table
`slopes` (component, slope, standard error, p-value), the number of
units `n_units`, the `flags` raised (a character vector, empty when
none) and the column names used.

## Details

- the exposure does not vary within units, so only the between-unit
  slope exists;

- only the between-unit slope is distinguishable from zero at `alpha`;

- the within-unit and between-unit slopes are both distinguishable from
  zero at `alpha` and disagree in sign.

The guard works on the rows given; subset to one arm beforehand when the
slope of interest is within that arm.

## See also

[`pollin_moderation()`](https://max578.github.io/pollinatoR/reference/pollin_moderation.md)
for the model-based within/between split.

## Examples

``` r
set.seed(1)
site <- rep(sprintf("s%02d", 1:10), each = 8)
site_bees <- rep(rnorm(10, 10, 3), each = 8)
bees <- site_bees + rnorm(80)
# Sites with more bees yield more, but plants within a site do not respond
pods <- 20 + 2 * site_bees + rnorm(80, 0, 2)
df <- data.frame(pods = pods, bees = bees, Site = site)
pollin_guard_ecological(df, "bees", "pods", unit = "Site")
#> Warning: Only the between-unit slope is distinguishable from zero (p = 3.4e-31 between, 0.0518 within); the association is between units, not within them.
#> <pollin_ecological>  pods ~ bees split by Site (10 units)
#>  component   slope std.error   p.value
#>     pooled  1.5788   0.13100 1.715e-19
#>     within -0.4887   0.24734 5.178e-02
#>    between  1.9108   0.09913 3.404e-31
#>   ! Only the between-unit slope is distinguishable from zero (p = 3.4e-31 between, 0.0518 within); the association is between units, not within them.
```

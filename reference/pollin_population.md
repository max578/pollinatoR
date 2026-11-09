# Declare a population of strata for an estimand

A population made of the levels (strata) of one column of the data, for
example distance bands, regions or varieties, each weighted by its area.
Passed to
[`pollin_estimand()`](https://max578.github.io/pollinatoR/reference/pollin_estimand.md)
as `population`, it carries the estimand from the trial to the
population: within each stratum the expected response under each arm
comes from the fitted model, and the population arm means are the
area-weighted sums of the stratum arm means:

## Usage

``` r
pollin_population(by, weights, label = NULL)
```

## Arguments

- by:

  Name of the column whose levels are the strata.

- weights:

  Named numeric vector of stratum areas, one positive value per level of
  `by`, named by level. Any unit of area serves: the weights are
  normalised to sum to 1 and the values given are kept as `area`.

- label:

  Optional short description of the population, for printing.

## Value

An object of class `pollin_population`: a list holding `by`, the
normalised `weights`, the `area` values as given, and `label`.

## Details

`high_pop = sum(w * high_stratum)` and `low_pop = sum(w * low_stratum)`,

with the weights `w` normalised to sum to 1. The estimand is formed from
the two population means, so on the ratio scale the population loss
share is `1 - sum(w * low_stratum) / sum(w * high_stratum)`: a ratio of
population means, not a weighted mean of the stratum loss shares. The
two agree when every stratum has the same expected response in the
`high` arm, or the same loss share; in general they differ.

The model must be fitted with the strata as its `by` column,
`pollin_fit(estimand, by = <by>)`, and the population estimand is read
with
[`pollin_draws()`](https://max578.github.io/pollinatoR/reference/pollin_draws.md).

## See also

[`pollin_estimand()`](https://max578.github.io/pollinatoR/reference/pollin_estimand.md),
[`pollin_draws()`](https://max578.github.io/pollinatoR/reference/pollin_draws.md).

## Examples

``` r
pollin_population(
  "Region",
  c(north = 1200, south = 800),
  label = "two growing regions"
)
#> <pollin_population>
#>   label   : two growing regions 
#>   strata  : 2 levels of Region 
#>   weights : north 0.6, south 0.4 (area, normalised)
```

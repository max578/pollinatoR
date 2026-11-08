# Summarise yield-component trade-offs (compensation)

Canola yield is the product of compensating components (pod number,
seeds per pod, seed size): a deficit in one is often offset by another,
which a component-by-component analysis can hide. This function returns
the correlation structure among the components across units and flags
negative associations as candidate compensation. It is a descriptive
companion to the per-component contribution models, not a substitute for
a joint multivariate model.

## Usage

``` r
pollin_compensation(components, method = "spearman")
```

## Arguments

- components:

  A numeric matrix or data frame, one column per yield component, one
  row per unit; or a
  [`pollin_data()`](https://max578.github.io/pollinatoR/reference/pollin_data.md)
  object, whose response columns are taken as the components (at least
  two).

- method:

  Correlation method passed to
  [`stats::cor()`](https://rdrr.io/r/stats/cor.html).

## Value

An object of class `pollin_compensation` with the correlation matrix and
a table of negatively-correlated component pairs.

## Examples

``` r
m <- data.frame(pods = c(10, 8, 6), seeds_per_pod = c(20, 24, 28))
pollin_compensation(m)
#> <pollin_compensation>  ( spearman correlations )
#>              a    b cor compensation
#>  seeds_per_pod pods  -1         TRUE
#>   1 negatively-correlated pair -> candidate compensation

# From a pollin_data object: its responses are the components
df <- data.frame(
  pods = c(10, 8, 6, 9, 7, 5),
  seeds_per_pod = c(20, 24, 28, 21, 25, 27),
  trt = rep(c("Cage", "Open"), 3),
  Site = rep(c("a", "b", "c"), each = 2)
)
pd <- pollin_data(df, c("pods", "seeds_per_pod"), "trt", unit = "Site")
pollin_compensation(pd)
#> <pollin_compensation>  ( spearman correlations )
#>              a    b        cor compensation
#>  seeds_per_pod pods -0.9428571         TRUE
#>   1 negatively-correlated pair -> candidate compensation
```

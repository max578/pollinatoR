# List the assumptions an estimand rests on

Returns the identification assumptions for the estimand's layout, family
and scale. Each is a statement about the experiment that the data cannot
confirm on their own.

## Usage

``` r
pollin_assumptions(x)
```

## Arguments

- x:

  A
  [`pollin_estimand()`](https://max578.github.io/pollinatoR/reference/pollin_estimand.md)
  object.

## Value

A data frame with columns `assumption` and `statement`.

## Examples

``` r
df <- data.frame(
  pods = rpois(48, 30),
  trt = rep(c("Cage", "Open"), 24),
  Site = rep(letters[1:4], each = 12)
)
pd <- pollin_data(df, "pods", "trt", unit = "Site")
ds <- pollin_design(pd, c(excluded = "Cage", open = "Open"), "multi_site")
est <- pollin_estimand(ds, "dependence", "ratio", "loss")
pollin_assumptions(est)
#>        assumption
#> 1     consistency
#> 2   excludability
#> 3 no_interference
#> 4      assignment
#> 5  positive_means
#> 6      population
#>                                                                                                                      statement
#> 1                                                  Each arm is one well-defined treatment, applied the same way on every unit.
#> 2 The exclusion device changes the response only by excluding pollinators (no cage effect on microclimate, pests or handling).
#> 3   One unit's treatment does not change another unit's response, for example by displacing pollinators to neighbouring units.
#> 4                                                    Arms were randomised within each site; the contrast is read within sites.
#> 5                                                     The expected response is positive in both arms, so the ratio is defined.
#> 6               The estimand refers to the units in the trial; carrying it to other fields or seasons is a further assumption.
```

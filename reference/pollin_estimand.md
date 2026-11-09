# Declare the estimand of a pollination experiment

Names the contrast, the scale it is reported on and the basis of any
share, for a declared
[`pollin_design()`](https://max578.github.io/pollinatoR/reference/pollin_design.md).
Two estimand families are available:

## Usage

``` r
pollin_estimand(design, family, scale, basis, population = "trial")
```

## Arguments

- design:

  A
  [`pollin_design()`](https://max578.github.io/pollinatoR/reference/pollin_design.md)
  object.

- family:

  `"deficit"` or `"dependence"`. There is no default.

- scale:

  One of `"difference"`, `"ratio"`, `"rii"` or `"lnrr"`. There is no
  default.

- basis:

  `"loss"` or `"gain"`. There is no default.

- population:

  The population the estimand refers to: `"trial"` (the default, the
  units in the experiment) or a
  [`pollin_population()`](https://max578.github.io/pollinatoR/reference/pollin_population.md)
  object, strata of one column weighted by area. A population estimand
  is fitted with the strata as `by` and read with
  [`pollin_draws()`](https://max578.github.io/pollinatoR/reference/pollin_draws.md).

## Value

An object of class `pollin_estimand`: a list holding `design`, `family`,
`scale`, `basis`, `population`, the `contrast` (the treatment levels of
the `high` and `low` arms) and the `assumptions` table.

## Details

- `"deficit"` (E1): supplemented versus open, the yield gap left by the
  pollination the crop received.

- `"dependence"` (E2): excluded versus the design's reference arm (open,
  or supplemented), the yield owed to pollinators.

Writing `high` for the expected response in the higher arm (open or
supplemented) and `low` for the lower one, the scales are:

|  |  |
|----|----|
| Scale | Value |
| `"difference"` | `high - low` |
| `"ratio"`, `basis = "loss"` | `1 - low / high`, the share of `high` lost without the contrast |
| `"ratio"`, `basis = "gain"` | `high / low - 1`, the gain over `low` |
| `"rii"` | `(high - low) / (high + low)`, equal to `D / (2 - D)` for the loss share `D` |
| `"lnrr"` | `log(high / low)` |

A loss share and a gain share answer different questions: an excluded
arm yielding 75% of the open arm is a 25% loss share and a 33% gain.
`basis` has no default; every analysis states it. It changes the value
only on the ratio scale and is recorded on every scale.

## See also

[`pollin_assumptions()`](https://max578.github.io/pollinatoR/reference/pollin_assumptions.md)
for the assumptions the estimand rests on.

## Examples

``` r
df <- data.frame(
  pods = rpois(48, 30),
  trt = rep(c("Cage", "Open"), 24),
  Site = rep(letters[1:4], each = 12)
)
pd <- pollin_data(df, "pods", "trt", unit = "Site")
ds <- pollin_design(pd, c(excluded = "Cage", open = "Open"), "multi_site")
pollin_estimand(ds, family = "dependence", scale = "ratio", basis = "loss")
#> <pollin_estimand>
#>   family     : dependence (E2) 
#>   contrast   : Open vs Cage 
#>   scale      : ratio 
#>   basis      : loss 
#>   population : trial 
#>   layout     : multi_site 
#>   assumptions: 6 ; see `pollin_assumptions()`
```

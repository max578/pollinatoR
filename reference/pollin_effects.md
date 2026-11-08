# Estimate the effect on the estimand's scale

Reads the estimand's contrast from each fitted model by G-computation:
every row of the model data is predicted under the `high` arm and under
the `low` arm, the predictions are averaged, and the effect is formed
from the two averaged predictions. A ratio-type scale is a ratio of
averaged predictions, never an average of unit ratios. With
`rr = high / low`:

## Usage

``` r
pollin_effects(
  fit,
  by = fit$by,
  interval = c("delta", "bootstrap"),
  n_boot = 999L,
  seed = NULL,
  conf_level = 0.95
)
```

## Arguments

- fit:

  A
  [`pollin_fit()`](https://max578.github.io/pollinatoR/reference/pollin_fit.md)
  object.

- by:

  The column the fit was made with, `pollin_fit(..., by = )`. Defaults
  to the fit's `by`, giving one contrast per level from the fixed
  treatment-by-level effects, each averaged over the rows at that level.
  `NULL` gives one contrast averaged over all rows. A column the fit was
  not made with is an error: refit with `pollin_fit(..., by = )`.

- interval:

  `"delta"` (the default) or `"bootstrap"`.

- n_boot:

  Number of parametric bootstrap refits, used when
  `interval = "bootstrap"`.

- seed:

  Optional integer seed for the bootstrap simulation. Each response is
  simulated from this seed, and the session's random number state is
  restored afterwards.

- conf_level:

  Confidence level of the interval.

## Value

An object of class `pollin_effects`: a list holding `table`, with one
row per response (and per `by` level) and the columns `response`,
`contrast`, any `by` column, `estimate`, `std.error`, `conf.low`,
`conf.high`, `p.value`, `family` (the estimand family), `scale`,
`basis`, `distribution` (the model family) and `interval`; and
`interval`, `n_boot`, `seed`, `conf_level`, `by` and the `estimand`.

## Details

|                |          |                                               |
|----------------|----------|-----------------------------------------------|
| Scale          | Basis    | Estimate                                      |
| `"difference"` | either   | `high - low`                                  |
| `"ratio"`      | `"loss"` | `D = 1 - 1 / rr`                              |
| `"ratio"`      | `"gain"` | `rr - 1`                                      |
| `"rii"`        | either   | `(rr - 1) / (rr + 1)`, equal to `D / (2 - D)` |
| `"lnrr"`       | either   | `log(rr)`                                     |

With `interval = "delta"`, the interval for a ratio-type scale is a
delta interval for `log(rr)`, transformed to the scale, so the loss
share stays below 1 and the RII inside (-1, 1); its p-value tests
`rr = 1`. The delta interval carries the uncertainty in the fixed
effects only. With `interval = "bootstrap"`, responses are simulated
from the fitted model (new random effects included), the model is
refitted to each simulated set, and the interval is the percentile
interval of the refitted estimates; `std.error` is their standard
deviation and `p.value` is `NA`.

## See also

[`pollin_estimand()`](https://max578.github.io/pollinatoR/reference/pollin_estimand.md)
for the scales and the basis.

## Examples

``` r
# \donttest{
set.seed(1)
d <- data.frame(
  trt = rep(rep(c("Cage", "Open"), each = 10L), 8L),
  Site = rep(paste0("s", seq_len(8L)), each = 20L)
)
d$pods <- rnbinom(nrow(d), mu = exp(3 + 0.3 * (d$trt == "Open")), size = 10)
pd <- pollin_data(d, "pods", "trt", unit = "Site")
ds <- pollin_design(pd, c(excluded = "Cage", open = "Open"), "rcbd")
est <- pollin_estimand(ds, "dependence", "ratio", "loss")
fit <- pollin_fit(est)
pollin_effects(fit)
#> <pollin_effects>
#>   family   : dependence (E2) 
#>   scale    : ratio 
#>   basis    : loss 
#>   estimate : 1 - Cage/Open (loss share) 
#>   interval : 95% delta 
#>   note     : delta intervals cover fixed-effect uncertainty only
#>  response     contrast  estimate std.error  conf.low conf.high      p.value
#>      pods Open vs Cage 0.2688808 0.0439437 0.1774744 0.3501293 1.882781e-07
#>      family scale basis distribution interval
#>  dependence ratio  loss      nbinom2    delta
pollin_effects(fit, interval = "bootstrap", n_boot = 50L, seed = 1L)
#> <pollin_effects>
#>   family   : dependence (E2) 
#>   scale    : ratio 
#>   basis    : loss 
#>   estimate : 1 - Cage/Open (loss share) 
#>   interval : 95% parametric bootstrap (50 refits) 
#>  response     contrast  estimate  std.error conf.low conf.high p.value
#>      pods Open vs Cage 0.2688808 0.03479058 0.210163 0.3186358      NA
#>      family scale basis distribution  interval
#>  dependence ratio  loss      nbinom2 bootstrap
# }
```

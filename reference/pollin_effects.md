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
  interval = c("delta", "delta_fixed", "bootstrap", "posterior"),
  n_boot = 999L,
  seed = NULL,
  conf_level = 0.95,
  n_draws = 4000L
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

  `"delta"` (the default for a glmmTMB fit), `"delta_fixed"`,
  `"bootstrap"`, or `"posterior"` (for a flexyBayes fit).

- n_boot:

  Number of parametric bootstrap refits, used when
  `interval = "bootstrap"`.

- seed:

  Optional integer seed for the bootstrap simulation or the posterior
  draws. Each response is drawn from this seed, and the session's random
  number state is restored afterwards.

- conf_level:

  Confidence level of the interval.

- n_draws:

  Number of posterior draws, used when `interval = "posterior"`.

## Value

An object of class `pollin_effects`: a list holding `table`, with one
row per response (and per `by` level) and the columns `response`,
`contrast`, any `by` column, `estimate`, `std.error`, `conf.low`,
`conf.high`, `p.value`, `family` (the estimand family), `scale`,
`basis`, `distribution` (the model family) and `interval`; and
`interval`, `n_boot`, `seed`, `conf_level`, `by`, the `estimand` and
`n_draws` (posterior intervals only).

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
`rr = 1`. The delta interval is taken over the fixed effects and the
random effects jointly, with their covariance from the inverse of the
joint precision of the Laplace approximation: the estimand averages over
the blocks and slopes of this trial, so the random effects enter at
their realised values, with their uncertainty and their correlation with
the fixed effects. A zero-inflated fit, or one whose joint precision is
singular (a variance at zero), is an error that points to the bootstrap.
`interval = "delta_fixed"` is the delta interval of pollinatoR 0.3.1 and
earlier, from `marginaleffects`: it varies the fixed effects and holds
the random effects at their modes, so it is wider than the joint
interval wherever the estimand depends on both, as on a design with a
random treatment slope or under a logit link. With
`interval = "bootstrap"`, responses are simulated from the fitted model
(new random effects included), the model is refitted to each simulated
set, and the interval is the percentile interval of the refitted
estimates; `std.error` is their standard deviation and `p.value` is
`NA`.

A fit with `engine = "flexyBayes"` has the `"posterior"` interval, its
default and only one: the estimate is the median of the posterior draws
of the estimand from
[`pollin_draws()`](https://max578.github.io/pollinatoR/reference/pollin_draws.md),
the interval their equal-tailed interval, `std.error` their standard
deviation, and `p.value` is `NA`. With the same `seed` and `n_draws`,
the row equals the summary of
`pollin_draws(fit, n = n_draws, seed = seed)`.

When the data declare a denominator
([`pollin_data()`](https://max578.github.io/pollinatoR/reference/pollin_data.md)),
each response is a count of successes out of that many trials, fitted as
binomial, and the arm means are expected proportions (successes /
trials), averaged over the rows with equal weight whatever their number
of trials; the contrast is between those averaged proportions.

The table reports the trial contrast, overall or per level of `by`. A
population declared with
[`pollin_population()`](https://max578.github.io/pollinatoR/reference/pollin_population.md)
is not used here; the population estimand is read with
[`pollin_draws()`](https://max578.github.io/pollinatoR/reference/pollin_draws.md).

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
#>   interval : 95% delta, fixed and random effects jointly 
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

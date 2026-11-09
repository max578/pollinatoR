# Draw the estimand from its sampling distribution

Returns draws of the two G-computation arm means of a fitted estimand
(the predictions under the `high` and the `low` arm, averaged over the
rows of the model data, as in
[`pollin_effects()`](https://max578.github.io/pollinatoR/reference/pollin_effects.md))
and the estimand formed from each pair of draws on the estimand's scale
and basis. A quantity computed from the estimand is then computed per
draw and summarised last.

## Usage

``` r
pollin_draws(
  x,
  n = 4000L,
  method = c("simulation", "bootstrap", "posterior"),
  response = NULL,
  seed = NULL
)
```

## Arguments

- x:

  A
  [`pollin_fit()`](https://max578.github.io/pollinatoR/reference/pollin_fit.md)
  object.

- n:

  Number of draws.

- method:

  `"simulation"` (the default for a glmmTMB fit), `"bootstrap"`, or
  `"posterior"` (for a flexyBayes fit).

- response:

  The response to draw. Required when the fit holds more than one
  response.

- seed:

  Optional seed for the draws. The session's random number state is
  restored afterwards. Posterior draws with a seed set both R's
  generator and INLA's, and are reproducible.

## Value

An object of class `pollin_draws`: a list holding

- `high`, `low`: numeric matrices of draws of the arm means, one row per
  draw and one column per level of the fit's `by` column, or a single
  column `"overall"` when the fit has no `by`;

- `estimate`: the matrix of the estimand on its scale and basis, formed
  from `high` and `low` draw by draw;

- `estimand`, `response`, `by`, `method`, `engine`, `n` (the number of
  draws), `seed`, and `population` (the
  [`pollin_population()`](https://max578.github.io/pollinatoR/reference/pollin_population.md)
  object, or `NULL` for the trial population).

## Details

Three methods are available:

- `"simulation"` (the default for a glmmTMB fit) draws the fixed effects
  and the random effects jointly from a multivariate normal centred on
  the estimates and the conditional modes, with the covariance of the
  delta interval of
  [`pollin_effects()`](https://max578.github.io/pollinatoR/reference/pollin_effects.md)
  (the inverse of the joint precision of the Laplace approximation), and
  recomputes the arm means for each draw. The draws carry the
  uncertainty in the realised block and slope effects of the trial and
  their correlation with the fixed effects.

- `"bootstrap"` keeps the arm means of the parametric bootstrap refits
  that
  [`pollin_effects()`](https://max578.github.io/pollinatoR/reference/pollin_effects.md)
  runs with `interval = "bootstrap"`: responses are simulated from the
  fitted model (new random effects included) and the model is refitted
  to each simulated set. With the same `seed` and `n = n_boot`, the
  percentile interval of the draws is the bootstrap interval of
  [`pollin_effects()`](https://max578.github.io/pollinatoR/reference/pollin_effects.md).
  Refits that fail are dropped, with a warning, so the number of draws
  can be smaller than `n`.

- `"posterior"`, the default and the only method for a fit with
  `engine = "flexyBayes"`, draws the fixed effects, the random effects
  and the variance parameters jointly from the posterior
  ([`flexyBayes::fb_inla_posterior_sample()`](https://aagi-aus.github.io/flexyBayes/reference/fb_extension.html))
  and recomputes the arm means for each draw with the prediction rule of
  [`pollin_effects()`](https://max578.github.io/pollinatoR/reference/pollin_effects.md):
  every row keeps its own random effects, here the drawn ones, and under
  each arm a row takes its unit's effect for that arm.

When the estimand refers to a population declared with
[`pollin_population()`](https://max578.github.io/pollinatoR/reference/pollin_population.md),
a further column `"population"` holds the population arm means, the
area-weighted sums of the stratum arm means, and the estimand formed
from them: a ratio of population means, not a weighted mean of the
stratum estimands.

## See also

[`pollin_effects()`](https://max578.github.io/pollinatoR/reference/pollin_effects.md)
for estimates and intervals;
[`pollin_population()`](https://max578.github.io/pollinatoR/reference/pollin_population.md)
for a population of strata.

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
dr <- pollin_draws(pollin_fit(est), n = 2000L, seed = 1L)
dr
#> <pollin_draws>
#>   family   : dependence (E2) 
#>   scale    : ratio 
#>   basis    : loss 
#>   estimate : 1 - Cage/Open (loss share) 
#>   response : pods 
#>   method   : simulation of the fixed and random effects 
#>   engine   : glmmTMB 
#>   draws    : 2000 
#>   columns  : overall 
summary(dr)
#>   response   level    median  conf.low conf.high     family scale basis
#> 1     pods overall 0.2691091 0.1741931 0.3524436 dependence ratio  loss
#>       method    n
#> 1 simulation 2000
# }
```

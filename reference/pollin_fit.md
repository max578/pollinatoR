# Fit the outcome models for an estimand

Fits one generalised linear mixed model per response in the estimand's
data, with the pollination treatment as the fixed effect and the random
terms of the design: `response ~ treatment + <design random terms>`.
With `by`, the fixed part is `treatment * by`, so the treatment effect
has its own fixed estimate at each level of `by`. The treatment is
releveled so the estimand's `low` arm is the reference level. Every
treatment level in the data is kept in the fit; the contrast is read
from the fitted model by
[`pollin_effects()`](https://max578.github.io/pollinatoR/reference/pollin_effects.md).

## Usage

``` r
pollin_fit(
  estimand,
  engine = "glmmTMB",
  family = NULL,
  by = NULL,
  pooled = NULL,
  ...
)
```

## Arguments

- estimand:

  A
  [`pollin_estimand()`](https://max578.github.io/pollinatoR/reference/pollin_estimand.md)
  object.

- engine:

  The fitting engine: `"glmmTMB"` (the default) or `"flexyBayes"`, a
  Bayesian fit on INLA through flexyBayes (see Details). The flexyBayes
  engine is experimental: its posterior intervals are calibrated on
  designs with unit intercepts and run slightly narrow on designs with a
  random treatment slope.

- family:

  `NULL` (the default) chooses a family per response from its values
  (count, proportion in (0, 1), positive continuous, or continuous with
  zeros). A single `glmmTMB` family object applies to every response; a
  named list gives one per response, and responses it does not name are
  chosen from their values.

- by:

  Optional name of a column of the data (for example variety or region).
  The treatment effect is fitted separately at each of its levels, as an
  interaction, and
  [`pollin_effects()`](https://max578.github.io/pollinatoR/reference/pollin_effects.md)
  then reports one contrast per level. The column is treated as a
  factor. An estimand whose population is a
  [`pollin_population()`](https://max578.github.io/pollinatoR/reference/pollin_population.md)
  must be fitted with `by` set to the population's `by` column.

- pooled:

  Optional name of a column of the data that is not a unit column (for
  example variety). Its levels get a random intercept and a random
  treatment slope, so the treatment effect at each level is pooled
  towards the overall effect;
  [`pollin_cate()`](https://max578.github.io/pollinatoR/reference/pollin_cate.md)
  reads the per-level effects. It must differ from `by`. A slope on a
  unit column is declared with `pollin_design(slope = )` instead. Rows
  with a missing value in the column are dropped.

- ...:

  Passed to
  [`glmmTMB::glmmTMB()`](https://rdrr.io/pkg/glmmTMB/man/glmmTMB.html),
  or to
  [`flexyBayes::flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.html)
  on the flexyBayes engine (for example `prior` or `seed`).

## Value

An object of class `pollin_fit`: a list holding the `estimand`, the
`engine`, the fitted models `fits` (named by response), the model `data`
for each response, the right-hand-side terms `rhs`, the `families` (as
glmmTMB family objects on both engines), the `by` column (or `NULL`),
and the further fitting arguments `args`.

## Details

Rows with a missing value in the response, the treatment, a unit column
or the `by` column are dropped for that response's model.

With `engine = "flexyBayes"` each response is fitted by
`flexyBayes::flexybayes(backend = "inla")`, which needs the flexyBayes
and INLA packages. The fixed part is that of the glmmTMB engine. Each
unit level is a random intercept; a nested level enters as one grouping
column `pollin_<outer>_<inner>` added to the model data. flexyBayes fits
a correlated random slope on brms only, so on INLA the treatment slope
of the design is a random intercept per unit plus an independent effect
per unit-by-arm cell (column `pollin_<unit>_<treatment>`): equal
variances across arms and a non-negative correlation between them, the
exchangeable form of the slope. A `pooled` column enters the same way.
The families map as `nbinom2()` to `"negative_binomial"`,
[`poisson()`](https://rdrr.io/r/stats/family.html) to `"poisson"`,
[`binomial()`](https://rdrr.io/r/stats/family.html) to `"binomial"`,
`Gamma(link = "log")` to `"gamma"` and `beta_family()` to `"beta"`; any
other family is refused. Priors are passed through `...` as
`flexybayes()` takes them; the seed of the posterior draws is set in
[`pollin_draws()`](https://max578.github.io/pollinatoR/reference/pollin_draws.md).
The backend is INLA only, since the posterior draws of a brms fit are
not in the flexyBayes extension API version 1. A fit whose convergence
verdict
([`flexyBayes::fb_convergence()`](https://aagi-aus.github.io/flexyBayes/reference/fb_extension.html))
is `FALSE` is an error, and one with no verdict a warning. Estimates and
intervals of a flexyBayes fit come from its posterior draws; see
[`pollin_draws()`](https://max578.github.io/pollinatoR/reference/pollin_draws.md).

## See also

[`pollin_effects()`](https://max578.github.io/pollinatoR/reference/pollin_effects.md)
for estimates and intervals on the estimand's scale.

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
pollin_fit(est)
#> <pollin_fit>
#>   engine    : glmmTMB 
#>   family    : dependence 
#>   contrast  : Open vs Cage 
#>   model     : <response> ~ trt + (1 | Site) 
#>   pods      : nbinom2, 160 rows
# }
```

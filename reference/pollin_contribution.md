# Estimate the pollinator contribution per outcome

Fits, for each outcome, a generalised linear mixed model with the
pollination treatment as the fixed effect of interest and the random
structure of the design
([`pollin_fit()`](https://max578.github.io/pollinatoR/reference/pollin_fit.md)),
then reports the contribution as a model-based average treatment effect,
the reference-minus-exclusion difference in expected response, with
propagated uncertainty, via G-computation
([`pollin_effects()`](https://max578.github.io/pollinatoR/reference/pollin_effects.md)).
A random treatment slope, where the design carries one, absorbs
unobserved between-unit heterogeneity in the contribution (for example
an unrecorded genotype). Family is chosen per outcome from its values
unless supplied.

## Usage

``` r
pollin_contribution(
  x,
  outcomes = NULL,
  by = NULL,
  families = list(),
  interval = c("delta", "delta_fixed"),
  ...
)
```

## Arguments

- x:

  A
  [`pollin_estimand()`](https://max578.github.io/pollinatoR/reference/pollin_estimand.md)
  object with `family = "dependence"` and `scale = "difference"`.

- outcomes:

  Outcomes to model. `NULL` (the default) models every response.

- by:

  Optional name of a column for conditional contrasts, one per level
  (for example a region).

- families:

  Optional named list of `glmmTMB` family objects, one per outcome;
  outcomes it does not name are chosen from their values.

- interval:

  `"delta"` (the default) or `"delta_fixed"`, the interval of
  [`pollin_effects()`](https://max578.github.io/pollinatoR/reference/pollin_effects.md);
  `"delta_fixed"` reproduces the intervals of pollinatoR 0.3.1 and
  earlier.

- ...:

  Passed to
  [`glmmTMB::glmmTMB()`](https://rdrr.io/pkg/glmmTMB/man/glmmTMB.html).

## Value

An object of class `pollin_contribution`: a list with `table` (a tidy
data frame of contrasts with FDR-adjusted p-values and the standard
deviation of the random treatment slope) and `fits` (the fitted models,
for diagnostics). From an estimand, `table` carries the columns of
[`pollin_effects()`](https://max578.github.io/pollinatoR/reference/pollin_effects.md),
and the list also holds the `estimand`.

## Details

The treatment effect is causal by design here – the exclusion contrast
is an experimental manipulation – so no covariate adjustment is applied;
see
[`pollin_moderation()`](https://max578.github.io/pollinatoR/reference/pollin_moderation.md)
for the observational "which insects" question.

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
est <- pollin_estimand(ds, "dependence", "difference", "loss")
pollin_contribution(est)
#> <pollin_contribution>  (open - exclusion, model-based ATE)
#>  response     contrast estimate std.error conf.low conf.high      p.value
#>      pods Open vs Cage   7.3875  1.436307  4.57239  10.20261 2.698125e-07
#>      family      scale basis distribution interval slope_sd        p.fdr
#>  dependence difference  loss      nbinom2    delta       NA 2.698125e-07
# }
```

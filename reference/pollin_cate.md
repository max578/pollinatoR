# Estimate the effect at each level of a grouping, with partial pooling

Reads the estimand's contrast at each level of a column that carries a
random treatment slope: the unit column of the design's slope
(`pollin_design(slope = )`, a slope on the outermost unit), or the
column declared with `pollin_fit(pooled = )`. Each level's effect is
formed by G-computation over the rows at that level, every row keeping
its random effects, so a level with little data is pulled towards the
overall effect.

## Usage

``` r
pollin_cate(fit, by = NULL, conf_level = 0.95, n = 4000L, seed = NULL)
```

## Arguments

- fit:

  A
  [`pollin_fit()`](https://max578.github.io/pollinatoR/reference/pollin_fit.md)
  object with a random treatment slope.

- by:

  The column whose levels are reported. `NULL` (the default) takes the
  fit's `pooled` column, or else the unit column of the design's slope.

- conf_level:

  Level of the interval.

- n:

  Number of draws for the interval.

- seed:

  Optional seed for the draws. The session's random number state is
  restored afterwards.

## Value

An object of class `pollin_cate`: a list holding `table`, with one row
per response and level and the columns `response`, `contrast`, `level`,
`estimate`, `std.error`, `conf.low`, `conf.high`, `pooled` (`TRUE`), `n`
(rows of the model data at the level), `family`, `scale`, `basis`,
`distribution` and `interval` (`"conditional"` on glmmTMB, `"posterior"`
on flexyBayes); and `by`, `engine`, `conf_level`, `n`, `seed` and the
`estimand`.

## Details

On the glmmTMB engine the estimate uses the fixed effects and the
conditional modes of the random effects. The interval is approximate:
the fixed and random effects are drawn jointly from a multivariate
normal with the covariance of the delta interval of
[`pollin_effects()`](https://max578.github.io/pollinatoR/reference/pollin_effects.md)
(the inverse of the joint precision of the Laplace approximation); it is
the percentile interval of the effects recomputed per draw, and
`std.error` their standard deviation. On the flexyBayes engine the
effects come from the posterior draws of
[`pollin_draws()`](https://max578.github.io/pollinatoR/reference/pollin_draws.md):
the estimate is their median and the interval their equal-tailed
interval.

The table holds pooled contrasts only. A fit that pools a column does
not also fit its levels as fixed effects, so the unpooled per-level
contrasts come from a second fit, `pollin_fit(estimand, by = )` without
the slope on that column, read with
[`pollin_effects()`](https://max578.github.io/pollinatoR/reference/pollin_effects.md);
the two side by side show what the pooling moved. The print names which
table is which.

## See also

[`pollin_effects()`](https://max578.github.io/pollinatoR/reference/pollin_effects.md)
for the trial contrast and the fixed per-level contrasts.

## Examples

``` r
# \donttest{
set.seed(3)
d <- data.frame(
  trt = rep(rep(c("Cage", "Open"), each = 8L), 12L),
  Site = rep(sprintf("s%02d", seq_len(12L)), each = 16L)
)
slope <- rnorm(12L, 0, 0.3)[as.integer(factor(d$Site))]
d$pods <- rnbinom(
  nrow(d),
  mu = exp(3 + (0.3 + slope) * (d$trt == "Open")),
  size = 20
)
pd <- pollin_data(d, "pods", "trt", unit = "Site")
ds <- pollin_design(pd, c(excluded = "Cage", open = "Open"), "multi_site")
est <- pollin_estimand(ds, "dependence", "ratio", "loss")
pollin_cate(pollin_fit(est), n = 1000L, seed = 1L)
#> <pollin_cate>
#>   by       : Site 
#>   scale    : ratio 
#>   basis    : loss 
#>   estimate : 1 - Cage/Open (loss share) 
#>   interval : 95% approximate, fixed and random effects jointly (1000 draws) 
#> Pooled: partial pooling through the random treatment slope on Site 
#>  response     contrast level    estimate  std.error    conf.low conf.high
#>      pods Open vs Cage   s01  0.13739779 0.11025817 -0.09685577 0.3214224
#>      pods Open vs Cage   s02  0.24772791 0.09322353  0.04517668 0.4191425
#>      pods Open vs Cage   s03  0.32379538 0.08455471  0.13054402 0.4647362
#>      pods Open vs Cage   s04  0.09358949 0.11281137 -0.14025235 0.2936593
#>      pods Open vs Cage   s05  0.34539898 0.07681851  0.18353729 0.4758282
#>      pods Open vs Cage   s06  0.35561185 0.07694023  0.17929007 0.4798230
#>      pods Open vs Cage   s07  0.31238229 0.08185953  0.13821557 0.4609481
#>      pods Open vs Cage   s08  0.32119972 0.08358701  0.13768183 0.4603141
#>      pods Open vs Cage   s09 -0.08410725 0.14108188 -0.37440616 0.1480293
#>      pods Open vs Cage   s10  0.59856163 0.05512882  0.47202134 0.6933275
#>      pods Open vs Cage   s11  0.01855484 0.13001571 -0.26382857 0.2379834
#>      pods Open vs Cage   s12  0.08055170 0.11286257 -0.18585308 0.2701437
#>  pooled  n     family scale basis distribution    interval
#>    TRUE 16 dependence ratio  loss      nbinom2 conditional
#>    TRUE 16 dependence ratio  loss      nbinom2 conditional
#>    TRUE 16 dependence ratio  loss      nbinom2 conditional
#>    TRUE 16 dependence ratio  loss      nbinom2 conditional
#>    TRUE 16 dependence ratio  loss      nbinom2 conditional
#>    TRUE 16 dependence ratio  loss      nbinom2 conditional
#>    TRUE 16 dependence ratio  loss      nbinom2 conditional
#>    TRUE 16 dependence ratio  loss      nbinom2 conditional
#>    TRUE 16 dependence ratio  loss      nbinom2 conditional
#>    TRUE 16 dependence ratio  loss      nbinom2 conditional
#>    TRUE 16 dependence ratio  loss      nbinom2 conditional
#>    TRUE 16 dependence ratio  loss      nbinom2 conditional
#> Unpooled: not carried by this fit, which pools Site; the fixed per-level contrasts come from `pollin_fit(by = "Site")` without the slope, read with pollin_effects().
if (requireNamespace("flexyBayes", quietly = TRUE) &&
  requireNamespace("INLA", quietly = TRUE)) {
  fit_b <- pollin_fit(est, engine = "flexyBayes")
  pollin_cate(fit_b, n = 1000L, seed = 1L)
}
#> flexyBayes: by default the residual and each random-effect standard deviation have a half_normal(scale = 3) prior, on the log scale of the mean, with the scale divided by sd(x) for a random slope on x. See ?fb_prior to set other priors; options(flexyBayes.silence_default_prior_note = TRUE) hides this note.
#> <pollin_cate>
#>   by       : Site 
#>   scale    : ratio 
#>   basis    : loss 
#>   estimate : 1 - Cage/Open (loss share) 
#>   interval : 95% posterior, equal-tailed (1000 draws) 
#> Pooled: partial pooling through the random treatment slope on Site 
#>  response     contrast level    estimate  std.error    conf.low conf.high
#>      pods Open vs Cage   s01  0.16674220 0.10434660 -0.05612786 0.3459806
#>      pods Open vs Cage   s02  0.28683257 0.08799859  0.09308358 0.4386692
#>      pods Open vs Cage   s03  0.33744471 0.07941388  0.15906565 0.4702732
#>      pods Open vs Cage   s04  0.10475245 0.11239027 -0.12921443 0.2943819
#>      pods Open vs Cage   s05  0.30099501 0.08193948  0.13078942 0.4419927
#>      pods Open vs Cage   s06  0.34364290 0.08211179  0.16630915 0.4997855
#>      pods Open vs Cage   s07  0.32363296 0.07942347  0.13784865 0.4615435
#>      pods Open vs Cage   s08  0.28340283 0.08173638  0.10509075 0.4255753
#>      pods Open vs Cage   s09 -0.06606611 0.13840611 -0.37845676 0.1603209
#>      pods Open vs Cage   s10  0.59125765 0.05446966  0.47019192 0.6832297
#>      pods Open vs Cage   s11  0.01719318 0.12336060 -0.24395693 0.2221635
#>      pods Open vs Cage   s12  0.12871296 0.10900000 -0.09798762 0.3039113
#>  pooled  n     family scale basis distribution  interval
#>    TRUE 16 dependence ratio  loss      nbinom2 posterior
#>    TRUE 16 dependence ratio  loss      nbinom2 posterior
#>    TRUE 16 dependence ratio  loss      nbinom2 posterior
#>    TRUE 16 dependence ratio  loss      nbinom2 posterior
#>    TRUE 16 dependence ratio  loss      nbinom2 posterior
#>    TRUE 16 dependence ratio  loss      nbinom2 posterior
#>    TRUE 16 dependence ratio  loss      nbinom2 posterior
#>    TRUE 16 dependence ratio  loss      nbinom2 posterior
#>    TRUE 16 dependence ratio  loss      nbinom2 posterior
#>    TRUE 16 dependence ratio  loss      nbinom2 posterior
#>    TRUE 16 dependence ratio  loss      nbinom2 posterior
#>    TRUE 16 dependence ratio  loss      nbinom2 posterior
#> Unpooled: not carried by this fit, which pools Site; the fixed per-level contrasts come from `pollin_fit(by = "Site")` without the slope, read with pollin_effects().
# }
```

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
the fixed effects are drawn from their sampling distribution and,
independently of them, each coefficient of the slope term from its
conditional variance, the other random terms held at their modes; it is
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
#>   interval : 95% approximate, fixed effects and conditional variance of the slope (1000 draws) 
#> Pooled: partial pooling through the random treatment slope on Site 
#>  response     contrast level    estimate  std.error    conf.low conf.high
#>      pods Open vs Cage   s01  0.13739779 0.15129199 -0.19507033 0.3858641
#>      pods Open vs Cage   s02  0.24772791 0.12959204 -0.04257898 0.4516558
#>      pods Open vs Cage   s03  0.32379538 0.11369481  0.06210526 0.5055297
#>      pods Open vs Cage   s04  0.09358949 0.15234740 -0.25038131 0.3281023
#>      pods Open vs Cage   s05  0.34539898 0.11526112  0.08358402 0.5402006
#>      pods Open vs Cage   s06  0.35561185 0.10955999  0.11940716 0.5374972
#>      pods Open vs Cage   s07  0.31238229 0.11667265  0.06816681 0.5102606
#>      pods Open vs Cage   s08  0.32119972 0.11858481  0.04391080 0.5070337
#>      pods Open vs Cage   s09 -0.08410725 0.19050122 -0.52560199 0.2147417
#>      pods Open vs Cage   s10  0.59856163 0.07065014  0.43360633 0.7165147
#>      pods Open vs Cage   s11  0.01855484 0.17185899 -0.37350365 0.2998889
#>      pods Open vs Cage   s12  0.08055170 0.16233668 -0.28683945 0.3441590
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
#>      pods Open vs Cage   s01  0.16674266 0.10434647 -0.05612691 0.3459809
#>      pods Open vs Cage   s02  0.28683263 0.08799850  0.09308358 0.4386692
#>      pods Open vs Cage   s03  0.33744459 0.07941382  0.15906553 0.4702729
#>      pods Open vs Cage   s04  0.10475349 0.11239011 -0.12921308 0.2943823
#>      pods Open vs Cage   s05  0.30099481 0.08193942  0.13078973 0.4419928
#>      pods Open vs Cage   s06  0.34364239 0.08211176  0.16630914 0.4997855
#>      pods Open vs Cage   s07  0.32363274 0.07942342  0.13784837 0.4615433
#>      pods Open vs Cage   s08  0.28340272 0.08173633  0.10509112 0.4255752
#>      pods Open vs Cage   s09 -0.06606523 0.13840585 -0.37845422 0.1603209
#>      pods Open vs Cage   s10  0.59125697 0.05446968  0.47019221 0.6832287
#>      pods Open vs Cage   s11  0.01719350 0.12336040 -0.24395613 0.2221641
#>      pods Open vs Cage   s12  0.12871303 0.10899986 -0.09798621 0.3039117
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

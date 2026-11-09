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
#>      pods Open vs Cage   s01  0.16674205 0.10434664 -0.05612837 0.3459803
#>      pods Open vs Cage   s02  0.28683263 0.08799862  0.09308378 0.4386693
#>      pods Open vs Cage   s03  0.33744486 0.07941389  0.15906581 0.4702735
#>      pods Open vs Cage   s04  0.10475197 0.11239033 -0.12921525 0.2943817
#>      pods Open vs Cage   s05  0.30099507 0.08193949  0.13078905 0.4419926
#>      pods Open vs Cage   s06  0.34364320 0.08211179  0.16630927 0.4997855
#>      pods Open vs Cage   s07  0.32363315 0.07942348  0.13784874 0.4615435
#>      pods Open vs Cage   s08  0.28340290 0.08173639  0.10509057 0.4255755
#>      pods Open vs Cage   s09 -0.06606694 0.13840626 -0.37845935 0.1603209
#>      pods Open vs Cage   s10  0.59125819 0.05446968  0.47019209 0.6832304
#>      pods Open vs Cage   s11  0.01719295 0.12336070 -0.24395705 0.2221631
#>      pods Open vs Cage   s12  0.12871292 0.10900005 -0.09798838 0.3039110
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

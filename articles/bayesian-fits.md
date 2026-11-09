# Bayesian fits with flexyBayes

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`pollinatoR`](https://github.com/max578/pollinatoR)`)`

## A second engine

[`pollin_fit()`](https://max578.github.io/pollinatoR/reference/pollin_fit.md)
fits with glmmTMB by default. With `engine = "flexyBayes"` it fits the
same model as a Bayesian mixed model, through the flexyBayes package on
INLA. Estimates and intervals then come from posterior draws, which
carry the uncertainty in the variance components as well as in the fixed
effects. The rest of the package takes those draws as it takes the
glmmTMB ones: effects, populations, values and guards work unchanged.

In planted-truth simulations the posterior 95% intervals cover close to
95% on designs with unit intercepts and on binomial counts. On designs
with a random treatment slope and 12 sites they cover about 91% to 93%.
An exact posterior computed by Stan gives the same coverage there, so
the shortfall comes from the interval itself with few sites. The package
NEWS gives the full results.

This vignette fits one simulated trial on both engines and compares the
results: the overall effect, its value, and the effect in each variety.

## Installing the engine

Neither package is on CRAN. flexyBayes installs from GitHub, and INLA
from its own repository.

\
`# install.packages("remotes")`\
`remotes``::`[`install_github`](https://rdrr.io/pkg/remotes/man/install_github.html)`(``"AAGI-AUS/flexyBayes"``)`\
[`install.packages`](https://rdrr.io/r/utils/install.packages.html)`(`\
`  ``"INLA"``,`\
`  repos ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(`[`getOption`](https://rdrr.io/r/base/options.html)`(``"repos"``)``, INLA ``=`` ``"https://inla.r-inla-download.org/R/stable"``)``,`\
`  dependencies ``=`` ``TRUE`\
`)`

The engine is checked when a fit is made. Without both packages,
`engine = "flexyBayes"` stops with a message naming the one that is
missing.

## A simulated trial

Eight canola varieties are grown in six blocks. In every block each
variety has bagged plants, closed to insects, and open plants. The
response is the number of pods per plant. Varieties differ in how much
they gain from insect visits: each has its own treatment effect on the
log scale, drawn around a mean of 0.4.

\
[`set.seed`](https://rdrr.io/r/base/Random.html)`(``2026``)`\
`varieties`` ``<-`` `[`sprintf`](https://rdrr.io/r/base/sprintf.html)`(``"V%d"``, ``1``:``8``)`\
`blocks`` ``<-`` `[`sprintf`](https://rdrr.io/r/base/sprintf.html)`(``"B%d"``, ``1``:``6``)`\
`trial`` ``<-`` `[`expand.grid`](https://rdrr.io/r/base/expand.grid.html)`(`\
`  plant ``=`` ``1``:``3``,`\
`  trt ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"Bagged"``, ``"Open"``)``,`\
`  block ``=`` ``blocks``,`\
`  variety ``=`` ``varieties``,`\
`  stringsAsFactors ``=`` ``FALSE`\
`)`\
`block_effect`` ``<-`` `[`setNames`](https://rdrr.io/r/stats/setNames.html)`(`[`rnorm`](https://rdrr.io/r/stats/Normal.html)`(``6``, ``0``, ``0.2``)``, ``blocks``)`\
`variety_effect`` ``<-`` `[`setNames`](https://rdrr.io/r/stats/setNames.html)`(`[`rnorm`](https://rdrr.io/r/stats/Normal.html)`(``8``, ``0``, ``0.2``)``, ``varieties``)`\
`variety_slope`` ``<-`` `[`setNames`](https://rdrr.io/r/stats/setNames.html)`(`[`rnorm`](https://rdrr.io/r/stats/Normal.html)`(``8``, ``0.4``, ``0.25``)``, ``varieties``)`\
`log_mu`` ``<-`` ``3`` ``+`` ``block_effect``[``trial``$``block``]`` ``+`` ``variety_effect``[``trial``$``variety``]`` ``+`\
`  ``variety_slope``[``trial``$``variety``]`` ``*`` ``(``trial``$``trt`` ``==`` ``"Open"``)`\
`trial``$``pods`` ``<-`` `[`rnbinom`](https://rdrr.io/r/stats/NegBinomial.html)`(`[`nrow`](https://rdrr.io/r/base/nrow.html)`(``trial``)``, mu ``=`` `[`exp`](https://rdrr.io/r/base/Log.html)`(``log_mu``)``, size ``=`` ``15``)`\
[`nrow`](https://rdrr.io/r/base/nrow.html)`(``trial``)`\
`#> [1] 288`\
\
`pd`` ``<-`` `[`pollin_data`](https://max578.github.io/pollinatoR/reference/pollin_data.md)`(``trial``, response ``=`` ``"pods"``, treatment ``=`` ``"trt"``, unit ``=`` ``"block"``)`\
`ds`` ``<-`` `[`pollin_design`](https://max578.github.io/pollinatoR/reference/pollin_design.md)`(``pd``, arms ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``excluded ``=`` ``"Bagged"``, open ``=`` ``"Open"``)``, layout ``=`` ``"rcbd"``)`\
`est`` ``<-`` `[`pollin_estimand`](https://max578.github.io/pollinatoR/reference/pollin_estimand.md)`(``ds``, ``"dependence"``, ``"ratio"``, ``"loss"``)`

The planted loss share of each variety is `1 - exp(-slope)`:

\
`planted`` ``<-`` ``1`` ``-`` `[`exp`](https://rdrr.io/r/base/Log.html)`(``-``variety_slope``)`\
[`round`](https://rdrr.io/r/base/Round.html)`(``planted``, ``2``)`\
`#>    V1    V2    V3    V4    V5    V6    V7    V8 `\
`#> -0.27  0.52  0.43  0.37  0.18  0.44  0.27  0.30`

## The same model on both engines

The fixed part and the random terms are the same on both engines: the
treatment as a fixed effect and a random intercept per block. A negative
binomial family is set explicitly, so both engines use it.

\
`# One INLA fit per model instead of the best of three.`\
[`options`](https://rdrr.io/r/base/options.html)`(``flexyBayes.inla_best_of ``=`` ``1L``, flexyBayes.silence_default_prior_note ``=`` ``TRUE``)`

\
`fit_tmb`` ``<-`` `[`pollin_fit`](https://max578.github.io/pollinatoR/reference/pollin_fit.md)`(``est``, family ``=`` ``glmmTMB``::`[`nbinom2`](https://rdrr.io/pkg/glmmTMB/man/nbinom2.html)`(``)``)`

\
`fit_fb`` ``<-`` `[`pollin_fit`](https://max578.github.io/pollinatoR/reference/pollin_fit.md)`(``est``, engine ``=`` ``"flexyBayes"``, family ``=`` ``glmmTMB``::`[`nbinom2`](https://rdrr.io/pkg/glmmTMB/man/nbinom2.html)`(``)``)`\
`fit_fb`\
`#> <pollin_fit>`\
`#>   engine    : flexyBayes `\
`#>   family    : dependence `\
`#>   contrast  : Open vs Bagged `\
`#>   model     : <response> ~ trt + block `\
`#>   pods      : nbinom2, 288 rows`

[`pollin_effects()`](https://max578.github.io/pollinatoR/reference/pollin_effects.md)
reads the loss share from each fit. On glmmTMB the interval is the delta
interval. On flexyBayes it is the posterior interval: the estimate is
the median of the posterior draws of the loss share, the interval their
equal-tailed 95% interval, and `p.value` is `NA`.

\
`eff_tmb`` ``<-`` `[`pollin_effects`](https://max578.github.io/pollinatoR/reference/pollin_effects.md)`(``fit_tmb``)`\
`eff_tmb``$``table``[``, `[`c`](https://rdrr.io/r/base/c.html)`(``"estimate"``, ``"std.error"``, ``"conf.low"``, ``"conf.high"``, ``"interval"``)``]`\
`#>    estimate  std.error  conf.low conf.high interval`\
`#> 1 0.3289669 0.03075109 0.2659063 0.3866105    delta`

\
`eff_fb`` ``<-`` `[`pollin_effects`](https://max578.github.io/pollinatoR/reference/pollin_effects.md)`(``fit_fb``, n_draws ``=`` ``1000L``, seed ``=`` ``1L``)`\
[`rbind`](https://rdrr.io/r/base/cbind.html)`(``eff_tmb``$``table``, ``eff_fb``$``table``)``[`\
`  , `[`c`](https://rdrr.io/r/base/c.html)`(``"estimate"``, ``"std.error"``, ``"conf.low"``, ``"conf.high"``, ``"interval"``)`\
`]`\
`#>    estimate  std.error  conf.low conf.high  interval`\
`#> 1 0.3289669 0.03075109 0.2659063 0.3866105     delta`\
`#> 2 0.3272635 0.03177409 0.2624175 0.3864073 posterior`

With a few hundred plants and a block intercept as the only random term,
the two engines agree closely. The posterior interval also covers the
uncertainty in the block variance, which the delta interval does not;
with few blocks this can widen it.

## From posterior draws to a value

On a flexyBayes fit,
[`pollin_draws()`](https://max578.github.io/pollinatoR/reference/pollin_draws.md)
uses `method = "posterior"`, its only method there. Each draw takes the
fixed effects, the block effects and the variance parameters jointly
from the posterior, and the two arm means are recomputed from them. A
seed sets both R’s random number generator and INLA’s, so the draws
repeat exactly.

\
`dr_fb`` ``<-`` `[`pollin_draws`](https://max578.github.io/pollinatoR/reference/pollin_draws.md)`(``fit_fb``, n ``=`` ``1000L``, method ``=`` ``"posterior"``, seed ``=`` ``1L``)`\
`dr_fb`\
`#> <pollin_draws>`\
`#>   family   : dependence (E2) `\
`#>   scale    : ratio `\
`#>   basis    : loss `\
`#>   estimate : 1 - Bagged/Open (loss share) `\
`#>   response : pods `\
`#>   method   : joint posterior draws (INLA) `\
`#>   engine   : flexyBayes `\
`#>   draws    : 1000 `\
`#>   columns  : overall`

The draws go into
[`pollin_value()`](https://max578.github.io/pollinatoR/reference/pollin_value.md)
exactly as simulation draws do. The market inputs below are
illustrative.

\
`canola`` ``<-`` `[`pollin_market`](https://max578.github.io/pollinatoR/reference/pollin_market.md)`(`\
`  quantity ``=`` ``3600``,`\
`  price ``=`` ``700``,`\
`  unit ``=`` ``"t"``,`\
`  currency ``=`` ``"AUD"``,`\
`  year ``=`` ``2025``,`\
`  source ``=`` ``"Illustrative inputs"``,`\
`  crop ``=`` ``"canola"`\
`)`\
`dr_tmb`` ``<-`` `[`pollin_draws`](https://max578.github.io/pollinatoR/reference/pollin_draws.md)`(``fit_tmb``, n ``=`` ``1000L``, seed ``=`` ``1L``)`\
[`rbind`](https://rdrr.io/r/base/cbind.html)`(`\
`  `[`summary`](https://rdrr.io/r/base/summary.html)`(`[`pollin_value`](https://max578.github.io/pollinatoR/reference/pollin_value.md)`(``dr_tmb``, ``canola``, channel ``=`` ``"dependence"``)``)``,`\
`  `[`summary`](https://rdrr.io/r/base/summary.html)`(`[`pollin_value`](https://max578.github.io/pollinatoR/reference/pollin_value.md)`(``dr_fb``, ``canola``, channel ``=`` ``"dependence"``)``)`\
`)``[``, `[`c`](https://rdrr.io/r/base/c.html)`(``"median"``, ``"conf.low"``, ``"conf.high"``, ``"currency"``, ``"method"``)``]`\
`#>     median conf.low conf.high currency     method`\
`#> 1 824536.6 662938.5  977626.5      AUD simulation`\
`#> 2 824704.1 661292.1  973746.5      AUD  posterior`

The vignette *From dependence to value* covers the channels, populations
and the double-count guard. Each of them takes posterior draws.

## The effect in each variety

The varieties were planted with different loss shares. There are two
ways to estimate them.

The unpooled way fits a separate treatment effect for each variety, as a
fixed interaction, with `pollin_fit(by = "variety")`, and reads them
with
[`pollin_effects()`](https://max578.github.io/pollinatoR/reference/pollin_effects.md).
Each variety’s estimate uses only that variety’s plants.

\
`fit_by`` ``<-`` `[`pollin_fit`](https://max578.github.io/pollinatoR/reference/pollin_fit.md)`(``est``, family ``=`` ``glmmTMB``::`[`nbinom2`](https://rdrr.io/pkg/glmmTMB/man/nbinom2.html)`(``)``, by ``=`` ``"variety"``)`\
`unpooled`` ``<-`` `[`pollin_effects`](https://max578.github.io/pollinatoR/reference/pollin_effects.md)`(``fit_by``)``$``table`

The pooled way gives the varieties a random intercept and a random
treatment slope, with `pollin_fit(pooled = "variety")`, and reads the
effect at each level with
[`pollin_cate()`](https://max578.github.io/pollinatoR/reference/pollin_cate.md).
A variety whose estimate is uncertain, or far from the rest, is pulled
towards the overall effect. A slope on a unit column, such as site in a
multi-site study, is declared in the design with
`pollin_design(slope = )` instead.

\
`fit_pool_tmb`` ``<-`` `[`pollin_fit`](https://max578.github.io/pollinatoR/reference/pollin_fit.md)`(``est``, family ``=`` ``glmmTMB``::`[`nbinom2`](https://rdrr.io/pkg/glmmTMB/man/nbinom2.html)`(``)``, pooled ``=`` ``"variety"``)`\
`cate_tmb`` ``<-`` `[`pollin_cate`](https://max578.github.io/pollinatoR/reference/pollin_cate.md)`(``fit_pool_tmb``, n ``=`` ``1000L``, seed ``=`` ``1L``)`

\
`fit_pool_fb`` ``<-`` `[`pollin_fit`](https://max578.github.io/pollinatoR/reference/pollin_fit.md)`(`\
`  ``est``,`\
`  engine ``=`` ``"flexyBayes"``,`\
`  family ``=`` ``glmmTMB``::`[`nbinom2`](https://rdrr.io/pkg/glmmTMB/man/nbinom2.html)`(``)``,`\
`  pooled ``=`` ``"variety"`\
`)`\
`cate_fb`` ``<-`` `[`pollin_cate`](https://max578.github.io/pollinatoR/reference/pollin_cate.md)`(``fit_pool_fb``, n ``=`` ``1000L``, seed ``=`` ``1L``)`\
`cate_fb`\
`#> <pollin_cate>`\
`#>   by       : variety `\
`#>   scale    : ratio `\
`#>   basis    : loss `\
`#>   estimate : 1 - Bagged/Open (loss share) `\
`#>   interval : 95% posterior, equal-tailed (1000 draws) `\
`#> Pooled: partial pooling through the random treatment slope on variety `\
`#>  response       contrast level    estimate  std.error    conf.low conf.high`\
`#>      pods Open vs Bagged    V1 -0.03764906 0.12817868 -0.32706420 0.1824306`\
`#>      pods Open vs Bagged    V2  0.52106736 0.04963085  0.42457569 0.6117435`\
`#>      pods Open vs Bagged    V3  0.41298404 0.05809910  0.28692986 0.5161863`\
`#>      pods Open vs Bagged    V4  0.30969132 0.06842307  0.17060543 0.4301072`\
`#>      pods Open vs Bagged    V5  0.15062597 0.09294575 -0.05051148 0.3107033`\
`#>      pods Open vs Bagged    V6  0.28531157 0.06768514  0.13552347 0.3989920`\
`#>      pods Open vs Bagged    V7  0.41258609 0.06275980  0.27923350 0.5295922`\
`#>      pods Open vs Bagged    V8  0.26678752 0.07837492  0.09901901 0.4024847`\
`#>  pooled  n     family scale basis distribution  interval`\
`#>    TRUE 36 dependence ratio  loss      nbinom2 posterior`\
`#>    TRUE 36 dependence ratio  loss      nbinom2 posterior`\
`#>    TRUE 36 dependence ratio  loss      nbinom2 posterior`\
`#>    TRUE 36 dependence ratio  loss      nbinom2 posterior`\
`#>    TRUE 36 dependence ratio  loss      nbinom2 posterior`\
`#>    TRUE 36 dependence ratio  loss      nbinom2 posterior`\
`#>    TRUE 36 dependence ratio  loss      nbinom2 posterior`\
`#>    TRUE 36 dependence ratio  loss      nbinom2 posterior`\
`` #> Unpooled: not carried by this fit, which pools variety; the fixed per-level contrasts come from `pollin_fit(by = "variety")` without the slope, read with pollin_effects(). ``

The print keeps the two kinds apart: the table holds the pooled
contrasts, and a line points to the unpooled ones. Set side by side with
the planted values:

\
`compare`` ``<-`` `[`data.frame`](https://rdrr.io/r/base/data.frame.html)`(`\
`  variety ``=`` ``varieties``,`\
`  planted ``=`` `[`round`](https://rdrr.io/r/base/Round.html)`(``planted``, ``2``)``,`\
`  unpooled ``=`` `[`round`](https://rdrr.io/r/base/Round.html)`(``unpooled``$``estimate``, ``2``)``,`\
`  pooled_glmmTMB ``=`` `[`round`](https://rdrr.io/r/base/Round.html)`(``cate_tmb``$``table``$``estimate``, ``2``)`\
`)`\
`compare`\
`#>    variety planted unpooled pooled_glmmTMB`\
`#> V1      V1   -0.27    -0.14          -0.10`\
`#> V2      V2    0.52     0.55           0.49`\
`#> V3      V3    0.43     0.43           0.45`\
`#> V4      V4    0.37     0.31           0.33`\
`#> V5      V5    0.18     0.12           0.17`\
`#> V6      V6    0.44     0.28           0.33`\
`#> V7      V7    0.27     0.43           0.37`\
`#> V8      V8    0.30     0.26           0.26`

\
`compare``$``pooled_flexyBayes`` ``<-`` `[`round`](https://rdrr.io/r/base/Round.html)`(``cate_fb``$``table``$``estimate``, ``2``)`\
`compare`\
`#>    variety planted unpooled pooled_glmmTMB pooled_flexyBayes`\
`#> V1      V1   -0.27    -0.14          -0.10             -0.04`\
`#> V2      V2    0.52     0.55           0.49              0.52`\
`#> V3      V3    0.43     0.43           0.45              0.41`\
`#> V4      V4    0.37     0.31           0.33              0.31`\
`#> V5      V5    0.18     0.12           0.17              0.15`\
`#> V6      V6    0.44     0.28           0.33              0.29`\
`#> V7      V7    0.27     0.43           0.37              0.41`\
`#> V8      V8    0.30     0.26           0.26              0.27`

Pooling pulls the extreme varieties inwards: V1 rises on both engines,
and V2 and V7 fall. This does not bring every variety closer to its
planted value. V1 was planted with a negative loss share and moves away
from it, while V7 moves towards it. Partial pooling lowers the error
over the varieties as a whole, not for each one.

The two engines differ in how they get the interval. On flexyBayes it
comes from the posterior draws, as for the overall effect. On glmmTMB it
is approximate: the fixed effects are drawn from their sampling
distribution and, independently of them, each variety’s slope from its
conditional variance. glmmTMB gives the conditional variance of each
coefficient on its own, not their joint covariance, so the approximation
cannot carry the correlation between them.

\
`cate_tmb``$``table``[``, `[`c`](https://rdrr.io/r/base/c.html)`(``"level"``, ``"estimate"``, ``"conf.low"``, ``"conf.high"``)``]`\
`#>   level   estimate     conf.low conf.high`\
`#> 1    V1 -0.0997176 -0.531228745 0.2059457`\
`#> 2    V2  0.4878510  0.288503895 0.6231699`\
`#> 3    V3  0.4502143  0.266566464 0.5886328`\
`#> 4    V4  0.3255154  0.102652785 0.4975262`\
`#> 5    V5  0.1719329 -0.136818373 0.3901447`\
`#> 6    V6  0.3327006  0.085463620 0.5134375`\
`#> 7    V7  0.3687068  0.155724480 0.5415963`\
`#> 8    V8  0.2630115  0.002286701 0.4538072`

In this trial they come out wider than the unpooled delta intervals,
although pooling should narrow them. The flexyBayes intervals are as
narrow as the unpooled ones or narrower.

When glmmTMB estimates the slope variance at zero, every variety’s
pooled effect equals the overall one, and the conditional variances
carry no information.
[`pollin_cate()`](https://max578.github.io/pollinatoR/reference/pollin_cate.md)
stops on such a fit rather than report intervals; the unpooled contrasts
or the flexyBayes engine remain.

## What the engine refuses

The flexyBayes engine refuses, by name, what it cannot fit faithfully.

flexyBayes can also fit with brms, but the posterior draws of a brms fit
are not part of version 1 of the flexyBayes extension API, which is what
`pollinatoR` reads. The brms backend is refused:

\
[`try`](https://rdrr.io/r/base/try.html)`(`[`pollin_fit`](https://max578.github.io/pollinatoR/reference/pollin_fit.md)`(``est``, engine ``=`` ``"flexyBayes"``, backend ``=`` ``"brms"``)``)`\
`#> Error : The flexyBayes engine fits on INLA only: posterior draws of a brms fit are not in the flexyBayes extension API version 1.`

Five families are mapped to their flexyBayes counterparts on INLA:
`nbinom2()` (as `negative_binomial`),
[`poisson()`](https://rdrr.io/r/stats/family.html),
[`binomial()`](https://rdrr.io/r/stats/family.html),
`Gamma(link = "log")` (as `gamma`) and `beta_family()` (as `beta`), each
with its log or logit link. Any other family, such as Tweedie, is
refused:

\
[`try`](https://rdrr.io/r/base/try.html)`(`[`pollin_fit`](https://max578.github.io/pollinatoR/reference/pollin_fit.md)`(``est``, engine ``=`` ``"flexyBayes"``, family ``=`` ``glmmTMB``::`[`tweedie`](https://rdrr.io/pkg/glmmTMB/man/nbinom2.html)`(``)``)``)`\
`` #> Error : Family 'tweedie' of 'pods' has no flexyBayes counterpart on INLA. The flexyBayes engine fits nbinom2, poisson, binomial, Gamma, beta; choose one with `family =`, or use `engine = "glmmTMB"`. ``

flexyBayes on INLA does not fit a correlated random treatment slope, the
unstructured slope that glmmTMB fits as `(trt | variety)`. The
flexyBayes route fits a variety intercept plus an independent effect for
each variety-by-arm cell instead. That is a restricted form of the
unstructured slope: the variance is the same in both arms and the
correlation between arms cannot be negative. The pooled effects on the
two engines therefore come from slightly different models, and where the
true arms vary unequally between varieties they can disagree for that
reason alone.

## Priors

Priors are passed to
[`flexyBayes::flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.html)
through `...` as written, in the flexyBayes prior language of
[`flexyBayes::fb_prior()`](https://aagi-aus.github.io/flexyBayes/reference/fb_prior.html).
`pollinatoR` adds no priors of its own. Without one, flexyBayes uses its
default, a half-normal prior on each random-effect standard deviation
and a nearly flat prior on the fixed effects. Here the treatment
coefficient, the log of the open-to-bagged ratio, gets a normal prior
with standard deviation 1:

\
`fit_prior`` ``<-`` `[`pollin_fit`](https://max578.github.io/pollinatoR/reference/pollin_fit.md)`(`\
`  ``est``,`\
`  engine ``=`` ``"flexyBayes"``,`\
`  family ``=`` ``glmmTMB``::`[`nbinom2`](https://rdrr.io/pkg/glmmTMB/man/nbinom2.html)`(``)``,`\
`  prior ``=`` ``flexyBayes``::`[`fb_prior`](https://aagi-aus.github.io/flexyBayes/reference/fb_prior.html)`(``b``(``"trtOpen"``)`` ``~`` ``normal``(``mean ``=`` ``0``, sd ``=`` ``1``)``)`\
`)`\
[`pollin_effects`](https://max578.github.io/pollinatoR/reference/pollin_effects.md)`(``fit_prior``, n_draws ``=`` ``1000L``, seed ``=`` ``1L``)``$``table``[`\
`  , `[`c`](https://rdrr.io/r/base/c.html)`(``"estimate"``, ``"conf.low"``, ``"conf.high"``)`\
`]`\
`#>    estimate conf.low conf.high`\
`#> 1 0.3267644 0.261867  0.385872`

[`flexyBayes::prior_summary()`](https://aagi-aus.github.io/flexyBayes/reference/prior_summary.html)
on a fitted model, `fit_prior$fits$pods`, lists each prior and the INLA
parameter it reached.

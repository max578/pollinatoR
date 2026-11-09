# Getting started with pollinatoR

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`pollinatoR`](https://github.com/max578/pollinatoR)`)`

## What the package does

`pollinatoR` estimates how much of a crop’s yield depends on insect
pollinators, from experiments that compare plants with and without
pollinator access: caged or bagged plants against open ones, and, where
a trial has it, a hand-supplemented arm. The methods apply to any crop
studied with such a design.

An analysis is built in layers. Each layer is one function that returns
an object you can print and check before going on.

| Layer | Function | What you declare or get |
|----|----|----|
| Data | [`pollin_data()`](https://max578.github.io/pollinatoR/reference/pollin_data.md) | Which columns are responses, treatment, units, groups and exposures |
| Design | [`pollin_design()`](https://max578.github.io/pollinatoR/reference/pollin_design.md) | Which treatment level is which arm, and the layout of the experiment |
| Estimand | [`pollin_estimand()`](https://max578.github.io/pollinatoR/reference/pollin_estimand.md) | The contrast, the scale and the basis of the number you want |
| Fit | [`pollin_fit()`](https://max578.github.io/pollinatoR/reference/pollin_fit.md) | One mixed model per response |
| Effects | [`pollin_effects()`](https://max578.github.io/pollinatoR/reference/pollin_effects.md) | The estimate with its interval, on the declared scale |
| Guards | `pollin_guard_*()` | Checks against common errors of analysis and reading |

## A simulated study

Twelve sites, sixteen plants per site, half caged and half open. Two
yield components are recorded per plant, along with counts of honey bees
and other insects. Open access raises pod count by an amount that varies
between sites, and slightly more on plants with more honey bees.

\
[`set.seed`](https://rdrr.io/r/base/Random.html)`(``1``)`\
`n_site`` ``<-`` ``12L`\
`n_plant`` ``<-`` ``16L`\
`site`` ``<-`` `[`rep`](https://rdrr.io/r/base/rep.html)`(`[`sprintf`](https://rdrr.io/r/base/sprintf.html)`(``"s%02d"``, `[`seq_len`](https://rdrr.io/r/base/seq.html)`(``n_site``)``)``, each ``=`` ``n_plant``)`\
`trt`` ``<-`` `[`rep`](https://rdrr.io/r/base/rep.html)`(`[`rep`](https://rdrr.io/r/base/rep.html)`(`[`c`](https://rdrr.io/r/base/c.html)`(``"Cage"``, ``"Open"``)``, each ``=`` ``n_plant`` ``/`` ``2L``)``, ``n_site``)`\
`honeybee`` ``<-`` `[`rep`](https://rdrr.io/r/base/rep.html)`(`[`rpois`](https://rdrr.io/r/stats/Poisson.html)`(``n_site``, ``8``)``, each ``=`` ``n_plant``)`` ``+`` `[`rpois`](https://rdrr.io/r/stats/Poisson.html)`(``n_site`` ``*`` ``n_plant``, ``1``)`\
`wild`` ``<-`` `[`rep`](https://rdrr.io/r/base/rep.html)`(`[`rpois`](https://rdrr.io/r/stats/Poisson.html)`(``n_site``, ``4``)``, each ``=`` ``n_plant``)`` ``+`` `[`rpois`](https://rdrr.io/r/stats/Poisson.html)`(``n_site`` ``*`` ``n_plant``, ``1``)`\
`idx`` ``<-`` `[`as.integer`](https://rdrr.io/r/base/integer.html)`(`[`factor`](https://rdrr.io/r/base/factor.html)`(``site``)``)`\
`site_effect`` ``<-`` `[`rnorm`](https://rdrr.io/r/stats/Normal.html)`(``n_site``, ``0``, ``0.3``)``[``idx``]`\
`site_slope`` ``<-`` `[`rnorm`](https://rdrr.io/r/stats/Normal.html)`(``n_site``, ``0``, ``0.25``)``[``idx``]`\
`open`` ``<-`` ``trt`` ``==`` ``"Open"`\
`pod_mu`` ``<-`` `[`exp`](https://rdrr.io/r/base/Log.html)`(``1.8`` ``+`` ``site_effect`` ``+`` ``(``0.35`` ``+`` ``site_slope`` ``+`` ``0.015`` ``*`` ``honeybee``)`` ``*`` ``open``)`\
`seed_mu`` ``<-`` `[`exp`](https://rdrr.io/r/base/Log.html)`(``2.6`` ``+`` ``site_effect`` ``+`` ``(``0.12`` ``+`` ``site_slope`` ``/`` ``2``)`` ``*`` ``open``)`\
`pods`` ``<-`` `[`rnbinom`](https://rdrr.io/r/stats/NegBinomial.html)`(``n_site`` ``*`` ``n_plant``, mu ``=`` ``pod_mu``, size ``=`` ``15``)`\
`seeds_per_pod`` ``<-`` `[`rnbinom`](https://rdrr.io/r/stats/NegBinomial.html)`(``n_site`` ``*`` ``n_plant``, mu ``=`` ``seed_mu``, size ``=`` ``15``)`\
`canola`` ``<-`` `[`data.frame`](https://rdrr.io/r/base/data.frame.html)`(``site``, ``trt``, ``honeybee``, ``wild``, ``pods``, ``seeds_per_pod``)`\
[`head`](https://rdrr.io/r/utils/head.html)`(``canola``)`\
`#>   site  trt honeybee wild pods seeds_per_pod`\
`#> 1  s01 Cage        7    3    3             9`\
`#> 2  s01 Cage        7    5    0             4`\
`#> 3  s01 Cage        8    5    4             1`\
`#> 4  s01 Cage        7    2    2             4`\
`#> 5  s01 Cage        7    2    4             6`\
`#> 6  s01 Cage       10    2    2             7`

## Data

[`pollin_data()`](https://max578.github.io/pollinatoR/reference/pollin_data.md)
pairs the table with a map of what each column is. The roles are checked
against the data, and the structural checks of
[`pollin_validate()`](https://max578.github.io/pollinatoR/reference/pollin_validate.md)
are run and stored.

\
`pd`` ``<-`` `[`pollin_data`](https://max578.github.io/pollinatoR/reference/pollin_data.md)`(`\
`  ``canola``,`\
`  response ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"pods"``, ``"seeds_per_pod"``)``,`\
`  treatment ``=`` ``"trt"``,`\
`  unit ``=`` ``"site"``,`\
`  exposures ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"honeybee"``, ``"wild"``)`\
`)`\
`pd`\
`#> <pollin_data>`\
`#>   rows        : 192 `\
`#>   response    : pods, seeds_per_pod `\
`#>   treatment   : trt ( Cage, Open )`\
`#>   unit        : site `\
`#>   exposures   : honeybee, wild `\
`` #>   checks      : 0 warnings; see `$checks` ``

## Design

[`pollin_design()`](https://max578.github.io/pollinatoR/reference/pollin_design.md)
says which treatment level is the excluded arm and which is the open
arm, and names the layout. The layout sets the random structure of the
models: here each site carries a random intercept and a random treatment
slope, so the pollinator effect may differ between sites.

\
`ds`` ``<-`` `[`pollin_design`](https://max578.github.io/pollinatoR/reference/pollin_design.md)`(`\
`  ``pd``,`\
`  arms ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``excluded ``=`` ``"Cage"``, open ``=`` ``"Open"``)``,`\
`  layout ``=`` ``"multi_site"`\
`)`\
`ds`\
`#> <pollin_design>`\
`#>   layout    : multi_site `\
`#>   arms      : excluded = Cage, open = Open `\
`#>   reference : open `\
`#>   units     : site `\
`#>   random    : (trt | site)`

The five layouts are `crd`, `rcbd`, `split_plot`, `matched_pairs` and
`multi_site`;
[`?pollin_design`](https://max578.github.io/pollinatoR/reference/pollin_design.md)
gives the random structure each one implies. A random treatment slope
needs enough sites to estimate its variance; with only a handful of
sites, `slope = FALSE` in
[`pollin_design()`](https://max578.github.io/pollinatoR/reference/pollin_design.md)
fits a random intercept alone.

## Estimand

[`pollin_estimand()`](https://max578.github.io/pollinatoR/reference/pollin_estimand.md)
names the number the analysis is for. The dependence (E2) compares
excluded plants with the reference arm; the deficit (E1), available when
a trial has a supplemented arm, compares supplemented with open plants.
The scale is one of `"difference"`, `"ratio"`, `"rii"` or `"lnrr"`.

The basis has no default. A loss share,
$`1 - E[Y_\mathrm{excl}]/E[Y_\mathrm{open}]`$, is the part of the
present yield that pollinators supply. A gain,
$`E[Y_\mathrm{open}]/E[Y_\mathrm{excl}] - 1`$, is the rise over a crop
without them. Excluded plants yielding 75% of open ones give a loss
share of 25% and a gain of 33%. They answer different questions, and the
package records which one every result is.

\
`est`` ``<-`` `[`pollin_estimand`](https://max578.github.io/pollinatoR/reference/pollin_estimand.md)`(``ds``, family ``=`` ``"dependence"``, scale ``=`` ``"ratio"``, basis ``=`` ``"loss"``)`\
`est`\
`#> <pollin_estimand>`\
`#>   family     : dependence (E2) `\
`#>   contrast   : Open vs Cage `\
`#>   scale      : ratio `\
`#>   basis      : loss `\
`#>   population : trial `\
`#>   layout     : multi_site `\
`` #>   assumptions: 6 ; see `pollin_assumptions()` ``\
`knitr``::`[`kable`](https://rdrr.io/pkg/knitr/man/kable.html)`(`[`pollin_assumptions`](https://max578.github.io/pollinatoR/reference/pollin_assumptions.md)`(``est``)``)`

| assumption | statement |
|:---|:---|
| consistency | Each arm is one well-defined treatment, applied the same way on every unit. |
| excludability | The exclusion device changes the response only by excluding pollinators (no cage effect on microclimate, pests or handling). |
| no_interference | One unit’s treatment does not change another unit’s response, for example by displacing pollinators to neighbouring units. |
| assignment | Arms were randomised within each site; the contrast is read within sites. |
| positive_means | The expected response is positive in both arms, so the ratio is defined. |
| population | The estimand refers to the units in the trial; carrying it to other fields or seasons is a further assumption. |

The assumptions are statements about the experiment that the data cannot
confirm. They belong in the methods section of any report.

## Fit and effects

[`pollin_fit()`](https://max578.github.io/pollinatoR/reference/pollin_fit.md)
fits one `glmmTMB` model per response, with the treatment as the fixed
effect and the design’s random terms. The family is chosen from each
response’s values unless supplied.

\
`fit`` ``<-`` `[`pollin_fit`](https://max578.github.io/pollinatoR/reference/pollin_fit.md)`(``est``)`\
`fit`\
`#> <pollin_fit>`\
`#>   engine    : glmmTMB `\
`#>   family    : dependence `\
`#>   contrast  : Open vs Cage `\
`#>   model     : <response> ~ trt + (trt | site) `\
`#>   pods      : nbinom2, 192 rows`\
`#>   seeds_per_pod: nbinom2, 192 rows`

[`pollin_effects()`](https://max578.github.io/pollinatoR/reference/pollin_effects.md)
reads the contrast from each model: every plant is predicted under each
arm, the predictions are averaged, and the scale is applied to the two
averages. A ratio-type estimate is a ratio of averaged predictions, not
an average of plant-level ratios.

\
`eff`` ``<-`` `[`pollin_effects`](https://max578.github.io/pollinatoR/reference/pollin_effects.md)`(``fit``)`\
`eff`\
`#> <pollin_effects>`\
`#>   family   : dependence (E2) `\
`#>   scale    : ratio `\
`#>   basis    : loss `\
`#>   estimate : 1 - Cage/Open (loss share) `\
`#>   interval : 95% delta, fixed and random effects jointly `\
`#>       response     contrast    estimate  std.error   conf.low conf.high`\
`#>           pods Open vs Cage 0.360842985 0.04175635  0.2735314 0.4376609`\
`#>  seeds_per_pod Open vs Cage 0.001949038 0.05577679 -0.1135833 0.1054951`\
`#>       p.value     family scale basis distribution interval`\
`#>  7.312430e-12 dependence ratio  loss      nbinom2    delta`\
`#>  9.721519e-01 dependence ratio  loss      nbinom2    delta`

A parametric bootstrap interval is available with
`interval = "bootstrap"`; *A three-arm trial: deficit and dependence*
shows it.

For the contribution on the difference scale, with a false discovery
rate adjustment across responses and the standard deviation of the
site-level treatment slope,
[`pollin_contribution()`](https://max578.github.io/pollinatoR/reference/pollin_contribution.md)
wraps the fit and the effects.

\
`est_diff`` ``<-`` `[`pollin_estimand`](https://max578.github.io/pollinatoR/reference/pollin_estimand.md)`(`\
`  ``ds``,`\
`  family ``=`` ``"dependence"``,`\
`  scale ``=`` ``"difference"``,`\
`  basis ``=`` ``"loss"`\
`)`\
`contribution`` ``<-`` `[`pollin_contribution`](https://max578.github.io/pollinatoR/reference/pollin_contribution.md)`(``est_diff``)`\
`contribution``$``table``[`\
`  `[`c`](https://rdrr.io/r/base/c.html)`(``"response"``, ``"estimate"``, ``"conf.low"``, ``"conf.high"``, ``"slope_sd"``, ``"p.fdr"``)`\
`]`\
`#>        response   estimate  conf.low conf.high  slope_sd        p.fdr`\
`#> 1          pods 3.70949151  2.637200  4.781783 0.2170980 2.398355e-11`\
`#> 2 seeds_per_pod 0.02988304 -1.647862  1.707628 0.1327942 9.721517e-01`

## Which insects drive the effect

[`pollin_moderation()`](https://max578.github.io/pollinatoR/reference/pollin_moderation.md)
fits treatment-by-exposure interactions, with each exposure split into a
between-site part (the site mean) and a within-site part (the deviation
from it). The within-site interaction is estimated from variation inside
a site, where anything constant at site level, such as an unrecorded
variety, is held fixed.

\
[`pollin_moderation`](https://max578.github.io/pollinatoR/reference/pollin_moderation.md)`(``ds``, outcome ``=`` ``"pods"``)`\
`#> <pollin_moderation>  outcome: pods  (treatment x exposure)`\
`#>                      term    estimate  std.error   p.value            component`\
`#>   trtOpen:honeybee_within -0.06756023 0.06933630 0.3298656      within (robust)`\
`#>  trtOpen:honeybee_between  0.02280788 0.09032579 0.8006493 between (vulnerable)`\
`#>       trtOpen:wild_within  0.01511041 0.06998086 0.8290485      within (robust)`\
`#>      trtOpen:wild_between  0.03201095 0.09144191 0.7262872 between (vulnerable)`\
`#>      p.fdr`\
`#>  0.6597312`\
`#>         NA`\
`#>  0.8290485`\
`#>         NA`

Within a site, honey bee counts here differ by only a bee or two between
plants, too little for this small trial to separate the planted effect
from zero. The within-site terms are the ones to read; the between-site
terms carry whatever differs between sites along with the insects.

The split itself is
[`pollin_within_between()`](https://max578.github.io/pollinatoR/reference/pollin_within_between.md).

\
`wb`` ``<-`` `[`pollin_within_between`](https://max578.github.io/pollinatoR/reference/pollin_within_between.md)`(``canola``$``honeybee``, group ``=`` ``canola``$``site``)`\
[`head`](https://rdrr.io/r/utils/head.html)`(`[`as.data.frame`](https://rdrr.io/r/base/as.data.frame.html)`(``wb``)``)`\
`#>    within between`\
`#> 1 -0.1875  7.1875`\
`#> 2 -0.1875  7.1875`\
`#> 3  0.8125  7.1875`\
`#> 4 -0.1875  7.1875`\
`#> 5 -0.1875  7.1875`\
`#> 6  2.8125  7.1875`

## Guards

Guards stop or warn on common analysis errors.
[`pollin_guard_spurious()`](https://max578.github.io/pollinatoR/reference/pollin_guard_spurious.md)
tests the regression of a dependence index on open-arm yield across
units. The open-arm value sits on both axes, so part of any slope is
built in. The guard permutes the excluded arm while holding the open arm
fixed and compares the observed slope with that null.

\
`open_mean`` ``<-`` `[`tapply`](https://rdrr.io/r/base/tapply.html)`(``canola``$``pods``[``open``]``, ``canola``$``site``[``open``]``, ``mean``)`\
`cage_mean`` ``<-`` `[`tapply`](https://rdrr.io/r/base/tapply.html)`(``canola``$``pods``[``!``open``]``, ``canola``$``site``[``!``open``]``, ``mean``)`\
[`pollin_guard_spurious`](https://max578.github.io/pollinatoR/reference/pollin_guard_spurious.md)`(``open_mean``, ``cage_mean``, n_perm ``=`` ``200``, seed ``=`` ``1``)`\
`#> <pollin_spurious>  index ~ open-yield slope guard`\
`#>   observed slope : 0.0137`\
`#>   coupling null  : mean 0.0528, sd 0.0164 (n_perm = 200)`\
`#>   permutation p  : 1.000  -> slope consistent with mechanical coupling (SPURIOUS risk)`

[`pollin_guard_bounded()`](https://max578.github.io/pollinatoR/reference/pollin_guard_bounded.md)
refuses to average site-level ratios unless asked, because the mean of
ratios is not the ratio of means.

\
[`try`](https://rdrr.io/r/base/try.html)`(`[`pollin_guard_bounded`](https://max578.github.io/pollinatoR/reference/pollin_guard_bounded.md)`(``ds``, index ``=`` ``"rii"``, response ``=`` ``"pods"``)``)`\
`` #> Error : Averaging bounded per-unit ratios is not the ratio of expected values; the mean of 12 unit-level RII values is not the trial's RII. The RII of the arm means is 0.2194. Set `allow = TRUE` to average them anyway. ``

The other two guards are
[`pollin_guard_ecological()`](https://max578.github.io/pollinatoR/reference/pollin_guard_ecological.md),
which warns when an exposure slope exists only between units, and
[`pollin_guard_posttreatment()`](https://max578.github.io/pollinatoR/reference/pollin_guard_posttreatment.md),
which refuses an adjustment covariate declared to be affected by the
treatment. Both are shown in *A matched-pair field study*.

## Descriptive indices

[`pollin_index()`](https://max578.github.io/pollinatoR/reference/pollin_index.md)
and
[`pollin_rii()`](https://max578.github.io/pollinatoR/reference/pollin_rii.md)
give the Relative Interaction Index, the log response ratio and the
dependence ratio from unit-level arm means, with intervals where one is
available. They are for comparison with published work; the model-based
estimate above is the primary one.

\
[`pollin_index`](https://max578.github.io/pollinatoR/reference/pollin_index.md)`(``ds``, response ``=`` ``"pods"``, n_boot ``=`` ``500``, seed ``=`` ``1``)`\
`#>              index  estimate  conf.low conf.high`\
`#> 1              RII 0.2194222 0.1003973 0.3440582`\
`#> 2             lnRR 0.4460981 0.1599727 0.7322235`\
`#> 3 dependence_ratio 0.3598790        NA        NA`\
`#>                                                              note`\
`#> 1                                                    bootstrap CI`\
`#> 2                                                  closed-form CI`\
`#> 3 point; do not regress on open yield (see pollin_guard_spurious)`

[`pollin_compensation()`](https://max578.github.io/pollinatoR/reference/pollin_compensation.md)
summarises how the yield components trade off against each other across
plants.

\
[`pollin_compensation`](https://max578.github.io/pollinatoR/reference/pollin_compensation.md)`(``pd``)`\
`#> <pollin_compensation>  ( spearman correlations )`\
`#>              a    b       cor compensation`\
`#>  seeds_per_pod pods 0.4296379        FALSE`\
`#>   0 negatively-correlated pairs -> candidate compensation`

## Checking the data

[`pollin_validate()`](https://max578.github.io/pollinatoR/reference/pollin_validate.md)
reports missing values in declared columns, expected factor levels and
any identities supplied, as counts only.

\
[`pollin_validate`](https://max578.github.io/pollinatoR/reference/pollin_validate.md)`(``pd``, expect_levels ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(``trt ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"Cage"``, ``"Open"``)``)``)`\
`#> <pollin_validate>`\
`#>       check status    detail`\
`#>        rows   INFO       192`\
`#>        cols   INFO         6`\
`#>  levels:trt   PASS Cage,Open`

## Moving from `pollin_spec()`

`pollin_spec()` was removed at version 0.2.0. It declared the data and
the design in one call; they are now two layers.

\
`# Before`\
`spec`` ``<-`` ``pollin_spec``(`\
`  ``canola``,`\
`  outcomes ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"pods"``, ``"seeds_per_pod"``)``,`\
`  treatment ``=`` ``"trt"``,`\
`  control ``=`` ``"Cage"``,`\
`  random ``=`` ``"site"``,`\
`  exposures ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"honeybee"``, ``"wild"``)``,`\
`  latent_unit ``=`` ``"site"`\
`)`\
[`pollin_contribution`](https://max578.github.io/pollinatoR/reference/pollin_contribution.md)`(``spec``)`\
[`pollin_moderation`](https://max578.github.io/pollinatoR/reference/pollin_moderation.md)`(``spec``, outcome ``=`` ``"pods"``)`

\
`# Now`\
`pd`` ``<-`` `[`pollin_data`](https://max578.github.io/pollinatoR/reference/pollin_data.md)`(`\
`  ``canola``,`\
`  response ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"pods"``, ``"seeds_per_pod"``)``,`\
`  treatment ``=`` ``"trt"``,`\
`  unit ``=`` ``"site"``,`\
`  exposures ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"honeybee"``, ``"wild"``)`\
`)`\
`ds`` ``<-`` `[`pollin_design`](https://max578.github.io/pollinatoR/reference/pollin_design.md)`(``pd``, `[`c`](https://rdrr.io/r/base/c.html)`(``excluded ``=`` ``"Cage"``, open ``=`` ``"Open"``)``, ``"multi_site"``)`\
[`pollin_contribution`](https://max578.github.io/pollinatoR/reference/pollin_contribution.md)`(`[`pollin_estimand`](https://max578.github.io/pollinatoR/reference/pollin_estimand.md)`(``ds``, ``"dependence"``, ``"difference"``, basis ``=`` ``"loss"``)``)`\
[`pollin_moderation`](https://max578.github.io/pollinatoR/reference/pollin_moderation.md)`(``ds``, outcome ``=`` ``"pods"``)`

`outcomes` becomes `response`, `random` becomes `unit`, `control`
becomes the `excluded` entry of `arms`, and `latent_unit` becomes the
`slope` argument of
[`pollin_design()`](https://max578.github.io/pollinatoR/reference/pollin_design.md)
(the `multi_site` layout puts the slope on the site by default).

## Further reading

- *A three-arm trial: deficit and dependence*
- *A split-plot trial: dependence by variety*
- *A matched-pair field study*
- *From dependence to value*

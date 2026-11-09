# A matched-pair field study

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`pollinatoR`](https://github.com/max578/pollinatoR)`)`

## The design

Field studies of crop pollination often work in commercial paddocks
rather than on a trial site. In each field, a caged quadrat is matched
with a nearby open quadrat, and several plants are harvested from each.
Fields differ in soil, sowing date and management, so the comparison is
read within each field, where those differences cancel.

The data below are simulated with the shape of such a study in canola:
24 fields in three regions, one caged and one open quadrat per field,
five plants per quadrat. Honey bee visits are counted on each
open-quadrat plant. Two yield components are recorded per plant: pod
count and seed weight.

\
[`set.seed`](https://rdrr.io/r/base/Random.html)`(``31``)`\
`n_field`` ``<-`` ``24L`\
`n_plant`` ``<-`` ``5L`\
`fields`` ``<-`` `[`sprintf`](https://rdrr.io/r/base/sprintf.html)`(``"F%02d"``, `[`seq_len`](https://rdrr.io/r/base/seq.html)`(``n_field``)``)`\
`fertility`` ``<-`` `[`setNames`](https://rdrr.io/r/stats/setNames.html)`(`[`rnorm`](https://rdrr.io/r/stats/Normal.html)`(``n_field``, ``0``, ``0.25``)``, ``fields``)`\
`region`` ``<-`` `[`setNames`](https://rdrr.io/r/stats/setNames.html)`(`[`rep`](https://rdrr.io/r/base/rep.html)`(`[`c`](https://rdrr.io/r/base/c.html)`(``"North"``, ``"Central"``, ``"South"``)``, each ``=`` ``8L``)``, ``fields``)`\
\
`study`` ``<-`` `[`expand.grid`](https://rdrr.io/r/base/expand.grid.html)`(`\
`  plant ``=`` `[`seq_len`](https://rdrr.io/r/base/seq.html)`(``n_plant``)``,`\
`  trt ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"Caged"``, ``"Open"``)``,`\
`  field ``=`` ``fields``,`\
`  stringsAsFactors ``=`` ``FALSE`\
`)`\
`study``$``region`` ``<-`` ``region``[``study``$``field``]`\
`f`` ``<-`` ``fertility``[``study``$``field``]`\
`open`` ``<-`` ``study``$``trt`` ``==`` ``"Open"`\
\
`# Busier fields are also more fertile; visits do not change yield within a field`\
`study``$``visits`` ``<-`` `[`ifelse`](https://rdrr.io/r/base/ifelse.html)`(``open``, `[`rpois`](https://rdrr.io/r/stats/Poisson.html)`(`[`nrow`](https://rdrr.io/r/base/nrow.html)`(``study``)``, `[`exp`](https://rdrr.io/r/base/Log.html)`(``1.5`` ``+`` ``1.2`` ``*`` ``f``)``)``, ``0``)`\
`pod_mu`` ``<-`` ``120`` ``*`` `[`exp`](https://rdrr.io/r/base/Log.html)`(``f``)`` ``*`` `[`ifelse`](https://rdrr.io/r/base/ifelse.html)`(``open``, ``1``, ``0.8``)`\
`seed_mu`` ``<-`` ``6`` ``*`` `[`exp`](https://rdrr.io/r/base/Log.html)`(``f``)`` ``*`` `[`ifelse`](https://rdrr.io/r/base/ifelse.html)`(``open``, ``1``, ``0.85``)`\
`study``$``pods`` ``<-`` `[`rnbinom`](https://rdrr.io/r/stats/NegBinomial.html)`(`[`nrow`](https://rdrr.io/r/base/nrow.html)`(``study``)``, mu ``=`` ``pod_mu``, size ``=`` ``25``)`\
`study``$``seed_g`` ``<-`` `[`rgamma`](https://rdrr.io/r/stats/GammaDist.html)`(`[`nrow`](https://rdrr.io/r/base/nrow.html)`(``study``)``, shape ``=`` ``20``, rate ``=`` ``20`` ``/`` ``seed_mu``)`\
`study``$``rainfall`` ``<-`` `[`round`](https://rdrr.io/r/base/Round.html)`(``300`` ``+`` ``200`` ``*`` ``f``)`\
`study``$``veg_biomass`` ``<-`` `[`round`](https://rdrr.io/r/base/Round.html)`(`[`rgamma`](https://rdrr.io/r/stats/GammaDist.html)`(`[`nrow`](https://rdrr.io/r/base/nrow.html)`(``study``)``, shape ``=`` ``15``, rate ``=`` ``15`` ``/`` ``(``90`` ``*`` `[`exp`](https://rdrr.io/r/base/Log.html)`(``f``)``)``)``, ``1``)`\
[`head`](https://rdrr.io/r/utils/head.html)`(``study``)`\
`#>   plant   trt field region visits pods   seed_g rainfall veg_biomass`\
`#> 1     1 Caged   F01  North      0   89 7.323259      303        99.2`\
`#> 2     2 Caged   F01  North      0  114 5.490651      303        93.1`\
`#> 3     3 Caged   F01  North      0   97 4.727265      303        95.3`\
`#> 4     4 Caged   F01  North      0   88 6.071478      303        83.0`\
`#> 5     5 Caged   F01  North      0  103 5.830595      303       139.6`\
`#> 6     1  Open   F01  North      5  123 6.222178      303        78.1`

Pollinator access was planted to add 25% to pod count over the caged
plants (a loss share of $`1 - 0.8 = 0.20`$) and about 18% to seed weight
(a loss share of 0.15).

## Data, design and estimand

Each field is one matched pair. In a `matched_pairs` layout the pair
carries a random intercept, and the design warns if any pair lacks an
arm.

\
`pd`` ``<-`` `[`pollin_data`](https://max578.github.io/pollinatoR/reference/pollin_data.md)`(`\
`  ``study``,`\
`  response ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"pods"``, ``"seed_g"``)``,`\
`  treatment ``=`` ``"trt"``,`\
`  unit ``=`` ``"field"``,`\
`  groups ``=`` ``"region"``,`\
`  exposures ``=`` ``"visits"``,`\
`  adjust ``=`` ``"rainfall"`\
`)`\
`pd`\
`#> <pollin_data>`\
`#>   rows        : 240 `\
`#>   response    : pods, seed_g `\
`#>   treatment   : trt ( Caged, Open )`\
`#>   unit        : field `\
`#>   groups      : region `\
`#>   exposures   : visits `\
`#>   adjust      : rainfall `\
`` #>   checks      : 0 warnings; see `$checks` ``\
`ds`` ``<-`` `[`pollin_design`](https://max578.github.io/pollinatoR/reference/pollin_design.md)`(`\
`  ``pd``,`\
`  arms ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``excluded ``=`` ``"Caged"``, open ``=`` ``"Open"``)``,`\
`  layout ``=`` ``"matched_pairs"`\
`)`\
`ds`\
`#> <pollin_design>`\
`#>   layout    : matched_pairs `\
`#>   arms      : excluded = Caged, open = Open `\
`#>   reference : open `\
`#>   units     : field `\
`#>   random    : (1 | field)`

The dependence compares caged plants with open ones. The assumptions for
this layout include one that a randomised plot trial does not need: the
two quadrats in a field must be alike apart from the cage.

\
`est_loss`` ``<-`` `[`pollin_estimand`](https://max578.github.io/pollinatoR/reference/pollin_estimand.md)`(``ds``, family ``=`` ``"dependence"``, scale ``=`` ``"ratio"``, basis ``=`` ``"loss"``)`\
`knitr``::`[`kable`](https://rdrr.io/pkg/knitr/man/kable.html)`(`[`pollin_assumptions`](https://max578.github.io/pollinatoR/reference/pollin_assumptions.md)`(``est_loss``)``)`

| assumption | statement |
|:---|:---|
| consistency | Each arm is one well-defined treatment, applied the same way on every unit. |
| excludability | The exclusion device changes the response only by excluding pollinators (no cage effect on microclimate, pests or handling). |
| no_interference | One unit’s treatment does not change another unit’s response, for example by displacing pollinators to neighbouring units. |
| assignment | Arms were randomised within each pair, and the pair members are alike apart from the arm. |
| positive_means | The expected response is positive in both arms, so the ratio is defined. |
| population | The estimand refers to the units in the trial; carrying it to other fields or seasons is a further assumption. |

## Loss share and gain, side by side

A loss share answers “what part of the yield now harvested do
pollinators supply?”. A gain answers “how much would a pollinator-free
crop rise if they returned?”. For an economic valuation of pollination
the loss share is the right number. Both come from the same pair of
expected values, so the gain is always the larger of the two: with a
loss share $`D`$, the gain is $`D/(1 - D)`$.

\
`est_gain`` ``<-`` `[`pollin_estimand`](https://max578.github.io/pollinatoR/reference/pollin_estimand.md)`(``ds``, family ``=`` ``"dependence"``, scale ``=`` ``"ratio"``, basis ``=`` ``"gain"``)`\
`loss`` ``<-`` `[`pollin_effects`](https://max578.github.io/pollinatoR/reference/pollin_effects.md)`(`[`pollin_fit`](https://max578.github.io/pollinatoR/reference/pollin_fit.md)`(``est_loss``)``)`\
`gain`` ``<-`` `[`pollin_effects`](https://max578.github.io/pollinatoR/reference/pollin_effects.md)`(`[`pollin_fit`](https://max578.github.io/pollinatoR/reference/pollin_fit.md)`(``est_gain``)``)`\
`loss`\
`#> <pollin_effects>`\
`#>   family   : dependence (E2) `\
`#>   scale    : ratio `\
`#>   basis    : loss `\
`#>   estimate : 1 - Caged/Open (loss share) `\
`#>   interval : 95% delta `\
`#>   note     : delta intervals cover fixed-effect uncertainty only`\
`#>  response      contrast  estimate  std.error  conf.low conf.high      p.value`\
`#>      pods Open vs Caged 0.2020465 0.02339658 0.1548468 0.2466103 1.384345e-14`\
`#>    seed_g Open vs Caged 0.1585004 0.02456385 0.1089523 0.2052934 3.383183e-09`\
`#>      family scale basis distribution interval`\
`#>  dependence ratio  loss      nbinom2    delta`\
`#>  dependence ratio  loss        Gamma    delta`\
[`rbind`](https://rdrr.io/r/base/cbind.html)`(``loss``$``table``, ``gain``$``table``)``[``, `[`c`](https://rdrr.io/r/base/c.html)`(``"response"``, ``"basis"``, ``"estimate"``, ``"conf.low"``, ``"conf.high"``)``]`\
`#>   response basis  estimate  conf.low conf.high`\
`#> 1     pods  loss 0.2020465 0.1548468 0.2466103`\
`#> 2   seed_g  loss 0.1585004 0.1089523 0.2052934`\
`#> 3     pods  gain 0.2532059 0.1832175 0.3273342`\
`#> 4   seed_g  gain 0.1883547 0.1222744 0.2583260`

Each response gets its own model and family: a negative binomial for the
pod counts and a Gamma with a log link for seed weight, chosen from the
values.

## Further scales

The difference scale keeps the response’s units (pods per plant, grams
per plant), which is what a yield budget needs. The RII and the log
response ratio are the scales most often used to compare studies.

\
`scale_rows`` ``<-`` `[`lapply`](https://rdrr.io/r/base/lapply.html)`(`[`c`](https://rdrr.io/r/base/c.html)`(``"difference"``, ``"rii"``, ``"lnrr"``)``, ``function``(``s``)`` ``{`\
`  ``est`` ``<-`` `[`pollin_estimand`](https://max578.github.io/pollinatoR/reference/pollin_estimand.md)`(``ds``, family ``=`` ``"dependence"``, scale ``=`` ``s``, basis ``=`` ``"loss"``)`\
`  `[`pollin_effects`](https://max578.github.io/pollinatoR/reference/pollin_effects.md)`(`[`pollin_fit`](https://max578.github.io/pollinatoR/reference/pollin_fit.md)`(``est``)``)``$``table`\
`}``)`\
[`do.call`](https://rdrr.io/r/base/do.call.html)`(``rbind``, ``scale_rows``)``[``, `[`c`](https://rdrr.io/r/base/c.html)`(``"response"``, ``"scale"``, ``"estimate"``, ``"conf.low"``, ``"conf.high"``)``]`\
`#>   response      scale   estimate    conf.low  conf.high`\
`#> 1     pods difference 24.3988155 17.79598211 31.0016488`\
`#> 2   seed_g difference  0.9644760  0.63345191  1.2955001`\
`#> 3     pods        rii  0.1123758  0.08392086  0.1406477`\
`#> 4   seed_g        rii  0.0860714  0.05761479  0.1143883`\
`#> 5     pods       lnrr  0.2257050  0.16823742  0.2831726`\
`#> 6   seed_g       lnrr  0.1725698  0.11535733  0.2297822`

For comparison with earlier work,
[`pollin_index()`](https://max578.github.io/pollinatoR/reference/pollin_index.md)
gives the descriptive indices from the field-level arm means, with
intervals.

\
[`pollin_index`](https://max578.github.io/pollinatoR/reference/pollin_index.md)`(``ds``, response ``=`` ``"pods"``, n_boot ``=`` ``500``, seed ``=`` ``1``)`\
`#>              index  estimate   conf.low conf.high`\
`#> 1              RII 0.1156905 0.04546906 0.1754487`\
`#> 2             lnRR 0.2324216 0.10074084 0.3641024`\
`#> 3 dependence_ratio 0.2073881         NA        NA`\
`#>                                                              note`\
`#> 1                                                    bootstrap CI`\
`#> 2                                                  closed-form CI`\
`#> 3 point; do not regress on open yield (see pollin_guard_spurious)`

## From the sampled fields to the growing area

The 24 fields are eight per region, but the regions grow different areas
of canola. Suppose the North grows 2,000 ha, the Central region 1,000 ha
and the South 500 ha.
[`pollin_population()`](https://max578.github.io/pollinatoR/reference/pollin_population.md)
declares the regions as strata weighted by area. The model is fitted
with the region as `by`, and
[`pollin_draws()`](https://max578.github.io/pollinatoR/reference/pollin_draws.md)
returns the loss share per region and for the population.

\
`regions`` ``<-`` `[`pollin_population`](https://max578.github.io/pollinatoR/reference/pollin_population.md)`(`\
`  ``"region"``,`\
`  `[`c`](https://rdrr.io/r/base/c.html)`(``North ``=`` ``2000``, Central ``=`` ``1000``, South ``=`` ``500``)``,`\
`  label ``=`` ``"canola area by region"`\
`)`\
`est_pop`` ``<-`` `[`pollin_estimand`](https://max578.github.io/pollinatoR/reference/pollin_estimand.md)`(``ds``, ``"dependence"``, ``"ratio"``, ``"loss"``, population ``=`` ``regions``)`\
`dr`` ``<-`` `[`pollin_draws`](https://max578.github.io/pollinatoR/reference/pollin_draws.md)`(`\
`  `[`pollin_fit`](https://max578.github.io/pollinatoR/reference/pollin_fit.md)`(``est_pop``, by ``=`` ``"region"``)``,`\
`  n ``=`` ``1000L``,`\
`  response ``=`` ``"seed_g"``,`\
`  seed ``=`` ``1L`\
`)`\
[`summary`](https://rdrr.io/r/base/summary.html)`(``dr``)``[``, `[`c`](https://rdrr.io/r/base/c.html)`(``"level"``, ``"median"``, ``"conf.low"``, ``"conf.high"``)``]`\
`#>        level     median     conf.low conf.high`\
`#> 1    Central 0.25858044  0.180760050 0.3318839`\
`#> 2      North 0.09242441 -0.003609891 0.1755858`\
`#> 3      South 0.11473639  0.022422333 0.1999879`\
`#> 4 population 0.13594660  0.075141990 0.1946263`

The population loss share is
$`1 - \sum w\,\mathrm{caged} / \sum w\,\mathrm{open}`$, a ratio of the
area-weighted arm means. It is not the area-weighted mean of the
regional loss shares, shown below for comparison.

\
`w`` ``<-`` ``regions``$``weights`\
`weighted_mean`` ``<-`` ``dr``$``estimate``[``, `[`names`](https://rdrr.io/r/base/names.html)`(``w``)``]`` `[`%*%`](https://rdrr.io/r/base/matmult.html)` ``w`\
[`quantile`](https://rdrr.io/r/stats/quantile.html)`(``weighted_mean``, `[`c`](https://rdrr.io/r/base/c.html)`(``0.5``, ``0.025``, ``0.975``)``)`\
`#>        50%       2.5%      97.5% `\
`#> 0.14228785 0.08327175 0.19839943`

The ratio of means weights each region by its area times its open yield,
that is, by its share of the seed harvested. A loss share is a share of
harvested seed, so that is the weighting it needs. In these estimates
the North yields most and loses least, so the weighted mean of the
regional shares comes out higher. (The regions were simulated with the
same loss share; their estimates differ by chance.) The gap widens as
regions differ more in yield and in dependence.

## Guard: a between-field slope is not a plant-level effect

The fields with more honey bee visits also yield more. Read across
fields, that looks like a pollinator effect on yield. In these data it
is not: the busier fields are the more fertile ones, and within a field
the plants with more visits yield no more than the others.
[`pollin_guard_ecological()`](https://max578.github.io/pollinatoR/reference/pollin_guard_ecological.md)
splits the exposure into its between-field and within-field parts and
warns when only the between-field slope is distinguishable from zero.

\
`open_plants`` ``<-`` ``study``[``study``$``trt`` ``==`` ``"Open"``, ``]`\
[`pollin_guard_ecological`](https://max578.github.io/pollinatoR/reference/pollin_guard_ecological.md)`(`\
`  ``open_plants``,`\
`  exposure ``=`` ``"visits"``,`\
`  response ``=`` ``"pods"``,`\
`  unit ``=`` ``"field"`\
`)`\
`#> Warning: Only the between-unit slope is distinguishable from zero (p = 1.44e-10`\
`#> between, 0.578 within); the association is between units, not within them.`\
`#> <pollin_ecological>  pods ~ visits split by field (24 units)`\
`#>  component   slope std.error   p.value`\
`#>     pooled  5.6100     1.454 1.873e-04`\
`#>     within  0.8665     1.552 5.776e-01`\
`#>    between 16.6666     2.369 1.443e-10`\
`#>   ! Only the between-unit slope is distinguishable from zero (p = 1.44e-10 between, 0.578 within); the association is between units, not within them.`

## Guard: adjusting for a downstream covariate

A cage does not change the rainfall a field receives, so adjusting for
rainfall is safe. Vegetative biomass is different: plants that set fewer
pods may put more into leaves and stems, so biomass sits downstream of
the treatment. Adjusting for it would remove part of the pollination
effect from the contrast. Whether a variable is downstream is a
judgement about the crop that the data cannot settle, so
[`pollin_guard_posttreatment()`](https://max578.github.io/pollinatoR/reference/pollin_guard_posttreatment.md)
takes it as a declaration and checks the declared adjustment covariates
against it.

\
[`pollin_guard_posttreatment`](https://max578.github.io/pollinatoR/reference/pollin_guard_posttreatment.md)`(``pd``, downstream ``=`` ``"veg_biomass"``)`\
\
`pd_bad`` ``<-`` `[`pollin_data`](https://max578.github.io/pollinatoR/reference/pollin_data.md)`(`\
`  ``study``,`\
`  response ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"pods"``, ``"seed_g"``)``,`\
`  treatment ``=`` ``"trt"``,`\
`  unit ``=`` ``"field"``,`\
`  adjust ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"rainfall"``, ``"veg_biomass"``)`\
`)`\
[`try`](https://rdrr.io/r/base/try.html)`(`[`pollin_guard_posttreatment`](https://max578.github.io/pollinatoR/reference/pollin_guard_posttreatment.md)`(``pd_bad``, downstream ``=`` ``"veg_biomass"``)``)`\
`` #> Error : Adjustment covariate `veg_biomass` is declared downstream of `trt`. Adjusting for a downstream covariate removes part of the treatment effect; drop it from `adjust`. ``

The first call returns silently because rainfall is not declared
downstream. The second stops before any analysis uses the biomass
adjustment.

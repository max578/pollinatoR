# A three-arm trial: deficit and dependence

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`pollinatoR`](https://github.com/max578/pollinatoR)`)`

## The design

A three-arm trial puts three pollination treatments side by side: plants
*excluded* from insects (caged or bagged), plants left *open* to
whatever pollinators visit, and plants *supplemented* by hand with extra
pollen. The three arms answer two different questions.

- The deficit (E1) compares supplemented with open plants. It asks how
  much yield was left on the table by the pollination the crop actually
  received.
- The dependence (E2) compares excluded plants with a reference arm. It
  asks how much yield is owed to pollinators. The reference is usually
  the open arm; in a three-arm trial it can also be the supplemented
  arm, which gives the dependence under full pollination.

This vignette simulates a small randomised complete block trial: eight
blocks, each holding all three arms, with six plants per arm in each
block. The response is pod count per plant.

\
[`set.seed`](https://rdrr.io/r/base/Random.html)`(``2024``)`\
`n_block`` ``<-`` ``8L`\
`n_plant`` ``<-`` ``6L`\
`arms`` ``<-`` `[`c`](https://rdrr.io/r/base/c.html)`(``"Bagged"``, ``"Open"``, ``"Hand"``)`\
`trial`` ``<-`` `[`expand.grid`](https://rdrr.io/r/base/expand.grid.html)`(`\
`  plant ``=`` `[`seq_len`](https://rdrr.io/r/base/seq.html)`(``n_plant``)``,`\
`  trt ``=`` ``arms``,`\
`  block ``=`` `[`sprintf`](https://rdrr.io/r/base/sprintf.html)`(``"B%d"``, `[`seq_len`](https://rdrr.io/r/base/seq.html)`(``n_block``)``)``,`\
`  stringsAsFactors ``=`` ``FALSE`\
`)`\
`block_effect`` ``<-`` `[`rnorm`](https://rdrr.io/r/stats/Normal.html)`(``n_block``, ``0``, ``0.15``)`\
`arm_mean`` ``<-`` `[`c`](https://rdrr.io/r/base/c.html)`(``Bagged ``=`` ``40``, Open ``=`` ``60``, Hand ``=`` ``70``)`\
[`names`](https://rdrr.io/r/base/names.html)`(``block_effect``)`` ``<-`` `[`sprintf`](https://rdrr.io/r/base/sprintf.html)`(``"B%d"``, `[`seq_len`](https://rdrr.io/r/base/seq.html)`(``n_block``)``)`\
`mu`` ``<-`` ``arm_mean``[``trial``$``trt``]`` ``*`` `[`exp`](https://rdrr.io/r/base/Log.html)`(``block_effect``[``trial``$``block``]``)`\
`trial``$``pods`` ``<-`` `[`rnbinom`](https://rdrr.io/r/stats/NegBinomial.html)`(`[`nrow`](https://rdrr.io/r/base/nrow.html)`(``trial``)``, mu ``=`` ``mu``, size ``=`` ``30``)`\
[`head`](https://rdrr.io/r/utils/head.html)`(``trial``)`\
`#>   plant    trt block pods`\
`#> 1     1 Bagged    B1   32`\
`#> 2     2 Bagged    B1   55`\
`#> 3     3 Bagged    B1   54`\
`#> 4     4 Bagged    B1   66`\
`#> 5     5 Bagged    B1   65`\
`#> 6     6 Bagged    B1   55`

The planted arm means are 40, 60 and 70 pods. On these values the
dependence on the open arm is a loss share of $`1 - 40/60 = 0.33`$, and
the deficit is a loss share of $`1 - 60/70 = 0.14`$.

## Data, design and estimand

The data step names the columns; the design step maps each treatment
level to its role and fixes the layout. In a randomised complete block
layout the block carries a random intercept.

\
`pd`` ``<-`` `[`pollin_data`](https://max578.github.io/pollinatoR/reference/pollin_data.md)`(``trial``, response ``=`` ``"pods"``, treatment ``=`` ``"trt"``, unit ``=`` ``"block"``)`\
`ds`` ``<-`` `[`pollin_design`](https://max578.github.io/pollinatoR/reference/pollin_design.md)`(`\
`  ``pd``,`\
`  arms ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``excluded ``=`` ``"Bagged"``, open ``=`` ``"Open"``, supplemented ``=`` ``"Hand"``)``,`\
`  layout ``=`` ``"rcbd"`\
`)`\
`ds`\
`#> <pollin_design>`\
`#>   layout    : rcbd `\
`#>   arms      : excluded = Bagged, open = Open, supplemented = Hand `\
`#>   reference : open `\
`#>   units     : block `\
`#>   random    : (1 | block)`

The estimand step names the contrast, the scale and the basis. Neither
the scale nor the basis has a default.

\
`deficit`` ``<-`` `[`pollin_estimand`](https://max578.github.io/pollinatoR/reference/pollin_estimand.md)`(``ds``, family ``=`` ``"deficit"``, scale ``=`` ``"ratio"``, basis ``=`` ``"loss"``)`\
`deficit`\
`#> <pollin_estimand>`\
`#>   family     : deficit (E1) `\
`#>   contrast   : Hand vs Open `\
`#>   scale      : ratio `\
`#>   basis      : loss `\
`#>   population : trial `\
`#>   layout     : rcbd `\
`` #>   assumptions: 6 ; see `pollin_assumptions()` ``\
`knitr``::`[`kable`](https://rdrr.io/pkg/knitr/man/kable.html)`(`[`pollin_assumptions`](https://max578.github.io/pollinatoR/reference/pollin_assumptions.md)`(``deficit``)``)`

| assumption | statement |
|:---|:---|
| consistency | Each arm is one well-defined treatment, applied the same way on every unit. |
| excludability | Hand supplementation changes the response only through the pollen added (no damage from handling). |
| no_interference | One unit’s treatment does not change another unit’s response, for example by displacing pollinators to neighbouring units. |
| assignment | Arms were randomised to units within each block. |
| positive_means | The expected response is positive in both arms, so the ratio is defined. |
| population | The estimand refers to the units in the trial; carrying it to other fields or seasons is a further assumption. |

The assumption table is the part of the analysis the data cannot check.
For the deficit, the key one is that hand supplementation adds pollen
and nothing else; for the dependence, it is that the bag or cage
excludes insects without changing the plant’s microclimate.

## Loss share or gain

The same pair of arm means can be reported in two ways. The loss share
is the fraction of the higher arm’s yield that is lost without the
contrast, $`1 - \mathrm{low}/\mathrm{high}`$. The gain is the increase
over the lower arm, $`\mathrm{high}/\mathrm{low} - 1`$. For the
dependence planted here, the loss share is 0.33 and the gain is 0.50.
The two answer different questions: the loss share is the part of
current yield that pollinators supply, which is what a valuation of
their service needs, and the gain is how much a crop without them would
rise if they returned. Reporting a gain where a loss share was meant
overstates the value by a factor of $`1/(1 - D)`$.

`basis` is therefore stated every time, and printed on every result.

\
`fit_effect`` ``<-`` ``function``(``design``, ``family``, ``scale``, ``basis``)`` ``{`\
`  ``est`` ``<-`` `[`pollin_estimand`](https://max578.github.io/pollinatoR/reference/pollin_estimand.md)`(``design``, family ``=`` ``family``, scale ``=`` ``scale``, basis ``=`` ``basis``)`\
`  `[`pollin_effects`](https://max578.github.io/pollinatoR/reference/pollin_effects.md)`(`[`pollin_fit`](https://max578.github.io/pollinatoR/reference/pollin_fit.md)`(``est``)``)``$``table`\
`}`\
`dep_loss`` ``<-`` ``fit_effect``(``ds``, ``"dependence"``, ``"ratio"``, ``"loss"``)`\
`dep_gain`` ``<-`` ``fit_effect``(``ds``, ``"dependence"``, ``"ratio"``, ``"gain"``)`\
[`rbind`](https://rdrr.io/r/base/cbind.html)`(``dep_loss``, ``dep_gain``)``[``, `[`c`](https://rdrr.io/r/base/c.html)`(``"contrast"``, ``"basis"``, ``"estimate"``, ``"conf.low"``, ``"conf.high"``)``]`\
`#>         contrast basis  estimate  conf.low conf.high`\
`#> 1 Open vs Bagged  loss 0.3059873 0.2411880 0.3652530`\
`#> 2 Open vs Bagged  gain 0.4408958 0.3178495 0.5754309`

The interval for a ratio-type scale is a delta interval for the log
ratio of the averaged predictions, transformed to the scale, so a loss
share stays below one.

Both intervals cover the planted values. A parametric bootstrap gives a
second interval: responses are simulated from the fitted model, the
model is refitted to each simulated set, and the interval is taken from
the percentiles of the refitted estimates. Use several hundred refits in
an analysis; the 100 here keep the example quick.

\
`est_dep`` ``<-`` `[`pollin_estimand`](https://max578.github.io/pollinatoR/reference/pollin_estimand.md)`(``ds``, family ``=`` ``"dependence"``, scale ``=`` ``"ratio"``, basis ``=`` ``"loss"``)`\
[`pollin_effects`](https://max578.github.io/pollinatoR/reference/pollin_effects.md)`(`[`pollin_fit`](https://max578.github.io/pollinatoR/reference/pollin_fit.md)`(``est_dep``)``, interval ``=`` ``"bootstrap"``, n_boot ``=`` ``100L``, seed ``=`` ``1L``)`\
`#> <pollin_effects>`\
`#>   family   : dependence (E2) `\
`#>   scale    : ratio `\
`#>   basis    : loss `\
`#>   estimate : 1 - Bagged/Open (loss share) `\
`#>   interval : 95% parametric bootstrap (100 refits) `\
`#>  response       contrast  estimate  std.error  conf.low conf.high p.value`\
`#>      pods Open vs Bagged 0.3059873 0.03188554 0.2424225 0.3655036      NA`\
`#>      family scale basis distribution  interval`\
`#>  dependence ratio  loss      nbinom2 bootstrap`

## Both estimands on several scales

The estimate is always formed from the model’s expected response in each
arm, averaged over the trial’s plants: a ratio of averaged predictions,
never an average of plant-level ratios. The scale only changes how that
pair of expected values is reported.

\
`grid`` ``<-`` `[`expand.grid`](https://rdrr.io/r/base/expand.grid.html)`(`\
`  family ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"deficit"``, ``"dependence"``)``,`\
`  scale ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"difference"``, ``"ratio"``, ``"rii"``, ``"lnrr"``)``,`\
`  stringsAsFactors ``=`` ``FALSE`\
`)`\
`rows`` ``<-`` `[`lapply`](https://rdrr.io/r/base/lapply.html)`(`[`seq_len`](https://rdrr.io/r/base/seq.html)`(`[`nrow`](https://rdrr.io/r/base/nrow.html)`(``grid``)``)``, ``function``(``i``)`` ``{`\
`  ``fit_effect``(``ds``, ``grid``$``family``[``i``]``, ``grid``$``scale``[``i``]``, basis ``=`` ``"loss"``)`\
`}``)`\
`effects`` ``<-`` `[`do.call`](https://rdrr.io/r/base/do.call.html)`(``rbind``, ``rows``)`\
`effects``[``, `[`c`](https://rdrr.io/r/base/c.html)`(``"family"``, ``"contrast"``, ``"scale"``, ``"estimate"``, ``"conf.low"``, ``"conf.high"``)``]`\
`#>       family       contrast      scale    estimate    conf.low  conf.high`\
`#> 1    deficit   Hand vs Open difference  8.81219559  2.83154691 14.7928443`\
`#> 2 dependence Open vs Bagged difference 19.94446792 14.75439295 25.1345429`\
`#> 3    deficit   Hand vs Open      ratio  0.11909515  0.04065030  0.1911257`\
`#> 4 dependence Open vs Bagged      ratio  0.30598730  0.24118799  0.3652530`\
`#> 5    deficit   Hand vs Open        rii  0.06331801  0.02074683  0.1056600`\
`#> 6 dependence Open vs Bagged        rii  0.18062869  0.13713119  0.2234309`\
`#> 7    deficit   Hand vs Open       lnrr  0.12680566  0.04149962  0.2121117`\
`#> 8 dependence Open vs Bagged       lnrr  0.36526501  0.27600121  0.4545288`

The RII and the loss share are linked by $`\mathrm{RII} = D/(2 - D)`$,
and the log response ratio is $`\log(\mathrm{high}/\mathrm{low})`$, so
the three ratio-type rows for each family carry the same information on
different scales. All rows here use `basis = "loss"`; the basis changes
the value only on the ratio scale. The difference scale is in pods per
plant.

## Dependence against the supplemented arm

With `reference = "supplemented"` the dependence compares excluded
plants with hand-pollinated ones. This is the dependence under full
pollination, an upper reference for the crop’s need for pollinators,
while the open reference gives the dependence realised in the field.

\
`ds_sup`` ``<-`` `[`pollin_design`](https://max578.github.io/pollinatoR/reference/pollin_design.md)`(`\
`  ``pd``,`\
`  arms ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``excluded ``=`` ``"Bagged"``, open ``=`` ``"Open"``, supplemented ``=`` ``"Hand"``)``,`\
`  layout ``=`` ``"rcbd"``,`\
`  reference ``=`` ``"supplemented"`\
`)`\
`est_sup`` ``<-`` `[`pollin_estimand`](https://max578.github.io/pollinatoR/reference/pollin_estimand.md)`(``ds_sup``, ``"dependence"``, scale ``=`` ``"ratio"``, basis ``=`` ``"loss"``)`\
[`pollin_effects`](https://max578.github.io/pollinatoR/reference/pollin_effects.md)`(`[`pollin_fit`](https://max578.github.io/pollinatoR/reference/pollin_fit.md)`(``est_sup``)``)`\
`#> <pollin_effects>`\
`#>   family   : dependence (E2) `\
`#>   scale    : ratio `\
`#>   basis    : loss `\
`#>   estimate : 1 - Bagged/Hand (loss share) `\
`#>   interval : 95% delta `\
`#>   note     : delta intervals cover fixed-effect uncertainty only`\
`#>  response       contrast  estimate  std.error  conf.low conf.high      p.value`\
`#>      pods Hand vs Bagged 0.3886408 0.02758782 0.3321065 0.4403898 1.096722e-27`\
`#>      family scale basis distribution interval`\
`#>  dependence ratio  loss      nbinom2    delta`

## A guard on plant-level ratios

A common shortcut computes the RII or the dependence ratio for each
block (or each plant) and averages them. The mean of bounded unit-level
ratios is not the ratio of the expected values, and
[`pollin_guard_bounded()`](https://max578.github.io/pollinatoR/reference/pollin_guard_bounded.md)
refuses it unless asked.

\
[`try`](https://rdrr.io/r/base/try.html)`(`[`pollin_guard_bounded`](https://max578.github.io/pollinatoR/reference/pollin_guard_bounded.md)`(``ds``, index ``=`` ``"dependence"``)``)`\
`` #> Error : Averaging bounded per-unit ratios is not the ratio of expected values; the mean of 8 unit-level dependence values is not the trial's dependence. The dependence of the arm means is 0.3075. Set `allow = TRUE` to average them anyway. ``

When the average is wanted for comparison with earlier work,
`allow = TRUE` returns it with the statement attached and the ratio of
the arm means beside it.

\
[`pollin_guard_bounded`](https://max578.github.io/pollinatoR/reference/pollin_guard_bounded.md)`(``ds``, index ``=`` ``"dependence"``, allow ``=`` ``TRUE``)`\
`#> <pollin_bounded>  average of unit-level dependence values`\
`#>   mean of 8 units : 0.3004`\
`#>   dependence of arm means : 0.3075`\
`#>   Averaging bounded per-unit ratios is not the ratio of expected values; the mean of 8 unit-level dependence values is not the trial's dependence.`

Here the block means differ only a little, so the two numbers are close.
They separate when units differ widely in yield, which is common across
sites and varieties.

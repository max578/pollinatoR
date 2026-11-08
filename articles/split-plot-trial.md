# A split-plot trial: dependence by variety

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`pollinatoR`](https://github.com/max578/pollinatoR)`)`

## The design

Variety trials are often laid out as split plots. Each block holds one
whole plot per variety, and each whole plot is split in two: one half
caged against insects, the other left open. The pollination treatment is
the sub-plot factor. Two sources of error follow: whole plots differ
from one another, and so do sub-plots within a whole plot. The
pollination contrast is read within whole plots, where the variety and
the patch of ground are shared.

The trial below is simulated: four blocks, eight canola varieties, and
three plants harvested from each sub-plot. Varieties differ both in
yield and in how much they depend on insects, and the two are unrelated.

\
[`set.seed`](https://rdrr.io/r/base/Random.html)`(``11``)`\
`n_block`` ``<-`` ``4L`\
`varieties`` ``<-`` `[`sprintf`](https://rdrr.io/r/base/sprintf.html)`(``"V%d"``, ``1``:``8``)`\
`open_pods`` ``<-`` `[`setNames`](https://rdrr.io/r/stats/setNames.html)`(`[`c`](https://rdrr.io/r/base/c.html)`(``70``, ``95``, ``80``, ``110``, ``60``, ``90``, ``100``, ``75``)``, ``varieties``)`\
`dependence`` ``<-`` `[`c`](https://rdrr.io/r/base/c.html)`(``0.10``, ``0.30``, ``0.20``, ``0.15``, ``0.35``, ``0.25``, ``0.05``, ``0.30``)`\
[`names`](https://rdrr.io/r/base/names.html)`(``dependence``)`` ``<-`` ``varieties`\
\
`trial`` ``<-`` `[`expand.grid`](https://rdrr.io/r/base/expand.grid.html)`(`\
`  plant ``=`` ``1``:``3``,`\
`  trt ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"Cage"``, ``"Open"``)``,`\
`  variety ``=`` ``varieties``,`\
`  block ``=`` `[`sprintf`](https://rdrr.io/r/base/sprintf.html)`(``"B%d"``, `[`seq_len`](https://rdrr.io/r/base/seq.html)`(``n_block``)``)``,`\
`  stringsAsFactors ``=`` ``FALSE`\
`)`\
`# Varieties are randomised to the whole plots of each block`\
`whole_plot`` ``<-`` `[`unlist`](https://rdrr.io/r/base/unlist.html)`(`[`lapply`](https://rdrr.io/r/base/lapply.html)`(`[`seq_len`](https://rdrr.io/r/base/seq.html)`(``n_block``)``, ``function``(``b``)`` `[`sample`](https://rdrr.io/r/base/sample.html)`(``8L``)``)``)`\
`blocks`` ``<-`` `[`sprintf`](https://rdrr.io/r/base/sprintf.html)`(``"B%d"``, `[`seq_len`](https://rdrr.io/r/base/seq.html)`(``n_block``)``)`\
[`names`](https://rdrr.io/r/base/names.html)`(``whole_plot``)`` ``<-`` `[`paste`](https://rdrr.io/r/base/paste.html)`(`[`rep`](https://rdrr.io/r/base/rep.html)`(``blocks``, each ``=`` ``8L``)``, ``varieties``)`\
`trial``$``plot`` ``<-`` `[`sprintf`](https://rdrr.io/r/base/sprintf.html)`(``"W%d"``, ``whole_plot``[`[`paste`](https://rdrr.io/r/base/paste.html)`(``trial``$``block``, ``trial``$``variety``)``]``)`\
\
`block_effect`` ``<-`` `[`setNames`](https://rdrr.io/r/stats/setNames.html)`(`[`rnorm`](https://rdrr.io/r/stats/Normal.html)`(``n_block``, ``0``, ``0.10``)``, ``blocks``)`\
`plot_effect`` ``<-`` `[`setNames`](https://rdrr.io/r/stats/setNames.html)`(`[`rnorm`](https://rdrr.io/r/stats/Normal.html)`(`[`length`](https://rdrr.io/r/base/length.html)`(``whole_plot``)``, ``0``, ``0.10``)``, `[`names`](https://rdrr.io/r/base/names.html)`(``whole_plot``)``)`\
`cage`` ``<-`` ``trial``$``trt`` ``==`` ``"Cage"`\
`mu`` ``<-`` ``open_pods``[``trial``$``variety``]`` ``*`\
`  `[`ifelse`](https://rdrr.io/r/base/ifelse.html)`(``cage``, ``1`` ``-`` ``dependence``[``trial``$``variety``]``, ``1``)`` ``*`\
`  `[`exp`](https://rdrr.io/r/base/Log.html)`(``block_effect``[``trial``$``block``]`` ``+`` ``plot_effect``[`[`paste`](https://rdrr.io/r/base/paste.html)`(``trial``$``block``, ``trial``$``variety``)``]``)`\
`trial``$``pods`` ``<-`` `[`rnbinom`](https://rdrr.io/r/stats/NegBinomial.html)`(`[`nrow`](https://rdrr.io/r/base/nrow.html)`(``trial``)``, mu ``=`` ``mu``, size ``=`` ``30``)`\
[`head`](https://rdrr.io/r/utils/head.html)`(``trial``)`\
`#>   plant  trt variety block plot pods`\
`#> 1     1 Cage      V1    B1   W2   60`\
`#> 2     2 Cage      V1    B1   W2   49`\
`#> 3     3 Cage      V1    B1   W2   46`\
`#> 4     1 Open      V1    B1   W2   80`\
`#> 5     2 Open      V1    B1   W2   79`\
`#> 6     3 Open      V1    B1   W2   78`

## Data, design and estimand

The unit nesting is block, then whole plot. Variety is a grouping
column: it labels the whole plots but is not a level of the nesting.

\
`pd`` ``<-`` `[`pollin_data`](https://max578.github.io/pollinatoR/reference/pollin_data.md)`(`\
`  ``trial``,`\
`  response ``=`` ``"pods"``,`\
`  treatment ``=`` ``"trt"``,`\
`  unit ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"block"``, ``"plot"``)``,`\
`  groups ``=`` ``"variety"`\
`)`\
`ds`` ``<-`` `[`pollin_design`](https://max578.github.io/pollinatoR/reference/pollin_design.md)`(`\
`  ``pd``,`\
`  arms ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``excluded ``=`` ``"Cage"``, open ``=`` ``"Open"``)``,`\
`  layout ``=`` ``"split_plot"`\
`)`\
`ds`\
`#> <pollin_design>`\
`#>   layout    : split_plot `\
`#>   arms      : excluded = Cage, open = Open `\
`#>   reference : open `\
`#>   units     : block / plot `\
`#>   random    : (1 | block) + (1 | block:plot)`

The `split_plot` layout gives each block and each whole plot its own
random intercept. The assumptions for this layout name the sub-plot
randomisation.

\
`est`` ``<-`` `[`pollin_estimand`](https://max578.github.io/pollinatoR/reference/pollin_estimand.md)`(``ds``, family ``=`` ``"dependence"``, scale ``=`` ``"ratio"``, basis ``=`` ``"loss"``)`\
`knitr``::`[`kable`](https://rdrr.io/pkg/knitr/man/kable.html)`(`[`pollin_assumptions`](https://max578.github.io/pollinatoR/reference/pollin_assumptions.md)`(``est``)``)`

| assumption | statement |
|:---|:---|
| consistency | Each arm is one well-defined treatment, applied the same way on every unit. |
| excludability | The exclusion device changes the response only by excluding pollinators (no cage effect on microclimate, pests or handling). |
| no_interference | One unit’s treatment does not change another unit’s response, for example by displacing pollinators to neighbouring units. |
| assignment | Arms were randomised to sub-plots within each whole plot. |
| positive_means | The expected response is positive in both arms, so the ratio is defined. |
| population | The estimand refers to the units in the trial; carrying it to other fields or seasons is a further assumption. |

## The trial as a whole

Over all eight varieties, the dependence on the loss basis is the share
of open-plot yield lost under the cage. The same contrast on the gain
basis is the rise from caged to open plots. Both are given here,
labelled, because they are easily confused: the gain is always the
larger of the two.

\
`loss`` ``<-`` `[`pollin_effects`](https://max578.github.io/pollinatoR/reference/pollin_effects.md)`(`[`pollin_fit`](https://max578.github.io/pollinatoR/reference/pollin_fit.md)`(``est``)``)`\
`loss`\
`#> <pollin_effects>`\
`#>   family   : dependence (E2) `\
`#>   scale    : ratio `\
`#>   basis    : loss `\
`#>   estimate : 1 - Cage/Open (loss share) `\
`#>   interval : 95% delta `\
`#>   note     : delta intervals cover fixed-effect uncertainty only`\
`#>  response     contrast  estimate  std.error  conf.low conf.high      p.value`\
`#>      pods Open vs Cage 0.2234674 0.02784709 0.1669243 0.2761728 1.754401e-12`\
`#>      family scale basis distribution interval`\
`#>  dependence ratio  loss      nbinom2    delta`\
`est_gain`` ``<-`` `[`pollin_estimand`](https://max578.github.io/pollinatoR/reference/pollin_estimand.md)`(``ds``, family ``=`` ``"dependence"``, scale ``=`` ``"ratio"``, basis ``=`` ``"gain"``)`\
`gain`` ``<-`` `[`pollin_effects`](https://max578.github.io/pollinatoR/reference/pollin_effects.md)`(`[`pollin_fit`](https://max578.github.io/pollinatoR/reference/pollin_fit.md)`(``est_gain``)``)`\
`gain``$``table``[`[`c`](https://rdrr.io/r/base/c.html)`(``"basis"``, ``"estimate"``, ``"conf.low"``, ``"conf.high"``)``]`\
`#>   basis  estimate  conf.low conf.high`\
`#> 1  gain 0.2877759 0.2003711 0.3815451`

## Dependence by variety

Breeders want the dependence of each variety. Giving `by = "variety"` to
[`pollin_fit()`](https://max578.github.io/pollinatoR/reference/pollin_fit.md)
fits the treatment-by-variety interaction, so each variety has its own
fixed treatment effect, and
[`pollin_effects()`](https://max578.github.io/pollinatoR/reference/pollin_effects.md)
then reports one contrast per variety.

\
`fit_v`` ``<-`` `[`pollin_fit`](https://max578.github.io/pollinatoR/reference/pollin_fit.md)`(``est``, by ``=`` ``"variety"``)`\
`fit_v`\
`#> <pollin_fit>`\
`#>   engine    : glmmTMB `\
`#>   family    : dependence `\
`#>   contrast  : Open vs Cage `\
`#>   model     : <response> ~ trt * variety + (1 | block) + (1 | block:plot) `\
`#>   by        : variety `\
`#>   pods      : nbinom2, 192 rows`\
`by_variety`` ``<-`` `[`pollin_effects`](https://max578.github.io/pollinatoR/reference/pollin_effects.md)`(``fit_v``)``$``table`\
`by_variety``$``planted`` ``<-`` ``dependence``[``by_variety``$``variety``]`\
`by_variety``[`[`c`](https://rdrr.io/r/base/c.html)`(``"variety"``, ``"estimate"``, ``"conf.low"``, ``"conf.high"``, ``"planted"``)``]`\
`#>   variety    estimate    conf.low  conf.high planted`\
`#> 1      V1  0.09209871 -0.08495258 0.24025734    0.10`\
`#> 2      V2  0.19019322  0.03605114 0.31968691    0.30`\
`#> 3      V3  0.24345723  0.09822532 0.36529937    0.20`\
`#> 4      V4  0.06684707 -0.10409827 0.21132529    0.15`\
`#> 5      V5  0.41342001  0.29394563 0.51267764    0.35`\
`#> 6      V6  0.37635114  0.25716216 0.47641615    0.25`\
`#> 7      V7 -0.07374589 -0.27114251 0.09299687    0.05`\
`#> 8      V8  0.41243944  0.29586199 0.50971626    0.30`

Each variety’s estimate rests on its own 24 plants, so the intervals are
wide. Seven of the eight intervals cover the planted value, about what
95% intervals should give. A variety-by-variety reading needs more
blocks than this.

On the difference scale the same fit gives pods per plant lost under the
cage, which mixes a variety’s dependence with its yield.

\
`est_diff`` ``<-`` `[`pollin_estimand`](https://max578.github.io/pollinatoR/reference/pollin_estimand.md)`(`\
`  ``ds``,`\
`  family ``=`` ``"dependence"``,`\
`  scale ``=`` ``"difference"``,`\
`  basis ``=`` ``"loss"`\
`)`\
[`pollin_effects`](https://max578.github.io/pollinatoR/reference/pollin_effects.md)`(`[`pollin_fit`](https://max578.github.io/pollinatoR/reference/pollin_fit.md)`(``est_diff``, by ``=`` ``"variety"``)``)``$``table``[`\
`  `[`c`](https://rdrr.io/r/base/c.html)`(``"variety"``, ``"estimate"``, ``"conf.low"``, ``"conf.high"``)`\
`]`\
`#>   variety  estimate   conf.low conf.high`\
`#> 1      V1  6.304846  -5.366142 17.975833`\
`#> 2      V2 16.148730   2.550347 29.747114`\
`#> 3      V3 20.279707   7.091321 33.468093`\
`#> 4      V4  6.859823  -9.857398 23.577044`\
`#> 5      V5 28.659101  17.663295 39.654907`\
`#> 6      V6 36.238474  21.561680 50.915267`\
`#> 7      V7 -6.848110 -23.129306  9.433086`\
`#> 8      V8 32.717201  20.366582 45.067819`

## Guard: dependence against yield across varieties

A tempting next step is to regress each variety’s dependence index on
its open yield, asking whether high-yielding varieties need insects
more. When the index is the RII and the yield axis is the open-plot
yield, the open value sits on both axes and part of any slope is built
in.
[`pollin_guard_spurious()`](https://max578.github.io/pollinatoR/reference/pollin_guard_spurious.md)
permutes the caged values across varieties, holding the open values
fixed, to show the slope that the shared term alone produces.

\
`open_mean`` ``<-`` `[`tapply`](https://rdrr.io/r/base/tapply.html)`(``trial``$``pods``[``!``cage``]``, ``trial``$``variety``[``!``cage``]``, ``mean``)`\
`cage_mean`` ``<-`` `[`tapply`](https://rdrr.io/r/base/tapply.html)`(``trial``$``pods``[``cage``]``, ``trial``$``variety``[``cage``]``, ``mean``)`\
[`pollin_guard_spurious`](https://max578.github.io/pollinatoR/reference/pollin_guard_spurious.md)`(``open_mean``, ``cage_mean``, n_perm ``=`` ``500``, seed ``=`` ``1``)`\
`#> <pollin_spurious>  index ~ open-yield slope guard`\
`#>   observed slope : -0.0031`\
`#>   coupling null  : mean 0.0055, sd 0.0047 (n_perm = 500)`\
`#>   permutation p  : 0.687  -> slope consistent with mechanical coupling (SPURIOUS risk)`

Eight varieties make a small panel, and the permutation null is wide.
The planted dependences were chosen independently of yield, and the
guard does not separate the observed slope from what the shared term
would give.

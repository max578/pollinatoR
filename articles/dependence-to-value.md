# From dependence to value

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`pollinatoR`](https://github.com/max578/pollinatoR)`)`

## From a ratio to a dollar figure

A dependence ratio says what share of a crop’s yield is owed to
pollinators. Multiplied by the crop’s production and price, it gives the
value of that pollination: the dependence-ratio method of Gallai et
al. (2009). This vignette takes a fitted estimand through to such a
value, with an interval. The value is computed for each draw of the
estimand, and the interval comes from the value draws.

Each printed value lists the estimand and its basis, the population, the
price year, the currency and the source of the market inputs.

## A simulated trial

The trial has three arms in twelve blocks: plants *bagged* against
insects, plants left *open*, and plants pollinated by *hand*. The blocks
lie in canola fields at three distances from native vegetation, four
blocks per distance band. Open plants yield less the further the field
is from vegetation; bagged and hand-pollinated plants do not depend on
distance. The response is seed yield per plant in grams.

\
[`set.seed`](https://rdrr.io/r/base/Random.html)`(``2025``)`\
`n_block`` ``<-`` ``12L`\
`n_plant`` ``<-`` ``6L`\
`blocks`` ``<-`` `[`sprintf`](https://rdrr.io/r/base/sprintf.html)`(``"B%02d"``, `[`seq_len`](https://rdrr.io/r/base/seq.html)`(``n_block``)``)`\
`bands`` ``<-`` `[`c`](https://rdrr.io/r/base/c.html)`(``"0-100 m"``, ``"100-300 m"``, ``"over 300 m"``)`\
`band`` ``<-`` `[`setNames`](https://rdrr.io/r/stats/setNames.html)`(`[`rep`](https://rdrr.io/r/base/rep.html)`(``bands``, each ``=`` ``4L``)``, ``blocks``)`\
`trial`` ``<-`` `[`expand.grid`](https://rdrr.io/r/base/expand.grid.html)`(`\
`  plant ``=`` `[`seq_len`](https://rdrr.io/r/base/seq.html)`(``n_plant``)``,`\
`  trt ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"Bagged"``, ``"Open"``, ``"Hand"``)``,`\
`  block ``=`` ``blocks``,`\
`  stringsAsFactors ``=`` ``FALSE`\
`)`\
`trial``$``band`` ``<-`` ``band``[``trial``$``block``]`\
`block_effect`` ``<-`` `[`setNames`](https://rdrr.io/r/stats/setNames.html)`(`[`rnorm`](https://rdrr.io/r/stats/Normal.html)`(``n_block``, ``0``, ``0.1``)``, ``blocks``)`\
`open_mean`` ``<-`` `[`setNames`](https://rdrr.io/r/stats/setNames.html)`(`[`c`](https://rdrr.io/r/base/c.html)`(``2.4``, ``2.1``, ``1.8``)``, ``bands``)`\
`arm_mean`` ``<-`` `[`ifelse`](https://rdrr.io/r/base/ifelse.html)`(`\
`  ``trial``$``trt`` ``==`` ``"Bagged"``,`\
`  ``1.4``,`\
`  `[`ifelse`](https://rdrr.io/r/base/ifelse.html)`(``trial``$``trt`` ``==`` ``"Hand"``, ``2.6``, ``open_mean``[``trial``$``band``]``)`\
`)`\
`mu`` ``<-`` ``arm_mean`` ``*`` `[`exp`](https://rdrr.io/r/base/Log.html)`(``block_effect``[``trial``$``block``]``)`\
`trial``$``yield`` ``<-`` `[`rgamma`](https://rdrr.io/r/stats/GammaDist.html)`(`[`nrow`](https://rdrr.io/r/base/nrow.html)`(``trial``)``, shape ``=`` ``25``, rate ``=`` ``25`` ``/`` ``mu``)`\
\
`pd`` ``<-`` `[`pollin_data`](https://max578.github.io/pollinatoR/reference/pollin_data.md)`(``trial``, response ``=`` ``"yield"``, treatment ``=`` ``"trt"``, unit ``=`` ``"block"``)`\
`ds`` ``<-`` `[`pollin_design`](https://max578.github.io/pollinatoR/reference/pollin_design.md)`(`\
`  ``pd``,`\
`  arms ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``excluded ``=`` ``"Bagged"``, open ``=`` ``"Open"``, supplemented ``=`` ``"Hand"``)``,`\
`  layout ``=`` ``"rcbd"`\
`)`

## A population of distance bands

The trial’s blocks are not a sample of the district’s canola. Suppose a
land use map shows 300 ha of canola within 100 m of native vegetation,
500 ha at 100 to 300 m and 1,200 ha beyond.
[`pollin_population()`](https://max578.github.io/pollinatoR/reference/pollin_population.md)
declares the bands as strata weighted by those areas.

\
`district`` ``<-`` `[`pollin_population`](https://max578.github.io/pollinatoR/reference/pollin_population.md)`(`\
`  ``"band"``,`\
`  `[`c`](https://rdrr.io/r/base/c.html)`(``"0-100 m"`` ``=`` ``300``, ``"100-300 m"`` ``=`` ``500``, ``"over 300 m"`` ``=`` ``1200``)``,`\
`  label ``=`` ``"district canola, by distance to vegetation"`\
`)`\
`district`\
`#> <pollin_population>`\
`#>   label   : district canola, by distance to vegetation `\
`#>   strata  : 3 levels of band `\
`#>   weights : 0-100 m 0.15, 100-300 m 0.25, over 300 m 0.60 (area, normalised)`

The estimand takes the population, and fitting with the bands as `by`
gives each band its own treatment effect.

\
`dep`` ``<-`` `[`pollin_estimand`](https://max578.github.io/pollinatoR/reference/pollin_estimand.md)`(``ds``, ``"dependence"``, ``"ratio"``, ``"loss"``, population ``=`` ``district``)`\
`fit_dep`` ``<-`` `[`pollin_fit`](https://max578.github.io/pollinatoR/reference/pollin_fit.md)`(``dep``, by ``=`` ``"band"``)`

## Draws of the estimand

[`pollin_draws()`](https://max578.github.io/pollinatoR/reference/pollin_draws.md)
draws the fixed effects from their sampling distribution and recomputes
the expected yield in each arm for every draw. The loss share is formed
from each pair of arm means. A further column holds the population: the
area-weighted arm means, and the loss share formed from them.

\
`dr_dep`` ``<-`` `[`pollin_draws`](https://max578.github.io/pollinatoR/reference/pollin_draws.md)`(``fit_dep``, n ``=`` ``1000L``, seed ``=`` ``1L``)`\
`dr_dep`\
`#> <pollin_draws>`\
`#>   family   : dependence (E2) `\
`#>   scale    : ratio `\
`#>   basis    : loss `\
`#>   estimate : 1 - Bagged/Open (loss share) `\
`#>   response : yield `\
`#>   method   : simulation of the fixed and random effects `\
`#>   engine   : glmmTMB `\
`#>   draws    : 1000 `\
`#>   columns  : 0-100 m, 100-300 m, over 300 m, population`\
[`summary`](https://rdrr.io/r/base/summary.html)`(``dr_dep``)``[``, `[`c`](https://rdrr.io/r/base/c.html)`(``"level"``, ``"median"``, ``"conf.low"``, ``"conf.high"``)``]`\
`#>        level    median   conf.low conf.high`\
`#> 1    0-100 m 0.4066923 0.33474198 0.4681777`\
`#> 2  100-300 m 0.3338761 0.24476716 0.4011106`\
`#> 3 over 300 m 0.1691020 0.07298917 0.2519253`\
`#> 4 population 0.2614497 0.20835286 0.3115552`

The planted loss shares are 0.42, 0.33 and 0.22 in the three bands. The
population loss share is a ratio of population means,
$`1 - \sum w\,\mathrm{low} / \sum w\,\mathrm{high}`$, and so weights
each band by its area and its open yield.

## Market inputs

[`pollin_market()`](https://max578.github.io/pollinatoR/reference/pollin_market.md)
records current production and its price. Current production means the
crop as it is pollinated now, the open arm in the trial’s terms. The
figures below are illustrative.

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
`canola`\
`#> <pollin_market>`\
`#>   crop     : canola `\
`#>   quantity : 3,600 t `\
`#>   price    : 700 AUD per t `\
`#>   year     : 2025 `\
`#>   source   : Illustrative inputs`

## The dependence channel

[`pollin_value()`](https://max578.github.io/pollinatoR/reference/pollin_value.md)
needs a channel; there is no default. The `dependence` channel values
the production that would be lost if pollinators were wholly excluded.

\
`val_dep`` ``<-`` `[`pollin_value`](https://max578.github.io/pollinatoR/reference/pollin_value.md)`(``dr_dep``, ``canola``, channel ``=`` ``"dependence"``)`\
`val_dep`\
`#> <pollin_value>`\
`#>   channel    : dependence - value lost under total exclusion of pollinators `\
`#>   value      : 658,900 AUD (95% interval 525,000 to 785,100) `\
`#>   basis      : 1 - Bagged/Open (loss share)`\
`#>   level      : population `\
`#>   population : district canola, by distance to vegetation: 3 strata of band, weighted by area `\
`#>   price year : 2025 `\
`#>   currency   : AUD `\
`#>   market     : 3,600 t at 700 AUD per t `\
`#>   source     : Illustrative inputs `\
`#>   attribution: agent all, product yield, crop canola `\
`#>   method     : simulation of the fixed and random effects (1000 draws) `\
`#>   assumptions:`\
`#>     - No price response: the price is held fixed whatever the change in production.`\
`#>     - Total-loss counterfactual: the value is the production lost if pollinators were wholly excluded, with no substitution by other pollinators or inputs.`\
`#>     - The estimand draws are independent of the market draws.`\
`` #>     - The population is the strata of `band` (0-100 m 0.15, 100-300 m 0.25, over 300 m 0.60), weighted by area. Within each stratum the expected response under each arm comes from the fitted model and is taken to hold across the stratum's area; the population estimand is formed from the area-weighted arm means. ``

The value is computed per draw as $`D \times P \times Q`$, with
$`D = 1 - \mathrm{excluded}/\mathrm{open}`$. Production without
pollinators would be $`Q \times \mathrm{excluded}/\mathrm{open}`$, a
loss of $`D \times Q`$.

Price and production are rarely known exactly. They can be given as
draws, one per estimand draw, matched by position.

\
[`set.seed`](https://rdrr.io/r/base/Random.html)`(``7``)`\
`canola_unc`` ``<-`` `[`pollin_market`](https://max578.github.io/pollinatoR/reference/pollin_market.md)`(`\
`  quantity ``=`` ``3600``,`\
`  price ``=`` ``700``,`\
`  unit ``=`` ``"t"``,`\
`  currency ``=`` ``"AUD"``,`\
`  year ``=`` ``2025``,`\
`  source ``=`` ``"Illustrative inputs with uncertainty"``,`\
`  quantity_draws ``=`` `[`rlnorm`](https://rdrr.io/r/stats/Lognormal.html)`(``1000L``, `[`log`](https://rdrr.io/r/base/Log.html)`(``3600``)``, ``0.10``)``,`\
`  price_draws ``=`` `[`rlnorm`](https://rdrr.io/r/stats/Lognormal.html)`(``1000L``, `[`log`](https://rdrr.io/r/base/Log.html)`(``700``)``, ``0.15``)``,`\
`  crop ``=`` ``"canola"`\
`)`\
[`summary`](https://rdrr.io/r/base/summary.html)`(`[`pollin_value`](https://max578.github.io/pollinatoR/reference/pollin_value.md)`(``dr_dep``, ``canola_unc``, channel ``=`` ``"dependence"``)``)``[`\
`  , `[`c`](https://rdrr.io/r/base/c.html)`(``"median"``, ``"conf.low"``, ``"conf.high"``, ``"currency"``)`\
`]`\
`#>     median conf.low conf.high currency`\
`#> 1 657314.2 433421.1  992872.4      AUD`

The median barely moves. The interval is wider because price and
production are now uncertain too.

## The deficit channel

The `deficit` channel values the yield lost because the crop’s
pollination falls short of hand pollination: the gain if every plant
were pollinated as well as by hand. It needs a supplemented arm, and a
deficit estimand with `basis = "gain"`.

\
`def`` ``<-`` `[`pollin_estimand`](https://max578.github.io/pollinatoR/reference/pollin_estimand.md)`(``ds``, ``"deficit"``, ``"ratio"``, ``"gain"``, population ``=`` ``district``)`\
`dr_def`` ``<-`` `[`pollin_draws`](https://max578.github.io/pollinatoR/reference/pollin_draws.md)`(`[`pollin_fit`](https://max578.github.io/pollinatoR/reference/pollin_fit.md)`(``def``, by ``=`` ``"band"``)``, n ``=`` ``1000L``, seed ``=`` ``1L``)`\
[`pollin_value`](https://max578.github.io/pollinatoR/reference/pollin_value.md)`(``dr_def``, ``canola``, channel ``=`` ``"deficit"``)`\
`#> <pollin_value>`\
`#>   channel    : deficit - value of closing the gap to hand pollination `\
`#>   value      : 837,600 AUD (95% interval 607,800 to 1,074,000) `\
`#>   basis      : Hand/Open - 1 (gain)`\
`#>   level      : population `\
`#>   population : district canola, by distance to vegetation: 3 strata of band, weighted by area `\
`#>   price year : 2025 `\
`#>   currency   : AUD `\
`#>   market     : 3,600 t at 700 AUD per t `\
`#>   source     : Illustrative inputs `\
`#>   attribution: agent all, product yield, crop canola `\
`#>   method     : simulation of the fixed and random effects (1000 draws) `\
`#>   assumptions:`\
`#>     - No price response: the price is held fixed whatever the change in production.`\
`#>     - Gap counterfactual: the value is the production gained if the whole crop received hand-supplemented pollination, at no change in costs.`\
`#>     - The estimand draws are independent of the market draws.`\
`` #>     - The population is the strata of `band` (0-100 m 0.15, 100-300 m 0.25, over 300 m 0.60), weighted by area. Within each stratum the expected response under each arm comes from the fitted model and is taken to hold across the stratum's area; the population estimand is formed from the area-weighted arm means. ``

Here the value per draw is
$`(\mathrm{hand}/\mathrm{open} - 1) \times P \times Q`$, the gain over
current production. In both channels the ratio is taken relative to the
open arm, because $`Q`$ is open-pollinated production. That is why each
channel accepts one ratio only: a loss share for `dependence` and a gain
for `deficit`.

## Why the dependence channel refuses a gain

The loss share $`D`$ and the gain
$`g = \mathrm{open}/\mathrm{excluded} - 1`$ come from the same two arm
means, and $`g = D/(1 - D)`$. A gain multiplied by current production
overstates the loss by the factor $`1/(1 - D)`$, so the dependence
channel refuses it.

\
`dep_gain`` ``<-`` `[`pollin_estimand`](https://max578.github.io/pollinatoR/reference/pollin_estimand.md)`(``ds``, ``"dependence"``, ``"ratio"``, ``"gain"``, population ``=`` ``district``)`\
`dr_gain`` ``<-`` `[`pollin_draws`](https://max578.github.io/pollinatoR/reference/pollin_draws.md)`(`[`pollin_fit`](https://max578.github.io/pollinatoR/reference/pollin_fit.md)`(``dep_gain``, by ``=`` ``"band"``)``, n ``=`` ``1000L``, seed ``=`` ``1L``)`\
[`try`](https://rdrr.io/r/base/try.html)`(`[`pollin_value`](https://max578.github.io/pollinatoR/reference/pollin_value.md)`(``dr_gain``, ``canola``, channel ``=`` ``"dependence"``)``)`\
`` #> Error : The `dependence` channel needs the share of current production lost without pollinators, D = 1 - low/high. A gain g = high/low - 1 relates to it as g = D / (1 - D), so valuing a gain overstates the loss by the factor 1/(1 - D). Refit the estimand with basis = "loss": `pollin_estimand(design, "dependence", "ratio", "loss")`. ``

A descriptive index from
[`pollin_index()`](https://max578.github.io/pollinatoR/reference/pollin_index.md)
or a plain number is refused as well, because neither records the
estimand it came from.

## Sensitivity: two populations

The area weights are an input, and a different map gives a different
value. To see how much, value the estimand under a second population and
compare. Here a second population puts most of the canola close to
vegetation.

\
`near`` ``<-`` `[`pollin_population`](https://max578.github.io/pollinatoR/reference/pollin_population.md)`(`\
`  ``"band"``,`\
`  `[`c`](https://rdrr.io/r/base/c.html)`(``"0-100 m"`` ``=`` ``1000``, ``"100-300 m"`` ``=`` ``600``, ``"over 300 m"`` ``=`` ``400``)``,`\
`  label ``=`` ``"canola mostly near vegetation"`\
`)`\
`dep_near`` ``<-`` `[`pollin_estimand`](https://max578.github.io/pollinatoR/reference/pollin_estimand.md)`(``ds``, ``"dependence"``, ``"ratio"``, ``"loss"``, population ``=`` ``near``)`\
`dr_near`` ``<-`` `[`pollin_draws`](https://max578.github.io/pollinatoR/reference/pollin_draws.md)`(`[`pollin_fit`](https://max578.github.io/pollinatoR/reference/pollin_fit.md)`(``dep_near``, by ``=`` ``"band"``)``, n ``=`` ``1000L``, seed ``=`` ``1L``)`\
[`rbind`](https://rdrr.io/r/base/cbind.html)`(`\
`  `[`summary`](https://rdrr.io/r/base/summary.html)`(`[`pollin_value`](https://max578.github.io/pollinatoR/reference/pollin_value.md)`(``dr_dep``, ``canola``, channel ``=`` ``"dependence"``)``)``,`\
`  `[`summary`](https://rdrr.io/r/base/summary.html)`(`[`pollin_value`](https://max578.github.io/pollinatoR/reference/pollin_value.md)`(``dr_near``, ``canola``, channel ``=`` ``"dependence"``)``)`\
`)``[``, `[`c`](https://rdrr.io/r/base/c.html)`(``"median"``, ``"conf.low"``, ``"conf.high"``, ``"currency"``)``]`\
`#>     median conf.low conf.high currency`\
`#> 1 658853.1 525049.2  785119.1      AUD`\
`#> 2 882757.4 762560.2 1004226.4      AUD`

The first row is the district map, the second the map with most canola
near vegetation. With production held at 3,600 t, the second value is
higher, because the crop near vegetation loses a larger share without
pollinators.

## Adding values without counting twice

Values can be added only when they count different things.
[`pollin_guard_doublecount()`](https://max578.github.io/pollinatoR/reference/pollin_guard_doublecount.md)
adds values draw by draw and stops when two of them overlap. It compares
the pollinators each exclusion removed (`agent`, all of them here), the
product valued and the crop.

The value of canola seed and the value of the oil pressed from it are
not separate contributions: the oil’s value includes the seed’s. The
guard refuses the sum.

\
`oil`` ``<-`` `[`pollin_market`](https://max578.github.io/pollinatoR/reference/pollin_market.md)`(`\
`  quantity ``=`` ``1500``,`\
`  price ``=`` ``1600``,`\
`  unit ``=`` ``"t"``,`\
`  currency ``=`` ``"AUD"``,`\
`  year ``=`` ``2025``,`\
`  source ``=`` ``"Illustrative inputs"``,`\
`  crop ``=`` ``"canola"`\
`)`\
[`try`](https://rdrr.io/r/base/try.html)`(`[`pollin_guard_doublecount`](https://max578.github.io/pollinatoR/reference/pollin_guard_doublecount.md)`(`\
`  seed ``=`` `[`pollin_value`](https://max578.github.io/pollinatoR/reference/pollin_value.md)`(``dr_dep``, ``canola``, channel ``=`` ``"dependence"``, product ``=`` ``"seed"``)``,`\
`  oil ``=`` `[`pollin_value`](https://max578.github.io/pollinatoR/reference/pollin_value.md)`(``dr_dep``, ``oil``, channel ``=`` ``"dependence"``, product ``=`` ``"oil"``)`\
`)``)`\
`` #> Error : The values `seed` and `oil` are for two products ('seed' and 'oil') attributed to 'all' for the same crop (canola). The value of the processed product already contains the value of the raw one, so adding them counts it again; value one product. ``

Values for two crops, each from its own trial, do not overlap. A second
small trial in faba bean gives a second value, and the guard adds the
two.

\
[`set.seed`](https://rdrr.io/r/base/Random.html)`(``11``)`\
`faba_trial`` ``<-`` `[`expand.grid`](https://rdrr.io/r/base/expand.grid.html)`(`\
`  plant ``=`` `[`seq_len`](https://rdrr.io/r/base/seq.html)`(``6L``)``,`\
`  trt ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"Bagged"``, ``"Open"``)``,`\
`  block ``=`` `[`sprintf`](https://rdrr.io/r/base/sprintf.html)`(``"F%d"``, `[`seq_len`](https://rdrr.io/r/base/seq.html)`(``8L``)``)``,`\
`  stringsAsFactors ``=`` ``FALSE`\
`)`\
`faba_mu`` ``<-`` `[`ifelse`](https://rdrr.io/r/base/ifelse.html)`(``faba_trial``$``trt`` ``==`` ``"Open"``, ``30``, ``24``)`\
`faba_trial``$``yield`` ``<-`` `[`rgamma`](https://rdrr.io/r/stats/GammaDist.html)`(`[`nrow`](https://rdrr.io/r/base/nrow.html)`(``faba_trial``)``, shape ``=`` ``20``, rate ``=`` ``20`` ``/`` ``faba_mu``)`\
`faba_ds`` ``<-`` `[`pollin_design`](https://max578.github.io/pollinatoR/reference/pollin_design.md)`(`\
`  `[`pollin_data`](https://max578.github.io/pollinatoR/reference/pollin_data.md)`(``faba_trial``, ``"yield"``, ``"trt"``, unit ``=`` ``"block"``)``,`\
`  `[`c`](https://rdrr.io/r/base/c.html)`(``excluded ``=`` ``"Bagged"``, open ``=`` ``"Open"``)``,`\
`  ``"rcbd"`\
`)`\
`faba_dr`` ``<-`` `[`pollin_draws`](https://max578.github.io/pollinatoR/reference/pollin_draws.md)`(`\
`  `[`pollin_fit`](https://max578.github.io/pollinatoR/reference/pollin_fit.md)`(`[`pollin_estimand`](https://max578.github.io/pollinatoR/reference/pollin_estimand.md)`(``faba_ds``, ``"dependence"``, ``"ratio"``, ``"loss"``)``)``,`\
`  n ``=`` ``1000L``,`\
`  seed ``=`` ``2L`\
`)`\
`faba`` ``<-`` `[`pollin_market`](https://max578.github.io/pollinatoR/reference/pollin_market.md)`(`\
`  quantity ``=`` ``900``,`\
`  price ``=`` ``450``,`\
`  unit ``=`` ``"t"``,`\
`  currency ``=`` ``"AUD"``,`\
`  year ``=`` ``2025``,`\
`  source ``=`` ``"Illustrative inputs"``,`\
`  crop ``=`` ``"faba bean"`\
`)`\
[`pollin_guard_doublecount`](https://max578.github.io/pollinatoR/reference/pollin_guard_doublecount.md)`(`\
`  canola ``=`` `[`pollin_value`](https://max578.github.io/pollinatoR/reference/pollin_value.md)`(``dr_dep``, ``canola``, channel ``=`` ``"dependence"``, product ``=`` ``"seed"``)``,`\
`  faba ``=`` `[`pollin_value`](https://max578.github.io/pollinatoR/reference/pollin_value.md)`(``faba_dr``, ``faba``, channel ``=`` ``"dependence"``, product ``=`` ``"seed"``)`\
`)`\
`#> <pollin_value_total>  sum of 2 values `\
`#>   total      : 730,400 AUD (95% interval 593,900 to 856,900) `\
`#>   price year : 2025 `\
`#>   draws      : 1000 (added draw by draw, paired by index)`\
`#>   components :`\
`#>     - canola: 658,900 AUD; dependence; agent all, product seed, crop canola; source: Illustrative inputs`\
`#>     - faba: 72,410 AUD; dependence; agent all, product seed, crop faba bean; source: Illustrative inputs`

## What the value assumes

Each printed value lists its assumptions:

- No price response. The price is held fixed, although losing a large
  share of a crop’s production would raise its price.
- A total-loss counterfactual, for the dependence channel. The value is
  the production lost if every pollinator were removed, with nothing in
  their place. It is not the loss from a partial decline.
- Independence of the estimand and the market draws. A season that is
  poor for pollination may also be a season of high prices; the draws
  are paired as if the two were unrelated.

Within each stratum of a population, the trial’s expected yields are
taken to hold across the stratum’s area, and the areas are taken as
exact.

The value is gross: the revenue from the production at stake, before any
costs. The default simulation draws take the fixed effects and the
trial’s own block effects jointly, so the interval covers the
uncertainty in both. It does not cover new blocks drawn afresh;
[`pollin_draws()`](https://max578.github.io/pollinatoR/reference/pollin_draws.md)
with `method = "bootstrap"` includes those, at the cost of refitting the
model.

## References

Gallai, N., Salles, J.-M., Settele, J., & Vaissière, B. E. (2009).
Economic valuation of the vulnerability of world agriculture confronted
with pollinator decline. *Ecological Economics*, *68*(3), 810–821.
<https://doi.org/10.1016/j.ecolecon.2008.06.014>

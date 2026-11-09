# Value a dependence or deficit estimand

Multiplies draws of a ratio estimand by price and current production,
draw by draw, and summarises the value draws last: the median and a
percentile interval. The ratio is formed from the arm-mean draws of
[`pollin_draws()`](https://max578.github.io/pollinatoR/reference/pollin_draws.md)
for each draw, so the uncertainty of the estimand (and of the market
inputs, when given as draws) carries through to the value.

## Usage

``` r
pollin_value(
  draws,
  market,
  channel,
  level = NULL,
  conf_level = 0.95,
  agent = "all",
  product = NULL
)
```

## Arguments

- draws:

  A
  [`pollin_draws()`](https://max578.github.io/pollinatoR/reference/pollin_draws.md)
  object.

- market:

  A
  [`pollin_market()`](https://max578.github.io/pollinatoR/reference/pollin_market.md)
  object.

- channel:

  `"dependence"` or `"deficit"`. There is no default.

- level:

  The column of the draws to value: a level of the fit's `by`,
  `"overall"` or `"population"`. The default is `"population"` when the
  draws hold it, else the only column; with several columns and no
  population, `level` must be given.

- conf_level:

  Level of the percentile interval.

- agent:

  The pollinators the estimand's exclusion removes, one string: `"all"`
  (the default) for a contrast that excludes every pollinator, or the
  label of a guild, for example `"managed honey bees"` or
  `"wild insects"`.

- product:

  The product valued, one string (for example `"seed"`, `"oil"` or
  `"fruit"`). The default is the name of the response.

## Value

An object of class `pollin_value`: a list holding

- `draws`: the numeric vector of value draws, in `currency`;

- `summary`: a one-row data frame with `channel`, `level`, `median`,
  `conf.low`, `conf.high`, `currency`, `year`, `unit`, `basis`,
  `population`, `method` and `n`;

- `channel`, `estimand`, `level`, `market`;

- `assumptions`: the statements the value rests on;

- `conf_level`: the level of the interval in `summary`;

- `attribution`: a list of `agent`, `product` and `crop` (`NULL` when
  the market names none).

## Details

The quantity in
[`pollin_market()`](https://max578.github.io/pollinatoR/reference/pollin_market.md)
is current production, the crop as it is pollinated now, so each channel
needs the ratio taken relative to it:

|  |  |  |
|----|----|----|
| Channel | Estimand | Value per draw |
| `"dependence"` | `family = "dependence"`, `scale = "ratio"`, `basis = "loss"` | `D x P x Q`, with `D = 1 - low/high`: the value lost if pollinators were wholly excluded |
| `"deficit"` | `family = "deficit"`, `scale = "ratio"`, `basis = "gain"` | `(E[Y_supp]/E[Y_open] - 1) x P x Q`: the value of closing the gap to hand pollination on current production |

Any other estimand is refused, with the call that declares the right
one. In particular a gain `g = high/low - 1` is refused by the
dependence channel: it relates to the loss share by `g = D / (1 - D)`,
so valuing a gain overstates the loss by the factor `1 / (1 - D)`.

A dependence estimand whose reference arm is `"supplemented"` (potential
dependence) is valued and labelled as the value at full pollination,
which is not the contribution of the pollination the crop receives now.

Draws of quantity or price in the market are paired with the estimand
draws by index; their number must equal the number of estimand draws.
Only a
[`pollin_draws()`](https://max578.github.io/pollinatoR/reference/pollin_draws.md)
object can be valued. A descriptive index from
[`pollin_index()`](https://max578.github.io/pollinatoR/reference/pollin_index.md)
or
[`pollin_rii()`](https://max578.github.io/pollinatoR/reference/pollin_rii.md),
or a plain number, carries no estimand.

Each value carries an attribution: the pollinators the exclusion removed
(`agent`), the product valued (`product`) and the crop (from
[`pollin_market()`](https://max578.github.io/pollinatoR/reference/pollin_market.md)).
[`pollin_guard_doublecount()`](https://max578.github.io/pollinatoR/reference/pollin_guard_doublecount.md)
reads it before values are added. For the deficit channel, `agent` names
the pollinators whose shortfall hand pollination closes.

## See also

[`pollin_draws()`](https://max578.github.io/pollinatoR/reference/pollin_draws.md),
[`pollin_market()`](https://max578.github.io/pollinatoR/reference/pollin_market.md).

## Examples

``` r
# \donttest{
set.seed(1)
d <- data.frame(
  trt = rep(rep(c("Cage", "Open"), each = 10L), 8L),
  Site = rep(paste0("s", seq_len(8L)), each = 20L)
)
d$seed <- rgamma(nrow(d), 20, 20 / exp(1 + 0.3 * (d$trt == "Open")))
pd <- pollin_data(d, "seed", "trt", unit = "Site")
ds <- pollin_design(pd, c(excluded = "Cage", open = "Open"), "rcbd")
est <- pollin_estimand(ds, "dependence", "ratio", "loss")
dr <- pollin_draws(pollin_fit(est), n = 2000L, seed = 1L)
mk <- pollin_market(12000, 650, "t", "AUD", 2024, "Illustrative inputs")
pollin_value(dr, mk, channel = "dependence")
#> <pollin_value>
#>   channel    : dependence - value lost under total exclusion of pollinators 
#>   value      : 2,003,000 AUD (95% interval 1,617,000 to 2,369,000) 
#>   basis      : 1 - Cage/Open (loss share)
#>   level      : overall 
#>   population : trial 
#>   price year : 2024 
#>   currency   : AUD 
#>   market     : 12,000 t at 650 AUD per t 
#>   source     : Illustrative inputs 
#>   attribution: agent all, product seed, crop not declared 
#>   method     : simulation of the fixed and random effects (2000 draws) 
#>   assumptions:
#>     - No price response: the price is held fixed whatever the change in production.
#>     - Total-loss counterfactual: the value is the production lost if pollinators were wholly excluded, with no substitution by other pollinators or inputs.
#>     - The estimand draws are independent of the market draws.
#>     - The estimand refers to the units in the trial; carrying it to other fields or seasons is a further assumption.
# }
```

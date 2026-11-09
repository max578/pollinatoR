# Guard against adding values that count one contribution twice

Adds two or more
[`pollin_value()`](https://max578.github.io/pollinatoR/reference/pollin_value.md)
results, draw by draw, after checking that their attributions do not
overlap. Each value carries the pollinators its exclusion removed
(`agent`), the product valued and the crop. For values of the same crop
the guard refuses

## Usage

``` r
pollin_guard_doublecount(..., nested = NULL, conf_level = 0.95)

# S3 method for class 'pollin_value_total'
summary(object, ...)
```

## Arguments

- ...:

  Two or more
  [`pollin_value()`](https://max578.github.io/pollinatoR/reference/pollin_value.md)
  objects, optionally named.

- nested:

  `NULL` (the default) or `TRUE` refuses two different named guilds of
  one crop; `FALSE` asserts that their exclusion contrasts are separate
  and not nested, and allows the sum.

- conf_level:

  Level of the percentile interval of the total.

- object:

  A `pollin_value_total` object.

## Value

An object of class `pollin_value_total`: a list holding `draws` (the
per-draw total), `summary` (a one-row data frame with `median`,
`conf.low`, `conf.high`, `currency`, `year`, `n` and `components`),
`components` (the values added), `assertion` (the statement made by
`nested = FALSE`, or `NULL`) and `conf_level`.

## Details

- a part plus its total: a value for `agent = "all"` with a value for
  any other agent;

- the same agent twice for the same product;

- two products of the same agent, for example seed plus oil: the value
  of the output already contains the value of the seed;

- two different named guilds, for example managed plus wild pollinators,
  unless `nested = FALSE` asserts that they come from separate,
  non-nested exclusion contrasts. The assertion is recorded and printed
  with the total.

A value whose market names no crop is taken to be of the same crop as
any other, since the guard cannot rule out the overlap; declare `crop`
in
[`pollin_market()`](https://max578.github.io/pollinatoR/reference/pollin_market.md).
Agents and products are compared ignoring case and surrounding spaces.

The values must share currency and price year and have the same number
of draws. They are added draw by draw, paired by index: values drawn
from separate fits are summed as independent, and values read from one
[`pollin_draws()`](https://max578.github.io/pollinatoR/reference/pollin_draws.md)
object keep their correlation.

## See also

[`pollin_value()`](https://max578.github.io/pollinatoR/reference/pollin_value.md).

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
canola <- pollin_market(
  12000, 650, "t", "AUD", 2024, "Illustrative inputs", crop = "canola"
)
faba <- pollin_market(
  3000, 400, "t", "AUD", 2024, "Illustrative inputs", crop = "faba bean"
)
# Two crops: no overlap
pollin_guard_doublecount(
  canola = pollin_value(dr, canola, channel = "dependence"),
  faba = pollin_value(dr, faba, channel = "dependence")
)
#> <pollin_value_total>  sum of 2 values 
#>   total      : 2,312,000 AUD (95% interval 1,866,000 to 2,733,000) 
#>   price year : 2024 
#>   draws      : 2000 (added draw by draw, paired by index)
#>   components :
#>     - canola: 2,003,000 AUD; dependence; agent all, product seed, crop canola; source: Illustrative inputs
#>     - faba: 308,200 AUD; dependence; agent all, product seed, crop faba bean; source: Illustrative inputs
# Seed plus oil of one crop is refused
try(pollin_guard_doublecount(
  pollin_value(dr, canola, channel = "dependence", product = "seed"),
  pollin_value(dr, canola, channel = "dependence", product = "oil")
))
#> Error : The values `value 1` and `value 2` are for two products ('seed' and 'oil') attributed to 'all' for the same crop (canola). The value of the processed product already contains the value of the raw one, so adding them counts it again; value one product.
# }
```

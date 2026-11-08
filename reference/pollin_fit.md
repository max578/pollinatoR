# Fit the outcome models for an estimand

Fits one generalised linear mixed model per response in the estimand's
data, with the pollination treatment as the fixed effect and the random
terms of the design: `response ~ treatment + <design random terms>`.
With `by`, the fixed part is `treatment * by`, so the treatment effect
has its own fixed estimate at each level of `by`. The treatment is
releveled so the estimand's `low` arm is the reference level. Every
treatment level in the data is kept in the fit; the contrast is read
from the fitted model by
[`pollin_effects()`](https://max578.github.io/pollinatoR/reference/pollin_effects.md).

## Usage

``` r
pollin_fit(estimand, engine = "glmmTMB", family = NULL, by = NULL, ...)
```

## Arguments

- estimand:

  A
  [`pollin_estimand()`](https://max578.github.io/pollinatoR/reference/pollin_estimand.md)
  object.

- engine:

  The fitting engine. Only `"glmmTMB"` is available; `"flexyBayes"`
  arrives in version 0.2.0.

- family:

  `NULL` (the default) chooses a family per response from its values
  (count, proportion in (0, 1), positive continuous, or continuous with
  zeros). A single `glmmTMB` family object applies to every response; a
  named list gives one per response, and responses it does not name are
  chosen from their values.

- by:

  Optional name of a column of the data (for example variety or region).
  The treatment effect is fitted separately at each of its levels, as an
  interaction, and
  [`pollin_effects()`](https://max578.github.io/pollinatoR/reference/pollin_effects.md)
  then reports one contrast per level. The column is treated as a
  factor.

- ...:

  Passed to
  [`glmmTMB::glmmTMB()`](https://rdrr.io/pkg/glmmTMB/man/glmmTMB.html).

## Value

An object of class `pollin_fit`: a list holding the `estimand`, the
`engine`, the fitted models `fits` (named by response), the model `data`
for each response, the right-hand-side terms `rhs`, the `families`, the
`by` column (or `NULL`), and the further fitting arguments `args`.

## Details

Rows with a missing value in the response, the treatment, a unit column
or the `by` column are dropped for that response's model.

## See also

[`pollin_effects()`](https://max578.github.io/pollinatoR/reference/pollin_effects.md)
for estimates and intervals on the estimand's scale.

## Examples

``` r
# \donttest{
set.seed(1)
d <- data.frame(
  trt = rep(rep(c("Cage", "Open"), each = 10L), 8L),
  Site = rep(paste0("s", seq_len(8L)), each = 20L)
)
d$pods <- rnbinom(nrow(d), mu = exp(3 + 0.3 * (d$trt == "Open")), size = 10)
pd <- pollin_data(d, "pods", "trt", unit = "Site")
ds <- pollin_design(pd, c(excluded = "Cage", open = "Open"), "rcbd")
est <- pollin_estimand(ds, "dependence", "ratio", "loss")
pollin_fit(est)
#> <pollin_fit>
#>   engine    : glmmTMB 
#>   family    : dependence 
#>   contrast  : Open vs Cage 
#>   model     : <response> ~ trt + (1 | Site) 
#>   pods      : nbinom2, 160 rows
# }
```

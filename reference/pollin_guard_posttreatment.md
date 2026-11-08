# Guard against adjusting for a covariate downstream of treatment

Adjusting a treatment contrast for a variable the treatment itself
changes (for example vegetative weight, or a pollinator visit count
under a cage) blocks part of the effect and can open a collider path.
Whether a covariate is downstream of treatment is a fact about the
system that the data cannot reveal. The declaration is the user's:
`downstream` names the variables the user holds to be affected by
treatment, and the guard stops if any adjustment covariate in `x` is
among them.

## Usage

``` r
pollin_guard_posttreatment(x, downstream)
```

## Arguments

- x:

  A
  [`pollin_data()`](https://max578.github.io/pollinatoR/reference/pollin_data.md),
  [`pollin_design()`](https://max578.github.io/pollinatoR/reference/pollin_design.md)
  or
  [`pollin_estimand()`](https://max578.github.io/pollinatoR/reference/pollin_estimand.md)
  object.

- downstream:

  Character vector of variables declared downstream of treatment.
  `character(0)` declares none. There is no default.

## Value

`x`, invisibly, when no adjustment covariate is declared downstream.

## Examples

``` r
df <- data.frame(
  pods = rpois(40, 30),
  trt = rep(c("Cage", "Open"), 20),
  Site = rep(letters[1:4], each = 10),
  rain = rep(c(250, 300, 320, 410), each = 10),
  veg_weight = rnorm(40, 100, 10)
)
pd <- pollin_data(df, "pods", "trt", unit = "Site", adjust = "rain")
pollin_guard_posttreatment(pd, downstream = "veg_weight")

pd2 <- pollin_data(df, "pods", "trt", unit = "Site",
                   adjust = c("rain", "veg_weight"))
try(pollin_guard_posttreatment(pd2, downstream = "veg_weight"))
#> Error : Adjustment covariate `veg_weight` is declared downstream of `trt`. Adjusting for a downstream covariate removes part of the treatment effect; drop it from `adjust`.
```

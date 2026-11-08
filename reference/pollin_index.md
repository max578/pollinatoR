# A menu of uncertainty-aware dependence indices

Returns several dependence indices for an open-versus-exclusion
contrast, each with an interval where one is available: the Relative
Interaction Index (bootstrap interval), the log response ratio
(closed-form interval, the meta-analytic workhorse), and the
proportional dependence ratio. Provided so an analysis can report a
defensible, interval-bearing index instead of a bare ratio.

## Usage

``` r
pollin_index(
  open,
  cage,
  n_boot = 2000L,
  conf = 0.95,
  seed = NULL,
  response = NULL,
  unit = NULL
)
```

## Arguments

- open:

  Numeric vector of the open-arm (pollinator-accessible) values, or a
  [`pollin_design()`](https://max578.github.io/pollinatoR/reference/pollin_design.md)
  with an `"excluded"` and an `"open"` arm.

- cage:

  Numeric vector of the exclusion-arm values. Leave it empty when `open`
  is a design.

- n_boot:

  Number of bootstrap resamples.

- conf:

  Confidence level.

- seed:

  Optional integer seed for reproducibility.

- response:

  With a design, the response column to use. It may be left `NULL` when
  the data declare one response.

- unit:

  With a design, the unit column the arm means are taken over, one of
  the design's unit nesting. `NULL` (the default) uses the outermost
  unit; an inner unit is taken within its outer units.

## Value

A data frame of indices with estimates and intervals where available.

## Details

Given a
[`pollin_design()`](https://max578.github.io/pollinatoR/reference/pollin_design.md),
the indices are computed from unit-level arm means, as described in
[`pollin_rii()`](https://max578.github.io/pollinatoR/reference/pollin_rii.md).

## Examples

``` r
set.seed(1)
pollin_index(rpois(50, 9), rpois(50, 6), n_boot = 500, seed = 1)
#>              index  estimate  conf.low conf.high
#> 1              RII 0.2125984 0.1599841 0.2691237
#> 2             lnRR 0.4317824 0.3064261 0.5571387
#> 3 dependence_ratio 0.3506494        NA        NA
#>                                                              note
#> 1                                                    bootstrap CI
#> 2                                                  closed-form CI
#> 3 point; do not regress on open yield (see pollin_guard_spurious)

# From a design with two responses: name the one to use
df <- data.frame(
  pods = rpois(48, 30),
  seed_wt = rgamma(48, 20, 5),
  trt = rep(c("Cage", "Open"), 24),
  Site = rep(letters[1:4], each = 12)
)
pd <- pollin_data(df, c("pods", "seed_wt"), "trt", unit = "Site")
ds <- pollin_design(pd, c(excluded = "Cage", open = "Open"), "multi_site")
pollin_index(ds, response = "pods", n_boot = 200, seed = 1)
#>              index    estimate    conf.low  conf.high
#> 1              RII -0.02251271 -0.08115501 0.02117663
#> 2             lnRR -0.04503303 -0.16107187 0.07100581
#> 3 dependence_ratio -0.04606241          NA         NA
#>                                                              note
#> 1                                                    bootstrap CI
#> 2                                                  closed-form CI
#> 3 point; do not regress on open yield (see pollin_guard_spurious)
```

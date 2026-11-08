# Estimate exposure moderation of the contribution

Asks which exposures (e.g. honey-bee versus other-insect abundance)
drive the open-minus-exclusion contribution, as treatment-by-exposure
interactions. Each exposure enters through a within-between (Mundlak)
decomposition: a between-unit part (the unit mean) and a within-unit
part (the deviation from it). The **within-unit** treatment-by-exposure
coefficient is robust to any unit-level confounder – observed or not,
including an unrecorded genotype – because it is identified from
variation that occurs within the unit, where the confounder is constant.
The between-unit part is reported but flagged as confounder-vulnerable.
Adjustment covariates enter as fixed effects.

## Usage

``` r
pollin_moderation(x, outcome, exposures = NULL, family = NULL, ...)
```

## Arguments

- x:

  A
  [`pollin_design()`](https://max578.github.io/pollinatoR/reference/pollin_design.md)
  object whose data declare `exposures`, or a deprecated
  [`pollin_spec()`](https://max578.github.io/pollinatoR/reference/pollin_spec.md)
  object with `exposures` set.

- outcome:

  Single outcome to model.

- exposures:

  Exposures to use. `NULL` (the default) uses every declared exposure.

- family:

  Optional `glmmTMB` family (chosen from the outcome's values if
  missing).

- ...:

  Passed to
  [`glmmTMB::glmmTMB()`](https://rdrr.io/pkg/glmmTMB/man/glmmTMB.html).

## Value

An object of class `pollin_moderation` with a tidy `table` of the
treatment-by-exposure terms (each tagged `within` or `between`),
FDR-adjusted across the within terms, and the fitted model in `fit`.

## Details

From a
[`pollin_design()`](https://max578.github.io/pollinatoR/reference/pollin_design.md),
the unit of the split is the design's slope unit (with the unit levels
outside it), or the outermost unit when the design carries no slope; the
treatment reference level is the excluded arm, or the open arm when the
design has no excluded arm. A deprecated
[`pollin_spec()`](https://max578.github.io/pollinatoR/reference/pollin_spec.md)
is still accepted and splits on its `latent_unit`.

## Examples

``` r
# \donttest{
set.seed(1)
site <- rep(paste0("s", seq_len(14L)), each = 20L)
trt <- rep(rep(c("Cage", "Open"), each = 10L), 14L)
bees <- rep(rpois(14L, 8), each = 20L) + rpois(280L, 1)
y <- rpois(length(trt), exp(1.6 + 0.3 * (trt == "Open") +
                            0.02 * bees * (trt == "Open")))
d <- data.frame(y = y, trt = trt, Site = site, bees = bees)
pd <- pollin_data(d, "y", "trt", unit = "Site", exposures = "bees")
ds <- pollin_design(pd, c(excluded = "Cage", open = "Open"), "multi_site")
pollin_moderation(ds, outcome = "y")
#> <pollin_moderation>  outcome: y  (treatment x exposure)
#>                  term    estimate  std.error    p.value            component
#>   trtOpen:bees_within -0.11348610 0.04948954 0.02184036      within (robust)
#>  trtOpen:bees_between  0.08148486 0.04938641 0.09895419 between (vulnerable)
#>       p.fdr
#>  0.02184036
#>          NA
# }
```

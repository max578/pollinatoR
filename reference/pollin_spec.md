# Declare a pollinator-exclusion study

Captures the structure of an open-versus-exclusion experiment in one
object so the downstream verbs share a single, validated description of
the data. The specification is deliberately general: outcomes, the
exclusion treatment, the nested random structure, optional exposures
(the pollinator community), optional adjustment covariates, an optional
fixed grouping, and an optional latent unit on which a random treatment
slope is carried to absorb unobserved between-unit heterogeneity (for
example an unrecorded genotype).

## Usage

``` r
pollin_spec(
  data,
  outcomes,
  treatment,
  control,
  random,
  exposures = NULL,
  adjust = NULL,
  group = NULL,
  latent_unit = NULL
)
```

## Arguments

- data:

  A data frame, one row per observational unit (e.g. one plant).

- outcomes:

  Character vector of outcome column names (yield components).

- treatment:

  Name of the treatment column (a factor with two or more levels, e.g.
  caged versus open).

- control:

  Treatment level that represents exclusion / the reference (e.g.
  `"Cage"`). The contribution is computed relative to this level.

- random:

  Character vector of grouping columns giving the nesting, outer to
  inner (e.g. `c("Site", "Distance", "Transect")`). Translated into
  nested random intercepts.

- exposures:

  Optional character vector of exposure columns (e.g. insect abundances)
  used by
  [`pollin_moderation()`](https://max578.github.io/pollinatoR/reference/pollin_moderation.md).

- adjust:

  Optional character vector of adjustment covariates (e.g. landscape,
  climate) – context, not the inferential target.

- group:

  Optional name of a fixed grouping column (e.g. region).

- latent_unit:

  Optional name of the unit carrying the random treatment slope
  (defaults to the outermost `random`). This absorbs unobserved
  unit-level modifiers of the contribution, such as genotype.

## Value

An object of class `pollin_spec`.

## Details

`pollin_spec()` is deprecated:
[`pollin_data()`](https://max578.github.io/pollinatoR/reference/pollin_data.md)
and
[`pollin_design()`](https://max578.github.io/pollinatoR/reference/pollin_design.md)
declare the same study in two layers. It warns once per session and is
removed at version 0.2.0.

## Examples

``` r
df <- data.frame(
  y = rpois(40, 5), trt = rep(c("Cage", "Open"), 20),
  Site = rep(letters[1:4], each = 10)
)
pollin_spec(df, outcomes = "y", treatment = "trt", control = "Cage",
            random = "Site")
#> Warning: `pollin_spec()` is deprecated; use `pollin_data()` + `pollin_design()`. It is removed at pollinatoR 0.2.0. This warning is shown once per session.
#> <pollin_spec>
#>   observations : 40 
#>   outcomes     : y 
#>   treatment    : trt (control = Cage ; 2 levels)
#>   random       : Site  [slope on: Site ]
```

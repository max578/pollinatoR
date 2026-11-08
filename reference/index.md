# Package index

## Data

Declare the unit-of-record table and check it.

- [`pollin_data()`](https://max578.github.io/pollinatoR/reference/pollin_data.md)
  : Declare the unit-of-record table of a pollination study
- [`pollin_validate()`](https://max578.github.io/pollinatoR/reference/pollin_validate.md)
  : Validate the structure and internal consistency of a study's data

## Design

Map treatment levels to arms and set the layout.

- [`pollin_design()`](https://max578.github.io/pollinatoR/reference/pollin_design.md)
  : Declare the design of a pollination experiment

## Estimand

Name the contrast, the scale and the basis.

- [`pollin_estimand()`](https://max578.github.io/pollinatoR/reference/pollin_estimand.md)
  : Declare the estimand of a pollination experiment
- [`pollin_assumptions()`](https://max578.github.io/pollinatoR/reference/pollin_assumptions.md)
  : List the assumptions an estimand rests on

## Fit

- [`pollin_fit()`](https://max578.github.io/pollinatoR/reference/pollin_fit.md)
  : Fit the outcome models for an estimand

## Effects

Estimates and intervals on the estimand’s scale.

- [`pollin_effects()`](https://max578.github.io/pollinatoR/reference/pollin_effects.md)
  : Estimate the effect on the estimand's scale
- [`pollin_contribution()`](https://max578.github.io/pollinatoR/reference/pollin_contribution.md)
  : Estimate the pollinator contribution per outcome
- [`pollin_moderation()`](https://max578.github.io/pollinatoR/reference/pollin_moderation.md)
  : Estimate exposure moderation of the contribution

## Guards

- [`pollin_guard_spurious()`](https://max578.github.io/pollinatoR/reference/pollin_guard_spurious.md)
  : Guard against the spurious index-versus-yield regression
- [`pollin_guard_bounded()`](https://max578.github.io/pollinatoR/reference/pollin_guard_bounded.md)
  : Guard against averaging bounded unit-level indices
- [`pollin_guard_ecological()`](https://max578.github.io/pollinatoR/reference/pollin_guard_ecological.md)
  : Guard against reading a between-unit slope as a within-unit one
- [`pollin_guard_posttreatment()`](https://max578.github.io/pollinatoR/reference/pollin_guard_posttreatment.md)
  : Guard against adjusting for a covariate downstream of treatment

## Indices and utilities

- [`pollin_index()`](https://max578.github.io/pollinatoR/reference/pollin_index.md)
  : A menu of uncertainty-aware dependence indices
- [`pollin_rii()`](https://max578.github.io/pollinatoR/reference/pollin_rii.md)
  : Relative Interaction Index with a bootstrap interval
- [`pollin_compensation()`](https://max578.github.io/pollinatoR/reference/pollin_compensation.md)
  : Summarise yield-component trade-offs (compensation)
- [`pollin_within_between()`](https://max578.github.io/pollinatoR/reference/pollin_within_between.md)
  : Within-between (Mundlak) decomposition of a numeric variable
- [`pollin_fdr()`](https://max578.github.io/pollinatoR/reference/pollin_fdr.md)
  : Benjamini-Hochberg false-discovery-rate adjustment

## Deprecated

- [`pollin_spec()`](https://max578.github.io/pollinatoR/reference/pollin_spec.md)
  : Declare a pollinator-exclusion study

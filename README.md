# pollinatoR

Uncertainty-aware pollinator-dependence analytics for pollinator-exclusion
experiments (open versus caged, and related contrasts). General to
open-versus-exclusion designs across crops; not specific to any one study.

## Why

The field largely reports a **bare ratio index** (e.g. the Relative Interaction
Index) regressed on yield, across small panels, with single univariate models.
`pollinatoR` replaces that with:

- **`pollin_contribution()`** — the contribution (open − exclusion) as a model-based
  average treatment effect with propagated uncertainty (`glmmTMB` +
  `marginaleffects`), carrying a random treatment slope on a latent unit to absorb
  unobserved between-unit heterogeneity (e.g. an unrecorded genotype).
- **`pollin_moderation()`** — which exposures drive the contribution, via a
  within-between (Mundlak) decomposition so the within-unit estimate is robust to
  unit-level confounders.
- **`pollin_index()` / `pollin_rii()`** — uncertainty-aware indices (RII with a
  bootstrap interval; the log response ratio with a closed-form interval).
- **`pollin_guard_spurious()`** — a permutation guard for the index-versus-yield
  regression, where index and yield share a term (Pearson 1897).
- **`pollin_validate()`** — structural / internal-consistency integrity checks.
- **`pollin_compensation()`** — yield-component trade-off summary.

## Dependencies

CRAN only. Core: `glmmTMB`, `marginaleffects`, `stats`. Optional (`Suggests`):
`DHARMa`, `emmeans`, `brms` (Bayesian, heavy — pulls a Stan toolchain), `boot`,
`fst`, `data.table`. No non-CRAN dependencies.

## Status

Development version (`0.0.0.9000`), released under the MIT licence. This is a
personal research-and-development package; the API may still change before a
`0.1.0` tag.

## Installation

From GitHub:

```r
# install.packages("pak")
pak::pak("max578/pollinatoR")

library(pollinatoR)
vignette("pollinatoR")
```

## Acknowledgements

Developed with support from the Grains Research and Development Corporation
(GRDC project USQ2401-002RTX) and Analytics for the Australian Grains Industry
(AAGI).

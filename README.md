# pollinatoR

`pollinatoR` estimates how much of a crop's yield depends on insect pollinators,
from experiments that compare plants with and without pollinator access: caged
or bagged plants against open ones, and hand-supplemented plants where a trial
has them. Estimates come from mixed models, with intervals, on a scale and basis
you declare. The methods apply to any crop studied with such a design.

## Installation

```r
# install.packages("remotes")
remotes::install_github("max578/pollinatoR")
```

## A minimal example

An analysis is declared in layers: the data, the design, the estimand. It is
then fitted, and the effect is read on the declared scale.

```r
library(pollinatoR)

set.seed(1)
d <- data.frame(
  trt = rep(rep(c("Cage", "Open"), each = 10L), 8L),
  Site = rep(paste0("s", seq_len(8L)), each = 20L)
)
d$pods <- rnbinom(nrow(d), mu = exp(3 + 0.3 * (d$trt == "Open")), size = 10)

pd <- pollin_data(d, response = "pods", treatment = "trt", unit = "Site")
ds <- pollin_design(pd, arms = c(excluded = "Cage", open = "Open"), layout = "rcbd")
est <- pollin_estimand(ds, family = "dependence", scale = "ratio", basis = "loss")
pollin_assumptions(est)
fit <- pollin_fit(est)
pollin_effects(fit)
```

`basis` has no default. A loss share (the part of present yield that pollinators
supply) and a gain (the rise over a crop without them) are different numbers:
excluded plants yielding 75% of open ones give a loss share of 25% and a gain of
33%.

## Functions by layer

| Layer | Function | Purpose |
|---|---|---|
| Data | `pollin_data()` | Declare responses, treatment, units, groups, exposures and adjustment columns |
| | `pollin_validate()` | Structural checks: missing values, expected levels, identities |
| Design | `pollin_design()` | Map treatment levels to the excluded, open and supplemented arms; set the layout |
| Estimand | `pollin_estimand()` | Deficit or dependence, on the difference, ratio, RII or lnRR scale, with a loss or gain basis |
| | `pollin_assumptions()` | The assumptions the estimand rests on |
| Fit | `pollin_fit()` | One `glmmTMB` model per response, optionally with the treatment effect fitted per level of a grouping column (`by`) |
| Effects | `pollin_effects()` | Estimates with delta or parametric bootstrap intervals, overall or per level of `by` |
| | `pollin_contribution()` | Dependence on the difference scale, with FDR-adjusted p-values and the treatment-slope SD |
| | `pollin_moderation()` | Treatment-by-exposure terms, split within and between units |
| Guards | `pollin_guard_spurious()` | Index-versus-yield slope against a permutation null |
| | `pollin_guard_bounded()` | Refuses to average unit-level ratios unless asked |
| | `pollin_guard_ecological()` | Warns when an exposure slope exists only between units |
| | `pollin_guard_posttreatment()` | Refuses an adjustment covariate declared downstream of treatment |
| Indices and utilities | `pollin_index()`, `pollin_rii()` | Descriptive indices with intervals, from arm means |
| | `pollin_compensation()` | Trade-offs among yield components |
| | `pollin_within_between()` | Between-unit and within-unit split of a variable |
| | `pollin_fdr()` | Benjamini-Hochberg adjustment |

`pollin_spec()` is deprecated and is removed at version 0.2.0; the
getting-started vignette shows the equivalent calls.

## Vignettes

```r
vignette("pollinatoR")
vignette("three-arm-trial", package = "pollinatoR")
vignette("split-plot-trial", package = "pollinatoR")
vignette("matched-pair-field-study", package = "pollinatoR")
```

## Status

Released under the MIT licence. This is a personal research and development
package, and the interface may still change before version 1.0.0.

## Acknowledgements

Developed with support from the Grains Research and Development Corporation
(GRDC project USQ2401-002RTX) and Analytics for the Australian Grains Industry
(AAGI).

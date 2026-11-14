# pollinatoR

`pollinatoR` estimates how much of a crop's yield depends on insect pollinators,
from experiments that compare plants with and without pollinator access: caged
or bagged plants against open ones, and hand-supplemented plants where a trial
has them. Estimates come from mixed models, with intervals, on a scale and basis
you declare. The estimand can be carried to a population of fields and on to
an economic value, with its uncertainty kept by working draw by draw. The
methods apply to any crop studied with such a design.

## Installation

```r
# install.packages("remotes")
remotes::install_github("max578/pollinatoR")
```

### Optional Bayesian engine

`pollin_fit(engine = "flexyBayes")` fits on INLA through the flexyBayes
package. Neither is on CRAN, and neither is needed for the rest of the package.

```r
remotes::install_github("AAGI-AUS/flexyBayes")
install.packages(
  "INLA",
  repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable"),
  dependencies = TRUE
)
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

## From dependence to value

Draws of the estimand, multiplied by price and current production, give the
value of pollination with an interval. The channel has no default: the
`dependence` channel values the production lost without pollinators and needs
a loss share.

```r
dr <- pollin_draws(fit, n = 1000L, seed = 1L)
mk <- pollin_market(quantity = 12000, price = 650, unit = "t",
                    currency = "AUD", year = 2024, source = "Illustrative inputs")
pollin_value(dr, mk, channel = "dependence")
```

## Functions by layer

| Layer | Function | Purpose |
|---|---|---|
| Data | `pollin_data()` | Declare responses, treatment, units, groups, exposures and adjustment columns |
| | `pollin_validate()` | Structural checks: missing values, expected levels, identities |
| Design | `pollin_design()` | Map treatment levels to the excluded, open, supplemented and sham arms; set the layout |
| Estimand | `pollin_estimand()` | Deficit, dependence or cage effect, on the difference, ratio, RII or lnRR scale, with a loss or gain basis; a total effect, or a direct effect with named mediators held fixed |
| | `pollin_assumptions()` | The assumptions the estimand rests on |
| Identification | `pollin_dag()` | Read a causal graph (`dagitty`): adjustment sets for the total and direct effects and for each exposure, and the class of each covariate |
| | `pollin_dag_check()` | Test the graph's implied independencies on the data |
| | `pollin_bound()` | Bound the loss share over a range of cage effects when there is no sham arm |
| | `pollin_sensitivity()` | E-value of a ratio of arm means; link-scale bias to the null for a moderation term |
| Fit | `pollin_fit()` | One model per response, on `glmmTMB` or, with `engine = "flexyBayes"`, a Bayesian fit on INLA; the treatment effect fitted per level of a column (`by`), or pooled across its levels (`pooled`) |
| Effects | `pollin_effects()` | Estimates with delta, parametric bootstrap or posterior intervals, overall or per level of `by` |
| | `pollin_cate()` | The effect at each level of a column with a random treatment slope, partially pooled |
| | `pollin_contribution()` | Dependence on the difference scale, with FDR-adjusted p-values and the treatment-slope SD |
| | `pollin_moderation()` | Treatment-by-exposure terms, split within and between units, with a graph's adjustment set and the exposure's slope in the excluded arm as a negative control |
| Population | `pollin_population()` | Strata of one column (distance bands, regions) weighted by area, passed to `pollin_estimand()` |
| Draws | `pollin_draws()` | Draws of the arm means and the estimand, by joint simulation of the fixed and random effects, parametric bootstrap, or from the posterior of a flexyBayes fit |
| Value | `pollin_market()` | Current production and price, as values or draws, with unit, currency, price year and source |
| | `pollin_value()` | Value per draw through the `dependence` or `deficit` channel, with the assumptions it rests on |
| Guards | `pollin_guard_spurious()` | Index-versus-yield slope against a permutation null |
| | `pollin_guard_bounded()` | Refuses to average unit-level ratios unless asked |
| | `pollin_guard_ecological()` | Warns when an exposure slope exists only between units |
| | `pollin_guard_posttreatment()` | Refuses an adjustment covariate downstream of treatment, declared or read from a graph |
| | `pollin_guard_doublecount()` | Adds values draw by draw and refuses overlapping attributions |
| Indices and utilities | `pollin_index()`, `pollin_rii()` | Descriptive indices with intervals, from arm means |
| | `pollin_compensation()` | Trade-offs among yield components |
| | `pollin_within_between()` | Between-unit and within-unit split of a variable |
| | `pollin_fdr()` | Benjamini-Hochberg adjustment |

`pollin_spec()` was removed at version 0.2.0; the getting-started vignette
shows the equivalent calls.

## Vignettes

```r
vignette("pollinatoR")
vignette("three-arm-trial", package = "pollinatoR")
vignette("split-plot-trial", package = "pollinatoR")
vignette("matched-pair-field-study", package = "pollinatoR")
vignette("identification", package = "pollinatoR")
vignette("dependence-to-value", package = "pollinatoR")
vignette("bayesian-fits", package = "pollinatoR")
```

## Status

Released under the MIT licence. This is a personal research and development
package, and the interface may still change before version 1.0.0.

## Acknowledgements

Developed with support from the Grains Research and Development Corporation
(GRDC project USQ2401-002RTX) and Analytics for the Australian Grains Industry
(AAGI).

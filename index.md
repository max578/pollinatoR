# pollinatoR

`pollinatoR` estimates how much of a crop’s yield depends on insect
pollinators, from experiments that compare plants with and without
pollinator access: caged or bagged plants against open ones, and
hand-supplemented plants where a trial has them. Estimates come from
mixed models, with intervals, on a scale and basis you declare. The
methods apply to any crop studied with such a design.

## Installation

\
`# install.packages("remotes")`\
`remotes``::`[`install_github`](https://rdrr.io/pkg/remotes/man/install_github.html)`(``"max578/pollinatoR"``)`

## A minimal example

An analysis is declared in layers: the data, the design, the estimand.
It is then fitted, and the effect is read on the declared scale.

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`pollinatoR`](https://github.com/max578/pollinatoR)`)`\
\
[`set.seed`](https://rdrr.io/r/base/Random.html)`(``1``)`\
`d`` ``<-`` `[`data.frame`](https://rdrr.io/r/base/data.frame.html)`(`\
`  trt ``=`` `[`rep`](https://rdrr.io/r/base/rep.html)`(`[`rep`](https://rdrr.io/r/base/rep.html)`(`[`c`](https://rdrr.io/r/base/c.html)`(``"Cage"``, ``"Open"``)``, each ``=`` ``10L``)``, ``8L``)``,`\
`  Site ``=`` `[`rep`](https://rdrr.io/r/base/rep.html)`(`[`paste0`](https://rdrr.io/r/base/paste.html)`(``"s"``, `[`seq_len`](https://rdrr.io/r/base/seq.html)`(``8L``)``)``, each ``=`` ``20L``)`\
`)`\
`d``$``pods`` ``<-`` `[`rnbinom`](https://rdrr.io/r/stats/NegBinomial.html)`(`[`nrow`](https://rdrr.io/r/base/nrow.html)`(``d``)``, mu ``=`` `[`exp`](https://rdrr.io/r/base/Log.html)`(``3`` ``+`` ``0.3`` ``*`` ``(``d``$``trt`` ``==`` ``"Open"``)``)``, size ``=`` ``10``)`\
\
`pd`` ``<-`` `[`pollin_data`](https://max578.github.io/pollinatoR/reference/pollin_data.md)`(``d``, response ``=`` ``"pods"``, treatment ``=`` ``"trt"``, unit ``=`` ``"Site"``)`\
`ds`` ``<-`` `[`pollin_design`](https://max578.github.io/pollinatoR/reference/pollin_design.md)`(``pd``, arms ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``excluded ``=`` ``"Cage"``, open ``=`` ``"Open"``)``, layout ``=`` ``"rcbd"``)`\
`est`` ``<-`` `[`pollin_estimand`](https://max578.github.io/pollinatoR/reference/pollin_estimand.md)`(``ds``, family ``=`` ``"dependence"``, scale ``=`` ``"ratio"``, basis ``=`` ``"loss"``)`\
[`pollin_assumptions`](https://max578.github.io/pollinatoR/reference/pollin_assumptions.md)`(``est``)`\
`fit`` ``<-`` `[`pollin_fit`](https://max578.github.io/pollinatoR/reference/pollin_fit.md)`(``est``)`\
[`pollin_effects`](https://max578.github.io/pollinatoR/reference/pollin_effects.md)`(``fit``)`

`basis` has no default. A loss share (the part of present yield that
pollinators supply) and a gain (the rise over a crop without them) are
different numbers: excluded plants yielding 75% of open ones give a loss
share of 25% and a gain of 33%.

## Functions by layer

| Layer | Function | Purpose |
|----|----|----|
| Data | [`pollin_data()`](https://max578.github.io/pollinatoR/reference/pollin_data.md) | Declare responses, treatment, units, groups, exposures and adjustment columns |
|  | [`pollin_validate()`](https://max578.github.io/pollinatoR/reference/pollin_validate.md) | Structural checks: missing values, expected levels, identities |
| Design | [`pollin_design()`](https://max578.github.io/pollinatoR/reference/pollin_design.md) | Map treatment levels to the excluded, open and supplemented arms; set the layout |
| Estimand | [`pollin_estimand()`](https://max578.github.io/pollinatoR/reference/pollin_estimand.md) | Deficit or dependence, on the difference, ratio, RII or lnRR scale, with a loss or gain basis |
|  | [`pollin_assumptions()`](https://max578.github.io/pollinatoR/reference/pollin_assumptions.md) | The assumptions the estimand rests on |
| Fit | [`pollin_fit()`](https://max578.github.io/pollinatoR/reference/pollin_fit.md) | One `glmmTMB` model per response, optionally with the treatment effect fitted per level of a grouping column (`by`) |
| Effects | [`pollin_effects()`](https://max578.github.io/pollinatoR/reference/pollin_effects.md) | Estimates with delta or parametric bootstrap intervals, overall or per level of `by` |
|  | [`pollin_contribution()`](https://max578.github.io/pollinatoR/reference/pollin_contribution.md) | Dependence on the difference scale, with FDR-adjusted p-values and the treatment-slope SD |
|  | [`pollin_moderation()`](https://max578.github.io/pollinatoR/reference/pollin_moderation.md) | Treatment-by-exposure terms, split within and between units |
| Guards | [`pollin_guard_spurious()`](https://max578.github.io/pollinatoR/reference/pollin_guard_spurious.md) | Index-versus-yield slope against a permutation null |
|  | [`pollin_guard_bounded()`](https://max578.github.io/pollinatoR/reference/pollin_guard_bounded.md) | Refuses to average unit-level ratios unless asked |
|  | [`pollin_guard_ecological()`](https://max578.github.io/pollinatoR/reference/pollin_guard_ecological.md) | Warns when an exposure slope exists only between units |
|  | [`pollin_guard_posttreatment()`](https://max578.github.io/pollinatoR/reference/pollin_guard_posttreatment.md) | Refuses an adjustment covariate declared downstream of treatment |
| Indices and utilities | [`pollin_index()`](https://max578.github.io/pollinatoR/reference/pollin_index.md), [`pollin_rii()`](https://max578.github.io/pollinatoR/reference/pollin_rii.md) | Descriptive indices with intervals, from arm means |
|  | [`pollin_compensation()`](https://max578.github.io/pollinatoR/reference/pollin_compensation.md) | Trade-offs among yield components |
|  | [`pollin_within_between()`](https://max578.github.io/pollinatoR/reference/pollin_within_between.md) | Between-unit and within-unit split of a variable |
|  | [`pollin_fdr()`](https://max578.github.io/pollinatoR/reference/pollin_fdr.md) | Benjamini-Hochberg adjustment |

[`pollin_spec()`](https://max578.github.io/pollinatoR/reference/pollin_spec.md)
is deprecated and is removed at version 0.2.0; the getting-started
vignette shows the equivalent calls.

## Vignettes

\
[`vignette`](https://rdrr.io/r/utils/vignette.html)`(``"pollinatoR"``)`\
[`vignette`](https://rdrr.io/r/utils/vignette.html)`(``"three-arm-trial"``, package ``=`` ``"pollinatoR"``)`\
[`vignette`](https://rdrr.io/r/utils/vignette.html)`(``"split-plot-trial"``, package ``=`` ``"pollinatoR"``)`\
[`vignette`](https://rdrr.io/r/utils/vignette.html)`(``"matched-pair-field-study"``, package ``=`` ``"pollinatoR"``)`

## Status

Released under the MIT licence. This is a personal research and
development package, and the interface may still change before version
1.0.0.

## Acknowledgements

Developed with support from the Grains Research and Development
Corporation (GRDC project USQ2401-002RTX) and Analytics for the
Australian Grains Industry (AAGI).

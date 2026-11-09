# pollinatoR 0.2.0 (development)

* Golden-value tests now also hold the 0.1.0 outputs: deficit and dependence
  on every scale, with delta and bootstrap intervals, overall and by level, and
  the layered contribution, moderation, validation and index results
  (`validation/golden/make_golden_v0.1.0.R`).
* `pollin_spec()` is removed, as announced in 0.1.0, with its print method
  and the paths that accepted it in `pollin_contribution()`,
  `pollin_moderation()` and `pollin_validate()`; these now take only the
  layered objects. The getting-started vignette shows the equivalent calls.
* `pollin_draws()` returns draws of the two arm means of a fitted estimand
  and of the estimand formed from them, per draw, on its scale and basis. The
  default `method = "simulation"` draws the fixed effects from their sampling
  distribution (fixed-effect uncertainty, as in the delta interval);
  `method = "bootstrap"` keeps the parametric bootstrap refits, whose
  percentile interval is the bootstrap interval of `pollin_effects()`.
  `summary()` gives the median and percentile interval per column.
* `pollin_population()` declares a population of strata (the levels of one
  column, such as distance bands or regions) weighted by area, and
  `pollin_estimand(population = )` accepts it besides `"trial"`. The estimand
  then states the strata and weights among its assumptions, the fit is made
  with the strata as `by`, and `pollin_draws()` adds a `"population"` column:
  the population arm means are the area-weighted sums of the stratum arm
  means, and the estimand is formed from them, a ratio of population means
  rather than a weighted mean of the stratum estimands. `pollin_effects()`
  reports the trial contrasts and does not use the population.
* `pollin_market()` declares the market inputs of a value: current
  production and price as scalars, each optionally with draws, with the unit,
  currency, price year and a source that prints with every value.
* `pollin_value()` values draws of an estimand, per draw, through one of two
  channels named by the user (there is no default): `"dependence"`, the
  loss share times price times production, the value lost under total
  exclusion; or `"deficit"`, the gain of hand pollination over open times
  price times production. Each channel accepts only its matching ratio
  estimand. A gain basis offered to the dependence channel is refused, since
  a gain `g = D / (1 - D)` overstates the loss by `1 / (1 - D)`; a
  descriptive index or a plain number is refused because it has no estimand.
  Market draws are paired with the estimand draws by index. Every print names
  the basis, population, price year, currency, source, method and the
  assumptions of the value.
* Each value carries an attribution: the pollinators its exclusion removed
  (`pollin_value(agent = )`, `"all"` by default), the product valued
  (`product`, the response name by default) and the crop
  (`pollin_market(crop = )`). `pollin_guard_doublecount()` adds two or more
  values draw by draw and refuses overlapping attributions within a crop: a
  part plus its total, the same agent twice, two products of one agent (seed
  plus oil), and two named guilds unless `nested = FALSE` asserts separate,
  non-nested contrasts, an assertion printed with the total. Values must
  share currency, price year and number of draws.

# pollinatoR 0.1.0

An analysis is now declared in layers: the data (`pollin_data()`), the design
(`pollin_design()`) and the estimand (`pollin_estimand()`), then fitted
(`pollin_fit()`) and read on the estimand's scale (`pollin_effects()`).

## Breaking changes

* The first argument of `pollin_validate()`, `pollin_contribution()` and
  `pollin_moderation()` is renamed from `spec` to `x`. Calls that pass it by
  position are unaffected.
* `pollin_contribution()` takes a `pollin_estimand` with
  `family = "dependence"` and `scale = "difference"`. Its table carries the
  columns of `pollin_effects()` plus `slope_sd` and `p.fdr`. With `by`, it fits
  the treatment-by-level interaction and reports each level's own contrast. A
  `pollin_spec` with `by` is fitted as before.
* `pollin_moderation()` takes a `pollin_design` whose data declare
  `exposures`. The default of `exposures` is now `NULL`, meaning every declared
  exposure.
* `pollin_spec()` is deprecated and is removed at version 0.2.0. It warns once
  per session. `pollin_contribution()`, `pollin_moderation()` and
  `pollin_validate()` still accept a `pollin_spec` until then.
* Removed six unused packages from Suggests (DHARMa, emmeans, brms, boot, fst,
  data.table).

## New features

* `pollin_data()` declares the unit-of-record table: response, treatment,
  unit nesting, and optional grouping, denominator, exposure and adjustment
  columns. Roles are checked against the data and the structural checks are
  stored with the object.
* `pollin_design()` declares the arms (excluded, open, supplemented), the
  layout (`crd`, `rcbd`, `split_plot`, `matched_pairs`, `multi_site`) and the
  reference arm. The layout sets the random structure; a random treatment
  slope sits on the site in `multi_site` and is opt-in elsewhere through
  `slope`. Units lacking an arm are counted in a warning.
* `pollin_estimand()` names the estimand: deficit (supplemented versus open)
  or dependence (excluded versus the reference arm), on the difference, ratio,
  RII or lnRR scale. `basis` (`"loss"` or `"gain"`) has no default, so a loss
  share cannot be read as a gain. `pollin_assumptions()` lists the assumptions
  the estimand rests on for its layout.
* `pollin_fit()` fits one `glmmTMB` model per response, with the treatment as
  the fixed effect and the design's random terms. The family is chosen from
  each response's values unless given, singly or per response. With `by`, the
  fixed part is `treatment * by`, so each level of `by` has its own treatment
  effect.
* `pollin_effects()` reports the estimand on its declared scale from averaged
  model predictions (a ratio of expected values, never a mean of unit
  ratios), with a delta interval computed on the log-ratio scale and
  transformed, so a loss share stays below 1 and the RII inside (-1, 1), or
  with a parametric bootstrap interval (`interval = "bootstrap"`, `n_boot`,
  `seed`). Every table records the family, scale and basis. Per-level
  contrasts come from a fit made with `by`; asking for a `by` the fit was not
  made with is an error. The printed result states that delta intervals cover
  fixed-effect uncertainty only.
* `pollin_index()` and `pollin_rii()` accept a `pollin_design` and work from
  unit-level arm means; `pollin_compensation()` accepts a `pollin_data` and
  uses its responses as the components. The vector and data-frame interfaces
  are unchanged and return the same numbers.
* `pollin_guard_bounded()` refuses to average unit-level RII or dependence
  values unless `allow = TRUE`, and reports the index of the arm means beside
  the average.
* `pollin_guard_ecological()` splits an exposure into between-unit and
  within-unit parts and warns when an exposure slope exists only between
  units, or when the two slopes disagree in sign.
* `pollin_guard_posttreatment()` stops when an adjustment covariate is among
  the variables the user declares to be downstream of the treatment.
* New vignettes: *A three-arm trial: deficit and dependence*, *A split-plot
  trial: dependence by variety* and *A matched-pair field study*. The getting-started vignette is rewritten around
  the layered workflow.

## Known limitations

* In the planted-truth recovery study (`validation/planted_truth_S1_report.md`),
  24 of 25 design-by-estimand cells meet the pre-set rule. The matched-pair
  dependence estimate on the difference scale is biased by -0.19 standard
  errors; its interval coverage (0.950) and the ratio-scale estimates pass.
  For matched-pair designs the ratio scale is the better-supported report.

## Moving from earlier calls

| Before | Now |
|---|---|
| `pollin_spec(data, outcomes, treatment, control, random)` | `pollin_data(data, response, treatment, unit)` then `pollin_design(x, arms = c(excluded = control, open = ...), layout)` |
| `pollin_spec(..., latent_unit = "Site")` | `pollin_design(..., layout = "multi_site")`, or `slope = "Site"` with another layout |
| `pollin_spec(..., exposures, adjust)` | `pollin_data(..., exposures, adjust)` |
| `pollin_contribution(spec)` | `pollin_contribution(pollin_estimand(design, "dependence", "difference", basis = "loss"))` |
| Any other scale or interval | `pollin_effects(pollin_fit(estimand))` |
| `pollin_contribution(spec, by = "Region")` | `pollin_effects(pollin_fit(estimand, by = "Region"))`, now one fixed treatment effect per region |
| `pollin_moderation(spec, outcome)` | `pollin_moderation(design, outcome)` |
| `pollin_validate(spec)` | `pollin_validate(pollin_data(...))`, also run when the object is built |
| `pollin_index(open_means, cage_means)` | Unchanged, or `pollin_index(design, response)` |

# pollinatoR 0.0.0.9000

* Initial development scaffold.
* `pollin_spec()` declares an exclusion-experiment study (outcomes, treatment,
  nested random structure, exposures, adjustment, latent grouping).
* `pollin_contribution()` estimates the per-outcome pollinator contribution
  (open minus exclusion) as a model-based average treatment effect with
  propagated uncertainty, carrying a random treatment slope on the latent unit
  to absorb unobserved between-unit (e.g. genotype) heterogeneity.
* `pollin_moderation()` estimates how the contribution varies with exposures
  (treatment-by-exposure conditional effects) using a within-between
  decomposition so the within-unit estimate is robust to site-level confounders.
* `pollin_index()` / `pollin_rii()` provide uncertainty-aware dependence indices,
  including the Relative Interaction Index with a bootstrap interval.
* `pollin_guard_spurious()` tests an index-versus-yield slope against a
  permutation null that isolates the mechanical (shared-term) coupling.
* `pollin_validate()` runs structural and internal-consistency integrity checks.
* `pollin_compensation()` summarises yield-component trade-offs.

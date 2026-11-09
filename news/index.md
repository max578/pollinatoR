# Changelog

## pollinatoR 0.3.0

This release adds a Bayesian engine and per-level effects with partial
pooling. A fit made through flexyBayes on INLA gives posterior draws,
and effects, populations, values and guards take them as they take the
glmmTMB draws. The vignette *Bayesian fits with flexyBayes* compares the
two engines. The flexyBayes engine is experimental in this release; the
validation below says where its intervals are and are not calibrated.

### New features

- `pollin_fit(engine = "flexyBayes")` fits each response with
  [`flexyBayes::flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.html)
  on INLA, with the fixed part of the glmmTMB engine and the design’s
  random terms. The families `nbinom2()`,
  [`poisson()`](https://rdrr.io/r/stats/family.html),
  [`binomial()`](https://rdrr.io/r/stats/family.html),
  `Gamma(link = "log")` and `beta_family()` are mapped to the flexyBayes
  families `negative_binomial`, `poisson`, `binomial`, `gamma` and
  `beta`, with their log or logit links. Priors pass through `...` as
  written. flexyBayes and INLA are suggested packages, and the extension
  API version of flexyBayes is checked when the fit is made.
- `pollin_draws(method = "posterior")`, the default and only method for
  a flexyBayes fit, draws the fixed effects, random effects and variance
  parameters jointly from the posterior and recomputes the arm means per
  draw by the G-computation rule of the glmmTMB engine. With a seed,
  which sets both R’s generator and INLA’s, the draws repeat exactly.
- `pollin_effects(interval = "posterior")`, the default for a flexyBayes
  fit, reports the posterior median, standard deviation and equal-tailed
  interval of the estimand, with `p.value` `NA`. New argument `n_draws`.
- [`pollin_cate()`](https://max578.github.io/pollinatoR/reference/pollin_cate.md)
  reports the contrast at each level of a column with a random treatment
  slope, partially pooled towards the overall effect. The column is the
  design’s slope unit (`pollin_design(slope = )`) or one declared with
  the new `pollin_fit(pooled = )`. On flexyBayes the estimate and
  interval come from the posterior draws. The print points to the
  unpooled per-level contrasts, from a fit with `by` and no slope, read
  with
  [`pollin_effects()`](https://max578.github.io/pollinatoR/reference/pollin_effects.md).

### Validation

- Planted-truth simulations, 200 trials per cell (100 for per-level
  effects), compare the engines and check posterior coverage
  (`validation/planted_truth_S3_report.md`). None of the three pre-set
  criteria is met in full:
  - The two engines agree on the estimate (within 0.25 glmmTMB standard
    errors) in 97% to 100% of trials in every cell. On designs with a
    random treatment slope, glmmTMB stops at a variance on its boundary
    in 4% to 12% of trials; flexyBayes fits every trial.
  - Posterior 95% intervals cover 0.945 to 0.960 on designs with unit
    intercepts (negative binomial, Poisson, deficit) and 0.932 to 0.939
    on the overall effect of pooled fits, but 0.909 to 0.925 on designs
    with a site slope, and 0.925 on a binary response. On slope designs
    the glmmTMB delta interval over-covers (0.989 to 1.000) at about 1.7
    times the posterior width.
  - [`pollin_cate()`](https://max578.github.io/pollinatoR/reference/pollin_cate.md)
    lowers the error of per-level effects by 12% to 15% on average
    against unpooled contrasts, but does better in only 69% to 79% of
    trials, against a 90% target. Posterior per-level coverage is 0.940
    to 0.941.

### Known limitations

- The flexyBayes engine is experimental: its intervals run slightly
  narrow on designs with a random treatment slope (see Validation).
- `pollin_data(denominator = )` does not reach the fit on either engine,
  so a count out of a known number of trials (seeds out of ovules)
  cannot yet be fitted as binomial; one 0/1 row per trial works.
- flexyBayes on INLA does not fit a correlated random treatment slope.
  The flexyBayes engine fits a unit intercept plus an independent
  unit-by-arm effect instead, a restricted form of glmmTMB’s
  unstructured slope: equal variances across arms and a non-negative
  correlation between them.
- The brms backend of flexyBayes is refused, because posterior draws of
  a brms fit are not in the flexyBayes extension API version 1. Families
  other than the five above, Tweedie among them, are refused by name.
- On glmmTMB,
  [`pollin_cate()`](https://max578.github.io/pollinatoR/reference/pollin_cate.md)
  intervals are approximate. The fixed effects and the slope
  coefficients are drawn independently, since glmmTMB gives the
  conditional variance of each coefficient on its own. When the slope
  variance is estimated at zero,
  [`pollin_cate()`](https://max578.github.io/pollinatoR/reference/pollin_cate.md)
  stops.

## pollinatoR 0.2.0

This release adds economic valuation. A fitted estimand can be scaled to
a population of fields and turned into a value with an interval,
computed draw by draw. Every value prints its basis, population, price
year, currency and source.

### Breaking changes

- `pollin_spec()` is removed, as announced in 0.1.0, with its print
  method.
  [`pollin_contribution()`](https://max578.github.io/pollinatoR/reference/pollin_contribution.md),
  [`pollin_moderation()`](https://max578.github.io/pollinatoR/reference/pollin_moderation.md)
  and
  [`pollin_validate()`](https://max578.github.io/pollinatoR/reference/pollin_validate.md)
  no longer accept a `pollin_spec` and take only the layered objects.
  The getting-started vignette shows the equivalent calls.

### New features

- [`pollin_draws()`](https://max578.github.io/pollinatoR/reference/pollin_draws.md)
  returns draws of the two arm means of a fitted estimand and of the
  estimand formed from them, on its scale and basis. The default
  `method = "simulation"` draws the fixed effects from their sampling
  distribution and covers fixed-effect uncertainty, as the delta
  interval does. `method = "bootstrap"` uses the parametric bootstrap
  refits of
  [`pollin_effects()`](https://max578.github.io/pollinatoR/reference/pollin_effects.md)
  and gives the same interval.
  [`summary()`](https://rdrr.io/r/base/summary.html) gives the median
  and percentile interval per column.
- [`pollin_population()`](https://max578.github.io/pollinatoR/reference/pollin_population.md)
  declares a population of strata, the levels of one column (distance
  bands, regions, varieties) weighted by area, and
  `pollin_estimand(population = )` accepts it. The fit is made with the
  strata as `by`, and
  [`pollin_draws()`](https://max578.github.io/pollinatoR/reference/pollin_draws.md)
  adds a `"population"` column. The population estimand is formed from
  the area-weighted arm means: a ratio of population means, not a
  weighted mean of the stratum estimands.
  [`pollin_effects()`](https://max578.github.io/pollinatoR/reference/pollin_effects.md)
  reports the trial contrasts and does not use the population.
- [`pollin_market()`](https://max578.github.io/pollinatoR/reference/pollin_market.md)
  declares current production and price, as scalars or with draws, with
  the unit, currency, price year and a source that prints with every
  value.
- [`pollin_value()`](https://max578.github.io/pollinatoR/reference/pollin_value.md)
  multiplies draws of an estimand by price and production. The channel
  has no default:
  - `"dependence"` takes a loss share and gives the value lost under
    total exclusion of pollinators;
  - `"deficit"` takes the gain of hand pollination over open and gives
    the value of closing that gap.

  Each channel accepts only its own estimand. A gain offered to the
  dependence channel is refused, since `g = D / (1 - D)` overstates the
  loss by `1 / (1 - D)`. A descriptive index or a plain number is
  refused because it has no estimand. Market draws are matched to the
  estimand draws by position.
- [`pollin_guard_doublecount()`](https://max578.github.io/pollinatoR/reference/pollin_guard_doublecount.md)
  adds values draw by draw and stops when two of them count the same
  pollination within a crop: a part plus its total, the same pollinators
  twice, two products of one crop (seed and oil), or two named guilds.
  For guilds, `nested = FALSE` confirms that they come from separate,
  non-nested exclusions, and the statement prints with the total. Each
  value records the pollinators removed (`pollin_value(agent = )`,
  `"all"` by default), the product and the crop
  (`pollin_market(crop = )`).
- New vignette *From dependence to value*. *A matched-pair field study*
  gains a section on scaling the estimand from the sampled fields to a
  growing area.

### Validation

- Golden-value tests check that the 0.1.0 outputs are unchanged: deficit
  and dependence on every scale, with delta and bootstrap intervals,
  overall and by level, and the contribution, moderation, validation and
  index results (`validation/golden/make_golden_v0.1.0.R`).
- Planted-truth calibration of value intervals over 1,000 simulated
  trials per cell (three-arm, split-plot and matched-pair designs, both
  channels, market inputs fixed or drawn, one two-stratum population):
  95% interval coverage between 0.937 and 0.962 in all 11 cells. The
  bias criterion (below 0.1 standard errors) holds in 9 of 11; see Known
  limitations (`validation/planted_truth_S2_report.md`).
- The published valuation of Yocgo et al. (2023, *Scientific Reports*
  13:20108) is reproduced from its supplementary data: totals to the
  dollar and all 20 published figures within their rounding
  (`validation/value_oracle/`).

### Known limitations

- `pollin_contribution(by = )` stops with “jacobian does not match the
  variance-covariance matrix” on nested multi-site data where the `by`
  variable is constant within the slope unit.
- The value summary is the median of the value draws. On the dependence
  channel with market draws it sits slightly low, by 0.13 to 0.14
  standard errors in two planted-truth cells, because the loss share is
  concave in the log ratio and the median of a product falls below the
  product of medians. Interval coverage is unaffected.
- Simulation draws cover fixed-effect uncertainty only. Use
  `method = "bootstrap"` when the random effects should vary as well.
- A value assumes no price response, a total-loss counterfactual for the
  dependence channel, and independence of the estimand and market draws.
  It is gross revenue, before costs.

## pollinatoR 0.1.0

An analysis is now declared in layers: the data
([`pollin_data()`](https://max578.github.io/pollinatoR/reference/pollin_data.md)),
the design
([`pollin_design()`](https://max578.github.io/pollinatoR/reference/pollin_design.md))
and the estimand
([`pollin_estimand()`](https://max578.github.io/pollinatoR/reference/pollin_estimand.md)),
then fitted
([`pollin_fit()`](https://max578.github.io/pollinatoR/reference/pollin_fit.md))
and read on the estimand’s scale
([`pollin_effects()`](https://max578.github.io/pollinatoR/reference/pollin_effects.md)).

### Breaking changes

- The first argument of
  [`pollin_validate()`](https://max578.github.io/pollinatoR/reference/pollin_validate.md),
  [`pollin_contribution()`](https://max578.github.io/pollinatoR/reference/pollin_contribution.md)
  and
  [`pollin_moderation()`](https://max578.github.io/pollinatoR/reference/pollin_moderation.md)
  is renamed from `spec` to `x`. Calls that pass it by position are
  unaffected.
- [`pollin_contribution()`](https://max578.github.io/pollinatoR/reference/pollin_contribution.md)
  takes a `pollin_estimand` with `family = "dependence"` and
  `scale = "difference"`. Its table carries the columns of
  [`pollin_effects()`](https://max578.github.io/pollinatoR/reference/pollin_effects.md)
  plus `slope_sd` and `p.fdr`. With `by`, it fits the treatment-by-level
  interaction and reports each level’s own contrast. A `pollin_spec`
  with `by` is fitted as before.
- [`pollin_moderation()`](https://max578.github.io/pollinatoR/reference/pollin_moderation.md)
  takes a `pollin_design` whose data declare `exposures`. The default of
  `exposures` is now `NULL`, meaning every declared exposure.
- `pollin_spec()` is deprecated and is removed at version 0.2.0. It
  warns once per session.
  [`pollin_contribution()`](https://max578.github.io/pollinatoR/reference/pollin_contribution.md),
  [`pollin_moderation()`](https://max578.github.io/pollinatoR/reference/pollin_moderation.md)
  and
  [`pollin_validate()`](https://max578.github.io/pollinatoR/reference/pollin_validate.md)
  still accept a `pollin_spec` until then.
- Removed six unused packages from Suggests (DHARMa, emmeans, brms,
  boot, fst, data.table).

### New features

- [`pollin_data()`](https://max578.github.io/pollinatoR/reference/pollin_data.md)
  declares the unit-of-record table: response, treatment, unit nesting,
  and optional grouping, denominator, exposure and adjustment columns.
  Roles are checked against the data and the structural checks are
  stored with the object.
- [`pollin_design()`](https://max578.github.io/pollinatoR/reference/pollin_design.md)
  declares the arms (excluded, open, supplemented), the layout (`crd`,
  `rcbd`, `split_plot`, `matched_pairs`, `multi_site`) and the reference
  arm. The layout sets the random structure; a random treatment slope
  sits on the site in `multi_site` and is opt-in elsewhere through
  `slope`. Units lacking an arm are counted in a warning.
- [`pollin_estimand()`](https://max578.github.io/pollinatoR/reference/pollin_estimand.md)
  names the estimand: deficit (supplemented versus open) or dependence
  (excluded versus the reference arm), on the difference, ratio, RII or
  lnRR scale. `basis` (`"loss"` or `"gain"`) has no default, so a loss
  share cannot be read as a gain.
  [`pollin_assumptions()`](https://max578.github.io/pollinatoR/reference/pollin_assumptions.md)
  lists the assumptions the estimand rests on for its layout.
- [`pollin_fit()`](https://max578.github.io/pollinatoR/reference/pollin_fit.md)
  fits one `glmmTMB` model per response, with the treatment as the fixed
  effect and the design’s random terms. The family is chosen from each
  response’s values unless given, singly or per response. With `by`, the
  fixed part is `treatment * by`, so each level of `by` has its own
  treatment effect.
- [`pollin_effects()`](https://max578.github.io/pollinatoR/reference/pollin_effects.md)
  reports the estimand on its declared scale from averaged model
  predictions (a ratio of expected values, never a mean of unit ratios),
  with a delta interval computed on the log-ratio scale and transformed,
  so a loss share stays below 1 and the RII inside (-1, 1), or with a
  parametric bootstrap interval (`interval = "bootstrap"`, `n_boot`,
  `seed`). Every table records the family, scale and basis. Per-level
  contrasts come from a fit made with `by`; asking for a `by` the fit
  was not made with is an error. The printed result states that delta
  intervals cover fixed-effect uncertainty only.
- [`pollin_index()`](https://max578.github.io/pollinatoR/reference/pollin_index.md)
  and
  [`pollin_rii()`](https://max578.github.io/pollinatoR/reference/pollin_rii.md)
  accept a `pollin_design` and work from unit-level arm means;
  [`pollin_compensation()`](https://max578.github.io/pollinatoR/reference/pollin_compensation.md)
  accepts a `pollin_data` and uses its responses as the components. The
  vector and data-frame interfaces are unchanged and return the same
  numbers.
- [`pollin_guard_bounded()`](https://max578.github.io/pollinatoR/reference/pollin_guard_bounded.md)
  refuses to average unit-level RII or dependence values unless
  `allow = TRUE`, and reports the index of the arm means beside the
  average.
- [`pollin_guard_ecological()`](https://max578.github.io/pollinatoR/reference/pollin_guard_ecological.md)
  splits an exposure into between-unit and within-unit parts and warns
  when an exposure slope exists only between units, or when the two
  slopes disagree in sign.
- [`pollin_guard_posttreatment()`](https://max578.github.io/pollinatoR/reference/pollin_guard_posttreatment.md)
  stops when an adjustment covariate is among the variables the user
  declares to be downstream of the treatment.
- New vignettes: *A three-arm trial: deficit and dependence*, *A
  split-plot trial: dependence by variety* and *A matched-pair field
  study*. The getting-started vignette is rewritten around the layered
  workflow.

### Known limitations

- In the planted-truth recovery study
  (`validation/planted_truth_S1_report.md`), 24 of 25 design-by-estimand
  cells meet the pre-set rule. The matched-pair dependence estimate on
  the difference scale is biased by -0.19 standard errors; its interval
  coverage (0.950) and the ratio-scale estimates pass. For matched-pair
  designs the ratio scale is the better-supported report.

### Moving from earlier calls

| Before | Now |
|----|----|
| `pollin_spec(data, outcomes, treatment, control, random)` | `pollin_data(data, response, treatment, unit)` then `pollin_design(x, arms = c(excluded = control, open = ...), layout)` |
| `pollin_spec(..., latent_unit = "Site")` | `pollin_design(..., layout = "multi_site")`, or `slope = "Site"` with another layout |
| `pollin_spec(..., exposures, adjust)` | `pollin_data(..., exposures, adjust)` |
| `pollin_contribution(spec)` | `pollin_contribution(pollin_estimand(design, "dependence", "difference", basis = "loss"))` |
| Any other scale or interval | `pollin_effects(pollin_fit(estimand))` |
| `pollin_contribution(spec, by = "Region")` | `pollin_effects(pollin_fit(estimand, by = "Region"))`, now one fixed treatment effect per region |
| `pollin_moderation(spec, outcome)` | `pollin_moderation(design, outcome)` |
| `pollin_validate(spec)` | `pollin_validate(pollin_data(...))`, also run when the object is built |
| `pollin_index(open_means, cage_means)` | Unchanged, or `pollin_index(design, response)` |

## pollinatoR 0.0.0.9000

- Initial development scaffold.
- `pollin_spec()` declares an exclusion-experiment study (outcomes,
  treatment, nested random structure, exposures, adjustment, latent
  grouping).
- [`pollin_contribution()`](https://max578.github.io/pollinatoR/reference/pollin_contribution.md)
  estimates the per-outcome pollinator contribution (open minus
  exclusion) as a model-based average treatment effect with propagated
  uncertainty, carrying a random treatment slope on the latent unit to
  absorb unobserved between-unit (e.g. genotype) heterogeneity.
- [`pollin_moderation()`](https://max578.github.io/pollinatoR/reference/pollin_moderation.md)
  estimates how the contribution varies with exposures
  (treatment-by-exposure conditional effects) using a within-between
  decomposition so the within-unit estimate is robust to site-level
  confounders.
- [`pollin_index()`](https://max578.github.io/pollinatoR/reference/pollin_index.md)
  /
  [`pollin_rii()`](https://max578.github.io/pollinatoR/reference/pollin_rii.md)
  provide uncertainty-aware dependence indices, including the Relative
  Interaction Index with a bootstrap interval.
- [`pollin_guard_spurious()`](https://max578.github.io/pollinatoR/reference/pollin_guard_spurious.md)
  tests an index-versus-yield slope against a permutation null that
  isolates the mechanical (shared-term) coupling.
- [`pollin_validate()`](https://max578.github.io/pollinatoR/reference/pollin_validate.md)
  runs structural and internal-consistency integrity checks.
- [`pollin_compensation()`](https://max578.github.io/pollinatoR/reference/pollin_compensation.md)
  summarises yield-component trade-offs.

# Planted-truth recovery rule, Stage 1 (criterion C2)

Fixed 2026-11-08, before any simulation run of the suite. This file is committed on
its own; the commit time is the evidence that the rule preceded the results. Nothing
here changes after the run. A changed rule is a new dated file next to this one,
with the reason, and both sets of results are reported.

## 1. What is recovered

Each cell recovers one estimand built with the public API:
`pollin_data()` -> `pollin_design()` -> `pollin_estimand(population = "trial")` ->
`pollin_fit(engine = "glmmTMB", family = glmmTMB::nbinom2())` ->
`pollin_effects(interval = "delta", conf_level = 0.95)`, with `by = NULL`.

**What `pollin_effects()` targets.** The delta path calls
`marginaleffects::avg_comparisons()` on the `glmmTMB` fit. Predictions use the
`glmmTMB` default `re.form = NULL`: every row of the model data is predicted under
the `high` arm and under the `low` arm *with that row's estimated random effects*,
and the predictions are averaged over the rows. (Checked before this rule was
written, on one simulated three-arm data set: the package's difference estimate,
11.866, equals the manual average of `predict(type = "response")` with
`re.form = NULL`, 11.866, and not the fixed-effects-only average, 10.682.) The
target is therefore the arms' expected responses **averaged over the realised rows
of the simulated trial, conditional on the realised random effects**, not the
population-marginal value over the random-effect distribution.

**Truth per seed.** For simulated data set `s` with rows `i = 1..n`, fixed
log-means `beta[arm]` and realised random effects `u_i` (the sum of the realised
random effects on row `i`'s units):

    H_s = mean_i exp(beta[high] + u_i)
    L_s = mean_i exp(beta[low]  + u_i)

and the truth on each scale is formed from `H_s` and `L_s` exactly as
`pollin_effects()` forms its estimate:

| Scale | Basis | Truth |
|---|---|---|
| difference | loss | `H_s - L_s` |
| ratio | loss | `1 - L_s / H_s` |
| ratio | gain | `H_s / L_s - 1` |
| rii | loss | `(H_s - L_s) / (H_s + L_s)` |
| lnrr | loss | `log(H_s / L_s)` |

All three designs carry random intercepts only (no random treatment slope), so on
the ratio-type scales `H_s / L_s = exp(beta[high] - beta[low])` for every seed: the
ratio-type truths are constants. The difference truth varies by seed with the
realised random effects. The mean of the negative binomial is `exp(eta)`, so `H_s`
and `L_s` are expected responses. `basis` changes only the ratio scale; the other
scales run with `basis = "loss"` (recorded, value unaffected).

## 2. Designs

All responses are counts simulated from NB2 (`glmmTMB::nbinom2()`, variance
`mu + mu^2 / theta`), and every fit uses `family = glmmTMB::nbinom2()` passed
explicitly (the fitted family matches the simulating family). Treatment levels are
`Cage` (excluded), `Open` (open) and `Hand` (supplemented). Every unit receives
every arm; no data are missing. The simulating model is exactly the model the
layout implies, so a miss is a property of the estimator and its interval, not of
misspecification.

| | D1 three-arm trial | D2 split-plot trial | D3 matched-pair field study |
|---|---|---|---|
| `layout` | `rcbd` | `split_plot` | `matched_pairs` |
| `unit` | `block` | `c("block", "wplot")` | `c("site", "pair")` |
| Structure | 8 blocks x 3 arms x 4 plants; arms applied to individual plants within blocks | 6 blocks x 4 whole plots x 3 arm sub-plots, one count per sub-plot | 20 sites x 5 pairs x 2 arms, one plant per arm per pair |
| Rows | 96 | 72 | 200 |
| Random terms fitted | `(1 \| block)` | `(1 \| block) + (1 \| block:wplot)` | `(1 \| site) + (1 \| site:pair)` |
| Random-effect SDs (simulated) | block 0.40 | block 0.30, whole plot 0.25 | site 0.50, pair 0.30 |
| NB2 `theta` | 8 | 15 | 5 |
| Arms | excluded, open, supplemented | excluded, open, supplemented | excluded, open |
| Mean, excluded | 20 | 50 | 150 |
| Mean, open | 30 (E2 ratio 1.5) | 70 (E2 ratio 1.4) | 187.5 (E2 ratio 1.25) |
| Mean, supplemented | 36 (E1 ratio 1.2) | 77 (E1 ratio 1.1) | none |
| E1 (deficit) | applies | applies | does not apply (no supplemented arm) |
| E2 (dependence) | applies, reference open | applies, reference open | applies, reference open |

"Mean" is `exp(beta[arm])`, the expected count at zero random effect. The whole-plot
factor of D2 (for example variety) is represented only through the whole-plot random
effect; the S1 model carries no fixed whole-plot term.

**Minimum design size.** The sizes above are the minimum sizes for which C2 is
claimed in S1: `rcbd` 8 blocks with 4 plants per arm per block; `split_plot` 6 blocks
x 4 whole plots; `matched_pairs` 20 sites x 5 pairs. No smaller design is run in S1.
A cell that misses at these sizes is a FAIL, reported as such; it is not rerun at a
larger size under this rule.

## 3. Grid

| Design | Families | Scale x basis | Cells |
|---|---|---|---|
| D1 | E1, E2 | difference, ratio-loss, ratio-gain, rii, lnrr | 10 |
| D2 | E1, E2 | the same five | 10 |
| D3 | E2 | the same five | 5 |

25 cells, **1,000 seeds per cell**. For a given design and seed, one simulated data
set feeds every cell of that design. Each cell builds its own estimand and its own
fit through the public API (one `pollin_fit()` per family x scale x basis on the
same data); no fit is shared between scales, and no internal function is called.

## 4. Seeds and parallelism

Seed `s` in `1..1000` for design `k` (D1 = 1, D2 = 2, D3 = 3) simulates its data
after `RNGkind("Mersenne-Twister", "Inversion", "Rejection")` and
`set.seed(100000 * k + s)`, set inside the task. Tasks run with
`parallel::mclapply()`; because every task sets its own seed, results do not depend
on the core count or the scheduling. Nothing random happens after simulation (the
delta interval is deterministic).

## 5. Per-seed outcomes and failures

For each cell and seed: the estimate, the reported `std.error`, `conf.low`,
`conf.high`, the truth, and whether the truth lies inside the closed interval
`[conf.low, conf.high]` (the **coverage event**).

A seed **fails** in a cell when `pollin_fit()` or `pollin_effects()` errors, or the
estimate, `std.error`, `conf.low` or `conf.high` is not finite. Failed seeds are
excluded from coverage and bias and counted (`n_failed`). Warnings are recorded and
counted (`n_warned`) but do not make a seed fail; the advisory `marginaleffects`
issues on every `glmmTMB` call ("only takes into account the uncertainty in
fixed-effect parameters") is excluded from that count, since it is raised on every
delta call by design. A cell with more than 2% failed seeds (more than 20) cannot
PASS; it is reported as FAIL (fits) whatever its coverage.

## 6. Pass rule (the plan's C2, made exact)

Over the `n_ok` successful seeds of a cell:

* **coverage** = share of seeds whose coverage event holds. Reported with its Wilson
  95% interval.
* **bias** = `mean(estimate - truth)` (truth per seed).
* **SE** for the bias criterion = the **Monte-Carlo standard deviation of the
  estimates** across seeds (the empirical SE). The mean reported `std.error` is
  reported alongside, not used in the rule.

A cell **PASSES** when all three hold:

1. coverage in `[0.93, 0.97]` (inclusive, the point estimate; the Wilson interval
   is reported, not used to pass);
2. `|bias| / SE < 0.1`;
3. `n_failed <= 20`.

Otherwise it FAILS, with the criterion it missed named.

Known in advance and stated now, not used to excuse a miss: the delta interval
carries the uncertainty in the fixed effects only; the estimated random effects
enter the averaged predictions as if known. On the ratio-type scales the random
intercepts cancel from `H / L`, so this matters only on the difference scale, where
it may move coverage in either direction. Coverage on ratio, rii and lnrr is the
same event within a family (each interval is a monotone transform of the `log(rr)`
interval and the truth transforms with it); their biases differ.

## 7. Outputs

`validation/planted_truth_S1.R` (the suite), `validation/planted_truth_S1_results.csv`
(one row per design x family x scale x basis: `n_ok`, `n_failed`, `n_warned`,
coverage, Wilson 95% interval, bias, empirical SE, bias/SE, mean reported SE, mean
interval width, verdict), `validation/planted_truth_S1_session.txt` (session info,
package versions, wall time) and `validation/planted_truth_S1_report.md` (each cell
PASS or FAIL against this rule, misses explained with the evidence).

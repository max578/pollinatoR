# Planted-truth recovery rule, Stage 3 (criteria C2, C3, C4)

Fixed 2026-11-09, before any simulation run of the suite whose results are reported.
This file is committed on its own; the commit time is the evidence that the rule
preceded the results. Nothing here changes after the run. A changed rule is a new
dated amendment next to this one, with the reason, written before the run it
governs, and both sets of results are reported.

**What was run before this rule was written.** The suite script was run on seed 1 of
every data design into a scratch directory, four times (the first three stopped on a script error,
the flexyBayes aggregation finding below, and an edit made mid-run), to confirm that the public
API calls run on both engines and to time them. Only the error, warning and timing
columns were read; no agreement, coverage, width or RMSE was computed from them or
looked at. Those checks found three things that shaped the rule:

1. On Poisson and binomial responses with factor predictors only, flexyBayes
   compresses the rows into cell sufficient statistics (`aggregate = "auto"`), and
   the posterior draws then fail with "subscript out of bounds" (Section 9). Every
   flexyBayes fit in the suite therefore passes `aggregate = FALSE` through
   `pollin_fit(...)`, a documented `flexybayes()` argument.
2. glmmTMB stopped with "singular convergence (7)" on seed 1 of the Gamma slope
   design and of the eight-level variety design, where the planted site-by-arm SD
   was 0.10 and 0.15. The slope-cell SD was set to 0.20 in every slope design
   before this rule was written, a size of between-site and between-variety
   variation in the pollination effect that is plausible for canola and keeps the
   slope variance away from zero. Boundary fits still count as failures
   (Section 6).
3. The C4 designs were reduced from four to three plants per block, variety and arm,
   so that the per-level contrasts are noisy enough for pooling to matter and the
   run fits the time budget.

## 1. What is recovered

Every cell builds its estimand with the public API: `pollin_data()` ->
`pollin_design()` -> `pollin_estimand()` -> `pollin_fit()` on each engine ->
`pollin_effects()` (overall or per level) or `pollin_cate()` (pooled per level).
No internal function is called.

**Truth.** As in S1 and S2 (`validation/recovery_rule_S1.md`, Section 1): both
engines form the arm means by G-computation, every row of the model data predicted
under the `high` arm and under the `low` arm with its own random effects, and the
predictions averaged over the rows. The target is the arms' expected responses
**averaged over the realised rows, conditional on the realised random effects**.
For simulated data set `s` with rows `i`, fixed linear predictors `beta[arm]`,
realised random effects `u_i` and, in slope designs, the realised unit-by-arm effect
`c[unit(i), arm]`, with inverse link `g`:

    H_s = mean_i g(beta[high] + u_i + c[unit(i), high])
    L_s = mean_i g(beta[low]  + u_i + c[unit(i), low])

| Scale | Basis | Truth |
|---|---|---|
| ratio | loss | `1 - L_s / H_s` |
| ratio | gain | `H_s / L_s - 1` |
| difference | loss | `H_s - L_s` |

The per-level truth of C4 is the same, with the means taken over the rows of the
level. Nothing about the truth depends on which engine is fitted.

**The two engines' random parts.** glmmTMB fits a random treatment slope as
`(trt | unit)`, an unstructured 2 x 2 covariance. The flexyBayes route on INLA fits
`unit + unit:trt`: a unit intercept plus an independent effect per unit-by-arm cell
(`validation/flexybayes_probe_S3.md`). Every slope design is therefore **planted
from the restricted form**: a unit intercept with SD `sd_unit` plus an independent
unit-by-arm effect with SD `sd_cell`. In the glmmTMB parameterisation that is an
intercept variance `sd_unit^2 + sd_cell^2`, a slope variance `2 sd_cell^2` and a
covariance `-sd_cell^2`, which the unstructured form contains. Both engines can
therefore represent the simulating model exactly, and a disagreement or a miss is a
property of the estimator and its interval, not of misspecification.

**The intervals compared.** glmmTMB: the delta interval of `pollin_effects()`
(fixed-effect uncertainty only, the conditional modes held as known), and the
approximate conditional interval of `pollin_cate()`. flexyBayes: the equal-tailed
95% posterior interval and the posterior median, from 1,000 posterior draws.

## 2. Designs

Treatment levels are `Cage` (excluded), `Open` (open) and `Hand` (supplemented).
Every unit receives every arm; no data are missing. "Mean" is `g(beta[arm])`, the
expected response at zero random effects.

| Data | Layout (`unit`) | Structure | Rows | Response, family | Random part (SDs, link scale) | Means |
|---|---|---|---|---|---|---|
| A1 | `multi_site` (`site`, slope on site) | 12 sites x 2 arms x 8 plants | 192 | count, NB2 `theta = 8` | site 0.40, site-by-arm 0.20 | Cage 20, Open 30 |
| A2 | `rcbd` (`block`) | 8 blocks x 2 arms x 6 plants | 96 | count, Poisson | block 0.40 | Cage 20, Open 30 |
| A3 | `multi_site` (`site`, slope on site) | 12 sites x 2 arms x 8 plants | 192 | positive continuous (seed weight per plant, g), Gamma, shape 10, log link | site 0.30, site-by-arm 0.20 | Cage 2.0, Open 2.6 |
| A4 | `rcbd` (`block`) | 8 blocks x 2 arms x 24 flowers | 384 | pod set per flower (0/1), binomial, logit link | block 0.40 | Cage 0.55, Open 0.75 |
| A5 | `rcbd` (`block`) | 8 blocks x 3 arms x 4 plants | 96 | count, NB2 `theta = 8` | block 0.40 | Cage 20, Open 30, Hand 36 |
| S8 | `rcbd` (`block`), `pooled = "variety"` | 4 blocks x 8 varieties x 2 arms x 3 plants | 192 | count, NB2 `theta = 8` | block 0.30, variety 0.30, variety-by-arm 0.20 | Cage 20, Open 30 |
| S12 | as S8 | 4 blocks x 12 varieties x 2 arms x 3 plants | 288 | as S8 | as S8 | as S8 |

Every fit passes its family explicitly: `glmmTMB::nbinom2()`, `stats::poisson()`,
`stats::Gamma(link = "log")` or `stats::binomial()`, the same object on both
engines.

**The binomial denominator.** `pollin_data(denominator = )` records what a response
is counted out of, but neither engine's fit uses it: `pollin_fit()` passes no
weights or trials to glmmTMB, and flexyBayes takes no trials argument (its
`weights` is refused on every family but the Gaussian). A seed-set response out of
`n` ovules with a size column therefore cannot be fitted as binomial through the
public API on either engine (Section 9). The binomial cell uses the form both
engines take: one row per flower, the response 0/1 (pod set or not), Bernoulli. Its
truth is the arms' mean pod-set probability.

## 3. C2 and C3 cells

| Cell | Data | Family | Channel | Scale, basis | Slope in the fit | Seeds |
|---|---|---|---|---|---|---|
| A1 | A1 | NB2 | dependence | ratio, loss | site | 200 |
| A2 | A2 | Poisson | dependence | ratio, loss | none | 200 |
| A3 | A3 | Gamma | dependence | ratio, loss | site | 200 |
| A4 | A4 | binomial | dependence | ratio, loss | none | 200 |
| A5 | A5 | NB2 | deficit | ratio, gain | none | 200 |
| A6 | A1 (same data sets as A1) | NB2 | dependence | difference, loss | site | 200 |
| P8 | S8 | NB2 | dependence | ratio, loss | `pooled = "variety"` | 100 |
| P12 | S12 | NB2 | dependence | ratio, loss | `pooled = "variety"` | 100 |

P8 and P12 read the **overall** estimand (`pollin_effects()`, `by = NULL`) from the
pooled fits that C4 makes on the S8 and S12 data, so they cost no extra fits and run
on the C4 seeds.

**Why the width ratio is not gated on slope cells.** A1, A3, A6, P8 and P12 carry a
random treatment slope. There the two engines target the same trial-conditional
estimand, but the glmmTMB delta interval varies the fixed treatment coefficient
while holding the conditional modes of the slope at their estimates. Under the
conditional estimand the mean of the unit-specific effects is what the data
determine; the fixed coefficient and the modes are negatively correlated given the
data, and holding the modes fixed drops that correlation, so the delta interval is
expected to be wider than the posterior interval, which keeps it. One fit seen
during development (not part of this suite) gave a loss share of 0.327 with delta
SE 0.063 and posterior SD 0.026, a width ratio near 0.4. A width ratio outside
`[0.8, 1.25]` on a slope cell would then reflect the delta interval's treatment of
the random effects, not a disagreement about the estimand. On slope cells the width
ratio is therefore **reported, not gated**, and **both engines' coverage of the
planted truth is reported** side by side, so the data show which interval is
calibrated for this estimand. The agreement criterion (posterior median within
0.25 glmmTMB SE) is gated on every cell. On cells without a slope (A2, A4, A5) the
width ratio is gated as the plan states. This narrows the plan's C2 wording on
slope cells; it is stated here, before any result, not read into one afterwards.

## 4. C4 cells

| Cell | Data | Levels | Pooled fit | Unpooled fit | Seeds |
|---|---|---|---|---|---|
| S8 | S8 | 8 varieties | `pollin_fit(pooled = "variety")`, read with `pollin_cate()` | `pollin_fit(by = "variety")` (fixed `trt * variety`, no slope), read with `pollin_effects()` | 100 |
| S12 | S12 | 12 varieties | as S8 | as S8 | 100 |

Both fits are made on both engines on the same data. On glmmTMB the unpooled
interval is the delta interval; on flexyBayes it is the posterior interval. The
per-level truth is the realised loss share of each variety (Section 1).

## 5. Thresholds (the plan's C2, C3, C4, made exact)

**C2, engine agreement.** Over the `n_ok` seeds of a cell on which both engines
succeed:

* **agreement share** = share of seeds with
  `|flexyBayes posterior median - glmmTMB estimate| <= 0.25 x glmmTMB SE`, the SE
  being the `std.error` of the delta row;
* **width ratio** = the median over seeds of
  `(flexyBayes interval width) / (glmmTMB delta interval width)`.

A cell PASSES C2 when the agreement share is `>= 0.95`, the width ratio lies in
`[0.8, 1.25]` (gated on A2, A4, A5 only; Section 3), and its failures are within the
limit of Section 6. **C2 is met** when all eight cells PASS.

**C3, Bayesian calibration.** Over the same seeds, coverage = share of seeds whose
planted truth lies inside the closed flexyBayes 95% posterior interval. A cell
PASSES C3 when coverage lies in `[0.93, 0.97]` and its failures are within the
limit. **C3 is met** when all eight cells PASS. The glmmTMB delta coverage of the
same truth is reported in every cell, and set beside the flexyBayes coverage on the
slope cells; it is not gated.

**C4, shrinkage recovery.** Per seed and engine, `RMSE = sqrt(mean over levels of
(estimate - truth)^2)` for the pooled estimates (`pollin_cate()`: the conditional
estimate on glmmTMB, the posterior median on flexyBayes) and for the unpooled ones
(`pollin_effects(by = "variety")`). A seed counts for an engine when both of that
engine's fits succeed. A cell PASSES C4 when

1. the pooled RMSE is strictly below the unpooled RMSE in `>= 90%` of seeds, on
   **each** engine separately;
2. the flexyBayes pooled per-level coverage, the share of (seed, level) pairs whose
   truth lies inside the posterior interval of `pollin_cate()`, lies in
   `[0.92, 0.98]`;
3. failures are within the limit on each of the four fits.

**C4 is met** when S8 and S12 both PASS. The glmmTMB pooled per-level coverage (its
interval is stated as approximate), and both engines' unpooled per-level coverage,
are reported, not gated.

**Monte Carlo error and borderline results.** The seed counts (Section 7) give a
coverage near 0.95 a Monte Carlo SE of about 0.015 at 200 seeds and 0.022 at 100;
the per-level coverage of C4 pools 800 and 1,200 (seed, level) pairs, a nominal SE
of about 0.008 and 0.006, which understates the true error because levels within a
seed share a fit. Each coverage is reported with its Monte Carlo SE and Wilson 95%
interval. Two times the SE at 200 seeds (0.031) is wider than the distance from 0.95
to either edge of the band (0.02), so a reading of "within 2 SE of an edge" that
covered values inside the band would leave no coverage able to pass. The rule is
therefore: a coverage **inside** the band is PASS; a coverage **outside** the band
by no more than 2 Monte Carlo SE is **BORDERLINE**, which is not a pass and is
reported as a miss of the criterion, with its size; anything further out is FAIL.
The same SE also means that a coverage inside the band is consistent with true
coverage some way either side of it; the Wilson interval is reported so that a
reader can see how far.

## 6. Failures and non-convergence

A seed **fails** for an engine in a cell when any call errors; when an estimate,
SE, interval end, or (for C4) any level's estimate or interval end is not finite;
when the glmmTMB fit reports a non-zero optimiser code or a Hessian that is not
positive definite (including the "singular convergence" stop that a slope variance
at zero produces); or when flexyBayes's convergence verdict is `FALSE` (an error in
`pollin_fit()`) or `NA` (no verdict, counted as not converged). Failed seeds are
counted per engine (`n_failed_*`) and excluded from every metric of that cell (in
C2 and C3, a seed counts only when both engines succeed). Warnings that do not stop
a fit are counted (`n_warned_*`) and do not fail the seed; the `marginaleffects`
advisory about fixed-effect uncertainty, raised on every delta call, is excluded.

The **ceiling** is 2% of the seeds of a cell per fit: at most 4 of 200 and 2 of 100.
A cell with more failures in any one of its fits is FAIL (fits), whatever its other
metrics; the metrics are still reported.

## 7. Seeds, draws and parallelism

**Counts.** 200 seeds per C2/C3 cell on A1-A5 data (A6 uses the A1 data sets), and
100 seeds per C4 cell (P8 and P12 use the same 100). The plan's S1 and S2 suites ran
1,000 seeds per cell on glmmTMB alone; an INLA fit with 1,000 posterior draws costs
about 15 to 40 seconds of one core here, against well under a second for glmmTMB, so
1,000 seeds per cell is beyond the compute budget of the run (about 90 minutes on
8 cores). The counts above are the budget's. If the timing of the smoke run
(below) shows the full run would exceed about 90 minutes, the counts are reduced by
a dated amendment written and committed before the full run starts, never after any
result of it is seen.

**Draws.** 1,000 posterior draws per flexyBayes call (`pollin_effects(n_draws = )`,
`pollin_cate(n = )`), and 1,000 draws for the glmmTMB interval of `pollin_cate()`.

**Seeds.** Seed `s` of data design `k` (A1 = 1, ..., A5 = 5, S8 = 7, S12 = 8)
simulates its data after `RNGkind("Mersenne-Twister", "Inversion", "Rejection")` and
`set.seed(3000000 + 100000 k + s)`, set inside the task (a range disjoint from S1 and
S2). Every draw call passes `seed = 23000000 + 100000 k + s`; the package sets R's
and INLA's generators from it and restores the session's state afterwards
(`validation/flexybayes_probe_S3.md`, Seeds). INLA runs with
`num.threads = "1:1"` (`INLA::inla.setOption()`) and
`options(flexyBayes.inla_best_of = 1L)` (one optimisation per fit, not the best of
three in fresh processes). The flexyBayes priors are the package defaults.

**Determinism.** Tasks run with `parallel::mclapply()`, one task per data set and
seed, each setting its own seeds. Before the full run, a smoke run of 3 seeds per
cell on 1 core and on 4 cores into scratch directories must give identical per-seed
outcomes; if it does not, the full run does not start and the cause is reported.
The smoke run's agreement, coverage and RMSE are not read.

## 8. Stated in advance

These are written now so that they cannot be read into a result afterwards. None is
used to excuse a miss.

1. **The delta interval on slope cells** is expected to be wider than the posterior
   (Section 3), so glmmTMB may over-cover there and the width ratio may lie well
   below 0.8. If, instead, the posterior under-covers, that is a C3 miss.
2. **The posterior depends on the priors.** With 8 to 12 groups the posterior of a
   random-effect variance is not dominated by the data, and the flexyBayes default
   priors shape the interval. Coverage of the planted truth is a frequentist
   property of a Bayesian interval; it is the plan's C3 and is checked as such.
3. **Coverage of a conditional estimand.** The posterior interval of the realised
   loss share carries the uncertainty in the realised random effects, which the
   delta interval does not; on intercept-only ratio cells the random intercepts
   cancel from `H / L` only under a log link, so A4 (logit link) is the
   intercept-only cell where the two intervals can differ most.
4. **Boundary fits.** A slope variance estimated at zero fails a glmmTMB seed
   (Section 6). With 8 to 12 groups this can happen on more than 2% of seeds; the
   cell is then FAIL (fits), reported as a property of glmmTMB on that design.
5. **The C4 shrinkage criterion is per seed.** Pooling lowers the expected RMSE
   whenever the slope variance is estimated well, but with 8 levels the per-seed
   comparison is noisy; 90% is a demanding share and S8 is the harder cell.
6. **Agreement tolerance.** On slope cells the glmmTMB SE is the wider of the two
   (item 1), so the agreement tolerance of 0.25 SE is looser there than on
   intercept-only cells.

**On a miss.** A cell that misses is reported as FAIL or BORDERLINE with the
criterion, the size of the miss and the evidence on its cause from the per-seed
outcomes. It is not rerun with more draws, more seeds, other priors, a larger design
or another method under this rule, and no package code is changed to move it.

## 9. Package findings recorded before the run

1. **flexyBayes aggregation breaks the posterior draws.** `R/fit_flexybayes.R`
   (`.flexy_new_fit()`) does not set `aggregate`, so `flexybayes()` takes its default
   `"auto"`, which on Poisson and binomial fits with factor predictors only
   compresses the rows into cells; the fit holds one `Predictor` row per cell, and
   `.flexy_arm_draws()` indexes `Predictor:1..n` by data row and stops with
   "subscript out of bounds". Reproducer: eight blocks x two arms x six plants,
   Poisson counts, `pollin_fit(est, engine = "flexyBayes", family = poisson())`,
   then `pollin_effects(fit)`; the posterior sample holds 16 `Predictor` rows for
   96 data rows. Passing `aggregate = FALSE` to `pollin_fit()` fits and draws.
2. **The binomial denominator is not used by any fit** (Section 2).

Neither is changed by this suite; both are reported.

## 10. Outputs

`validation/planted_truth_S3.R` (the suite), `validation/planted_truth_S3_results.csv`
(one row per cell: counts of seeds, failures and warnings per engine; for C2/C3
cells the agreement share, the median `|difference| / SE`, the width ratio and
whether it is gated, both verdicts, the flexyBayes coverage with Monte Carlo SE and
Wilson interval, the glmmTMB coverage, both engines' bias and median seconds per
seed; for C4 cells the win share and mean RMSEs per engine, the per-level coverages
and the verdicts), `validation/planted_truth_S3_seeds.csv.gz` (per seed, and per
level for C4), `validation/planted_truth_S3_session.txt` (session, versions, wall
time, error and warning tables) and `validation/planted_truth_S3_report.md` (a
table per criterion, each cell's verdict, misses explained with the evidence).

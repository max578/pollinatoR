# Planted-truth recovery rule, Stage 3.1 (criteria C2 and C3 for 0.3.1)

Fixed 2026-11-09, before any simulation run of the suite whose results are reported.
This file is committed on its own; the commit time is the evidence that the rule
preceded the results. Nothing here changes after the run. A changed rule is a new
dated amendment next to this one, with the reason, written before the run it
governs, and both sets of results are reported.

It re-tests what changed between pollinatoR 0.3.0 and 0.3.1 (development, branch
`s3-followup`): `pollin_data(denominator = )` now reaches the fit as a binomial
response of successes out of trials, on both engines. Everything not restated here
is as in `validation/recovery_rule_S3.md`, including its Amendment 1.

**What was run before this rule was written.**

1. The exploratory probe of the S3 slope-coverage shortfall,
   `validation/probe_S3.1_slope_coverage.md`, on its own seeds (Section 5).
2. The suite script `validation/planted_truth_S3.1.R` on seed 1 of each of B1, A2
   and A5, on 3 cores, into a scratch directory, to confirm that the public API runs
   on both engines with a denominator and to time the fits. Only the error, warning
   and timing columns were read: no fit failed, none warned, and the flexyBayes fits
   took about 10 s (B1, 2,793 trial rows), 9 s (A2) and 21 s (A5) of one core. No
   agreement, coverage or width was computed or looked at.

## 1. What is recovered, and what is not re-run

**Recovered.** As in S3, Section 1: every cell builds its estimand with the public
API (`pollin_data()` -> `pollin_design()` -> `pollin_estimand()` -> `pollin_fit()` on
each engine -> `pollin_effects()`), and the truth is the arms' expected responses
averaged over the realised rows, conditional on the realised random effects. With a
denominator, the arm mean is the expected proportion (successes / trials) averaged
over the rows with equal weight, whatever each row's number of trials, which is how
both engines form it (`R/fit.R`, `R/fit_flexybayes.R`). For data set `s` with plants
`i`, block effect `u_i` and fixed logit-scale predictors `beta[arm]`:

    H_s = mean_i plogis(beta[Open] + u_i),   L_s = mean_i plogis(beta[Cage] + u_i)

and the dependence loss share is `1 - L_s / H_s`. The ovule counts do not enter the
truth; they set how much the data say about it.

**The slope cells are not re-run.** A1, A3, A6, P8 and P12 (and the C4 cells S8 and
S12) are not part of this rule. The S3.1 probe found no setting to change on the
site-slope design: integration strategy, hyperparameters at the mode and three
variance-component priors left every coverage event unchanged on 100 paired seeds,
and an exact posterior (NUTS with the hyperparameters fixed at their planted values)
covered as INLA does (0.930 against 0.930, width ratio 0.992). The only package
change since S3 (`v0.3.0`) is the denominator route; a fit without a denominator
takes the same path as before (A2 and A5 check that), so a re-run of the slope
cells would test the same code with the same settings on new seeds. The label decision for slope designs (whether the
flexyBayes engine stays experimental there) is Max's, made on the S3 numbers plus
the probe, and is outside this rule.

**A4 is replaced by B1.** In S3 the binomial cell had to be fitted as one 0/1 row per
flower, because a denominator reached no fit. B1 is the same kind of design fitted
the way 0.3.1 now supports: seed counts out of a declared number of ovules.

## 2. Designs

| Data | Layout (`unit`) | Structure | Rows | Response, family | Random part (SD, link scale) | Means |
|---|---|---|---|---|---|---|
| B1 | `rcbd` (`block`) | 8 blocks x 2 arms x 6 plants | 96 | seeds out of ovules per plant, binomial, logit link; `pollin_data(denominator = "n")`; ovules drawn uniformly on 20 to 40 per plant | block 0.40 | Cage 0.55, Open 0.75 (probability a given ovule sets seed) |
| A2 | as S3 | 8 blocks x 2 arms x 6 plants | 96 | count, Poisson | block 0.40 | Cage 20, Open 30 |
| A5 | as S3 | 8 blocks x 3 arms x 4 plants | 96 | count, NB2 `theta = 8` | block 0.40 | Cage 20, Open 30, Hand 36 |

B1 follows A4's design (8 blocks, block SD 0.40, means 0.55 and 0.75 on the logit
scale, no slope), with plants as rows and a seed-set count out of each plant's
ovules in place of a 0/1 row per flower. Within a data set, the block effects are
drawn first, then the ovule counts, then the seed counts. Every fit passes its
family explicitly, the same object on both engines: `stats::binomial()`,
`stats::poisson()` or `glmmTMB::nbinom2()`. On glmmTMB, B1 is fitted as
`cbind(y, n - y)`; on flexyBayes as one 0/1 row per ovule (about 2,900 rows), with
the arm means rebuilt on the 96 plant rows. Every flexyBayes fit passes
`aggregate = FALSE`, as in S3.

## 3. Cells

| Cell | Data | Family | Channel | Scale, basis | Slope | Seeds | Width ratio gated |
|---|---|---|---|---|---|---|---|
| B1 | B1 | binomial (with denominator) | dependence | ratio, loss | none | 200 | yes |
| A2 | A2 | Poisson | dependence | ratio, loss | none | 100 | yes |
| A5 | A5 | NB2 | deficit | ratio, gain | none | 100 | yes |

A2 and A5 are controls: the same designs and cells as in S3, on new seeds, to show
that the 0.3.1 changes did not move a result that passed.

## 4. Thresholds (unchanged from S3)

**C2, engine agreement.** Over the `n_ok` seeds of a cell on which both engines
succeed: the **agreement share** (share of seeds with
`|flexyBayes posterior median - glmmTMB estimate| <= 0.25 x glmmTMB SE`, the SE
being the `std.error` of the delta row) must be `>= 0.95`, and the **width ratio**
(median over seeds of flexyBayes interval width / glmmTMB delta interval width) must
lie in `[0.8, 1.25]`, on every cell here (none has a slope). Failures must be within
Section 6's limit. A cell meeting all three PASSES C2.

**C3, Bayesian calibration.** Coverage = share of the `n_ok` seeds whose planted
truth lies inside the closed flexyBayes 95% posterior interval. A cell PASSES C3 when
coverage lies in `[0.93, 0.97]` and its failures are within the limit. The glmmTMB
delta coverage of the same truth is reported beside it, with its Monte Carlo SE, and
not gated.

**Borderline.** As S3: inside the band is PASS; outside it by no more than 2 Monte
Carlo SE is BORDERLINE, which is a miss reported with its size, not a pass; further
out is FAIL. Monte Carlo SE near 0.95 is about 0.015 at 200 seeds and 0.022 at 100.
Each coverage is reported with its Monte Carlo SE and Wilson 95% interval.

**Also reported, not gated:** misses below and above the interval; sd(z) and mean(z)
with z = (posterior median - truth) / posterior SD (the probe's finer measure of
width and location); bias on both engines; the interquartile range of the width
ratio; median seconds per fit; the median number of trials of B1.

**Verdict for 0.3.1.** B1 passing both C2 and C3, and A2 and A5 each passing both,
is the evidence the plan asks for the binomial-denominator route. A control cell that
fails or is BORDERLINE where it passed in S3 is reported as a regression, with the
evidence on its cause.

## 5. Seeds, draws and parallelism

**Seeds.** Seed `s` of data design `k` (A2 = 2, A5 = 5, B1 = 9) simulates its data
after `RNGkind("Mersenne-Twister", "Inversion", "Rejection")` and
`set.seed(41000000 + 100000 k + s)`, set inside the task. Every posterior draw call
passes `seed = 61000000 + 100000 k + s`. s runs from 1 to 200 (B1) or 1 to 100 (A2,
A5). These ranges are disjoint from every earlier seed:

| Suite | Data seeds | Draw seeds |
|---|---|---|
| S1 | 100,000 x id + s (id 1 to 3) | |
| S2 | 2,000,000 + 100,000 x id + s (id 1 to 4) | 1,000,000 x id + s; 5,000,000 + s |
| S3 | 3,000,000 + 100,000 k + s (k 1 to 8) | 23,000,000 + 100,000 k + s |
| S3.1 probe | 900,000 + s; 950,000 + s | 990,000 + s |
| **S3.1** | **41,000,000 + 100,000 k + s** | **61,000,000 + 100,000 k + s** |

**Draws.** 1,000 posterior draws per flexyBayes call
(`pollin_effects(n_draws = 1000)`). INLA runs with `num.threads = "1:1"` and
`options(flexyBayes.inla_best_of = 1)`. The flexyBayes priors are the package
defaults.

**Budget.** About 45 minutes. The single-seed timings put the full run (400 tasks)
at about 11 minutes of 8 cores without other load. If the smoke run below shows the
full run would exceed 45 minutes, the seed counts are cut by a dated amendment
written and committed before the full run, never after any result of it is seen.

**Determinism (as S3 Amendment 1).** Tasks run with `parallel::mclapply()`, one task
per data set and seed, each setting its own seeds. Before the full run, a smoke run
of 3 seeds per cell on 1 core and on 4 cores into scratch directories is compared.
The simulated data, planted truths, every glmmTMB result and every error column
must be bit-identical between the two; if they are not, the full run does not start
and the cause is reported. The flexyBayes results are required to agree only to the
run-to-run variation of the INLA fit, which is not bit-reproducible across processes
(S3 Amendment 1); the largest difference in posterior medians and interval ends is
reported in posterior SD, with whether any coverage or agreement event changed. Only
those differences are read from the smoke run, not its agreement, coverage or width.
In the full run, seeds whose coverage or agreement event lies within 0.15 posterior
SD of the deciding interval end or tolerance are counted and reported.

## 6. Failures and non-convergence (as S3, Section 6)

A seed fails for an engine when any call errors; when an estimate, SE or interval end
is not finite; when the glmmTMB fit reports a non-zero optimiser code or a Hessian
that is not positive definite; or when the flexyBayes convergence verdict is `FALSE`
(an error in `pollin_fit()`) or `NA`. Failed seeds are counted per engine and
excluded from the metrics (a seed counts only when both engines succeed). Warnings
that do not stop a fit are counted and do not fail the seed; the `marginaleffects`
advisory about fixed-effect uncertainty is excluded. The ceiling is 2% of a cell's
seeds per engine (4 of 200, 2 of 100); a cell over it is FAIL (fits), whatever its
other metrics, which are still reported. No change to how comparator failures count
is made here; none of these cells has a slope, where the S3 comparator failures
arose.

## 7. Stated in advance

1. **Trials give the binomial cell much more data than A4 had.** About 2,900 ovules
   over 96 plants against 384 flowers, so the arm-level proportions are tightly
   known and the interval is driven by the 8 block effects, which do not cancel from
   `H / L` under a logit link (S3, Section 8, item 3). The posterior interval carries
   the block-effect uncertainty of the realised blocks; the delta interval holds the
   conditional modes as known. The two can therefore differ in width on B1 more than
   on A2 and A5, and the width ratio is gated as the plan states.
2. **No overdispersion is planted.** Seeds are binomial given the plant's
   probability, and the fit is binomial; the design tests the denominator route, not
   robustness to extra-binomial variation.
3. **INLA run-to-run variation** can flip events near a boundary (Section 5); the
   count of such seeds is reported, and a verdict that rests on fewer events than
   that count is said to be not robust to a rerun.

**On a miss.** A cell that misses is reported as FAIL or BORDERLINE with the
criterion, the size of the miss and the evidence on its cause from the per-seed
outcomes. It is not rerun with more draws, more seeds, other priors, a larger design
or another method under this rule, and no package code is changed to move it.

## 8. Outputs

`validation/planted_truth_S3.1.R` (the suite), `validation/planted_truth_S3.1_results.csv`
(one row per cell: seed, failure and warning counts per engine; agreement share,
median `|difference| / SE`, width ratio and its interquartile range, C2 verdict;
flexyBayes coverage with Monte Carlo SE and Wilson interval, misses below and above,
sd(z), mean(z), C3 verdict on coverage alone and with the fits limit; glmmTMB
coverage and its Monte Carlo SE; both engines' bias; counts of seeds near a deciding
boundary; median seconds per fit; median trials),
`validation/planted_truth_S3.1_seeds.csv.gz` (per seed),
`validation/planted_truth_S3.1_session.txt` (session, versions, git HEAD, wall time,
error and warning tables) and `validation/planted_truth_S3.1_report.md` (tables,
verdicts, misses explained with the evidence).

# Planted-truth recovery, Stage 3 (C2 engine agreement, C3 Bayesian calibration, C4 shrinkage): report

Rule: `validation/recovery_rule_S3.md`, committed alone at `efa6aa7`
(2026-11-09 18:21:40 +0800), with Amendment 1 committed at `1153cde`
(2026-11-09 18:42:16 +0800). Both commits came before the full run. The amendment
replaced only the determinism condition of the smoke run (INLA fits are not
bit-reproducible across processes; see Determinism below). Run: 2026-11-09 18:44:03 to
20:45:38 (AWST), **121.6 min wall time** on 8 cores, 1,200 data set x seed tasks.
That is over both the rule's estimate (about 75 min) and the 90-minute budget, with
other jobs running on the machine during the run. The seed counts were not changed.
Suite: `validation/planted_truth_S3.R`. Per-seed outcomes:
`validation/planted_truth_S3_seeds.csv.gz`. Per-cell summary:
`validation/planted_truth_S3_results.csv`. Session:
`validation/planted_truth_S3_session.txt`.

## Verdict

| Criterion | Verdict | In one line |
|---|---|---|
| C2 engine agreement | **Not met.** 4 of 8 cells PASS | Every cell meets the agreement share (0.973 to 1.000), and the gated width ratios are inside the band. A1, A3, A6 and P8 FAIL on fits: glmmTMB stops at a slope variance on the boundary in 4% to 12% of seeds |
| C3 Bayesian calibration | **Not met.** 3 of 8 cells PASS, 1 BORDERLINE, 4 FAIL | The posterior under-covers on the site-slope cells (0.909 to 0.925) and on the binomial cell (0.925). The four FAILs are fits failures (glmmTMB), and on coverage alone three of them would be BORDERLINE and one (P8) a PASS |
| C4 shrinkage recovery | **Not met.** S8 and S12 both FAIL | Pooling lowers the mean RMSE on both engines by 12% to 15%, but the pooled RMSE beats the unpooled in only 69% to 79% of seeds, against the 90% gate. flexyBayes per-level coverage PASSES (0.941, 0.940) |

No flexyBayes fit failed or lost its convergence verdict in any cell. Every failure
is a glmmTMB fit stopped by the rule's convergence check: 14 seeds of the A1 data
(failing both A1 and A6, which share the fit's data), 8 of A3, 12 of S8 and 1 of
S12 (the pooled fit, failing both the P cell and the C4 cell). Over the 49 failed
C2/C3 seed-cells the messages were "singular convergence (7)" in 39,
"non-positive-definite Hessian" in 9 and a Hessian fallback in 1.

## C2: engine agreement

Agreement share = the share of seeds whose flexyBayes posterior median lies within
0.25 glmmTMB delta SE of the glmmTMB estimate (gate `>= 0.95`). Width ratio = the
median over seeds of the flexyBayes interval width divided by the glmmTMB delta
width (gate `[0.8, 1.25]`, applied to A2, A4 and A5 only). Failure limit: 4 of 200
seeds, or 2 of 100.

| Cell | Family | Scale | Slope | Seeds ok | glmmTMB failed | flexyBayes failed | Agreement share | Median abs diff / SE | Width ratio (gated) | Verdict |
|---|---|---|---|---|---|---|---|---|---|---|
| A1 | NB2 | ratio, loss | site | 186 / 200 | **14** | 0 | 0.989 | 0.057 | 0.628 (no) | **FAIL (fits)** |
| A2 | Poisson | ratio, loss | none | 200 / 200 | 0 | 0 | 1.000 | 0.027 | 0.996 (yes) | PASS |
| A3 | Gamma | ratio, loss | site | 192 / 200 | **8** | 0 | 0.984 | 0.055 | 0.556 (no) | **FAIL (fits)** |
| A4 | binomial | ratio, loss | none | 200 / 200 | 0 | 0 | 1.000 | 0.036 | 0.934 (yes) | PASS |
| A5 | NB2 | ratio, gain (deficit) | none | 200 / 200 | 0 | 0 | 1.000 | 0.027 | 1.005 (yes) | PASS |
| A6 | NB2 | difference | site | 186 / 200 | **14** | 0 | 0.973 | 0.058 | 0.593 (no) | **FAIL (fits)** |
| P8 | NB2 | ratio, loss | `pooled = "variety"`, 8 levels | 88 / 100 | **12** | 0 | 1.000 | 0.036 | 0.578 (no) | **FAIL (fits)** |
| P12 | NB2 | ratio, loss | `pooled = "variety"`, 12 levels | 99 / 100 | 1 | 0 | 0.990 | 0.049 | 0.532 (no) | PASS |

On every cell the engines agree on the estimate. The median gap is 0.03 to 0.06
glmmTMB SE, and between 1 and 32 seeds per cell lie within the INLA run-to-run
variation of the 0.25 SE tolerance (`n_near_agreement` in the results), against
4 to 10 seeds allowed to miss at the 95% gate. On cells without a slope the
intervals agree as well: width ratios of 0.93 to 1.01, interquartile 0.89 to 1.02.
On every slope cell the posterior interval is 0.53 to 0.63 of the delta width
(interquartile 0.48 to 0.72), and the delta SE is 1.6 to 1.9 times the posterior SD.
The rule anticipated this (Section 3) and did not gate it. The C3 table shows which
interval is the calibrated one.

The four C2 fails come from glmmTMB alone. On the site-slope designs (12 sites,
site-by-arm SD 0.20, so a slope SD of 0.28) the unstructured `(trt | site)` fit
lands on the boundary in 7% of seeds (A1/A6) and 4% (A3); with 8 varieties (P8) it
does so in 12%. flexyBayes fitted all of those seeds. Its posterior on the same
seeds covers the truth in 13 of 14 (A1), 7 of 8 (A3), 14 of 14 (A6) and 12 of 12
(P8), so excluding them does not flatter the posterior.

## C3: Bayesian calibration (both engines' coverage reported)

Coverage of the planted trial-conditional truth by the closed 95% interval, over
the seeds on which both engines succeeded. Band `[0.93, 0.97]`; BORDERLINE = outside
the band by no more than 2 Monte Carlo SE; a fits failure makes the cell FAIL
whatever its coverage (rule Section 6). glmmTMB coverage is reported, not gated.

| Cell | Slope | flexyBayes coverage (MC SE) | Wilson 95% | Coverage alone | Verdict | flexyBayes coverage, all 200 / 100 seeds | Misses below / above | glmmTMB delta coverage | Seeds near an interval end |
|---|---|---|---|---|---|---|---|---|---|
| A1 | site | **0.909** (0.021) | 0.859 to 0.942 | BORDERLINE | **FAIL (fits)** | 0.910 | 8 / 10 | **0.989** | 6 |
| A2 | none | 0.960 (0.014) | 0.923 to 0.980 | PASS | PASS | 0.960 | 5 / 3 | 0.955 | 7 |
| A3 | site | **0.922** (0.019) | 0.875 to 0.952 | BORDERLINE | **FAIL (fits)** | 0.920 | 11 / 5 | **1.000** | 3 |
| A4 | none | **0.925** (0.019) | 0.880 to 0.954 | BORDERLINE | **BORDERLINE** | 0.925 | 11 / 4 | 0.965 | 12 |
| A5 | none | 0.945 (0.016) | 0.904 to 0.969 | PASS | PASS | 0.945 | 6 / 5 | 0.940 | 6 |
| A6 | site | **0.925** (0.019) | 0.878 to 0.955 | BORDERLINE | **FAIL (fits)** | 0.930 | 6 / 8 | **0.995** | 2 |
| P8 | pooled, 8 | 0.932 (0.027) | 0.859 to 0.968 | PASS | **FAIL (fits)** | 0.940 | 4 / 2 | **1.000** | 3 |
| P12 | pooled, 12 | 0.939 (0.024) | 0.874 to 0.972 | PASS | PASS | 0.940 | 3 / 3 | **1.000** | 4 |

"Misses below" counts seeds whose truth lies below the interval.

**Answer to the question the slope cells were added for.** On the trial-conditional
estimand with a random treatment slope, the glmmTMB delta interval over-covers:
0.989 to 1.000 on all five slope cells. That is consistent with the reading in the
rule: it varies the fixed treatment coefficient while holding the slope modes fixed,
and so it states uncertainty the conditional estimand does not have. The posterior
interval is close to nominal but somewhat narrow on the site-slope designs (0.909
to 0.925, about 2 to 4 points under 0.95, each within 2 MC SE of the band edge) and
nominal on the pooled-variety designs (0.932 to 0.940 on the ok seeds, 0.940 on all
seeds). Neither engine's interval is calibrated on the site-slope cells. The
posterior misses by a few points on the narrow side; the delta interval misses by
4 to 5 points on the wide side, at about 1.7 times the posterior width.

**Reading of the posterior misses.** On A1 and A6 the misses fall on both sides
(8/10 and 6/8) and the bias is negligible (0.005 and 0.026 posterior SD). The
shortfall there is width, not location. On A3 (Gamma) and A4 (binomial) the misses
lean low (11/5, 11/4) while the mean bias stays small (0.10 and 0.01 SD), which
points to the shape of the interval more than its centre. Two causes are consistent
with this, and the run cannot tell them apart. (1) With 12 sites, the default prior
on the site-by-arm SD and INLA's integration over that hyperparameter may
understate how far the realised site-by-arm effects can sit from zero. (2) INLA's
Gaussian approximation to the latent field may be too light in the tails of a
non-Gaussian, non-log-linear functional (the binomial loss share is a ratio of
logit-scale means). Testing either needs a new dated rule (other priors, or a
sampler as reference). Under this rule nothing was rerun or changed.

Between 2 and 12 seeds per cell have their truth within 0.15 posterior SD of an
interval end, so an INLA rerun could flip them (Amendment 1). On A4 that is 12
seeds against a shortfall of 1 coverage event (185 of 200 against 186 for 0.93);
A4's BORDERLINE reading is therefore not robust to a rerun in either direction.

## C4: shrinkage recovery

Per seed and engine, RMSE over the levels of the pooled `pollin_cate()` estimates
against that of the unpooled `pollin_fit(by = "variety")` + `pollin_effects()`
estimates. Gate: pooled RMSE strictly lower in `>= 90%` of seeds on each engine;
flexyBayes pooled per-level coverage in `[0.92, 0.98]`; at most 2 of 100 failures
per fit.

| Cell | Levels | Engine | Seeds ok | Failed (pooled / unpooled) | Pooled wins | Mean RMSE pooled | Mean RMSE unpooled | Verdict |
|---|---|---|---|---|---|---|---|---|
| S8 | 8 | glmmTMB | 88 | **12** / 0 | **0.693** | 0.104 | 0.118 | **FAIL (fits, shrinkage)** |
| S8 | 8 | flexyBayes | 100 | 0 / 0 | **0.700** | 0.098 | 0.115 | **FAIL (shrinkage)** |
| S12 | 12 | glmmTMB | 99 | 1 / 0 | **0.717** | 0.105 | 0.121 | **FAIL (shrinkage)** |
| S12 | 12 | flexyBayes | 100 | 0 / 0 | **0.790** | 0.103 | 0.121 | **FAIL (shrinkage)** |

| Cell | flexyBayes pooled per-level coverage (nominal MC SE) | Wilson 95% | Verdict | Pairs near an end | glmmTMB pooled (approximate, not gated) | glmmTMB unpooled delta | flexyBayes unpooled |
|---|---|---|---|---|---|---|---|
| S8 | 0.941 (0.008), 800 pairs | 0.923 to 0.956 | PASS | 34 | 0.972 | 0.948 | 0.959 |
| S12 | 0.940 (0.007), 1,200 pairs | 0.925 to 0.952 | PASS | 42 | 0.970 | 0.931 | 0.939 |

**Reading of the shrinkage miss.** Pooling works on average: the pooled estimates
have 12% to 15% lower mean RMSE on both engines, and the median per-seed ratio of
pooled to unpooled RMSE is 0.82 to 0.86. The per-seed criterion fails because in a
fifth to a third of seeds pooling is worse. Its upper decile ratio is 1.10 to 1.29.
The design gave pooling little room: the realised SD of the true per-variety loss
shares (0.18 to 0.20) is larger than the per-level posterior SD of the unpooled
estimates (about 0.12), so the optimal shrinkage is mild and a misestimated slope
variance easily costs more than it gains. The flexyBayes win share is 0.80 (S8) and
0.86 (S12) in seeds where the realised spread of the level effects is below its
median, and 0.60 and 0.72 above it. The site-by-arm and variety-by-arm SD of 0.20
was raised from 0.10 to 0.15 before the rule, to keep glmmTMB off the boundary
(rule, preamble item 2). That choice made the levels more distinct and pooling's
advantage per seed smaller. It is recorded as the evidence on the cause, not as a
reason to discount the miss. The glmmTMB pooled fit also fails the fits limit on
S8 (12 of 100 boundary fits with 8 varieties).

## Agreement between engines on the pooled per-level effects

Both engines pool to a similar degree. Mean pooled RMSE is 0.104 and 0.098 on S8,
and 0.105 and 0.103 on S12. The glmmTMB conditional interval of `pollin_cate()`
over-covers (0.970 to 0.972), as its approximate construction suggests.

## Determinism (Amendment 1)

Data, planted truths, every glmmTMB result and every error column are bit-identical
across the four smoke runs on 1, 2 and 4 cores. The flexyBayes results are
reproducible to the INLA run-to-run variation. The fit itself differs across
processes in the fifth decimal of the log marginal likelihood, even single-threaded
with BLAS and OpenMP pinned. Across the smoke runs that variation moved posterior
medians by at most 0.011 SD and interval ends by at most 0.13 SD, and it changed no
coverage, agreement or win event. The seed counts near a deciding boundary are in
the tables above.

## Package findings (recorded in the rule, Section 9)

1. **`R/fit_flexybayes.R`, `.flexy_new_fit()` / `.flexy_fit_one()`, does not set
   `aggregate`.** `flexybayes()` then compresses Poisson and binomial rows with
   factor-only predictors into cells. `.flexy_arm_draws()` (the
   `latent[paste0("Predictor:", seq_len(nrow(d))), ]` line, about line 398) then
   stops with "subscript out of bounds". Reproducer: 8 blocks x 2 arms x 6 plants,
   Poisson counts, `pollin_fit(est, engine = "flexyBayes", family = poisson())`,
   `pollin_effects(fit)`: 16 `Predictor` rows for 96 data rows. The suite passes
   `aggregate = FALSE`.
2. **`pollin_data(denominator = )` reaches no fit**, so a binomial seed-set response
   out of `n` ovules cannot be fitted on either engine. The binomial cell used a
   Bernoulli response.

## Rule compliance

The rule and its amendment were committed before the full run. The single-seed
checks and the smoke runs read only errors, timings and between-run differences.
No cell was rerun, no seed count, design, prior or package code was changed after
the run, and the verdicts above are the rule's verdicts as computed by the suite.
Two departures are reported here: the wall time exceeded the stated budget, and
`results.csv` gives the C3 verdict of a fits-failed cell as `FAIL`. The table above
adds what its coverage alone would give.

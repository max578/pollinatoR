# Planted-truth recovery, Stage 3.1 (binomial denominator and controls): report

Rule: `validation/recovery_rule_S3.1.md`, committed alone at `c3a22bd`
(2026-11-09 22:10:06 +0800), before the smoke run and the full run. The probe that
informed it is `validation/probe_S3.1_slope_coverage.md` (`9af555e`). Run:
2026-11-09 22:12:28 to 22:34:47 (AWST), **22.3 min wall time** on 8 cores, 400 tasks,
inside the 45-minute budget; no amendment was needed. Suite:
`validation/planted_truth_S3.1.R`. Per-seed outcomes:
`validation/planted_truth_S3.1_seeds.csv.gz`. Per-cell summary:
`validation/planted_truth_S3.1_results.csv`. Session:
`validation/planted_truth_S3.1_session.txt`.

## Verdict

| Cell | Seeds ok | C2 agreement share | C2 width ratio | C2 | flexyBayes coverage (MC SE) | C3 | glmmTMB coverage (MC SE) |
|---|---|---|---|---|---|---|---|
| B1 (binomial, denominator) | 200 / 200 | 1.000 | **0.756** | **FAIL (width ratio)** | 0.960 (0.014) | PASS | 0.985 (0.009) |
| A2 (control, Poisson) | 100 / 100 | 1.000 | 0.989 | PASS | **0.980** (0.014) | **BORDERLINE** | 0.980 (0.014) |
| A5 (control, NB2, deficit) | 100 / 100 | 1.000 | 1.003 | PASS | 0.940 (0.024) | PASS | 0.950 (0.022) |

* **C2 is not met.** B1 misses the width-ratio band: the flexyBayes posterior
  interval is a median 0.756 of the glmmTMB delta width (interquartile 0.698 to
  0.809; only 29% of seeds reach 0.8). The engines agree on the estimate in every
  seed of every cell.
* **C3 is not met.** B1 PASSES (0.960). A2 is BORDERLINE on the high side: 98 of 100
  covered, one event above the band's 97. Under the rule this is a control
  regression against S3 (A2 covered 0.960 there), reported with the evidence on its
  cause below.
* No fit failed and no fit warned, on either engine, in any cell.

## B1: the binomial-denominator route

The denominator route works as specified on both engines: 200 of 200 seeds fitted,
2,690 to 3,012 trials per data set (median 2,871) expanded to one row per ovule on
flexyBayes, at a median 20 s per flexyBayes fit (90th percentile 21 s) against
0.3 s on glmmTMB. The posterior median lies within 0.25 glmmTMB SE of the glmmTMB
estimate in every seed (median gap 0.023 SE).

| Measure | flexyBayes posterior | glmmTMB delta |
|---|---|---|
| Coverage of the planted truth (MC SE) | 0.960 (0.014), Wilson 0.923 to 0.980 | 0.985 (0.009) |
| Misses below / above | 7 / 1 | 2 / 1 |
| sd(z), z = (estimate - truth) / SE | 0.975 | 0.748 |
| mean(z) | 0.196 | 0.143 |
| Bias (loss share) | 0.0037 | 0.0042 |
| SE, median ratio glmmTMB / flexyBayes | | 1.31 |

**Reading of the width-ratio miss.** The miss is the delta interval's, not the
posterior's. The posterior is calibrated for the planted trial-conditional estimand
(coverage 0.960, sd(z) 0.975); the delta interval is about 31% wider and
over-covers (0.985, sd(z) 0.748, so its SE is about a third larger than the spread
of its own estimates about the truth). This is the mechanism the S3 rule described
for the slope cells (Section 3) and anticipated for a logit-link intercept cell
(Section 8, item 3): the delta interval varies the fixed coefficients while holding
the conditional modes of the block effects as known, so the uncertainty in the
intercept, which reflects how far the 8 blocks' mean could sit from zero, enters
the interval although the trial-conditional estimand does not carry it (given the
data, the intercept and the block modes are negatively correlated, and the sum that
the estimand depends on is well determined). Under a log link the intercept cancels
from the ratio `H / L`, and on A2 and A5 the width ratios are 0.989 and 1.003. Under
the logit link it does not cancel. With 30 ovules per plant the plant-level
proportions are pinned down, so the intercept's share of the delta variance is
large; in S3, A4 (one 0/1 row per flower, less information per block) gave a width
ratio of 0.934 by the same mechanism, smaller in size. The rule gated the width
ratio on B1 because it has no slope, and under the rule the cell FAILS C2. Nothing
was rerun or regated.

**Location.** The flexyBayes misses lean low (7 below, 1 above), and mean(z) is
0.196 (about 2.8 times its SE of 0.07 at 200 seeds). glmmTMB shows the same lean
(mean(z) 0.143 in its own SE units; bias 0.0042 against 0.0037), so it is shared by
both engines: a small upward bias of the estimated loss share, about 0.004 on a
truth near 0.26, which is not specific to the posterior. The run cannot say more
about its source.

Six B1 seeds have their truth within 0.15 posterior SD of an interval end (4 covered,
2 not), and none lies that close to the agreement tolerance. With 192 of 200 covered,
flipping all of them one way gives 188 (0.940) or 194 (0.970), both inside the band,
so the C3 PASS is robust to the INLA run-to-run variation.

## Controls: A2 and A5

| Cell | S3 (200 seeds) coverage | S3.1 (100 new seeds) coverage | S3.1 glmmTMB | sd(z), flexyBayes / glmmTMB | Width ratio S3 / S3.1 | Seeds near an end |
|---|---|---|---|---|---|---|
| A2 | 0.960 | **0.980** (BORDERLINE) | 0.980 | 0.878 / 0.885 | 0.996 / 0.989 | 4 |
| A5 | 0.945 | 0.940 | 0.950 | 1.031 / 1.049 | 1.005 / 1.003 | 6 |

**Reading of the A2 BORDERLINE.** A2 over-covers by one event (98 of 100 against a
band edge of 97), 0.01 above the band and within 1 Monte Carlo SE of it. The
evidence points to Monte Carlo variation on these 100 seeds, not to a change in the
flexyBayes route:

* the glmmTMB delta interval, an independent method, covers exactly as often (0.980)
  and has the same sd(z) (0.885 against 0.878), so on these seeds the realised
  errors are small relative to both engines' SEs. A shared sd(z) of 0.88 is 1.7 SE
  (about 0.07 at 100 seeds) below 1;
* the width ratio is unchanged from S3 (0.989 against 0.996), and the engines agree
  in every seed;
* a fit without a denominator takes the same path in 0.3.1 as in 0.3.0
  (`R/fit.R`, `R/fit_flexybayes.R`; the denominator code runs only when one is
  declared).

Under the rule this is still a BORDERLINE miss of C3 on a control and is reported as
a regression against S3; it is not excused. A5 passes both criteria, matching S3.

## Determinism

The smoke run (3 seeds per cell, on 1 and on 4 cores) gave bit-identical simulated
data, trial counts, planted truths, glmmTMB results, error and warning columns. The
flexyBayes posterior medians differed by at most 1.7e-8 posterior SD and the
interval ends by at most 3.3e-6 SD; no coverage or agreement event changed (0 of 18).
The INLA variation here is far smaller than S3 measured on the slope designs (S3
Amendment 1), as S3 found for intercept-only Poisson and binomial fits. Seeds near a
deciding boundary in the full run: coverage 6 (B1), 4 (A2), 6 (A5); agreement 0, 0
and 2.

## Package findings

No package bug was found. The denominator route fitted every seed on both engines,
with no warning. One property of an existing interval is documented by this run:

* **The glmmTMB delta interval of `pollin_effects()` over-covers the
  trial-conditional estimand on a logit-link fit with a random intercept**, by about
  a third in SE on B1 (coverage 0.985, sd(z) 0.75), because it holds the conditional
  modes as known while varying the intercept (`interval = "delta"` states
  fixed-effect uncertainty only). The flexyBayes posterior is calibrated on the same
  seeds. This is the cause of the B1 C2 miss.

INLA writes `sh: /bin/kstat: No such file or directory` to standard error on every
fit in this environment (a CPU-probe call that the sandbox does not provide); it
does not affect the fits and is not recorded as a warning.

## Rule compliance

The rule was committed alone before the smoke run and the full run. Before it, only
the S3.1 probe (its own seeds) and a single-seed check of each cell (errors,
warnings and timings read) were run. The smoke run's comparison read only
between-run differences. No seed count, design, prior, threshold or package code was
changed after the rule, no cell was rerun, and the verdicts above are the suite's.
The slope cells (A1, A3, A6, P8, P12) and the C4 cells (S8, S12) were not re-run, as
the rule states; their S3 results and the probe stand as the evidence for the label
decision on slope designs, which is Max's.

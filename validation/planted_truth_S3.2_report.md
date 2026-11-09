# Planted-truth recovery, Stage 3.2 (joint delta interval, 0.4.0): report

Rule: `validation/recovery_rule_S3.2.md`, committed alone at `4f8e39f`
(2026-11-10 01:42 +0800), before the smoke run and the full run. Suite:
`validation/planted_truth_S3.2.R` (`b43f98f`, the HEAD recorded in the session file).
Per-seed outcomes: `planted_truth_S3.2_seeds.csv.gz`; per-cell summary:
`planted_truth_S3.2_results.csv`; session: `planted_truth_S3.2_session.txt`.
Full run: 800 tasks on 8 cores, **44.3 min** wall time, inside the 45-minute budget;
no amendment. The smoke run (16 tasks, 4 cores, scratch directory) was read for
errors, warnings and timings only.

## Verdict

All six cells PASS their gates.

| Cell | Gate | Seeds ok / run | Joint coverage (MC SE) | Width ratio, posterior / joint (IQR) | Verdict |
|---|---|---|---|---|---|
| B1 binomial, denominator | G3 | 200 / 200 | 0.935 (0.017) | 0.992 (0.975 to 1.011) | PASS |
| A2 Poisson (control) | G3 | 100 / 100 | 0.930 (0.026) | 0.984 (0.972 to 1.008) | PASS |
| A5 NB2, deficit (control) | G3 | 100 / 100 | 0.930 (0.026) | 1.004 (0.987 to 1.024) | PASS |
| A1 NB2, slope, 12 sites | G4 width | 94 / 100 | 0.957 (0.021) | 0.981 (0.960 to 0.998) | PASS |
| A3 Gamma, slope, 12 sites | G4 width | 99 / 100 | 0.960 (0.020) | 0.997 (0.973 to 1.015) | PASS |
| A1W NB2, slope, 24 sites | G4 coverage | 198 / 200 | 0.970 (0.012) | 0.990 (0.969 to 1.013) | PASS |

A2 and A5 sit on the lower edge of the coverage band (0.930 is inside the closed band
`[0.93, 0.97]`), and A1W on the upper edge; B1 is 0.005 above the lower edge.

## Before and after on the same seeds

| Cell | Coverage: joint / delta_fixed / posterior | sd(z): joint / delta_fixed / posterior | SE ratio, delta_fixed / joint (median) |
|---|---|---|---|
| B1 | 0.935 / 0.980 / 0.930 | 1.062 / 0.820 / 1.073 | 1.31 |
| A2 | 0.930 / 0.930 / 0.920 | 0.965 / 0.965 / 0.975 | 1.00 |
| A5 | 0.930 / 0.930 / 0.940 | 1.076 / 1.076 / 1.058 | 1.00 |
| A1 | 0.957 / 1.000 / 0.968 | 1.043 / 0.702 / 1.076 | 1.57 |
| A3 | 0.960 / 1.000 / 0.970 | 0.984 / 0.535 / 0.950 | 1.90 |
| A1W | 0.970 / 0.990 / 0.975 | 1.006 / 0.659 / 1.015 | 1.57 |

Where the estimand depends on the random effects (logit link, treatment slope), the
0.3.1 interval was 1.3 to 1.9 times too wide (sd(z) 0.54 to 0.82) and over-covered;
the joint interval has sd(z) between 0.98 and 1.06 and matches the posterior width
(median ratio 0.98 to 1.00). Under a log link without a slope the two delta intervals
are the same (ratio 1.00), as expected: the intercept cancels from the ratio.

The joint and fixed-only estimates differ by at most 3.5e-5 (A1W); the joint estimate
is the package's G-computation, and the fixed-only one comes from `marginaleffects`
(rule, amendment to G2).

## Fits and warnings

flexyBayes fitted every seed. glmmTMB did not converge on 9 seeds, all on slope cells
(A1 6, A3 1, A1W 2; non-positive-definite Hessian or singular convergence), which the
rule reports and does not gate; every glmmTMB warning in the run is on one of those 9
seeds. No counted seed warned, and no seed used the fallback that holds the variance
parameters at their estimates.

Misses (joint): B1 3 below / 10 above, A2 5 / 2, A5 4 / 3, A1 1 / 3, A3 0 / 4, A1W
5 / 1. Seeds with the truth within 0.15 SE of an interval end: 10 (B1), 4, 4, 6, 7, 6.
On A2 and A5, four such seeds can move coverage by 0.04 either way, so a PASS on the
band edge there is not robust to Monte Carlo variation. Their delta intervals are the
same under 0.3.1 and 0.4.0 (SE ratio 1.00), so this is not a change made by 0.4.0.

## Not settled by this run

The posterior covers 0.968 (A1) and 0.970 (A3) on these 12-site slope seeds, where
S3 and probe S3.1 measured 0.909 to 0.930 on other seeds of the same designs. The
difference is about 2 to 3 Monte Carlo SE at 100 seeds; this run cannot say whether
it is chance. It does not bear on the gates, which compare the joint interval to the
posterior on the same seeds.

# Planted-truth recovery, Stage 1 (criterion C2): report

Rule: `validation/recovery_rule_S1.md`, committed alone at `6f20304`
(2026-11-08 22:28:39 +0800) before the first suite run. The rule was not changed.
Run: 2026-11-08 22:30:15 to 22:47:10 (AWST), 16.9 min wall time on 8 cores,
3,000 design x seed tasks, 25,000 fits. Suite: `validation/planted_truth_S1.R`.
Per-seed outcomes: `validation/planted_truth_S1_seeds.csv.gz`; per-cell summary:
`validation/planted_truth_S1_results.csv`; session: `validation/planted_truth_S1_session.txt`.

## Verdict

**24 of 25 cells PASS; 1 FAILS** (D3, E2, difference scale, on bias). Every cell
ran its full 1,000 seeds; no seed failed.

The coverage event counted is: the per-seed truth lies inside the closed delta
interval `[conf.low, conf.high]`. Bias/SE uses the Monte-Carlo SD of the estimates.

| Design | Family | Scale | Basis | n_ok | Failed | Coverage (Wilson 95%) | Bias/SE | Verdict |
|---|---|---|---|---|---|---|---|---|
| D1 three-arm | E1 deficit | difference | loss | 1000 | 0 | 0.946 (0.930, 0.958) | -0.056 | PASS |
| D1 three-arm | E1 deficit | ratio | gain | 1000 | 0 | 0.941 (0.925, 0.954) | 0.016 | PASS |
| D1 three-arm | E1 deficit | ratio | loss | 1000 | 0 | 0.941 (0.925, 0.954) | -0.083 | PASS |
| D1 three-arm | E1 deficit | rii | loss | 1000 | 0 | 0.941 (0.925, 0.954) | -0.039 | PASS |
| D1 three-arm | E1 deficit | lnrr | loss | 1000 | 0 | 0.941 (0.925, 0.954) | -0.034 | PASS |
| D1 three-arm | E2 dependence | difference | loss | 1000 | 0 | 0.958 (0.944, 0.969) | -0.011 | PASS |
| D1 three-arm | E2 dependence | ratio | gain | 1000 | 0 | 0.936 (0.919, 0.950) | 0.071 | PASS |
| D1 three-arm | E2 dependence | ratio | loss | 1000 | 0 | 0.936 (0.919, 0.950) | -0.036 | PASS |
| D1 three-arm | E2 dependence | rii | loss | 1000 | 0 | 0.936 (0.919, 0.950) | 0.007 | PASS |
| D1 three-arm | E2 dependence | lnrr | loss | 1000 | 0 | 0.936 (0.919, 0.950) | 0.018 | PASS |
| D2 split-plot | E1 deficit | difference | loss | 1000 | 0 | 0.945 (0.929, 0.958) | -0.021 | PASS |
| D2 split-plot | E1 deficit | ratio | gain | 1000 | 0 | 0.937 (0.920, 0.951) | 0.030 | PASS |
| D2 split-plot | E1 deficit | ratio | loss | 1000 | 0 | 0.937 (0.920, 0.951) | -0.051 | PASS |
| D2 split-plot | E1 deficit | rii | loss | 1000 | 0 | 0.937 (0.920, 0.951) | -0.013 | PASS |
| D2 split-plot | E1 deficit | lnrr | loss | 1000 | 0 | 0.937 (0.920, 0.951) | -0.011 | PASS |
| D2 split-plot | E2 dependence | difference | loss | 1000 | 0 | 0.958 (0.944, 0.969) | -0.076 | PASS |
| D2 split-plot | E2 dependence | ratio | gain | 1000 | 0 | 0.938 (0.921, 0.951) | 0.006 | PASS |
| D2 split-plot | E2 dependence | ratio | loss | 1000 | 0 | 0.938 (0.921, 0.951) | -0.079 | PASS |
| D2 split-plot | E2 dependence | rii | loss | 1000 | 0 | 0.938 (0.921, 0.951) | -0.044 | PASS |
| D2 split-plot | E2 dependence | lnrr | loss | 1000 | 0 | 0.938 (0.921, 0.951) | -0.037 | PASS |
| D3 matched-pair | E2 dependence | difference | loss | 1000 | 0 | 0.950 (0.935, 0.962) | **-0.194** | **FAIL (bias)** |
| D3 matched-pair | E2 dependence | ratio | gain | 1000 | 0 | 0.939 (0.922, 0.952) | -0.004 | PASS |
| D3 matched-pair | E2 dependence | ratio | loss | 1000 | 0 | 0.939 (0.922, 0.952) | -0.070 | PASS |
| D3 matched-pair | E2 dependence | rii | loss | 1000 | 0 | 0.939 (0.922, 0.952) | -0.041 | PASS |
| D3 matched-pair | E2 dependence | lnrr | loss | 1000 | 0 | 0.939 (0.922, 0.952) | -0.037 | PASS |

E1 is not run on D3 (no supplemented arm), as the rule states.

## Observations, stated without changing any verdict

**Ratio-type coverage sits low in the window.** All 12 ratio, rii and lnrr
coverage values lie between 0.936 and 0.941: inside `[0.93, 0.97]`, and every
Wilson interval contains 0.95 except D1 E2 (upper limit 0.950). The mean reported
SE is 2 to 6% below the Monte-Carlo SD in every ratio-type cell (for example D1 E2
lnrr: 0.102 reported, 0.108 empirical), consistent with Wald intervals from a
small number of blocks. These cells pass; the margin is narrow.

**The one miss: D3 difference, bias -0.194 SE.** The estimate is on average 2.66
seeds per plant below a truth of 44.3 (6%), against an empirical SE of 13.7.
Coverage is 0.950, so the interval is not too narrow; the estimate is shifted.

Evidence that the shift comes from the level of the averaged predictions, not
from the treatment contrast:

1. On the same D3 data sets, the ratio-type scales pass with bias/SE between
   -0.07 and 0.00. The random intercepts cancel from `H / L`, so the ratio
   depends on the fixed treatment effect alone, and that is recovered.
2. A diagnostic outside the rule (300 seeds per design, scratch only, not part of
   any verdict) compared the fitted excluded-arm average, `mean(predict(type =
   "response"))` with `re.form = NULL`, with the realised truth `L_s`. The ratio
   fitted / true averaged 0.954 in D3, against 0.990 in D1 and 0.988 in D2; the
   fitted D3 treatment ratio averaged 0.998 of its true value.
3. The difference scale multiplies the contrast by that level. Under a log link
   the average of `exp(estimated random effect)` falls below the average of
   `exp(realised random effect)` when the estimates are shrunk towards zero
   (Jensen's inequality). D3 has the largest random-effect variance (site SD
   0.50, pair SD 0.30) and the least-informed level (two plants per pair), so its
   pair effects are shrunk hardest. D1 and D2 show the same sign on the
   difference scale (bias/SE -0.011 to -0.076), smaller and inside the rule.

So the difference-scale estimate under `population = "trial"` is biased downward
by the shrinkage of the random effects it averages over, and the bias grows with
the random-effect variance and with how few observations inform each unit. This is
a property of the plug-in G-computation estimator with estimated random effects,
not a coding error in `R/` (the delta path reproduces the manual average of the
conditional predictions). The delta interval carries fixed-effect uncertainty
only, as stated in the rule; in D3 that interval still covers at 0.950, but it does
not remove the shift.

**Warnings.** One D3 data set (seed 947) produced a `glmmTMB` non-positive-definite
Hessian warning in all five D3 cells. Its estimates and intervals were finite, so
by the rule it counts as successful (`n_warned = 1` in each D3 cell). For that fit
`marginaleffects` printed that it replaced the variance matrix with the nearest
positive-definite matrix before forming the delta interval (five messages in the
run log, one per cell). Its truth lay inside its interval in all five cells.
Dropping it would not change any verdict.

## Reproducibility

Each task sets `RNGkind("Mersenne-Twister", "Inversion", "Rejection")` and
`set.seed(100000 * k + s)`; a 3-seed run on 1 core and on 4 cores gave
byte-identical per-seed files. Versions: R 4.6.1, glmmTMB 1.1.14, TMB 1.9.25,
marginaleffects 0.32.0, Matrix 1.7.6. The session file records HEAD `8c4dcb2`
(read when the run ended); commits made in parallel during the run touched only
`validation/croppol/`, `vignettes/`, `NEWS.md`, `README.md` and `_pkgdown.yml`, so
the `R/` code exercised is that of `6f20304`.

## What C2 can claim

E1 and E2 recover a planted contrast on the ratio (loss and gain), RII and lnRR
scales in all three designs at the stated minimum sizes, and on the difference
scale in the three-arm and split-plot designs. On the difference scale in the
matched-pair field design (20 sites x 5 pairs, site SD 0.50, pair SD 0.30) the
estimate is biased low by 0.19 empirical SE, which fails the pre-registered bias
criterion. Any change to that, such as a different difference-scale target or a
larger design, needs a new dated rule and a new run, reported beside this one.

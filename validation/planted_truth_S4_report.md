# Planted-truth recovery, Stage 4 (identification, 0.5.0): report

Rule: `validation/recovery_rule_S4.md`, committed alone at `0c00a49`
(2026-11-14 13:07 +1030), before the smoke run and the full run. Suite:
`validation/planted_truth_S4.R` (`f0a8459`); the session file records HEAD `14c5b5c`,
which changes documentation only. Per-seed outcomes: `planted_truth_S4_seeds.csv.gz`;
per-gate summary: `planted_truth_S4_results.csv`; session: `planted_truth_S4_session.txt`.
Full run: 7,200 tasks on 8 cores, **27.3 min** wall time, inside the 60-minute budget;
no amendment. The smoke run went to a scratch directory and was read for errors,
warnings and timings only.

**Reproduced.** The full run was made twice from the same HEAD, in two processes
(13:13 to 13:48, 35.6 min; 13:41 to 14:09, 27.3 min). The two gate tables are
byte-identical (md5 `ac27b13548f4b29172b64edea37f3ea9`). The per-seed files differ
only in timings and in the flexyBayes X1B estimates (109 of 200 seeds, largest
difference 8e-7), which INLA does not make bit-reproducible across processes. The
files committed here are from the second run.

Planted values computed by the suite: M1 total 0.252861, M1 direct 0.139292.

## Verdict

All 17 gate rows PASS. No row is BORDERLINE.

| Gate | Cell | Quantity | Value (MC SE) | Requirement | Seeds ok / failed | Verdict |
|---|---|---|---|---|---|---|
| G-M1 | M1a | total coverage | 0.937 (0.008) | [0.93, 0.97] | 1000 / 0 | PASS |
| G-M1 | M1a | direct coverage | 0.945 (0.007) | [0.93, 0.97] | 1000 / 0 | PASS |
| G-M1 | M1a | labels as fitted | 1 | all seeds | 1000 / 0 | PASS |
| G-M1 | M1b | total coverage | 0.937 (0.008) | [0.93, 0.97] | 1000 / 0 | PASS |
| G-M1 | M1b | direct coverage | 0.932 (0.008) | [0.93, 0.97] | 1000 / 0 | PASS |
| G-M1 | M1b | labels as fitted | 1 | all seeds | 1000 / 0 | PASS |
| G-X1 | X1 | cage factor coverage | 0.957 (0.006) | [0.93, 0.97] | 1000 / 0 | PASS |
| G-X1 | X1 | net dependence coverage | 0.933 (0.008) | [0.93, 0.97] | 1000 / 0 | PASS |
| G-F1 | F1 | false-call rate, within | 0.059 (0.008) | <= 0.07 | 867 / 133 | PASS |
| G-F1 | F1 | false-call rate, between | 0.050 (0.007) | <= 0.07 | 867 / 133 | PASS |
| G-F1 | F1 | negative control rejects (between) | 0.999 (0.001) | >= 0.80 | 867 / 133 | PASS |
| G-F2a | F2 | moderation coverage, graph's set | 0.953 (0.008) | [0.93, 0.97] | 634 / 366 | PASS |
| G-F2b | F2 | unadjusted bias minus planted bias | -0.0024 (0.0013) | within 2 MC SE of 0 | 801 / 199 | PASS |
| G-F3 | F3 | between-only term coverage | 0.938 (0.010) | [0.93, 0.97] | 642 / 358 | PASS |
| G-X1b | X1B | flexyBayes cage rows returned | 1 | >= 0.98 | 200 / 0 | PASS |
| G-X2a | X2 | algebraic bound contains 0.25 | 1 | every inside seed | 661 / 0 | PASS |
| G-X2b | X2 | outer interval coverage, inside seeds | 0.995 (0.003) | >= 0.95 | 661 / 0 | PASS |

M1b direct (0.932) and X1 net dependence (0.933) sit 0.002 and 0.003 above the lower
edge of the band, inside one Monte Carlo SE of it.

## Total and direct effects (M1)

| Cell | Total: mean, bias, median width | Direct: mean, bias, median width | Misses below / above (total; direct) |
|---|---|---|---|
| M1a, 24 fields | 0.2506, -0.0023, 0.146 | 0.1368, -0.0025, 0.175 | 30 / 33; 29 / 26 |
| M1b, 44 fields | 0.2509, -0.0020, 0.108 | 0.1386, -0.0007, 0.129 | 33 / 30; 42 / 26 |

The two estimands differ by a factor of about 1.8 on the same data, and every row
carries the label it was fitted under. Adjusting for the mediator is a different
question, not a correction: the direct interval is wider than the total one by about
a fifth at both sizes.

## Moderation and identification (F1 to F3)

- **F1.** With an equal-shift confounder and no planted moderation, neither
  treatment-by-bees term rejects more often than 0.06. The negative control (the
  `bees_between` main term, which picks up the unmodelled fertility shift) rejects in
  0.999 of seeds between sites and 0.036 within sites, as it should: the confounder
  acts between sites only.
- **F2.** With the graph's adjustment set the moderation term covers 0.10 in 0.953
  of seeds (bias 0.0002). Without the set the estimate is biased by 0.153 on average
  (MC SE 0.0016), and that bias matches the planted bias `0.20 b_z` (mean 0.155; mean
  `b_z` 0.777) to within 0.0024, under 2 MC SE.
- **F3.** With an exposure constant within each site, every converged seed returns
  the between term with an estimate and standard error, names `bees` in
  `between_only`, and covers the planted 0.10 in 0.938 of seeds (bias 0.0000). Before
  0.5.0 this shape returned an all-missing row.

## Cage and bounds (X1, X2)

- **X1.** The cage-factor interval covers 1.10 in 0.957 of seeds (bias -0.0016,
  median width 0.214) and the dependence net of the cage, with `reference = "sham"`,
  covers 0.25 in 0.933 (bias -0.0029, median width 0.148).
- **X1B.** flexyBayes returned the cage row on all 200 seeds; posterior coverage
  0.940, median width 1.005 times the glmmTMB interval (reported, not gated).
- **X2.** On every inside seed (k in [1, 1.2]; 661 seeds) the algebraic bound
  contains 0.25, and the outer interval covers 0.25 in 0.995. On the 339 outside
  seeds (positive control) the outer interval covers 0.920, and the median bounds
  contain 0.25 in 0.316 of seeds against 0.787 inside: the bound loses the truth when
  the user's cage-factor range is wrong, as it should.

## Fits and warnings

No call errored, and no seed warned on M1, X1, X1B or X2. glmmTMB did not converge
on 1,056 fits, all in the `multi_site` cells with a random treatment slope, which the
rule reports and does not gate: F1 133, F2 unadjusted 199, F3 358 (690 in all), and
F2 with the graph's set 366. By first warning: singular convergence 914 (572 + 342),
non-positive-definite Hessian 142 (118 + 24). F3's planted slope SD (0.10) is small
against 44 sites, so about a third of its fits put the slope variance at the
boundary; the gates read on converged seeds, as in S3.2.

Median seconds per fit step: 0.33 to 0.59.

## Not settled by this run

Coverage on F2 and F3 is read on the 63 to 64% of seeds that converged. If
non-convergence is related to the size of the estimate on a seed, coverage over all
seeds could differ from the figures above; this run does not test that.

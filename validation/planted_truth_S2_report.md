# Planted-truth recovery, Stage 2 (criterion C2: value intervals): report

Rule: `validation/recovery_rule_S2.md`, committed alone at `2c8d01b`
(2026-11-09 17:10:43 +0800) before the first suite run. The rule was not changed.
Run: 2026-11-09 17:12:35 to 17:20:33 (AWST), 8.0 min wall time on 8 cores, 4,200
design x seed tasks (4,000 simulation, 200 bootstrap). Suite:
`validation/planted_truth_S2.R`. Per-seed outcomes:
`validation/planted_truth_S2_seeds.csv.gz`; per-cell summary:
`validation/planted_truth_S2_results.csv`; session:
`validation/planted_truth_S2_session.txt`.

## Verdict

**C2 is not met. 9 of the 11 C2 cells PASS; 2 FAIL on bias** (V06 and V10, both
dependence channel with market uncertainty on, bias/SE -0.140 and -0.126 against the
limit 0.1). **Coverage is inside `[0.93, 0.97]` in every cell**, 0.937 to 0.962.
Every cell ran its full seeds; no seed failed and none warned. The supplementary
bootstrap cell B01 passes.

The coverage event is: the per-seed true value lies inside the closed 95% percentile
interval of the value draws. Bias is the mean of (value median - true value); SE is
the mean over seeds of the SD of the value draws. Values are in AUD; width is the
median over seeds of the interval width.

| Cell | Design | Channel | Market | Method | Seeds | Coverage (MC SE) | Bias/SE | Median width | Ratio coverage | Verdict |
|---|---|---|---|---|---|---|---|---|---|---|
| V01 | D1 three-arm | dependence | off | simulation | 1000 | 0.952 (0.007) | -0.039 | 2,072,906 | 0.952 | PASS |
| V02 | D1 three-arm | dependence | on | simulation | 1000 | 0.956 (0.006) | -0.092 | 2,412,044 | 0.952 | PASS |
| V03 | D1 three-arm | deficit | off | simulation | 1000 | 0.940 (0.008) | 0.016 | 3,612,278 | 0.940 | PASS |
| V04 | D1 three-arm | deficit | on | simulation | 1000 | 0.937 (0.008) | -0.003 | 3,743,083 | 0.940 | PASS |
| V05 | D2 split-plot | dependence | off | simulation | 1000 | 0.944 (0.007) | -0.091 | 1,797,546 | 0.944 | PASS |
| V06 | D2 split-plot | dependence | on | simulation | 1000 | 0.947 (0.007) | **-0.140** | 2,082,676 | 0.944 | **FAIL (bias)** |
| V07 | D2 split-plot | deficit | off | simulation | 1000 | 0.946 (0.007) | 0.079 | 2,696,712 | 0.946 | PASS |
| V08 | D2 split-plot | deficit | on | simulation | 1000 | 0.945 (0.007) | 0.061 | 2,782,004 | 0.946 | PASS |
| V09 | D3 matched-pair | dependence | off | simulation | 1000 | 0.946 (0.007) | -0.082 | 1,603,055 | 0.946 | PASS |
| V10 | D3 matched-pair | dependence | on | simulation | 1000 | 0.952 (0.007) | **-0.126** | 1,781,740 | 0.946 | **FAIL (bias)** |
| V11 | D4 two-stratum population | dependence | off | simulation | 1000 | 0.962 (0.006) | -0.044 | 2,303,125 | 0.962 | PASS |
| B01 | D1 three-arm (supplementary) | dependence | off | bootstrap, 199 refits | 200 | 0.955 (0.015) | 0.015 | 2,013,958 | 0.955 | PASS |

Truth means (AUD): V01, V02 2.60e6 and 2.62e6; V03, V04 1.56e6 and 1.57e6; V05, V06
2.23e6 and 2.25e6; V07, V08 0.78e6 and 0.79e6; V09, V10 1.56e6 and 1.57e6; V11 1.80e6
(population loss share 0.231 on average across seeds). The market-on truth means
sit above the market-off ones by the lognormal mean factor, as the rule states.

## The two misses

Both misses are dependence-channel cells with market on, and in both the
market-off cell on the same data and the same estimand draws passes with bias/SE
close to the limit (V05 -0.091, V09 -0.082). Market on and market off share the fit
and the estimand draws, so the difference between a pair isolates the market
propagation. Decomposed in units of the market-on SE:

| Cell | Estimand bias (market-off bias, rescaled) | Shift of the value median | Shift of the truth | Bias/SE on | MC SE of bias/SE |
|---|---|---|---|---|---|
| V02 | -0.034 | -0.030 | +0.028 | -0.092 | 0.030 |
| V04 | +0.015 | -0.008 | +0.011 | -0.003 | 0.031 |
| V06 | -0.079 | -0.026 | +0.035 | -0.140 | 0.031 |
| V08 | +0.076 | -0.006 | +0.009 | +0.061 | 0.031 |
| V10 | -0.074 | -0.022 | +0.030 | -0.126 | 0.030 |

(bias on = estimand bias + median shift - truth shift.)

Two parts, both of the same sign on the dependence channel:

1. **The estimand part: the plug-in loss share is biased low.** The median of the
   simulation draws of `D` is the plug-in estimate `1 - exp(-delta_hat)`, and `D` is
   a concave function of the log ratio `delta`, so its plug-in is biased downward by
   Jensen's inequality; the relative bias of the ratio median is -0.8% (D1), -1.9%
   (D2) and -2.2% (D3), larger where the log ratio is less precisely estimated. The
   same sign appeared in S1 on the ratio-loss scale (bias/SE -0.036 to -0.083). On
   the deficit channel the gain `exp(delta) - 1` is convex and the bias is upward
   (V07 +0.076), which is why the deficit cells move the other way and pass.
2. **The market part: a median set against a lognormal truth.** With market on, the
   value median is 0.993 to 0.995 of the market-off median, while the true value is
   on average 1.0066 to 1.0087 of its market-off counterpart (expected 1.0063; the
   excess is the Monte Carlo spread of the realised true price and quantity). The
   two together add -0.04 to -0.06 SE. The rule (Section 8, item 3) anticipated this
   offset at about 0.03 SE; it was larger, because the median of `D x P x Q` falls
   below `median(D) x median(P x Q)` as well, which the advance estimate left out.

Neither part is a coding error. The ratio diagnostic shows the estimand draws
recover `D` with the expected coverage (0.944 and 0.946) and the small Jensen bias
of a plug-in ratio; the value draws are `D x P x Q` per draw, as `R/value.R`
specifies, and the median is the point value `pollin_value()` reports. The misses
are 1.3 (V06) and 0.9 (V10) Monte Carlo SEs beyond the limit, so the criterion is
missed by the rule but not by a margin the run can separate from 0.1; the cause,
though, is systematic and would persist with more seeds.

What the miss means for a user: on the dependence channel, the reported median value
sits a little below the expected value, by about 1 to 2% of the value in these
designs, about a tenth of a standard error. The interval is calibrated. Under the
rule this is documented, not tuned away; a change of point summary (for example the
mean of the value draws) or of the criterion is a new dated rule and a new run.

## Observations, stated without changing any verdict

**Coverage.** All twelve cells lie between 0.937 and 0.962. The S1 risk that
simulation draws under-cover on ratio scales did not materialise on these seeds: the
lowest is V04 at 0.937 (deficit, D1), and every Wilson 95% interval includes 0.95
(V04 0.920 to 0.950 at its upper limit; all intervals are in the results CSV).
With market on, coverage moved by -0.003 to +0.006 against market off.

**Population (V11).** The two-stratum population value covers at 0.962, the highest
of the suite, inside the window, with bias/SE -0.044. Its heaviest stratum has 5
blocks for 75% of the area.

**Bootstrap (B01).** 200 seeds, 199 refits each: coverage 0.955 (MC SE 0.015),
bias/SE 0.015, a median width 3% narrower than the simulation cell V01 on the same
data sets. Indicative only, as the rule states.

## Reproducibility

Each task sets `RNGkind("Mersenne-Twister", "Inversion", "Rejection")` and
`set.seed(2000000 + 100000 * k + s)`; estimand draws use `pollin_draws(seed =)`. A
3-seed smoke run (2 bootstrap seeds) on 1 core and on 4 cores gave byte-identical
per-seed files. Versions: R 4.6.1, glmmTMB 1.1.14, TMB 1.9.25, marginaleffects
0.32.0, Matrix 1.7.6. The session file records HEAD `4f146fb` (read when the run
ended). Commits made in parallel during the run touched `validation/value_oracle/`,
`vignettes/`, `NEWS.md`, `README.md`, `_pkgdown.yml` and, at 17:16, the print method
of `R/value.R` (`4f146fb`, output text only); the code was loaded at 17:12, so the
`R/` code exercised is that of `2c8d01b`, and no computed value depends on the print
change.

## What C2 can claim

Value intervals from simulation draws of the estimand are calibrated: 95% coverage
lies in `[0.93, 0.97]` for both channels, with and without market uncertainty, in the
three-arm, split-plot and matched-pair designs at the S1 minimum sizes and for a
two-stratum population value. The value median is unbiased within 0.1 SE on the
deficit channel and on the dependence channel without market uncertainty; on the
dependence channel with lognormal market uncertainty it is biased low by 0.13 to
0.14 SE in the split-plot and matched-pair designs, which fails the pre-registered
bias criterion.

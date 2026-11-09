# Planted-truth recovery rule, Stage 2 (criterion C2: value intervals)

Fixed 2026-11-09, before any simulation run of the suite. This file is committed on
its own; the commit time is the evidence that the rule preceded the results. Nothing
here changes after the run. A changed rule is a new dated file next to this one,
with the reason, and both sets of results are reported.

Before this rule was written the public API was run once on one simulated data set
of design D4 (below), to confirm that a population estimand fits, draws and values
without error. No coverage, bias or width was computed from it.

## 1. What is recovered

Each cell recovers one economic value built with the public API:
`pollin_data()` -> `pollin_design()` -> `pollin_estimand(scale = "ratio")` ->
`pollin_fit(engine = "glmmTMB", family = glmmTMB::nbinom2())` ->
`pollin_draws(n = 2000L, method = "simulation")` -> `pollin_market()` ->
`pollin_value(conf_level = 0.95)`.

**What the draws target.** `pollin_draws()` forms the arm means as
`pollin_effects()` does: every row of the model data is predicted under the `high`
arm and under the `low` arm with that row's estimated random effects
(`re.form = NULL`), and the predictions are averaged over the rows. The target is
therefore the arms' expected responses **averaged over the realised rows of the
simulated trial, conditional on the realised random effects**, the S1 truth
(`validation/recovery_rule_S1.md`, Section 1).

**Truth per seed.** For simulated data set `s` with rows `i = 1..n`, fixed log-means
`beta[arm]` and realised random effects `u_i`:

    H_s = mean_i exp(beta[high] + u_i)
    L_s = mean_i exp(beta[low]  + u_i)

| Channel | Estimand declared | Ratio truth `R_s` |
|---|---|---|
| `dependence` | `family = "dependence"`, `scale = "ratio"`, `basis = "loss"`; `high` = open, `low` = excluded | `D_s = 1 - L_s / H_s` |
| `deficit` | `family = "deficit"`, `scale = "ratio"`, `basis = "gain"`; `high` = supplemented, `low` = open | `G_s = H_s / L_s - 1` |

The true value is `V_s = R_s x P_s x Q_s`, where `P_s` and `Q_s` are the true price
and quantity of seed `s` (Section 3). For a trial population with random intercepts
only, `R_s` is a constant across seeds (the random intercepts cancel from `H / L`);
for the population cell it is not (Section 2, D4).

**Population truth (D4).** With strata `k`, area weights `w_k` normalised to sum to
1, and `H_sk`, `L_sk` the realised means above taken over the rows of stratum `k`:

    D_s = 1 - sum_k(w_k L_sk) / sum_k(w_k H_sk)

a ratio of population means, as `pollin_population()` defines it; not a weighted
mean of stratum loss shares.

## 2. Designs

Responses are counts simulated from NB2 (variance `mu + mu^2 / theta`), and every
fit passes `family = glmmTMB::nbinom2()` explicitly. Treatment levels are `Cage`
(excluded), `Open` (open) and `Hand` (supplemented). Every unit receives every arm;
no data are missing. The simulating model is the model the layout implies.

D1, D2 and D3 are the S1 designs unchanged (sizes, random-effect SDs, `theta` and
arm means as in `recovery_rule_S1.md`, Section 2). D4 is new.

| | D1 three-arm trial | D2 split-plot trial | D3 matched-pair field study | D4 two-stratum population |
|---|---|---|---|---|
| `layout` | `rcbd` | `split_plot` | `matched_pairs` | `rcbd` |
| `unit` | `block` | `c("block", "wplot")` | `c("site", "pair")` | `block` (block names unique across strata) |
| Structure | 8 blocks x 3 arms x 4 plants | 6 blocks x 4 whole plots x 3 arm sub-plots | 20 sites x 5 pairs x 2 arms | stratum `near`: 10 blocks; stratum `far`: 5 blocks; each block x 2 arms x 4 plants |
| Rows | 96 | 72 | 200 | 120 (80 near, 40 far) |
| Random-effect SDs | block 0.40 | block 0.30, whole plot 0.25 | site 0.50, pair 0.30 | block 0.40 |
| NB2 `theta` | 8 | 15 | 5 | 8 |
| Mean, excluded | 20 | 50 | 150 | near 20, far 20 |
| Mean, open | 30 | 70 | 187.5 | near 32, far 24 |
| Mean, supplemented | 36 | 77 | none | none |
| Channels | dependence, deficit | dependence, deficit | dependence | dependence |
| Ratio truth | D = 1/3; G = 0.2 | D = 2/7; G = 0.1 | D = 0.2 | per seed (stratum D: near 0.375, far 1/6) |

**D4 population.** `pollin_population("band", c(near = 1, far = 3))` (normalised
weights 0.25, 0.75), passed as `population =` to `pollin_estimand()`, fitted with
`pollin_fit(by = "band")` (fixed part `trt * band`, random `(1 | block)`). The design
is unbalanced in two ways: the stratum with three quarters of the area has a third of
the rows, and the strata differ in their loss share, so the ratio of population means
differs from the weighted mean of stratum loss shares. The value is read at
`level = "population"`.

## 3. Market inputs

Point values: `quantity = 12000` (t), `price = 650` (AUD per t), year 2024.

* **Market off.** Scalar `P = 650`, `Q = 12000`; `P_s = P`, `Q_s = Q` for every
  seed, so `V_s = R_s x 7.8e6`.
* **Market on.** The market uncertainty is a stated distribution:
  `P ~ lognormal(log 650, sdlog 0.10)` and `Q ~ lognormal(log 12000, sdlog 0.05)`,
  independent of each other and of the estimand. For each seed the **true** `P_s`
  and `Q_s` are themselves drawn once from that distribution, and
  `price_draws` and `quantity_draws` (2,000 each, one per estimand draw) are drawn
  independently from the same distribution. The truth `V_s = R_s x P_s x Q_s` is
  therefore a draw from the distribution the interval claims to describe, and a
  calibrated 95% interval covers it in 95% of seeds. This is the framing under which
  market uncertainty can be checked by coverage: the market draws express belief
  about an unknown price and quantity, and the planted truth is one realisation of
  that belief.

Market on and market off use the same fit and the same estimand draws of a seed;
only the market differs.

## 4. Grid

| Cell | Design | Channel | Market | Draw method | Seeds |
|---|---|---|---|---|---|
| V01 | D1 | dependence | off | simulation, n = 2000 | 1,000 |
| V02 | D1 | dependence | on | simulation, n = 2000 | 1,000 |
| V03 | D1 | deficit | off | simulation, n = 2000 | 1,000 |
| V04 | D1 | deficit | on | simulation, n = 2000 | 1,000 |
| V05 | D2 | dependence | off | simulation, n = 2000 | 1,000 |
| V06 | D2 | dependence | on | simulation, n = 2000 | 1,000 |
| V07 | D2 | deficit | off | simulation, n = 2000 | 1,000 |
| V08 | D2 | deficit | on | simulation, n = 2000 | 1,000 |
| V09 | D3 | dependence | off | simulation, n = 2000 | 1,000 |
| V10 | D3 | dependence | on | simulation, n = 2000 | 1,000 |
| V11 | D4 (population) | dependence | off | simulation, n = 2000 | 1,000 |
| B01 | D1 | dependence | off | bootstrap, `n = 199` refits | 200 |

V01-V11 are the C2 cells. **B01 is supplementary**: it runs the bootstrap fallback
the plan names, on 200 seeds (the first 200 of D1, the same data sets as V01), and
is reported against the same thresholds but does not enter the C2 verdict, since
C2 asks for 1,000 seeds. At 200 seeds the Monte Carlo SE of a coverage of 0.95 is
0.015, so its coverage is indicative only.

## 5. Seeds and parallelism

Seed `s` in `1..1000` for design `k` (D1 = 1, D2 = 2, D3 = 3, D4 = 4) simulates its
data after `RNGkind("Mersenne-Twister", "Inversion", "Rejection")` and
`set.seed(2000000 + 100000 * k + s)`, set inside the task (a seed range disjoint
from S1). In the same stream, after the data, the task draws `P_s`, `Q_s`, then the
2,000 price draws, then the 2,000 quantity draws. The estimand draws use
`pollin_draws(seed = 1000000 * k + s)`, which restores the stream afterwards. B01
calls `pollin_draws(method = "bootstrap", n = 199L, seed = 5000000 + s)`. Tasks run
with `parallel::mclapply()`; every task sets its own seeds, so results do not depend
on core count or scheduling.

## 6. Per-seed outcomes and failures

For each cell and seed: the value median, `conf.low` and `conf.high` from
`summary(pollin_value())`, the SD of the value draws, the truth `V_s`, and whether
`V_s` lies inside the closed interval `[conf.low, conf.high]` (the **coverage
event**). The ratio median and interval (from the same draws, at the same level) and
whether they cover `R_s` are recorded as a diagnostic.

A seed **fails** in a cell when any call errors or the median, `conf.low`,
`conf.high` or SD is not finite. Failed seeds are excluded and counted (`n_failed`).
Warnings are recorded and counted (`n_warned`), do not make a seed fail. A C2 cell
with more than 20 failed seeds (2%) cannot PASS; B01 cannot PASS with more than 4
(2% of 200).

## 7. Pass rule (the plan's C2, made exact)

Over the `n_ok` successful seeds of a cell:

* **coverage** = share of seeds whose coverage event holds, with its Monte Carlo SE
  `sqrt(c (1 - c) / n_ok)` and Wilson 95% interval;
* **bias** = `mean(median - V_s)`, the median being the point value
  `pollin_value()` reports;
* **SE** for the bias criterion = the **mean over seeds of the SD of the value
  draws**, the uncertainty the package states for its value. The Monte Carlo SD of
  the medians across seeds is reported alongside, not used in the rule;
* **width** = the median over seeds of `conf.high - conf.low`.

A cell **PASSES** when all three hold:

1. coverage in `[0.93, 0.97]` (inclusive, the point estimate; the Monte Carlo SE
   and Wilson interval are reported, not used to pass);
2. `|bias| / SE < 0.1`;
3. `n_failed` within the limit of Section 6.

Otherwise it FAILS, with the criterion it missed named. **C2 is met** when all of
V01-V11 PASS.

## 8. Stated in advance

These are written now so that they cannot be read into a result afterwards. None is
used to excuse a miss.

1. **Simulation draws carry fixed-effect uncertainty only** (as the delta interval
   does). With market off the value is `R x 7.8e6`, a positive multiple of the
   ratio, so its coverage event is the ratio's. S1 found ratio coverage of the delta
   interval between 0.936 and 0.941 on D1-D3; the plan records the risk that the
   simulation draws under-cover on ratio scales. The margin to 0.93 is thin.
2. **Market on widens the interval** with uncertainty the draws state exactly, so
   coverage is expected to move towards 0.95 as the market share of the variance
   grows.
3. **Median against a lognormal truth.** With market on, the median of
   `R x P x Q` lies below the mean of `P_s Q_s` by the factor
   `exp(-(0.10^2 + 0.05^2) / 2)`, about 0.6% of the value. Against an SE of roughly
   20% of the value this is a bias of about 0.03 SE, inside the criterion; it is
   counted as bias like any other, not removed.
4. **D4 has few blocks in its heaviest stratum** (5 blocks carry 75% of the area).
   It is the smallest design here for which C2 is claimed for a population value.
5. **The bootstrap (B01) simulates new random effects** for each refit, so its
   draws include random-effect variation the conditional truth does not; on the
   ratio scale with random intercepts only that variation cancels from `H / L`.

**On a miss.** A cell that misses is reported as FAIL with the criterion, the size
of the miss, and the evidence on its cause (the ratio diagnostic of Section 6 shows
whether the miss is in the estimand draws or in the market propagation). It is not
rerun with more draws, more seeds, a larger design or another method under this
rule, and no package code is changed to move it. A change is a new dated rule and a
new run reported beside this one.

## 9. Outputs

`validation/planted_truth_S2.R` (the suite), `validation/planted_truth_S2_results.csv`
(one row per cell: `n_ok`, `n_failed`, `n_warned`, truth mean, coverage, Monte Carlo
SE, Wilson 95% interval, bias, SE, bias/SE, Monte Carlo SD of the medians, median
width, ratio coverage, verdict), `validation/planted_truth_S2_seeds.csv.gz`
(per-seed outcomes), `validation/planted_truth_S2_session.txt` (session info,
package versions, wall time) and `validation/planted_truth_S2_report.md` (each cell
PASS or FAIL against this rule, misses explained with the evidence).

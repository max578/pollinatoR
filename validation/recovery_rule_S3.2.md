# Planted-truth recovery rule, Stage 3.2 (the joint delta interval, 0.4.0)

Fixed 2026-11-10, before any simulation run of the suite whose results are reported.
This file is committed on its own; the commit time is the evidence that the rule
preceded the results. Nothing here changes after the run. A changed rule is a new
dated amendment next to this one, with the reason, written before the run it
governs, and both sets of results are reported.

It tests what changed between pollinatoR 0.3.1 and 0.4.0 (branch `joint-interval`):
the glmmTMB delta interval of `pollin_effects()` is taken over the fixed and random
effects jointly, with their covariance from the inverse of the joint precision of the
Laplace approximation (`R/joint.R`). The 0.3.1 interval stays available as
`interval = "delta_fixed"`. Everything not restated here is as in
`validation/recovery_rule_S3.md` and `validation/recovery_rule_S3.1.md`.

**What was run before this rule was written.**

1. Exact checks, in `tests/testthat/test-joint.R`: arm designs reproduce
   `predict(type = "link")` to 1e-10; the joint delta estimate equals the
   G-computation of `.arm_means()` to 1e-10; the analytic gradient equals a central
   finite difference to 1e-6; with the variances held fixed, the joint covariance of
   a Gaussian mixed model equals the inverse of Henderson's mixed-model equations
   computed by hand, to 1e-6. All pass.
2. One feasibility probe on one Poisson data set with 12 sites and a random treatment
   slope (plan `PLAN_0.4.0_joint_interval_v0.1.md`): joint SE 0.0275 against
   fixed-only SE 0.0745. No coverage or engine comparison was computed.

**Known before the run, from S3 and S3.1** (on other seeds): on the slope cells the
fixed-only delta interval covers 0.989 to 1.000 and the flexyBayes posterior is 0.53
to 0.63 of its width; on B1 the posterior is 0.756 of its width; the exact posterior
on the 12-site slope design covers about 0.91 to 0.93 (probe S3.1, against NUTS), and
about nominal at 24 sites.

**Plan amendment recorded here (G2, point estimates).** The plan required every
golden estimate to stay identical. It cannot: the 0.3.1 delta estimate came from
`marginaleffects`, which differs from the package's own G-computation by up to
1.5e-5 relative; the joint delta estimate equals the G-computation to 1e-10. The
amended G2, already in the test suite before this rule: the 0.1.0 delta outputs are
reproduced to 1e-10 under `interval = "delta_fixed"`, and the joint delta estimate
equals `.arm_means()` to 1e-10.

## 1. Cells

| Cell | Data (as) | Family | Scale, basis | Slope | Sites / blocks | Seeds |
|---|---|---|---|---|---|---|
| B1 | S3.1 B1: 8 blocks x 6 plants x 2 arms, seeds out of 20 to 40 ovules | binomial | ratio, loss | no | 8 | 200 |
| A2 | S3.1 A2 (control) | Poisson | ratio, loss | no | 8 | 100 |
| A5 | S3.1 A5 (control, three arms) | NB2 | ratio, gain (deficit) | no | 8 | 100 |
| A1 | S3 A1: multi-site, site + site-by-arm SD 0.40 / 0.20, 8 per arm | NB2 | ratio, loss | yes | 12 | 100 |
| A3 | S3 A3: multi-site, Gamma, SD 0.30 / 0.20 | Gamma (log) | ratio, loss | yes | 12 | 100 |
| A1W | A1 with 24 sites | NB2 | ratio, loss | yes | 24 | 200 |

Each cell is fitted on both engines exactly as in S3.1 (`pollin_fit(engine =
"glmmTMB")`, `pollin_fit(engine = "flexyBayes", aggregate = FALSE)`), and the glmmTMB
fit is read three ways: `interval = "delta"` (joint, under test), `"delta_fixed"`
(0.3.1, reported for comparison), and the flexyBayes posterior. The truth is as in
S3 / S3.1: arm means at the realised random effects, averaged over the rows.

## 2. Gates

**G3, no-slope cells (B1, A2, A5).** Coverage of the joint delta interval in
`[0.93, 0.97]`, and the width ratio (median over seeds of flexyBayes interval width /
joint delta width) in `[0.8, 1.25]`.

**G4, slope cells.** At 12 sites (A1, A3): width ratio in `[0.8, 1.25]`; joint delta
coverage reported, not gated. At 24 sites (A1W): joint delta coverage in
`[0.93, 0.97]`; width ratio reported.

**Borderline** as S3: outside a band by no more than 2 Monte Carlo SE is BORDERLINE, a
miss reported with its size; further out is FAIL.

**Fits.** A seed counts when glmmTMB converged (optimiser code 0 and a positive
definite Hessian, as S3.1) and both engines returned a row. On the no-slope cells the
S3.1 failure limit applies (at most `floor(0.02 n)` failed seeds). On the slope cells
glmmTMB stopped at a slope variance on the boundary in 4% to 12% of S3 seeds, a
property of the optimiser that 0.4.0 does not change; there the failure share is
reported and not gated, and the gates are read on the converged seeds. Seeds whose
joint precision was singular and fell back to fixed variance parameters (the warning
of `R/joint.R`) are counted, kept, and reported.

**Also reported, not gated:** coverage and width of `delta_fixed`; misses below and
above; sd(z) and mean(z) with z = (estimate - truth) / SE for each interval; the
interquartile range of the width ratio; median seconds per fit.

**Verdict for 0.4.0.** All six cells passing their gates is the evidence the plan
asks for. A miss is reported with its numbers and the release decision goes to Max,
as at 0.3.0.

## 3. Seeds, draws and parallelism

Seed `s` of cell index `k` (B1 = 1, A2 = 2, A5 = 3, A1 = 4, A3 = 5, A1W = 6) simulates
its data after `RNGkind("Mersenne-Twister", "Inversion", "Rejection")` and
`set.seed(51000000 + 100000 k + s)`; posterior draws use
`seed = 71000000 + 100000 k + s`. Disjoint from every earlier range:

| Suite | Data seeds | Draw seeds |
|---|---|---|
| S1 | 100,000 x id + s | |
| S2 | 2,000,000 + 100,000 x id + s | 1,000,000 x id + s; 5,000,000 + s |
| S3 | 3,000,000 + 100,000 k + s | 23,000,000 + 100,000 k + s |
| S3.1 probe | 900,000 + s; 950,000 + s | 990,000 + s |
| S3.1 | 41,000,000 + 100,000 k + s | 61,000,000 + 100,000 k + s |
| **S3.2** | **51,000,000 + 100,000 k + s** | **71,000,000 + 100,000 k + s** |

1,000 posterior draws per flexyBayes call; INLA `num.threads = "1:1"`,
`flexyBayes.inla_best_of = 1`, package-default priors. Budget about 45 minutes on 8
cores; if the smoke run shows the full run would exceed it, seed counts are cut by a
dated amendment before the full run.

## 4. Outputs

`validation/planted_truth_S3.2.R` (suite), `planted_truth_S3.2_seeds.csv.gz`
(per seed), `planted_truth_S3.2_results.csv` (per cell),
`planted_truth_S3.2_session.txt`, and `planted_truth_S3.2_report.md`.

# Planted-truth recovery rule, Stage 4 (identification, 0.5.0)

Fixed 2026-11-14, before any simulation run of the suite whose results are reported.
This file is committed on its own; the commit time is the evidence that the rule
preceded the results. Nothing here changes after the run. A changed rule is a new
dated amendment next to this one, with the reason, written before the run it
governs, and both sets of results are reported.

It tests what Stage 4 adds on branch `stage4-causal` (plan
`PLAN_stage4_causal_field_v0.1.md`, criteria C3, C4 and C5): total and direct effects
with their labels, moderation with a graph's adjustment set and the negative control,
moderation on an exposure that is constant within sites, the `cage` family with a
sham arm, and `pollin_bound()`. Criteria C1, C2 and C6 are exact tests in the test
suite and are not restated here; C7 (the published oracle) is Step 7.

**What was run before this rule was written.** The Step 0 probes
(`validation/causal_probe_S4.R`, `causal_probe_S4.md`) and the unit tests of Steps
1-5 (suite 1135 pass, 0 fail at d22dc3f). Unit-test fixtures (`helper-field.R`) use
the same model shapes as cells M1 and F1 on fixed seeds and a smaller size; no
coverage, false-call rate or bias over seeds has been computed for any cell below.

## 1. Cells

All responses are simulated on a log link; every random effect is Gaussian.
`NB(mu, 10)` is a negative binomial with mean `mu` and size 10;
`Gamma(mu, 10)` has mean `mu` and shape 10.

| Cell | Layout, size | Response | Plants | Seeds |
|---|---|---|---|---|
| M1a | matched pairs, 24 fields, 5 plants per arm per field | NB | mediator | 1,000 |
| M1b | as M1a with 44 fields | NB | mediator | 1,000 |
| F1 | multi-site, 44 sites, 6 plants per arm per site | NB | equal-shift confounder, no moderation | 1,000 |
| F2 | as F1 | NB | contrast-modifying confounder, moderation | 1,000 |
| F3 | as F1, exposure constant within site | NB | moderation, no confounding | 1,000 |
| X1 | RCBD, 10 blocks, 8 plants per arm per block, three arms | Gamma | sham arm, cage factor | 1,000 (glmmTMB); 200 (flexyBayes) |
| X2 | RCBD, 10 blocks, 8 plants per arm per block, two arms | Gamma | cage factor drawn per seed | 1,000 |

**M1 (mediator; criterion C3).** Field intercept `u ~ N(0, 0.3)`. Vegetative weight
`veg ~ Gamma(m_a, 15)` with `m_a = 40` caged and `40 exp(0.25)` open, independent of
the field. Yield `y ~ NB(exp(3 + u + 0.15 open + 0.012 veg), 10)`. Graph: `trt -> veg
-> yield`, `trt -> yield`. Two estimands on the ratio scale, loss basis: *total*
(data declare no `adjust`; `effect = "total"`) and *direct* (data declare
`adjust = "veg"`; `effect = "direct"`, `mediators = "veg"`). Planted values:
direct `1 - exp(-0.15) = 0.139292`; total `1 - 1/R`,
`R = exp(0.15) ((1 - 0.012 * 40 / 15) / (1 - 0.012 * 40 exp(0.25) / 15))^15`,
the ratio of the arm means over the population of `veg` (the Gamma moment generating
function; the field effect cancels). The script computes `R` from this formula and
reports it to six decimals.

**F1 (equal-shift confounder; C4).** Site fertility `f ~ N(0, 1)`; site bee level
`b_s ~ Poisson(exp(2 + 0.5 f))`; each plant's count `bees = b_s + Poisson(2)`, so the
exposure varies within and between sites. Site intercept `u ~ N(0, 0.3)`; site
contrast `tau_s = 0.3 + N(0, 0.15)`; yield `NB(exp(3 + u + 0.4 f + tau_s open), 10)`.
Fertility shifts both arms alike and is not given to the model. No moderation is
planted: the true within and between treatment-by-bees terms are 0. The negative
control's planted strength is the fertility shift, 0.4 on the log scale per SD of `f`.
Fit: `pollin_moderation()` with no adjustment.

**F2 (contrast-modifying confounder; C4).** Site variable `z ~ N(0, 1)`; site bee
level `b_s ~ Poisson(exp(2 + 0.5 z))`, plant counts as F1. Let `zb_s` be the
standardised between-unit component of `bees` exactly as `pollin_moderation()`
forms it (the within-between split by site, then centred and scaled over rows). Site
contrast `tau_s = 0.3 + 0.10 zb_s + 0.20 z + N(0, 0.10)`; yield
`NB(exp(3 + u + 0.3 z + tau_s open), 10)`, `u ~ N(0, 0.3)`. Planted moderation (the
`trt:bees_between` coefficient) is 0.10; the within term is 0. Graph:
`z -> bees`, `z -> yield`, `bees -> yield`, `trt -> yield`; the adjustment set of
`bees` is `{z}`. Two fits: with `dag` (set `{z}` as main effect and treatment
interaction) and with no adjustment. The planted bias of the second, per seed, is
`0.20 b_z`, where `b_z` is the coefficient of `zb` in the site-level least-squares
regression of `z` on `zb` (one row per site, with intercept).

**F3 (site-level exposure; C4).** Site bee level `b_s ~ Poisson(exp(2))`, and
`bees = b_s` on every plant of the site (no within variation). `zb_s` as in F2.
Site contrast `tau_s = 0.3 + 0.10 zb_s + N(0, 0.10)`; yield
`NB(exp(3 + u + tau_s open), 10)`, `u ~ N(0, 0.3)`. Planted between term 0.10. Fit:
`pollin_moderation()` with no adjustment; it must report `bees` as between units only.

**X1 (sham arm; C5).** Block intercept `u ~ N(0, 0.3)`; open-arm mean `exp(1 + u)`;
sham `1.10 exp(1 + u)`; excluded `0.75 * 1.10 exp(1 + u)`. The cage multiplies the
yield by `k = 1.10` with and without pollinators; the dependence net of the cage is
0.25. Two estimands: the `cage` family on the ratio scale (planted `k = 1.10`), and
the dependence family, loss basis, with `reference = "sham"` (planted 0.25). Both on
glmmTMB (1,000 seeds); the cage estimand also on flexyBayes (200 seeds), to show the
engine route carries the family.

**X2 (bounds; C5).** As X1 without the sham arm: open `exp(1 + u)`, excluded
`0.75 k exp(1 + u)`, with `k ~ Uniform(0.95, 1.25)` drawn per seed. The observed loss
share is then `1 - 0.75 k`; the share net of the cage is 0.25. Fit the dependence
estimand (ratio, loss), take 1,000 draws, and call `pollin_bound(cage_factor = c(1,
1.2))`. Seeds with `k` in `[1, 1.2]` are *inside*; the others are *outside* and serve
as the positive control (the bound should then miss more often).

The truth in every cell is a constant (the multiplicative block, field and site
effects cancel in each ratio), except F2's per-seed planted bias.

## 2. Gates

Coverage is of the package's default 95% interval: the joint delta interval of
`pollin_effects()` for M1 and X1 on glmmTMB; the posterior interval on flexyBayes;
`estimate +/- qnorm(0.975) std.error` for a moderation term; the outer interval of
`pollin_bound()` for X2.

| Gate | Cell | Requirement |
|---|---|---|
| G-M1 | M1a, M1b | total estimate covers the planted total, and direct covers the planted direct, each in `[0.93, 0.97]`; every row is labelled `total` or `direct` as fitted (all seeds) |
| G-F1 | F1 | false-call rate (p < 0.05) of the `trt:bees_within` term and of the `trt:bees_between` term, each at most 0.07; the negative control's `bees_between` term rejects (p < 0.05) in at least 0.80 of seeds |
| G-F2a | F2 | with the graph's set, the `trt:bees_between` interval covers 0.10 in `[0.93, 0.97]` |
| G-F2b | F2 | without the set, `mean(est - 0.10) - mean(0.20 b_z)` is within 2 Monte Carlo SE of 0, the SE being `sd(est - 0.10 - 0.20 b_z) / sqrt(n)`; and `mean(est - 0.10)` is more than 2 Monte Carlo SE from 0 (the bias is there to detect) |
| G-F3 | F3 | no seed returns a missing estimate or standard error for `trt:bees_between`, and every seed names `bees` in `between_only`; the term's interval covers 0.10 in `[0.93, 0.97]` |
| G-X1 | X1, glmmTMB | the cage-factor interval covers 1.10 and the net dependence interval covers 0.25, each in `[0.93, 0.97]` |
| G-X1b | X1, flexyBayes | the cage estimand returns a row in at least 0.98 of seeds; coverage and the width ratio against glmmTMB are reported, not gated (S3 found flexyBayes intervals at 0.53 to 0.76 of the delta width) |
| G-X2a | X2 | algebra: in every inside seed, `pollin_bound()` run on that seed's draws object with `low` replaced by the planted `0.75 k` and `high` by 1 returns `lower <= 0.25 <= upper` (100%) |
| G-X2b | X2 | among inside seeds, the outer interval covers 0.25 in at least 0.95 |

**Borderline.** As S3: outside a band by no more than 2 Monte Carlo SE is BORDERLINE,
a miss reported with its size; further out is FAIL. For a one-sided requirement
(at most 0.07, at least 0.80, at least 0.95) the same 2 SE margin applies on its one
side.

**Fits.** A seed counts when glmmTMB converged (optimiser code 0 and a
positive-definite Hessian) and the call returned its row. On M1, X1 and X2 (no random
slope) at most `floor(0.02 n)` seeds may fail. On F1-F3 (`multi_site`, random
treatment slope) the failure share is reported and not gated, and the gates are read
on converged seeds, as in S3.2. Warnings are counted by message and reported.

**Also reported, not gated:** bias and median width of every gated interval; misses
below and above; the total and direct estimates side by side in M1 (mean over seeds);
F1's negative-control rejection for the within term; F2's mean `b_z`; X2's coverage
and the share of seeds whose median bounds contain 0.25, inside and outside; median
seconds per fit.

**Verdict for 0.5.0.** All gates passing is the evidence the plan asks for. A miss is
reported with its numbers and the release decision goes to Max, as at 0.3.0.

## 3. Seeds, draws and parallelism

Seed `s` of cell index `k` (M1a = 1, M1b = 2, F1 = 3, F2 = 4, F3 = 5, X1 = 6, X2 = 7;
the X1 flexyBayes run reuses X1's data seeds) simulates its data after
`RNGkind("Mersenne-Twister", "Inversion", "Rejection")` and
`set.seed(81000000 + 100000 k + s)`; draws (X2's `pollin_draws()` and X1's posterior)
use `seed = 91000000 + 100000 k + s`. Disjoint from every earlier range:

| Suite | Data seeds | Draw seeds |
|---|---|---|
| S1 | 100,000 x id + s | |
| S2 | 2,000,000 + 100,000 x id + s | 1,000,000 x id + s; 5,000,000 + s |
| S3 | 3,000,000 + 100,000 k + s | 23,000,000 + 100,000 k + s |
| S3.1 probe | 900,000 + s; 950,000 + s | 990,000 + s |
| S3.1 | 41,000,000 + 100,000 k + s | 61,000,000 + 100,000 k + s |
| S3.2 | 51,000,000 + 100,000 k + s | 71,000,000 + 100,000 k + s |
| **S4** | **81,000,000 + 100,000 k + s** | **91,000,000 + 100,000 k + s** |

1,000 draws per `pollin_draws()` call and per flexyBayes posterior; INLA
`num.threads = "1:1"`, `flexyBayes.inla_best_of = 1`, package-default priors.
Budget about 60 minutes on 8 cores. A smoke run at a small fraction of the seeds goes
to a scratch directory and is not reported; if it shows the full run would exceed the
budget, seed counts are cut by a dated amendment before the full run.

## 4. Outputs

`validation/planted_truth_S4.R` (suite), `planted_truth_S4_seeds.csv.gz` (per seed),
`planted_truth_S4_results.csv` (per gate), `planted_truth_S4_session.txt`, and
`planted_truth_S4_report.md`.

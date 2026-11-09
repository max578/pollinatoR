# Probe S3.1: flexyBayes posterior coverage on site-slope designs

2026-11-09. Run on pollinatoR 0.3.0 (the `v0.3.0` tag), flexyBayes 0.19.0, INLA
26.08.07, cmdstan 2.38.0. Synthetic data only. This probe is exploratory: it was not
pre-registered, its seeds are not evaluation seeds, and it is not a recovery rule.
It informs `validation/recovery_rule_S3.1.md`, which is.

Scripts, run from the package root:

* `validation/probe_S3.1/probe.R`: the public route under each setting;
* `validation/probe_S3.1/fixed_hyper.R`: INLA with every hyperparameter fixed at its
  planted value;
* `validation/probe_S3.1/fixed_hyper_nuts.R`: the same fixed-hyperparameter model
  sampled by NUTS (cmdstanr), the exact reference.

`probe.R` loads the package from the working tree; to reproduce the numbers below,
run it on a checkout of `v0.3.0`. The per-seed outputs are not kept in the repository;
they are summarised below.

## Question

S3 (`validation/planted_truth_S3_report.md`) found the flexyBayes 95% posterior
interval under-covering the planted trial-conditional loss share on the site-slope
cells (0.909 to 0.925) and on the binomial cell A4 (0.925). The probe asked which of
three causes is responsible, and whether one default setting on the flexyBayes route
would remove it:

* **H1.** INLA's integration over the hyperparameters is too coarse, or the
  hyperparameters are held at their mode.
* **H2.** The default prior on the variance components is too tight.
* **H3.** The posterior is correct, and the shortfall is the frequentist coverage of
  a correct Bayesian interval for this estimand on a 12-site design.

## Settings and seeds

Data as `validation/planted_truth_S3.R`. A1: 12 sites x 2 arms x 8 plants, NB2 size 8,
site SD 0.40, site-by-arm SD 0.20, means 20 (Cage) and 30 (Open). A4: 8 blocks x 2 arms
x 24 flowers, Bernoulli 0.55 and 0.75, block SD 0.40. Route:
`pollin_fit(engine = "flexyBayes", family, aggregate = FALSE, ...)` then
`pollin_effects(interval = "posterior", n_draws = 1000, conf_level = 0.95)`, with
`options(flexyBayes.inla_best_of = 1)` and INLA `num.threads = "1:1"`.

**Seeds** (disjoint from S3, whose data seeds are 3,000,000 + 100,000 k + s and draw
seeds 23,000,000 + 100,000 k + s): data seed 900,000 + s for every A1-type data set
and 950,000 + s for A4; draw seed 990,000 + s; s = 1 to 100 for the route settings
and 1 to 200 for the two fixed-hyperparameter checks. Every setting uses the same
data and draw seeds, so the settings are paired seed by seed.

Besides coverage, the probe reports sd(z), with z = (posterior median - truth) /
posterior SD. For a calibrated interval it is near 1, and at 100 seeds its SE is
about 0.07, which resolves width more finely than coverage does (Monte Carlo SE
0.022 to 0.029).

**How settings reached the fit.** `pollin_fit()` passes its `...` to `flexybayes()`,
so `prior =` and `prior_vc_sd =` are public routes. `flexybayes()` has no argument
that reaches INLA's `control.inla`, so the H1 settings were injected by replacing
flexyBayes's internal INLA call site (`.inla_call`) inside the probe process
(`utils::assignInNamespace()`, taken from `asNamespace("flexyBayes")`). `probe.R` also
reads the INLA object inside the fit (`fit$fits[[1]]$inla`) to record the number of
integration points and the posterior of the cell SD. Both reach into flexyBayes
internals; validation scripts are not shipped with the package, and neither is
available to a user. A prior on the site-by-arm term has to name pollinatoR's internal
column `pollin_site_trt`.

## Results

A1 design, 100 seeds unless stated. "Discordant" counts seeds whose coverage event
differs from the default's.

| Setting | Route | Seeds | Coverage | MC SE | Median width | Width / default | sd(z) | Misses below / above | Discordant |
|---|---|---|---|---|---|---|---|---|---|
| Default | public | 100 | 0.910 | 0.029 | 0.1566 | 1 | 1.12 | 5 / 4 | - |
| H1: `int.strategy = "eb"` (hyperparameters at the mode) | injected | 100 | 0.910 | 0.029 | 0.1561 | 0.990 | 1.14 | 5 / 4 | 0 |
| H1: finer grid (`dz = 0.25`, `diff.logdens = 10`) | injected | 100 | 0.910 | 0.029 | 0.1566 | 1.000 | 1.12 | 5 / 4 | 0 |
| H2: `half_normal(scale = 1)` on both SDs | `prior =` | 100 | 0.910 | 0.029 | 0.1559 | 1.000 | 1.12 | 5 / 4 | 0 |
| H2: `pc(upper = 1, prob = 0.01)` on both SDs | `prior =` | 100 | 0.910 | 0.029 | 0.1556 | 0.996 | 1.12 | 5 / 4 | 0 |
| H2: `lognormal(0, 1)` on the SD (boundary-avoiding) | `prior_vc_sd = 1` | 100 | 0.910 | 0.029 | 0.1558 | 0.999 | 1.13 | 5 / 4 | 0 |
| H3: default, 24 sites | public | 100 | 0.910 | 0.029 | 0.1133 | - | 1.00 | 4 / 5 | - |
| Hyperparameters fixed at planted values, INLA | direct INLA | 200 | 0.930 | 0.018 | 0.1571 | 0.997 (seeds 1 to 100) | 1.11 | 7 / 7 | - |
| Same, sampler skewness correction off | direct INLA | 200 | 0.930 | 0.018 | 0.1571 | - | 1.11 | 7 / 7 | - |
| Same model by NUTS (exact reference; 2 x 1,000 draws, R-hat at most 1.007, no divergences) | Stan | 200 | 0.930 | 0.018 | 0.1593 | INLA / NUTS 0.992 | 1.10 | 8 / 6 | 4 of 200 against fixed INLA |

"Width / default" is the paired median ratio. The posterior median of the
site-by-arm SD was 0.208 under the default prior, 0.191 under the PC prior and 0.222
under the lognormal, against 0.20 planted. On every 12-site setting, INLA used 15
integration points (`auto` on INLA 26.08, identical to `grid` here).

A4 design, 100 seeds, default route: coverage **0.960** (MC SE 0.020), median width
0.2183, sd(z) 1.00, misses 2 below and 2 above.

## Readings

* **H1: rejected.** The posterior sampler draws the hyperparameters across the
  15 stored integration points, so they are not fixed at the mode. Fixing them at the
  mode, or refining the grid, changed no coverage event and moved the width by at
  most 1%.
* **H2: rejected.** The default variance-component prior, `half_normal(scale = 3)` on
  each random-effect SD for log links, is wide. A tighter half-normal, a PC prior and
  a boundary-avoiding lognormal changed no coverage event (widths within 0.4%).
* **H3: supported, in a corrected form.** With every hyperparameter known, the exact
  NUTS posterior covers 0.930 with sd(z) 1.10, and INLA gives the same intervals
  (width ratio 0.992, coverage 0.930). The flexyBayes/INLA posterior is therefore
  the correct posterior, and no prior or integration setting can raise coverage on
  this design: the shortfall lies in the exact posterior itself. It is a property of
  the Bayesian interval for the trial-conditional estimand on a 12-site design, with
  the fixed effects held fixed across seeds rather than drawn from their prior, where
  Bayesian calibration is not guaranteed. The pooled evidence puts the route's
  coverage on A1 at about 0.91 to 0.93 (S3 0.909 on 186 seeds, probe 0.910 on 100,
  exact reference 0.930 on 200), with intervals about 5% to 10% short of calibrated
  width (sd(z) 1.06 to 1.12). The form "coverage rises with the number of sites" was
  not shown cleanly: at 24 sites coverage was again 0.910, but sd(z) was 1.00, which
  reads that 0.910 as Monte Carlo noise (1.4 MC SE).
* **A4.** The S3 BORDERLINE (0.925) did not reproduce on independent seeds (0.960,
  sd(z) 1.00), in line with the S3 report's note that 12 of its A4 seeds lay within
  the INLA run-to-run variation of an interval end.

**Implication for the default route.** No setting is passed by default: no public
setting changes the result, and the one that changes runtime (`eb`, about three
times faster) cannot reach INLA through `flexybayes()`. The C3 band `[0.93, 0.97]` is
tighter than the exact posterior attains on this design (0.930, MC SE 0.018).

Re-reading the S3 per-seed file the same way, over every seed on which flexyBayes
succeeded (200 or 100 per cell), gives sd(z) of 1.06 (A1), 1.04 (A3),
1.03 (A6), 0.98 (P8), 0.99 (P12) and 1.12 (A4), against 0.96 to 0.99 on A2 and A5. A3
also shows a mean z of 0.23, a location lean on the Gamma cell that the probe did
not test.

## Limits

* **48-site run cut.** The 48-site setting was stopped before it wrote any output, to
  keep total compute near the 60-minute budget (about 50 minutes, with other jobs on
  the machine). The trend of coverage with the number of sites therefore rests on
  12 and 24 sites only.
* **A3 and A6 not probed.** The Gamma slope cell (A3) and the difference-scale cell
  on A1 data (A6) were not run. The A3 location lean is a separate question.
* The H1 settings were injected through a flexyBayes internal, which a user cannot
  do; a change to `eb` would need a flexyBayes change.
* 100 seeds resolve a single coverage only to about 0.06 (95%). The between-setting
  comparison is conclusive because the settings are paired (0 discordant seeds out of
  100 across all 12-site settings), not because of the coverage levels.
* The plain `int.strategy = "grid"` setting was dropped after the smoke run because
  on INLA 26.08 it is identical to `auto`; the finer grid replaced it.

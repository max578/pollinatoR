# flexyBayes probe for Stage 3 (Step 0)

2026-11-09. flexyBayes 0.19.0, INLA 26.08.07, extension API version 1. Script:
`validation/flexybayes_probe_S3.R` (synthetic data; eight sites, two arms, ten plants
per arm and site). Only `flexybayes()` and the extension API v1 functions are called.

## Fit and timing

`flexybayes(pods ~ trt, random = ~ Site, family = "negative_binomial", backend = "inla")`
fits; `fb_convergence()` returns `engine = "inla"`, `ok = TRUE`, no reasons, with
evidence `mode_status`, `mlik` and `numerical_confirm`. The sampler fields
(`max_rhat`, `min_ess_bulk`, ...) are `NA` on INLA.

| Setting | Wall time per fit (160 rows) |
|---|---|
| Default, best of three fits (one in session, two in fresh processes) | 4.2 to 4.7 s |
| `options(flexyBayes.inla_best_of = 1L)` | 2.1 to 2.7 s |

The best-of-three fits agreed to 5.7e-5 in log marginal likelihood on this model.

## Latent rows of `fb_inla_posterior_sample()`

Each draw holds `latent` (a one-column matrix) and `hyperpar` (a named vector).

| Term | Latent rows | Level map |
|---|---|---|
| Linear predictor | `Predictor:1` to `Predictor:<n>`, one per data row, in data order | row order |
| Fixed effect | `<column>:1`, the column names of `fb_fixed_design(~ fixed, data)$X`: `(Intercept):1`, `trtOpen:1`, `Regionsouth:1`, `trtOpen:Regionsouth:1` | column name |
| Simple random term `~ Site` | `Site:1` to `Site:8` | `fb_model_terms(fit)$random_terms[[i]]$var_levels` |
| Nested term `~ Site:trt` | `Site_trt_id:1` to `Site_trt_id:16` | none: `var_levels` is `NULL` |
| `~ diag(trt):Site` | `Site_at_trt_Cage:k`, `Site_at_trt_Open:k` | `inner_levels` |

Hyperparameters on the negative binomial fit: `size for the nbinomial observations
(1/overdispersion)` and `Precision for <term>`.

A nested term carries no level map in the extension API, so the engine route enters
every multi-column grouping (a nested unit, or a unit-by-arm cell) as one combined
factor column of its own, a simple term whose `var_levels` give the map. With the
fixed part `trt * Region` and simple terms for site, site-by-arm cell and
site-by-plot, on unbalanced data (seven rows dropped), `X beta + Z u` rebuilt from the
latent rows equals the `Predictor` rows to 8.9e-16 over 50 draws. The rebuild is
therefore checked at run time on every draw.

## Random treatment slope by site

| ASReml syntax | INLA |
|---|---|
| `~ us(trt):Site` (the correlated slope that glmmTMB's `(trt \| Site)` fits) | refused: `correlated_covariance_requires_brms`, "flexyBayes cannot fit this model on INLA." |
| `(trt \| Site)` in the fixed formula | refused: "The bar grammar has no lowering for a random slope on a factor." |
| `~ diag(trt):Site` | fits; one independent site effect per arm, no site intercept |
| `~ Site + Site:trt` | fits; a site intercept plus an independent effect per site-arm cell |

The slope on INLA is therefore entered as `Site + Site:trt`: the treatment effect at
site `s` departs from the mean by `v[s, high] - v[s, low]`, with variance twice the
cell variance, and the shared site effect cancels from it. It is the exchangeable
form of the unstructured slope (equal arm variances, non-negative correlation
between arms). `diag(trt):Site` sets the correlation between arms to zero, so the
site main effect would enter each per-site contrast; it is not used.

## Seeds

| Call | Two calls identical |
|---|---|
| `fb_inla_posterior_sample(fit, 20, seed = 42, num.threads = "1:1")` twice | no |
| `set.seed(9)` then the same call, twice | yes, bit for bit |

The sample also moves R's random number state: INLA picks the hyperparameter
configuration of each draw with R's generator and samples the latent field with its
own. Reproducible draws need both `set.seed()` and `seed =`; the engine route sets
both from the user's seed and restores the session's state afterwards.

In the probe environment `INLA::inla.getOption("num.threads")` reads `"NA:1"` and a
sample without `num.threads` stops with "missing value where TRUE/FALSE needed"; the
engine route always passes `num.threads = "1:1"`.

## Family map

| glmmTMB family in `pollin_fit()` | flexyBayes string | Link (from `fb_model_terms()`) |
|---|---|---|
| `nbinom2()` | `negative_binomial` | log |
| `poisson()` | `poisson` | log |
| `binomial()` | `binomial` | logit |
| `Gamma(link = "log")` | `gamma` | log |
| `beta_family()` | `beta` | logit |

All five fit and converge on INLA. The glmmTMB family's link must equal the
flexyBayes link; another link is refused. `tweedie` is refused by flexyBayes
("Unsupported family \"tweedie\"") and by the engine route before fitting, as is
every other glmmTMB family. The glmmTMB spelling `nbinom2` is not a flexyBayes family
name.

# Draws probe, Stage 2 (Step 0)

2026-11-09. Question: can `marginaleffects::get_draws()` supply the simulation draws for
`pollin_draws(method = "simulation")`? It is used only if it reproduces a hand-written
multivariate-normal draw of the fixed effects on the same `glmmTMB` fit.

Versions: R 4.6.1, `marginaleffects` 0.32.0, `glmmTMB` 1.1.14. Fit: a simulated
randomised-block trial (6 blocks x 2 arms x 4 plants, Poisson counts with a block effect),
`pollin_fit()` with the chosen family (negative binomial), dependence on the lnRR scale.

## Verdict

**`get_draws()` cannot serve the simulation path.** Its formals are `(x, shape)`: it only
extracts draws already held by a `marginaleffects` object. The draws must come from
`inferences(method = "simulation")`, and on a `glmmTMB` fit that call fails
(`mean and sigma have non-conforming size`): `get_coef()` returns the 3 conditional fixed
effects while `get_vcov()` returns a 4 x 4 matrix that includes the dispersion parameter.
`inferences(method = "boot")` also fails on this fit ("Failed to refit the model").
`marginaleffects` simulation would further depend on `mvtnorm`, which is not an import.

The simulation path is therefore written by hand: fixed effects drawn as
`beta_hat + t(L) z`, with `L = chol(vcov(fit, full = FALSE)$cond)` and `z` from
`stats::rnorm()`. Each draw is pushed through the same prediction rule as `.arm_means()`
(`predict(type = "response")` with the `glmmTMB` default `re.form = NULL`, so every row keeps
its estimated random effects): the random-effect part of the linear predictor,
`predict(type = "link") - X beta_hat`, is held fixed and only `X beta` changes per draw.
Checked: at `beta_hat` this rule reproduces `.arm_means()` (`all.equal` TRUE), and the
random-effect offset is the same under both arms.

Side observation, recorded and not acted on: `avg_predictions()` gives arm means that differ
from `.arm_means()` by about 2e-5 (relative 1e-6). The delta interval in `pollin_effects()`
is built on `marginaleffects` and the bootstrap and draws paths on `.arm_means()`; the gap
is far below any reported precision.

With 4,000 hand draws, the 95% percentile interval of `log(rr)` is (0.25457, 0.48304)
against the delta interval (0.2569, 0.48154). The formal agreement gate on large `n` is
fixed in `tests/testthat/test-draws.R`.

## Commands

Run from the package root: `Rscript probe.R`, where `probe.R` is:

```r
pkgload::load_all(".", quiet = TRUE)
set.seed(5)
d <- expand.grid(plant = 1:4, trt = c("Cage", "Open"), Block = sprintf("b%d", 1:6), stringsAsFactors = FALSE)
blk <- as.integer(factor(d$Block))
d$pods <- rpois(nrow(d), exp(3 + 0.4 * (d$trt == "Open") + rnorm(6, 0, 0.2)[blk]))
d$plant <- NULL
pd <- pollin_data(d, "pods", "trt", unit = "Block")
ds <- pollin_design(pd, c(excluded = "Cage", open = "Open"), "rcbd")
fit <- pollin_fit(pollin_estimand(ds, "dependence", "lnrr", "loss"))
m <- fit$fits$pods; dd <- fit$data$pods
cat("get_draws formals:", paste(names(formals(marginaleffects::get_draws)), collapse = ", "), "\n")
cat("inferences formals:", paste(names(formals(marginaleffects::inferences)), collapse = ", "), "\n")
# Point arm means, package rule
am <- .arm_means(m, dd, "trt", "Open", "Cage", NULL)
cat("arm means (.arm_means):", am$high, am$low, "\n")
# marginaleffects simulation draws of the arm means
p <- suppressWarnings(marginaleffects::avg_predictions(m, variables = list(trt = c("Cage", "Open")), by = "trt"))
print(as.data.frame(p)[, c("trt", "estimate")])
R <- 4000L
set.seed(1)
t0 <- Sys.time()
ps <- tryCatch(suppressWarnings(marginaleffects::inferences(p, method = "simulation", R = R)), error = function(e) e)
if (inherits(ps, "error")) cat("inferences(simulation) ERROR:", conditionMessage(ps), "\n")
cat("length(coef(m)):", length(stats::coef(m)), " length(get_coef):", length(marginaleffects::get_coef(m)), " dim(get_vcov):", dim(marginaleffects::get_vcov(m, vcov = TRUE)), "\n")
gdp <- tryCatch(marginaleffects::get_draws(p), error = function(e) e)
cat("get_draws on a delta object:", if (inherits(gdp, "error")) conditionMessage(gdp) else paste(class(gdp), collapse = "/"), "\n")
pb <- tryCatch(suppressWarnings(marginaleffects::inferences(p, method = "boot", R = 50)), error = function(e) e)
cat("inferences(boot):", if (inherits(pb, "error")) conditionMessage(pb) else "ran", "\n")
cat("inferences(simulation) time:", format(Sys.time() - t0), "\n")
gd <- if (inherits(ps, "error")) NULL else marginaleffects::get_draws(ps, shape = "DxP")
cat("arm means avg_predictions minus .arm_means:", signif(p$estimate[order(p$trt)] - c(am$low, am$high), 3), "\n")

# Hand MVN draw: Cholesky of vcov(full = FALSE)$cond, stats::rnorm only
V <- vcov(m, full = FALSE)$cond
b <- glmmTMB::fixef(m)$cond
L <- chol(V)
set.seed(1)
Z <- matrix(rnorm(R * length(b)), R)
B <- sweep(Z %*% L, 2L, b, "+")
tt <- delete.response(terms(m))
nd_h <- dd; nd_h$trt <- factor("Open", levels = levels(dd$trt))
nd_l <- dd; nd_l$trt <- factor("Cage", levels = levels(dd$trt))
X_h <- model.matrix(tt, nd_h); X_l <- model.matrix(tt, nd_l)
re_h <- predict(m, newdata = nd_h, type = "link") - drop(X_h %*% b)
re_l <- predict(m, newdata = nd_l, type = "link") - drop(X_l %*% b)
linv <- family(m)$linkinv
hand_h <- colMeans(linv(X_h %*% t(B) + re_h)); hand_l <- colMeans(linv(X_l %*% t(B) + re_l))
cat("hand rule at beta_hat reproduces .arm_means:",
    all.equal(c(mean(linv(drop(X_h %*% b) + re_h)), mean(linv(drop(X_l %*% b) + re_l))), c(am$high[[1]], am$low[[1]])), "\n")
q <- c(0.025, 0.5, 0.975)

cat("log(rr) quantiles, hand     :", round(quantile(log(hand_h / hand_l), q), 5), "\n")
cat("arm-mean SDs hand      (Cage, Open):", round(c(sd(hand_l), sd(hand_h)), 4), "\n")
dl <- suppressWarnings(pollin_effects(fit))$table
cat("delta 95% for log(rr):", round(c(dl$conf.low, dl$conf.high), 5), "\n")
cat("hand log(rr) median, 95% (n=4000):", round(quantile(log(hand_h / hand_l), q), 5), "\n")
cat("chol check t(L) %*% L == V:", isTRUE(all.equal(crossprod(L), V)), "\n")
cat("predict(link) - X b constant across arms (RE offset):", isTRUE(all.equal(unname(re_h), unname(re_l))), "\n")
```

## Output

```
get_draws formals: x, shape 
inferences formals: x, method, R, conf_type, data_train, data_test, data_calib, conformal_score, estimator, ... 
arm means (.arm_means): 31.03659 21.45482 
   trt estimate
1 Cage 21.45480
2 Open 31.03656
inferences(simulation) ERROR: mean and sigma have non-conforming size 
length(coef(m)): 3  length(get_coef): 3  dim(get_vcov): 4 4 
This object does not include draws. The `posterior_draws` function only supports bayesian models produced by the `marginaleffects` or `predictions` functions of the `marginaleffects` package. 
get_draws on a delta object: NULL 
inferences(boot): Failed to refit the model. 
arm means avg_predictions minus .arm_means: -1.87e-05 -2.71e-05 
hand rule at beta_hat reproduces .arm_means: TRUE 
log(rr) quantiles, hand     : 0.25457 0.36797 0.48304 
arm-mean SDs hand      (Cage, Open): 2.4929 3.4729 
delta 95% for log(rr): 0.2569 0.48154 
hand log(rr) median, 95% (n=4000): 0.25457 0.36797 0.48304 
chol check t(L) %*% L == V: TRUE 
predict(link) - X b constant across arms (RE offset): TRUE 
```

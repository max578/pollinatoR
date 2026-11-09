# flexybayes_probe_S3.R -- Stage 3 Step 0 probe of the flexyBayes extension API.
#
# Run from the package root:
#   Rscript --vanilla validation/flexybayes_probe_S3.R
# Reads only flexyBayes's exported extension functions and `flexybayes()`.
# Findings are written up in validation/flexybayes_probe_S3.md.

suppressPackageStartupMessages(library(flexyBayes))
options(flexyBayes.silence_default_prior_note = TRUE)
cat(
  "flexyBayes", as.character(utils::packageVersion("flexyBayes")),
  "| INLA", as.character(utils::packageVersion("INLA")),
  "| extension API", fb_extension_api(), "\n"
)

set.seed(1)
d <- data.frame(
  trt = rep(rep(c("Cage", "Open"), each = 10L), 8L),
  Site = rep(paste0("s", 1:8), each = 20L),
  Plot = rep(rep(c("p1", "p2"), each = 5L), 16L)
)
d$pods <- stats::rnbinom(
  nrow(d),
  mu = exp(3 + 0.3 * (d$trt == "Open") + stats::rnorm(8L, 0, 0.3)[
    as.integer(factor(d$Site))
  ]),
  size = 10
)
d$Region <- ifelse(d$Site %in% paste0("s", 1:4), "north", "south")
d$trt <- factor(d$trt)

fit_inla <- function(fixed, random, data = d, family = "negative_binomial") {
  t0 <- Sys.time()
  f <- tryCatch(
    flexybayes(
      fixed,
      random = random,
      data = data,
      family = family,
      backend = "inla",
      verbose = FALSE
    ),
    error = function(e) e
  )
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(
    "\n--", deparse(fixed), "| random", deparse(random), "|", family,
    sprintf("| %.1f s\n", secs)
  )
  if (inherits(f, "error")) {
    cat("refused:", strsplit(conditionMessage(f), "\n")[[1L]][1L], sep = "\n  ")
    return(invisible(NULL))
  }
  f
}

latent_terms <- function(f) {
  set.seed(1)
  s <- fb_inla_posterior_sample(f, 1L, seed = 1L, num.threads = "1:1")
  rn <- rownames(s[[1L]]$latent)
  print(table(sub(":[0-9]+$", "", rn)))
  cat("hyperparameters:", names(s[[1L]]$hyperpar), sep = "\n  ")
}

# 1. Two-arm negative binomial with a site effect, timing with best-of-3 and 1 -
f3 <- fit_inla(pods ~ trt, ~Site)
print(fb_convergence(f3)[c("engine", "ok", "reasons")])
options(flexyBayes.inla_best_of = 1L)
f1 <- fit_inla(pods ~ trt, ~Site)
latent_terms(f1)

# 2. Seed behaviour -----------------------------------------------------------
same <- function(a, b) {
  identical(lapply(a, `[[`, "latent"), lapply(b, `[[`, "latent"))
}
a <- fb_inla_posterior_sample(f1, 20L, seed = 42L, num.threads = "1:1")
b <- fb_inla_posterior_sample(f1, 20L, seed = 42L, num.threads = "1:1")
cat("\nseed = 42 alone, identical:", same(a, b), "\n")
set.seed(9)
a <- fb_inla_posterior_sample(f1, 20L, seed = 42L, num.threads = "1:1")
set.seed(9)
b <- fb_inla_posterior_sample(f1, 20L, seed = 42L, num.threads = "1:1")
cat("set.seed(9) + seed = 42, identical:", same(a, b), "\n")
cat("R RNG state moved by a sample:", {
  set.seed(9)
  r0 <- .Random.seed
  invisible(fb_inla_posterior_sample(f1, 2L, seed = 42L, num.threads = "1:1"))
  !identical(r0, .Random.seed)
}, "\n")
cat("INLA num.threads option here:", format(INLA::inla.getOption("num.threads")), "\n")

# 3. Treatment slope by site in ASReml syntax ----------------------------------
fit_inla(pods ~ trt, ~ us(trt):Site)
f <- fit_inla(pods ~ trt, ~ diag(trt):Site)
latent_terms(f)
f <- fit_inla(pods ~ trt, ~ Site + Site:trt)
latent_terms(f)
str(lapply(fb_model_terms(f)$random_terms, `[[`, "var_levels"))

# 4. Combined grouping columns: simple terms with their levels ------------------
d2 <- d
d2$pollin_Site_trt <- interaction(d2$Site, d2$trt, drop = TRUE, sep = ":")
d2$pollin_Site_Plot <- interaction(d2$Site, d2$Plot, drop = TRUE, sep = ":")
d2 <- d2[-(1:7), ]
f <- fit_inla(
  pods ~ trt * Region,
  ~ Site + pollin_Site_trt + pollin_Site_Plot,
  data = d2
)
lv <- lapply(fb_model_terms(f)$random_terms, `[[`, "var_levels")
names(lv) <- vapply(fb_model_terms(f)$random_terms, `[[`, "", "var")
X <- fb_fixed_design(~ trt * Region, d2)$X
cat("fixed design columns:", colnames(X), "\n")
set.seed(3)
s <- fb_inla_posterior_sample(f, 50L, seed = 3L, num.threads = "1:1")
gap <- max(vapply(
  s,
  function(z) {
    lat <- z$latent[, 1L]
    eta <- X %*% lat[paste0(colnames(X), ":1")]
    for (v in names(lv)) {
      eta <- eta + lat[paste0(v, ":", match(as.character(d2[[v]]), lv[[v]]))]
    }
    max(abs(eta - lat[grep("^Predictor:", names(lat))]))
  },
  numeric(1)
))
cat("max |X beta + Z u - Predictor| over 50 draws:", format(gap), "\n")

# 5. Families --------------------------------------------------------------------
set.seed(2)
d$y01 <- stats::rbinom(nrow(d), 1L, 0.4)
d$cnt <- stats::rpois(nrow(d), 20)
d$pos <- stats::rgamma(nrow(d), 5, 1)
d$prop <- stats::rbeta(nrow(d), 4, 6)
for (z in list(
  c("y01", "binomial"), c("cnt", "poisson"), c("pods", "negative_binomial"),
  c("pos", "gamma"), c("prop", "beta"), c("pos", "tweedie"),
  c("pods", "nbinom2")
)) {
  f <- fit_inla(stats::as.formula(paste(z[1L], "~ trt")), ~Site, family = z[2L])
  if (!is.null(f)) {
    cat("  link:", fb_model_terms(f)$link, "| converged:", fb_convergence(f)$ok, "\n")
  }
}

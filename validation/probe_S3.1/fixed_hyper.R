# Oracle-hyperparameter check on A1-like data: INLA with the site SD, cell SD
# and NB size fixed at their planted values. With the hyperparameters known,
# the exact posterior of the trial-conditional loss share is calibrated over
# the generating distribution, so any shortfall here comes from the latent
# approximation / sampler, not from priors or hyperparameter integration.
# Usage (from the package root): Rscript validation/probe_S3.1/fixed_hyper.R <n_seeds> <cores> <out_csv> [mode]
#   mode = "gauss" (default inla.posterior.sample, skew.corr = TRUE)
#        | "noskew" (skew.corr = FALSE)
#        | "laplace" (control.inla strategy = "laplace", default sampler)
args <- commandArgs(trailingOnly = TRUE)
n_seeds <- as.integer(args[1L]); cores <- as.integer(args[2L]); out_csv <- args[3L]
mode <- if (length(args) >= 4L) args[4L] else "gauss"
suppressMessages(library(INLA))
inla.setOption(num.threads = "1:1")

one <- function(s) {
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(900000L + s)   # same data seeds as the probe's A1 settings
  n_site <- 12L; arms <- c("Cage", "Open"); beta <- log(c(Cage = 20, Open = 30))
  sites <- sprintf("s%02d", seq_len(n_site))
  g <- expand.grid(rep = 1:8, trt = arms, site = sites, stringsAsFactors = FALSE)
  u_site <- rnorm(n_site, 0, 0.40); names(u_site) <- sites
  u_cell <- matrix(rnorm(n_site * 2, 0, 0.20), n_site, dimnames = list(sites, arms))
  eta <- sapply(arms, function(a) beta[[a]] + u_site[g$site] + u_cell[cbind(g$site, a)])
  own <- eta[cbind(seq_len(nrow(g)), match(g$trt, arms))]
  g$y <- rnbinom(nrow(g), mu = exp(own), size = 8)
  truth <- 1 - mean(exp(eta[, "Cage"])) / mean(exp(eta[, "Open"]))

  g$site_i <- match(g$site, sites)
  g$cell_i <- (g$site_i - 1L) * 2L + match(g$trt, arms)
  g$open <- as.integer(g$trt == "Open")
  fx <- function(sd) list(prec = list(initial = log(1 / sd^2), fixed = TRUE))
  f <- y ~ 1 + open + f(site_i, model = "iid", hyper = fx(0.40)) +
    f(cell_i, model = "iid", hyper = fx(0.20))
  ci <- if (mode == "laplace") list(strategy = "laplace") else list()
  m <- inla(f, family = "nbinomial", data = g,
            control.family = list(hyper = list(size = list(initial = log(8), fixed = TRUE))),
            control.fixed = list(prec.intercept = 1e-4, prec = 1e-4),
            control.inla = ci, control.compute = list(config = TRUE))
  set.seed(990000L + s)
  smp <- inla.posterior.sample(1000L, m, seed = 990000L + s, num.threads = "1:1",
                               skew.corr = mode != "noskew")
  L <- do.call(cbind, lapply(smp, `[[`, "latent"))
  rn <- rownames(smp[[1]]$latent)
  b0 <- L["(Intercept):1", ]; b1 <- L["open:1", ]
  us <- L[paste0("site_i:", 1:n_site), , drop = FALSE]
  uc <- L[paste0("cell_i:", 1:(2 * n_site)), , drop = FALSE]
  si <- g$site_i
  e_c <- matrix(b0, nrow(g), 1000, byrow = TRUE) + us[si, ] + uc[(si - 1L) * 2L + 1L, ]
  e_o <- matrix(b0 + b1, nrow(g), 1000, byrow = TRUE) + us[si, ] + uc[(si - 1L) * 2L + 2L, ]
  r <- 1 - colMeans(exp(e_c)) / colMeans(exp(e_o))
  q <- quantile(r, c(0.025, 0.975), names = FALSE)
  # Second moments of the cell effects: sampler vs INLA marginal summaries
  sd_ratio <- median(apply(uc, 1, sd) / m$summary.random$cell_i$sd)
  data.frame(seed = s, truth = truth, est = median(r), sd = sd(r), low = q[1], high = q[2],
             sd_ratio_cell = sd_ratio)
}
res <- do.call(rbind, parallel::mclapply(seq_len(n_seeds), one, mc.cores = cores,
                                         mc.preschedule = FALSE))
write.csv(res, out_csv, row.names = FALSE)
cov <- mean(res$truth >= res$low & res$truth <= res$high)
z <- (res$est - res$truth) / res$sd
message(sprintf("fixed-hyper %s: n=%d cov=%.3f (MCSE %.3f) sd(z)=%.3f mean(z)=%.3f medwidth=%.4f sd_ratio_cell=%.3f",
                mode, nrow(res), cov, sqrt(cov * (1 - cov) / nrow(res)), sd(z), mean(z),
                median(res$high - res$low), median(res$sd_ratio_cell)))

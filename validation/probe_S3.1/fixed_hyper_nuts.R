# Exact reference for fixed_hyper.R: the same model (site SD 0.40, cell SD
# 0.20, NB size 8 fixed at the planted values; flat-ish fixed effects) sampled
# by NUTS (cmdstanr), on the same data seeds 900000 + s.
# Usage (from the package root): Rscript validation/probe_S3.1/fixed_hyper_nuts.R <n_seeds> <cores> <out_csv>
args <- commandArgs(trailingOnly = TRUE)
n_seeds <- as.integer(args[1L]); cores <- as.integer(args[2L]); out_csv <- args[3L]
code <- "
data { int N; int S; array[N] int y; vector[N] open; array[N] int site; array[N] int cell; }
parameters { real b0; real b1; vector[S] zs; vector[2 * S] zc; }
model {
  b0 ~ normal(0, 100); b1 ~ normal(0, 100); zs ~ std_normal(); zc ~ std_normal();
  y ~ neg_binomial_2_log(b0 + b1 * open + 0.40 * zs[site] + 0.20 * zc[cell], 8);
}"
f <- cmdstanr::write_stan_file(code, dir = tempdir())
mod <- cmdstanr::cmdstan_model(f, dir = tempdir(), quiet = TRUE)

one <- function(s) {
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(900000L + s)
  n_site <- 12L; arms <- c("Cage", "Open"); beta <- log(c(Cage = 20, Open = 30))
  sites <- sprintf("s%02d", seq_len(n_site))
  g <- expand.grid(rep = 1:8, trt = arms, site = sites, stringsAsFactors = FALSE)
  u_site <- rnorm(n_site, 0, 0.40); names(u_site) <- sites
  u_cell <- matrix(rnorm(n_site * 2, 0, 0.20), n_site, dimnames = list(sites, arms))
  eta <- sapply(arms, function(a) beta[[a]] + u_site[g$site] + u_cell[cbind(g$site, a)])
  own <- eta[cbind(seq_len(nrow(g)), match(g$trt, arms))]
  g$y <- rnbinom(nrow(g), mu = exp(own), size = 8)
  truth <- 1 - mean(exp(eta[, "Cage"])) / mean(exp(eta[, "Open"]))
  si <- match(g$site, sites); ci <- (si - 1L) * 2L + match(g$trt, arms)
  fit <- mod$sample(data = list(N = nrow(g), S = n_site, y = g$y, open = as.numeric(g$trt == "Open"),
                                site = si, cell = ci),
                    chains = 2, parallel_chains = 1, iter_warmup = 1000, iter_sampling = 1000,
                    seed = 990000L + s, refresh = 0, show_messages = FALSE, show_exceptions = FALSE)
  d <- fit$draws(format = "draws_matrix")
  b0 <- d[, "b0"]; b1 <- d[, "b1"]
  zs <- d[, paste0("zs[", 1:n_site, "]")]; zc <- d[, paste0("zc[", 1:(2 * n_site), "]")]
  us <- t(0.40 * zs); uc <- t(0.20 * zc)
  e_c <- matrix(b0, nrow(g), nrow(d), byrow = TRUE) + us[si, ] + uc[(si - 1L) * 2L + 1L, ]
  e_o <- matrix(b0 + b1, nrow(g), nrow(d), byrow = TRUE) + us[si, ] + uc[(si - 1L) * 2L + 2L, ]
  r <- 1 - colMeans(exp(e_c)) / colMeans(exp(e_o))
  q <- quantile(r, c(0.025, 0.975), names = FALSE)
  sm <- fit$summary(c("b0", "b1"))
  data.frame(seed = s, truth = truth, est = median(r), sd = sd(r), low = q[1], high = q[2],
             max_rhat = max(sm$rhat), n_div = sum(fit$diagnostic_summary(quiet = TRUE)$num_divergent))
}
res <- do.call(rbind, parallel::mclapply(seq_len(n_seeds), one, mc.cores = cores,
                                         mc.preschedule = FALSE))
write.csv(res, out_csv, row.names = FALSE)
cov <- mean(res$truth >= res$low & res$truth <= res$high)
z <- (res$est - res$truth) / res$sd
message(sprintf("fixed-hyper NUTS: n=%d cov=%.3f (MCSE %.3f) sd(z)=%.3f mean(z)=%.3f medwidth=%.4f max_rhat=%.3f div=%d",
                nrow(res), cov, sqrt(cov * (1 - cov) / nrow(res)), sd(z), mean(z),
                median(res$high - res$low), max(res$max_rhat), sum(res$n_div)))

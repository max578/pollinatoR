# causal_probe_S4.R -- Stage 4 Step 0 probes (a) dagitty and (b) moderation.
#
# Run from the package root:
#   Rscript --vanilla validation/causal_probe_S4.R
# Synthetic data only. Findings are written up in validation/causal_probe_S4.md.

suppressPackageStartupMessages({
  library(dagitty)
  pkgload::load_all(".", quiet = TRUE)
})
cat(
  "dagitty", as.character(utils::packageVersion("dagitty")),
  "| glmmTMB", as.character(utils::packageVersion("glmmTMB")),
  "| pollinatoR", as.character(utils::packageVersion("pollinatoR")), "\n"
)

# (a1) Minimal veg-weight graph (plan Step 0 command) -------------------------
g0 <- dagitty(
  "dag { trt -> yield; trt -> veg; veg -> yield; rain -> yield; rain -> veg }"
)
exposures(g0) <- "trt"
outcomes(g0) <- "yield"
cat("\n(a1) minimal graph\n total: ")
print(adjustmentSets(g0, effect = "total"))
cat(" direct: ")
print(adjustmentSets(g0, effect = "direct"))

# (a2) Field-study graph, branch H (veg measured at harvest) ------------------
# trt is randomised within site; site-level nodes: rain, forest, bees,
# fertility (latent). Interaction (bees modifies the trt contrast) cannot be
# drawn on a DAG; bees -> yield stands for it.
field_h <- "dag {
  fertility [latent]
  trt -> veg -> yield
  trt -> yield
  rain -> veg
  rain -> yield
  fertility -> veg
  fertility -> yield
  forest -> bees -> yield
  forest -> fertility
}"
gh <- dagitty(field_h)
exposures(gh) <- "trt"
outcomes(gh) <- "yield"
cat("\n(a2) field graph, branch H\n total: ")
print(adjustmentSets(gh, effect = "total"))
cat(" direct: ")
print(adjustmentSets(gh, effect = "direct"))
cat(" descendants of trt: ", paste(descendants(gh, "trt"), collapse = ", "), "\n")

# Exposure 'bees' (the who-delivers driver), total effect on yield
gb <- gh
exposures(gb) <- "bees"
cat(" bees -> yield, total: ")
print(adjustmentSets(gb, effect = "total"))

# (a3) Branch P (veg measured before flowering: no trt -> veg edge) ----------
gp <- dagitty(sub("trt -> veg -> yield", "veg -> yield", field_h, fixed = TRUE))
exposures(gp) <- "trt"
outcomes(gp) <- "yield"
cat("\n(a3) field graph, branch P\n total: ")
print(adjustmentSets(gp, effect = "total"))
cat(" descendants of trt: ", paste(descendants(gp, "trt"), collapse = ", "), "\n")

# (a4) localTests on a simulated site table (44 sites, site-level nodes) ------
# Only the site-level part of the graph is testable at site level.
set.seed(404)
n_site <- 44L
forest <- stats::rnorm(n_site)
fertility <- 0.5 * forest + stats::rnorm(n_site)
rain <- stats::rnorm(n_site)
bees <- 0.6 * forest + stats::rnorm(n_site)
veg <- 0.5 * rain + 0.5 * fertility + stats::rnorm(n_site)
site_tab <- data.frame(forest, rain, bees, veg)
g_site <- dagitty("dag {
  fertility [latent]
  rain -> veg
  fertility -> veg
  forest -> bees
  forest -> fertility
}")
cat("\n(a4) implied independencies (site level):\n")
print(impliedConditionalIndependencies(g_site))
lt <- localTests(g_site, data = site_tab, type = "cis")
print(round(lt, 3))
cat(" n per test:", n_site, "\n")

# (b) pollin_moderation() on a site-level-only exposure -----------------------
set.seed(1)
site <- rep(paste0("s", seq_len(14L)), each = 20L)
trt <- rep(rep(c("Cage", "Open"), each = 10L), 14L)
bees_site <- rep(stats::rpois(14L, 8), each = 20L)
y <- stats::rpois(
  length(trt),
  exp(1.6 + 0.3 * (trt == "Open") + 0.02 * bees_site * (trt == "Open"))
)
d <- data.frame(y = y, trt = trt, Site = site, bees = bees_site)
pd <- pollin_data(d, "y", "trt", unit = "Site", exposures = "bees")
ds <- pollin_design(pd, c(excluded = "Cage", open = "Open"), "multi_site")

cat("\n(b1) current behaviour, site-constant exposure:\n")
w <- character(0)
res <- withCallingHandlers(
  tryCatch(pollin_moderation(ds, outcome = "y"), error = function(e) e),
  warning = function(cnd) {
    w <<- c(w, conditionMessage(cnd))
    invokeRestart("muffleWarning")
  }
)
if (inherits(res, "error")) {
  cat(" ERROR:", conditionMessage(res), "\n")
} else {
  print(res$table)
}
cat(" warnings (", length(w), "):\n", paste0("  - ", unique(w), "\n"), sep = "")

# (b2) Fix sketch: drop a zero-variance within component, label between-only --
inputs <- .design_inputs(ds, ds$arms[["excluded"]])
d2 <- inputs$data
wb <- pollin_within_between(d2$bees, inputs$wb_group)
cat(
  "\n(b2) within variance:", signif(stats::var(wb$within), 3),
  "| between variance:", signif(stats::var(wb$between), 3), "\n"
)
d2$bees_between <- .z(wb$between)
fit2 <- .fit_glmm(
  "y",
  c(sprintf("%s * bees_between", inputs$treatment), inputs$random),
  d2,
  .pick_family(d2$y)
)
co <- summary(fit2)$coefficients$cond
print(signif(co[grepl(":", rownames(co)), , drop = FALSE], 4))
cat(" pdHess:", fit2$sdr$pdHess, "\n")
# Planted interaction on the link scale: 0.02 per bee; per SD of the site
# means that is 0.02 * sd(site means).
cat(
  " planted per-SD interaction:",
  signif(0.02 * stats::sd(wb$between), 4), "\n"
)

# make_golden_v0.1.0.R -- Capture the outputs of the layered API at 0.1.0.
#
# Run once, from the package root, on the 0.1.0 code (commit fdf7007):
#   Rscript --vanilla validation/golden/make_golden_v0.1.0.R
# Writes tests/testthat/_golden/golden_v0.1.0.rds. The file is the oracle for
# Stage 2 (criterion C1) and is never regenerated from later code.

pkgload::load_all(".", quiet = TRUE)
source("tests/testthat/helper-golden.R")
d <- golden_inputs_v010()
g <- golden_inputs()

pd <- pollin_data(
  d,
  c("pods", "seed_wt"),
  "trt",
  unit = "Block",
  groups = "Region",
  exposures = "bees"
)
ds <- pollin_design(pd, golden_arms_v010, "rcbd")
grid <- golden_grid_v010()
n_boot <- 19L

# E1 and E2 on every scale, delta and bootstrap intervals ---------------------

effects <- list()
for (fam in c("deficit", "dependence")) {
  for (i in seq_len(nrow(grid))) {
    key <- paste(fam, grid$scale[i], grid$basis[i], sep = "_")
    est <- pollin_estimand(ds, fam, grid$scale[i], grid$basis[i])
    fit <- suppressWarnings(pollin_fit(est))
    by_fit <- suppressWarnings(pollin_fit(est, by = "Region"))
    effects[[key]] <- list(
      estimand = unclass(est)[c("family", "scale", "basis", "population", "contrast")],
      assumptions = pollin_assumptions(est),
      delta = golden_effects(suppressWarnings(pollin_effects(fit))),
      delta_by = golden_effects(suppressWarnings(pollin_effects(by_fit))),
      boot = golden_effects(suppressWarnings(pollin_effects(
        fit,
        interval = "bootstrap",
        n_boot = n_boot,
        seed = 31L
      ))),
      boot_by = golden_effects(suppressWarnings(pollin_effects(
        by_fit,
        interval = "bootstrap",
        n_boot = n_boot,
        seed = 32L
      )))
    )
  }
}

# Layered replacements for the 0.0.0.9000 specification paths -----------------

pd_one <- pollin_data(g$plants, "pods", "trt", unit = "Site")
ds_one <- pollin_design(pd_one, c(excluded = "Cage", open = "Open"), "multi_site")
pd_nest <- pollin_data(
  g$plants,
  c("pods", "seed_wt"),
  "trt",
  unit = c("Site", "Plot"),
  groups = "Region",
  exposures = c("bees", "wild"),
  adjust = "rain"
)
ds_nest <- pollin_design(pd_nest, c(excluded = "Cage", open = "Open"), "multi_site")
e_one <- pollin_estimand(ds_one, "dependence", "difference", "loss")
e_nest <- pollin_estimand(ds_nest, "dependence", "difference", "loss")
strip_contribution <- function(x) unclass(x)[["table"]]

out <- list(
  effects = effects,
  n_boot = n_boot,
  random_terms = list(one = ds_one$random, nest = ds_nest$random),
  contribution_one = strip_contribution(suppressWarnings(pollin_contribution(e_one))),
  contribution_nest = strip_contribution(suppressWarnings(pollin_contribution(e_nest))),
  contribution_by = strip_contribution(suppressWarnings(pollin_contribution(
    pollin_estimand(ds, "dependence", "difference", "loss"),
    by = "Region"
  ))),
  moderation = unclass(suppressWarnings(pollin_moderation(
    ds_nest,
    outcome = "pods"
  )))[["table"]],
  validate = unclass(pollin_validate(
    pd_nest,
    expect_levels = list(trt = c("Cage", "Open")),
    identities = list(pods_positive = function(d) d$pods >= 0)
  )),
  index_design = pollin_index(ds, response = "pods", n_boot = 199L, seed = 13L),
  compensation_data = unclass(pollin_compensation(pd)),
  within_between = pollin_within_between(d$bees, d$Block)
)

attr(out, "provenance") <- list(
  package_version = as.character(utils::packageVersion("pollinatoR")),
  git_head = system("git rev-parse HEAD", intern = TRUE),
  glmmTMB = as.character(utils::packageVersion("glmmTMB")),
  TMB = as.character(utils::packageVersion("TMB")),
  marginaleffects = as.character(utils::packageVersion("marginaleffects")),
  R = R.version.string,
  created = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
)

dest <- file.path("tests", "testthat", "_golden", "golden_v0.1.0.rds")
if (file.exists(dest)) {
  stop("`", dest, "` exists; the golden file is never overwritten.", call. = FALSE)
}
saveRDS(out, dest)
message("Wrote ", dest)

# make_golden.R -- Capture the outputs of the ten exports at 0.0.0.9000.
#
# Run once, from the package root, on the 0.0.0.9000 code:
#   Rscript --vanilla validation/golden/make_golden.R
# Writes tests/testthat/_golden/golden_v0.0.0.9000.rds. The file is the oracle
# for the Stage 1 rewrite and is never regenerated from rewritten code.

pkgload::load_all(".", quiet = TRUE)
source("tests/testthat/helper-golden.R")
g <- golden_inputs()

# Specifications ---------------------------------------------------------------

sp_one <- pollin_spec(g$plants, "pods", "trt", "Cage", "Site")
sp_nest <- pollin_spec(
  g$plants,
  outcomes = c("pods", "seed_wt"),
  treatment = "trt",
  control = "Cage",
  random = c("Site", "Plot"),
  exposures = c("bees", "wild"),
  adjust = "rain",
  group = "Region"
)

# Model-based outputs ----------------------------------------------------------

contrib_one <- suppressWarnings(pollin_contribution(sp_one))
contrib_nest <- suppressWarnings(pollin_contribution(sp_nest))
contrib_by <- suppressWarnings(pollin_contribution(sp_nest, by = "Region"))
moder <- suppressWarnings(pollin_moderation(sp_nest, outcome = "pods"))

# Index, guard and utility outputs ---------------------------------------------

out <- list(
  random_terms = list(
    one = .random_terms(sp_one),
    nest = .random_terms(sp_nest)
  ),
  spec_nest = golden_strip(sp_nest),
  contribution_one = golden_strip(contrib_one),
  contribution_nest = golden_strip(contrib_nest),
  contribution_by = golden_strip(contrib_by),
  moderation = golden_strip(moder),
  rii = unclass(pollin_rii(g$open, g$cage, n_boot = 999L, seed = 11L)),
  index = pollin_index(g$open, g$cage, n_boot = 999L, seed = 11L),
  spurious = unclass(pollin_guard_spurious(
    g$open,
    g$cage,
    n_perm = 499L,
    decouple = g$vigour,
    seed = 12L
  )),
  compensation = unclass(pollin_compensation(g$components)),
  validate = unclass(pollin_validate(
    sp_nest,
    expect_levels = list(trt = c("Cage", "Open")),
    identities = list(pods_positive = function(d) d$pods >= 0)
  )),
  within_between = pollin_within_between(g$plants$bees, g$plants$Site),
  fdr = pollin_fdr(g$p)
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

dest <- file.path("tests", "testthat", "_golden", "golden_v0.0.0.9000.rds")
if (file.exists(dest)) {
  stop("`", dest, "` exists; the golden file is never overwritten.", call. = FALSE)
}
saveRDS(out, dest)
message("Wrote ", dest)

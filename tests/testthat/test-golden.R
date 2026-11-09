# Golden values: numbers carried from 0.0.0.9000 ------------------------------

gold <- readRDS(golden_path())
g <- golden_inputs()

# The 0.0.0.9000 study declarations, in the layered API.
golden_designs <- function(g) {
  arms <- c(excluded = "Cage", open = "Open")
  pd_one <- pollin_data(g$plants, "pods", "trt", unit = "Site")
  pd_nest <- pollin_data(
    g$plants,
    c("pods", "seed_wt"),
    "trt",
    unit = c("Site", "Plot"),
    groups = "Region",
    exposures = c("bees", "wild"),
    adjust = "rain"
  )
  list(
    one = pollin_design(pd_one, arms, "multi_site"),
    nest = pollin_design(pd_nest, arms, "multi_site")
  )
}

test_that("index, guard and utility outputs match 0.0.0.9000 exactly", {
  expect_identical(
    unclass(pollin_rii(g$open, g$cage, n_boot = 999L, seed = 11L)),
    gold$rii
  )
  expect_identical(
    pollin_index(g$open, g$cage, n_boot = 999L, seed = 11L),
    gold$index
  )
  expect_identical(
    unclass(pollin_guard_spurious(
      g$open,
      g$cage,
      n_perm = 499L,
      decouple = g$vigour,
      seed = 12L
    )),
    gold$spurious
  )
  expect_identical(
    unclass(pollin_compensation(g$components)),
    gold$compensation
  )
  expect_identical(
    pollin_within_between(g$plants$bees, g$plants$Site),
    gold$within_between
  )
  expect_identical(pollin_fdr(g$p), gold$fdr)
})

test_that("random terms and validation match 0.0.0.9000", {
  ds <- golden_designs(g)
  expect_identical(ds$one$random, gold$random_terms$one)
  expect_identical(ds$nest$random, gold$random_terms$nest)
  expect_identical(
    unclass(pollin_validate(
      ds$nest$data,
      expect_levels = list(trt = c("Cage", "Open")),
      identities = list(pods_positive = function(d) d$pods >= 0)
    )),
    gold$validate
  )
})

# The by-level contribution of 0.0.0.9000 averaged one treatment effect over
# the levels; since 0.1.0 each level has its own fixed effect, so that table
# has no counterpart and is guarded by the 0.1.0 golden values instead.
test_that("model-based contributions match 0.0.0.9000", {
  skip_on_cran()
  ds <- golden_designs(g)
  tol <- 1e-10
  num <- c("estimate", "std.error", "p.value", "conf.low", "conf.high")
  for (k in c("one", "nest")) {
    est <- pollin_estimand(ds[[k]], "dependence", "difference", "loss")
    tab <- suppressWarnings(
      pollin_contribution(est, interval = "delta_fixed")
    )$table
    old <- gold[[paste0("contribution_", k)]]$table
    expect_identical(tab$response, old$outcome, info = k)
    expect_identical(tab$distribution, old$family, info = k)
    expect_equal(tab[num], old[num], tolerance = tol, info = k)
    expect_equal(tab$slope_sd, old$slope_sd, tolerance = tol, info = k)
    expect_equal(tab$p.fdr, old$p.fdr, tolerance = tol, info = k)
  }
  mod <- suppressWarnings(pollin_moderation(ds$nest, outcome = "pods"))
  expect_equal(mod$table, gold$moderation$table, tolerance = tol)
  expect_identical(mod$outcome, gold$moderation$outcome)
})

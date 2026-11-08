# Golden values: numbers carried from 0.0.0.9000 ------------------------------

gold <- readRDS(golden_path())
g <- golden_inputs()

golden_specs <- function(g) {
  list(
    one = suppressWarnings(pollin_spec(
      g$plants,
      "pods",
      "trt",
      "Cage",
      "Site"
    )),
    nest = suppressWarnings(pollin_spec(
      g$plants,
      outcomes = c("pods", "seed_wt"),
      treatment = "trt",
      control = "Cage",
      random = c("Site", "Plot"),
      exposures = c("bees", "wild"),
      adjust = "rain",
      group = "Region"
    ))
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

test_that("specification and validation match 0.0.0.9000", {
  sp <- golden_specs(g)
  expect_identical(.random_terms(sp$one), gold$random_terms$one)
  expect_identical(.random_terms(sp$nest), gold$random_terms$nest)
  expect_identical(golden_strip(sp$nest), gold$spec_nest)
  expect_identical(
    unclass(pollin_validate(
      sp$nest,
      expect_levels = list(trt = c("Cage", "Open")),
      identities = list(pods_positive = function(d) d$pods >= 0)
    )),
    gold$validate
  )
})

test_that("model-based contributions match 0.0.0.9000", {
  skip_on_cran()
  sp <- golden_specs(g)
  tol <- 1e-10
  expect_equal(
    golden_strip(suppressWarnings(pollin_contribution(sp$one))),
    gold$contribution_one,
    tolerance = tol
  )
  expect_equal(
    golden_strip(suppressWarnings(pollin_contribution(sp$nest))),
    gold$contribution_nest,
    tolerance = tol
  )
  expect_equal(
    golden_strip(suppressWarnings(pollin_contribution(sp$nest, by = "Region"))),
    gold$contribution_by,
    tolerance = tol
  )
  expect_equal(
    golden_strip(suppressWarnings(pollin_moderation(
      sp$nest,
      outcome = "pods"
    ))),
    gold$moderation,
    tolerance = tol
  )
})

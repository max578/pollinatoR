# Golden values: numbers carried from 0.1.0 -----------------------------------
#
# The oracle is `_golden/golden_v0.1.0.rds`, written once from the 0.1.0 code
# by `validation/golden/make_golden_v0.1.0.R`. Tolerances are those of the
# Stage 1 golden tests: exact identity for the outputs that involve no model
# fit, and a relative tolerance of 1e-10 for the model-based outputs.

gold_v010 <- readRDS(golden_path_v010())
d_v010 <- golden_inputs_v010()
tol_v010 <- 1e-10

golden_design_v010 <- function(d) {
  pd <- pollin_data(
    d,
    c("pods", "seed_wt"),
    "trt",
    unit = "Block",
    groups = "Region",
    exposures = "bees"
  )
  pollin_design(pd, golden_arms_v010, "rcbd")
}

test_that("estimands and assumptions match 0.1.0 exactly", {
  ds <- golden_design_v010(d_v010)
  grid <- golden_grid_v010()
  for (fam in c("deficit", "dependence")) {
    for (i in seq_len(nrow(grid))) {
      key <- paste(fam, grid$scale[i], grid$basis[i], sep = "_")
      est <- pollin_estimand(ds, fam, grid$scale[i], grid$basis[i])
      gold <- gold_v010$effects[[key]]
      expect_identical(
        unclass(est)[c("family", "scale", "basis", "population", "contrast")],
        gold$estimand,
        info = key
      )
      expect_identical(pollin_assumptions(est), gold$assumptions, info = key)
    }
  }
})

test_that("delta_fixed effects on every scale match the 0.1.0 delta", {
  skip_on_cran()
  ds <- golden_design_v010(d_v010)
  grid <- golden_grid_v010()
  for (fam in c("deficit", "dependence")) {
    for (i in seq_len(nrow(grid))) {
      key <- paste(fam, grid$scale[i], grid$basis[i], sep = "_")
      est <- pollin_estimand(ds, fam, grid$scale[i], grid$basis[i])
      gold <- gold_v010$effects[[key]]
      fit <- suppressWarnings(pollin_fit(est))
      by_fit <- suppressWarnings(pollin_fit(est, by = "Region"))
      expect_equal(
        golden_fixed_as_delta(golden_effects(suppressWarnings(
          pollin_effects(fit, interval = "delta_fixed")
        ))),
        gold$delta,
        tolerance = tol_v010,
        info = key
      )
      expect_equal(
        golden_fixed_as_delta(golden_effects(suppressWarnings(
          pollin_effects(by_fit, interval = "delta_fixed")
        ))),
        gold$delta_by,
        tolerance = tol_v010,
        info = key
      )
    }
  }
})

test_that("bootstrap effects on every scale match 0.1.0", {
  skip_on_cran()
  ds <- golden_design_v010(d_v010)
  grid <- golden_grid_v010()
  n_boot <- gold_v010$n_boot
  for (fam in c("deficit", "dependence")) {
    for (i in seq_len(nrow(grid))) {
      key <- paste(fam, grid$scale[i], grid$basis[i], sep = "_")
      est <- pollin_estimand(ds, fam, grid$scale[i], grid$basis[i])
      gold <- gold_v010$effects[[key]]
      fit <- suppressWarnings(pollin_fit(est))
      by_fit <- suppressWarnings(pollin_fit(est, by = "Region"))
      boot <- suppressWarnings(pollin_effects(
        fit,
        interval = "bootstrap",
        n_boot = n_boot,
        seed = 31L
      ))
      boot_by <- suppressWarnings(pollin_effects(
        by_fit,
        interval = "bootstrap",
        n_boot = n_boot,
        seed = 32L
      ))
      expect_equal(golden_effects(boot), gold$boot, tolerance = tol_v010, info = key)
      expect_equal(
        golden_effects(boot_by),
        gold$boot_by,
        tolerance = tol_v010,
        info = key
      )
    }
  }
})

test_that("layered contribution, moderation and validation match 0.1.0", {
  skip_on_cran()
  g <- golden_inputs()
  pd_one <- pollin_data(g$plants, "pods", "trt", unit = "Site")
  ds_one <- pollin_design(
    pd_one,
    c(excluded = "Cage", open = "Open"),
    "multi_site"
  )
  pd_nest <- pollin_data(
    g$plants,
    c("pods", "seed_wt"),
    "trt",
    unit = c("Site", "Plot"),
    groups = "Region",
    exposures = c("bees", "wild"),
    adjust = "rain"
  )
  ds_nest <- pollin_design(
    pd_nest,
    c(excluded = "Cage", open = "Open"),
    "multi_site"
  )
  expect_identical(
    list(one = ds_one$random, nest = ds_nest$random),
    gold_v010$random_terms
  )
  e_one <- pollin_estimand(ds_one, "dependence", "difference", "loss")
  e_nest <- pollin_estimand(ds_nest, "dependence", "difference", "loss")
  e_by <- pollin_estimand(
    golden_design_v010(d_v010),
    "dependence",
    "difference",
    "loss"
  )
  expect_equal(
    golden_fixed_as_delta(suppressWarnings(
      pollin_contribution(e_one, interval = "delta_fixed")
    )$table),
    gold_v010$contribution_one,
    tolerance = tol_v010
  )
  expect_equal(
    golden_fixed_as_delta(suppressWarnings(
      pollin_contribution(e_nest, interval = "delta_fixed")
    )$table),
    gold_v010$contribution_nest,
    tolerance = tol_v010
  )
  expect_equal(
    golden_fixed_as_delta(suppressWarnings(
      pollin_contribution(e_by, by = "Region", interval = "delta_fixed")
    )$table),
    gold_v010$contribution_by,
    tolerance = tol_v010
  )
  expect_equal(
    suppressWarnings(pollin_moderation(ds_nest, outcome = "pods"))$table,
    gold_v010$moderation,
    tolerance = tol_v010
  )
  expect_identical(
    unclass(pollin_validate(
      pd_nest,
      expect_levels = list(trt = c("Cage", "Open")),
      identities = list(pods_positive = function(d) d$pods >= 0)
    )),
    gold_v010$validate
  )
})

test_that("indices and utilities from the layered objects match 0.1.0", {
  ds <- golden_design_v010(d_v010)
  expect_identical(
    pollin_index(ds, response = "pods", n_boot = 199L, seed = 13L),
    gold_v010$index_design
  )
  expect_identical(
    unclass(pollin_compensation(ds$data)),
    gold_v010$compensation_data
  )
  expect_identical(
    pollin_within_between(d_v010$bees, d_v010$Block),
    gold_v010$within_between
  )
})

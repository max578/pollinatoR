# Layer L5: G4 guards ---------------------------------------------------------

# Four sites with unequal yields; arm means chosen by hand.
guard_design <- function() {
  df <- data.frame(
    pods = c(10, 20, 10, 20, 30, 60, 30, 60, 50, 50, 50, 50, 5, 40, 5, 40),
    trt = rep(c("Cage", "Open"), 8),
    Site = rep(c("a", "b", "c", "d"), each = 4)
  )
  pd <- pollin_data(df, "pods", "trt", unit = "Site")
  pollin_design(pd, c(excluded = "Cage", open = "Open"), "multi_site")
}

test_that("bounded guard refuses to average unit-level ratios by default", {
  ds <- guard_design()
  expect_error(pollin_guard_bounded(ds, index = "rii"), "not the ratio")
  expect_error(pollin_guard_bounded(ds, index = "rii"), "allow = TRUE")
  expect_error(pollin_guard_bounded(c(0.1, 0.4), "dependence"), "not the ratio")
  expect_error(pollin_guard_bounded(ds), "`index` must be one of")
})

test_that("bounded guard averages when allowed and says so", {
  ds <- guard_design()
  b <- pollin_guard_bounded(ds, index = "rii", allow = TRUE)
  per_unit <- c(10 / 30, 30 / 90, 0, 35 / 45)
  expect_equal(b$per_unit$value, per_unit)
  expect_equal(b$mean_of_units, mean(per_unit))
  expect_equal(b$ratio_of_means, (42.5 - 23.75) / (42.5 + 23.75))
  expect_output(print(b), "not the ratio of expected values")
  d <- pollin_guard_bounded(ds, index = "dependence", allow = TRUE)
  expect_equal(d$mean_of_units, mean(c(0.5, 0.5, 0, 35 / 40)))
  expect_equal(d$ratio_of_means, 1 - 23.75 / 42.5)
  expect_warning(
    v <- pollin_guard_bounded(c(0.2, 0.4, NA), "rii", allow = TRUE),
    "1 unit with an undefined index"
  )
  expect_equal(v$mean_of_units, 0.3)
  expect_true(is.na(v$ratio_of_means))
  expect_output(print(v), "not available")
})

test_that("bounded guard checks its input", {
  expect_error(pollin_guard_bounded(c(0.2, 1.4), "rii"), "outside the bounds")
  expect_error(pollin_guard_bounded("a", "rii"), "numeric vector")
  expect_error(pollin_guard_bounded(0.2, "rii", allow = NA), "TRUE or FALSE")
  expect_warning(
    pollin_guard_bounded(c(0.2, Inf), "rii", allow = TRUE),
    "1 unit with an undefined index"
  )
})

# Ten sites; the exposure varies within sites.
eco_df <- function(within_slope) {
  set.seed(3L)
  site <- rep(sprintf("s%02d", 1:10), each = 8L)
  site_bees <- rep(seq(4, 22, by = 2), each = 8L)
  bees <- site_bees + stats::rnorm(80L)
  pods <- 20 +
    2 * site_bees +
    within_slope * (bees - site_bees) +
    stats::rnorm(80L, 0, 0.5)
  data.frame(pods = pods, bees = bees, Site = site, Plot = rep(1:2, 40L))
}

test_that("ecological guard warns when only the between slope carries it", {
  expect_warning(
    g <- pollin_guard_ecological(eco_df(0), "bees", "pods", unit = "Site"),
    "Only the between-unit slope"
  )
  expect_equal(g$slopes$component, c("pooled", "within", "between"))
  expect_length(g$flags, 1L)
  expect_warning(
    pollin_guard_ecological(eco_df(-2), "bees", "pods", unit = "Site"),
    "disagree in sign"
  )
  flat <- eco_df(0)
  flat$bees <- rep(seq(4, 22, by = 2), each = 8L)
  expect_warning(
    g <- pollin_guard_ecological(flat, "bees", "pods", unit = "Site"),
    "does not vary within units"
  )
  expect_true(is.na(g$slopes$slope[2L]))
  expect_output(print(g), "! `bees`")
})

test_that("ecological guard passes when within and between agree", {
  expect_no_warning(
    g <- pollin_guard_ecological(eco_df(2), "bees", "pods", unit = "Site")
  )
  expect_length(g$flags, 0L)
  expect_equal(g$slopes$slope[2L], 2, tolerance = 0.1)
  expect_equal(g$slopes$slope[3L], 2, tolerance = 0.1)
  expect_output(print(g), "slopes agree")
})

test_that("ecological guard takes a pollin_data object and checks input", {
  d <- eco_df(2)
  d$trt <- rep(c("Cage", "Open"), 40L)
  pd <- pollin_data(
    d,
    "pods",
    "trt",
    unit = c("Site", "Plot"),
    exposures = "bees"
  )
  g <- pollin_guard_ecological(pd, "bees")
  expect_identical(g$unit, "Site")
  g2 <- suppressWarnings(pollin_guard_ecological(pd, "bees", unit = "Plot"))
  expect_identical(g2$unit, c("Site", "Plot"))
  expect_error(pollin_guard_ecological(d, "bees", "pods"), "unit")
  expect_error(
    pollin_guard_ecological(d, "bees", "nope", unit = "Site"),
    "not found"
  )
  expect_error(
    pollin_guard_ecological(d, "trt", "pods", unit = "Site"),
    "must be numeric"
  )
  expect_error(
    pollin_guard_ecological(
      d[d$Site %in% c("s01", "s02"), ],
      "bees",
      "pods",
      unit = "Site"
    ),
    "at least three units"
  )
  expect_error(
    pollin_guard_ecological(d, "bees", "pods", unit = "Site", alpha = 2),
    "alpha"
  )
})

test_that("posttreatment guard refuses a declared downstream covariate", {
  df <- data.frame(
    pods = 1:8,
    trt = rep(c("Cage", "Open"), 4L),
    Site = rep(c("a", "b"), each = 4L),
    rain = rep(c(250, 300), each = 4L),
    veg_weight = 11:18
  )
  pd <- pollin_data(df, "pods", "trt", "Site", adjust = c("rain", "veg_weight"))
  expect_error(
    pollin_guard_posttreatment(pd, downstream = "veg_weight"),
    "`veg_weight` is declared downstream of `trt`"
  )
  ds <- pollin_design(pd, c(excluded = "Cage", open = "Open"), "multi_site")
  expect_error(
    pollin_guard_posttreatment(ds, c("veg_weight", "rain")),
    "covariates `rain`, `veg_weight` are"
  )
  es <- pollin_estimand(ds, "dependence", "ratio", "loss")
  expect_error(pollin_guard_posttreatment(es, "rain"), "downstream")
})

test_that("posttreatment guard passes when nothing adjusted is downstream", {
  df <- data.frame(
    pods = 1:8,
    trt = rep(c("Cage", "Open"), 4L),
    Site = rep(c("a", "b"), each = 4L),
    rain = rep(c(250, 300), each = 4L),
    veg_weight = 11:18
  )
  pd <- pollin_data(df, "pods", "trt", "Site", adjust = "rain")
  expect_identical(pollin_guard_posttreatment(pd, "veg_weight"), pd)
  expect_invisible(pollin_guard_posttreatment(pd, character(0)))
  expect_error(pollin_guard_posttreatment(pd), "character vector")
  expect_error(pollin_guard_posttreatment(df, "rain"), "`pollin_data`")
})

test_that("posttreatment guard refuses a graph descendant of treatment", {
  skip_if_not_installed("dagitty")
  df <- data.frame(
    pods = 1:8,
    trt = rep(c("Cage", "Open"), 4L),
    Site = rep(c("a", "b"), each = 4L),
    rain = rep(c(250, 300), each = 4L),
    veg_weight = 11:18
  )
  pd <- pollin_data(df, "pods", "trt", "Site", adjust = c("rain", "veg_weight"))
  ds <- pollin_design(pd, c(excluded = "Cage", open = "Open"), "multi_site")
  # Branch H: veg_weight is a mediator, so the graph alone refuses it
  g_h <- pollin_dag(
    ds,
    "dag { rain -> pods; trt -> veg_weight -> pods; trt -> pods }",
    "pods"
  )
  expect_error(
    pollin_guard_posttreatment(ds, dag = g_h),
    "`veg_weight` is in the graph downstream of `trt`"
  )
  # Both sources name it
  expect_error(
    pollin_guard_posttreatment(ds, "veg_weight", dag = g_h),
    "is declared and in the graph downstream"
  )
  # A direct effect holding veg_weight fixed passes the same graph
  es_d <- pollin_estimand(
    ds,
    "dependence",
    "ratio",
    "loss",
    effect = "direct",
    mediators = "veg_weight",
    dag = g_h
  )
  expect_identical(pollin_guard_posttreatment(es_d, dag = g_h), es_d)
  # A graph from another design is refused
  pd2 <- pollin_data(df, "pods", "trt", "Site", adjust = "rain")
  expect_error(pollin_guard_posttreatment(pd2, dag = g_h), "different design")
})

test_that("posttreatment guard passes a graph with no treatment descendant", {
  skip_if_not_installed("dagitty")
  df <- data.frame(
    pods = 1:8,
    trt = rep(c("Cage", "Open"), 4L),
    Site = rep(c("a", "b"), each = 4L),
    rain = rep(c(250, 300), each = 4L),
    veg_weight = 11:18
  )
  pd <- pollin_data(df, "pods", "trt", "Site", adjust = c("rain", "veg_weight"))
  ds <- pollin_design(pd, c(excluded = "Cage", open = "Open"), "multi_site")
  # Branch P: veg_weight measured before flowering, no trt -> veg_weight edge
  g_p <- pollin_dag(
    ds,
    "dag { rain -> pods; veg_weight -> pods; trt -> pods }",
    "pods"
  )
  expect_identical(pollin_guard_posttreatment(ds, dag = g_p), ds)
})

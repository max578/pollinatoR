# Layer L2: pollin_estimand() -------------------------------------------------

l2_design <- function(
  arms = c(excluded = "Cage", open = "Open"),
  reference = "open"
) {
  d <- data.frame(
    pods = rep(c(20, 30, 36), 8),
    trt = rep(c("Cage", "Open", "Hand"), 8),
    Site = rep(letters[1:4], each = 6)
  )
  d <- d[d$trt %in% arms, ]
  pd <- pollin_data(d, "pods", "trt", unit = "Site")
  pollin_design(pd, arms, "multi_site", reference = reference)
}

test_that("a missing basis, scale or family is an error (no defaults)", {
  ds <- l2_design()
  expect_error(pollin_estimand(ds, "dependence", "ratio"), "basis")
  expect_error(pollin_estimand(ds, "dependence", basis = "loss"), "scale")
  expect_error(pollin_estimand(ds, scale = "ratio", basis = "loss"), "family")
  expect_error(
    pollin_estimand(ds, "dependence", "ratio", basis = "share"),
    "`basis` must be one of"
  )
  expect_error(
    pollin_estimand(ds, "dependence", "ratio", "loss", population = "field"),
    "population"
  )
  expect_error(
    pollin_estimand(ds$data, "dependence", "ratio", "loss"),
    "`design`"
  )
})

test_that("dependence contrasts the excluded arm with the reference arm", {
  est <- pollin_estimand(l2_design(), "dependence", "ratio", "loss")
  expect_s3_class(est, "pollin_estimand")
  expect_identical(est$contrast, c(high = "Open", low = "Cage"))
  three <- c(excluded = "Cage", open = "Open", supplemented = "Hand")
  est_s <- pollin_estimand(
    l2_design(three, reference = "supplemented"),
    "dependence",
    "lnrr",
    "gain"
  )
  expect_identical(est_s$contrast, c(high = "Hand", low = "Cage"))
  expect_identical(est_s$basis, "gain")
})

test_that("deficit needs a supplemented arm and contrasts it with open", {
  expect_error(
    pollin_estimand(l2_design(), "deficit", "ratio", "loss"),
    "needs the supplemented arm"
  )
  ds <- l2_design(c(open = "Open", supplemented = "Hand"))
  expect_error(
    pollin_estimand(ds, "dependence", "ratio", "loss"),
    "needs the excluded arm"
  )
  est <- pollin_estimand(ds, "deficit", "difference", "loss")
  expect_identical(est$contrast, c(high = "Hand", low = "Open"))
})

test_that("a sham reference nets out the cage and says so", {
  d <- data.frame(
    pods = rep(c(20, 30, 29), 8),
    trt = rep(c("Cage", "Open", "Sham"), 8),
    Site = rep(letters[1:4], each = 6)
  )
  pd <- pollin_data(d, "pods", "trt", unit = "Site")
  arms <- c(excluded = "Cage", open = "Open", sham = "Sham")
  excl <- function(est) {
    a <- pollin_assumptions(est)
    a$statement[a$assumption == "excludability"]
  }
  ds_sham <- pollin_design(pd, arms, "multi_site", reference = "sham")
  est <- pollin_estimand(ds_sham, "dependence", "ratio", "loss")
  expect_identical(est$contrast, c(high = "Sham", low = "Cage"))
  expect_match(excl(est), "net of the cage effect")
  ds_open <- pollin_design(pd, arms, "multi_site")
  est_o <- pollin_estimand(ds_open, "dependence", "ratio", "loss")
  expect_identical(est_o$contrast, c(high = "Open", low = "Cage"))
  expect_match(excl(est_o), "sham versus open, and this contrast includes it")
  expect_match(
    excl(pollin_estimand(l2_design(), "dependence", "ratio", "loss")),
    "no cage effect on microclimate"
  )
})

test_that("assumptions follow the layout, family and scale", {
  a <- pollin_assumptions(pollin_estimand(
    l2_design(),
    "dependence",
    "ratio",
    "loss"
  ))
  expect_identical(
    a$assumption,
    c(
      "consistency",
      "excludability",
      "no_interference",
      "assignment",
      "positive_means",
      "population",
      "effect"
    )
  )
  expect_match(a$statement[a$assumption == "assignment"], "within each site")
  expect_match(a$statement[a$assumption == "excludability"], "exclusion device")
  ds <- l2_design(c(open = "Open", supplemented = "Hand"))
  a2 <- pollin_assumptions(pollin_estimand(ds, "deficit", "difference", "gain"))
  expect_false("positive_means" %in% a2$assumption)
  expect_match(
    a2$statement[a2$assumption == "excludability"],
    "supplementation"
  )
  expect_error(pollin_assumptions(ds), "pollin_estimand")
})

test_that("the estimand prints its basis", {
  out <- capture.output(print(
    pollin_estimand(l2_design(), "dependence", "ratio", "gain")
  ))
  expect_match(out[1], "<pollin_estimand>")
  expect_true(any(grepl("basis +: gain", out)))
  expect_true(any(grepl("Open vs Cage", out)))
})

test_that("the cage family contrasts sham with open and takes no basis", {
  d <- data.frame(
    pods = rep(c(20, 30, 33), 8),
    trt = rep(c("Cage", "Open", "Sham"), 8),
    Site = rep(letters[1:4], each = 6)
  )
  pd <- pollin_data(d, "pods", "trt", unit = "Site")
  ds <- pollin_design(
    pd,
    c(excluded = "Cage", open = "Open", sham = "Sham"),
    "multi_site"
  )
  est <- pollin_estimand(ds, "cage", "ratio")
  expect_identical(est$contrast, c(high = "Sham", low = "Open"))
  expect_identical(est$basis, NA_character_)
  expect_identical(pollin_estimand(ds, "cage", "lnrr", NULL)$basis, NA_character_)
  a <- pollin_assumptions(est)
  expect_match(a$statement[a$assumption == "excludability"], "sham device")
  expect_error(pollin_estimand(ds, "cage", "ratio", "loss"), "no basis")
  expect_error(pollin_estimand(ds, "cage", "rii"), "scale")
  expect_error(
    pollin_estimand(ds, "cage", "ratio", effect = "direct", mediators = "x"),
    "total effect"
  )
  out <- capture.output(print(est))
  expect_true(any(grepl("cage effect \\(sham versus open\\)", out)))
  expect_true(any(grepl("none \\(cage factor\\)", out)))
  expect_error(
    pollin_estimand(l2_design(), "cage", "ratio"),
    "needs the sham arm"
  )
})

test_that("a direct effect is never the default and says what it assumes", {
  d <- data.frame(
    pods = rep(c(20, 30), 12),
    trt = rep(c("Cage", "Open"), 12),
    Site = rep(letters[1:4], each = 6),
    veg = seq_len(24)
  )
  pd <- pollin_data(d, "pods", "trt", unit = "Site", adjust = "veg")
  ds <- pollin_design(pd, c(excluded = "Cage", open = "Open"), "multi_site")
  tot <- pollin_estimand(ds, "dependence", "ratio", "loss")
  expect_identical(tot$effect, "total")
  expect_null(tot$adjustment)
  expect_false("adjustment" %in% pollin_assumptions(tot)$assumption)
  expect_error(
    pollin_estimand(ds, "dependence", "ratio", "loss", mediators = "veg"),
    "only for a direct effect"
  )
  expect_error(
    pollin_estimand(ds, "dependence", "ratio", "loss", effect = "direct"),
    "needs `mediators`"
  )
  expect_error(
    pollin_estimand(ds, "dependence", "ratio", "loss",
      effect = "direct", mediators = "rain"),
    "declared `adjust` columns.*rain"
  )
  dir <- pollin_estimand(ds, "dependence", "ratio", "loss",
    effect = "direct", mediators = "veg")
  expect_identical(dir$adjustment, "veg")
  a <- pollin_assumptions(dir)
  expect_match(
    a$statement[a$assumption == "effect"],
    "no unmeasured variable causes both veg and the response"
  )
  expect_match(a$statement[a$assumption == "adjustment"], "analyst's declaration")
  expect_true(any(grepl("direct, holding veg fixed", capture.output(print(dir)))))
})

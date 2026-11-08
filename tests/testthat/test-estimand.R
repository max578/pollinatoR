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
      "population"
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

# Layer L0: pollin_data() -----------------------------------------------------

l0_df <- function() {
  data.frame(
    pods = c(30, 25, NA, 40, 22, 35),
    ovules = c(60, 50, 55, 70, 45, 65),
    trt = c("Cage", "Open", "Cage", "Open", "Cage", "Open"),
    Site = c("a", "a", "b", "b", "c", "c"),
    Plot = c("p1", "p1", "p1", "p1", "p2", "p2"),
    Region = c("n", "n", "n", "s", "s", "s"),
    bees = c(3, 5, 2, 6, 1, 4),
    rain = c(300, 300, 280, 280, 310, 310),
    stringsAsFactors = FALSE
  )
}

test_that("pollin_data stores the data, the role map and the checks", {
  pd <- pollin_data(
    l0_df(),
    response = "pods",
    treatment = "trt",
    unit = c("Site", "Plot"),
    groups = "Region",
    denominator = "ovules",
    exposures = "bees",
    adjust = "rain"
  )
  expect_s3_class(pd, "pollin_data")
  expect_identical(pd$data, l0_df())
  expect_identical(pd$n, 6L)
  expect_identical(
    names(pd$roles),
    c(
      "response",
      "treatment",
      "unit",
      "groups",
      "denominator",
      "exposures",
      "adjust"
    )
  )
  expect_identical(pd$roles$unit, c("Site", "Plot"))
  expect_null(pollin_data(l0_df(), "pods", "trt", "Site")$roles$groups)
  expect_identical(pd$checks, pollin_validate(pd)$table)
  expect_identical(pd$checks$check, c("rows", "cols", "NA:pods"))
  expect_identical(pd$checks$detail[3], "1")
})

test_that("pollin_data refuses roles that do not fit the data", {
  d <- l0_df()
  expect_error(pollin_data(as.list(d), "pods", "trt", "Site"), "data frame")
  expect_error(pollin_data(d, "pods", "trt", "Field"), "not found.*Field")
  expect_error(pollin_data(d, "pods", c("trt", "Site"), "Plot"), "one column")
  expect_error(pollin_data(d, character(0), "trt", "Site"), "`response`")
  expect_error(pollin_data(d, "pods", "trt", 1), "`unit`")
  expect_error(
    pollin_data(d, "pods", "trt", "Site", groups = "Site"),
    "one role only.*Site"
  )
  expect_error(pollin_data(d, "Region", "trt", "Site"), "numeric: Region")
  d1 <- d
  d1$trt <- c("Cage", NA, "Cage", "Cage", "Cage", "Cage")
  expect_error(pollin_data(d1, "pods", "trt", "Site"), "two distinct.*found 1")
  d2 <- d
  d2$ovules[2] <- -1
  expect_error(
    pollin_data(d2, "pods", "trt", "Site", denominator = "ovules"),
    "non-negative"
  )
  expect_error(
    pollin_data(d, "pods", "trt", "Site", denominator = c("ovules", "bees")),
    "one column"
  )
})

test_that("pollin_validate takes a pollin_data with levels and identities", {
  pd <- pollin_data(l0_df(), "pods", "trt", "Site", denominator = "ovules")
  v <- pollin_validate(
    pd,
    expect_levels = list(trt = c("Cage", "Open")),
    identities = list(pods_le_ovules = function(d) d$pods <= d$ovules)
  )
  expect_s3_class(v, "pollin_validate")
  tab <- v$table
  expect_identical(tab$status[tab$check == "levels:trt"], "PASS")
  expect_identical(tab$status[tab$check == "identity:pods_le_ovules"], "PASS")
  expect_error(pollin_validate(l0_df()), "pollin_data")
})

test_that("pollin_data prints its roles and the warning count", {
  pd <- pollin_data(l0_df(), "pods", "trt", "Site", exposures = "bees")
  out <- capture.output(print(pd))
  expect_match(out[1], "<pollin_data>")
  expect_true(any(grepl("exposures +: bees", out)))
  expect_true(any(grepl("1 warning", out)))
})

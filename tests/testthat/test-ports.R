# Step 7: index, RII and compensation accept the L0/L1 objects ----------------

port_design <- function(response = c("pods", "seed_wt")) {
  g <- golden_inputs()
  pd <- pollin_data(g$plants, response, "trt", unit = c("Site", "Plot"))
  pollin_design(pd, c(excluded = "Cage", open = "Open"), "multi_site")
}

test_that("a design gives the golden index and RII from site-level means", {
  gold <- readRDS(golden_path())
  ds <- port_design()
  expect_identical(
    pollin_index(ds, response = "pods", n_boot = 999L, seed = 11L),
    gold$index
  )
  expect_identical(
    unclass(pollin_rii(ds, response = "pods", n_boot = 999L, seed = 11L)),
    gold$rii
  )
})

test_that("arm means are taken at the requested unit level", {
  g <- golden_inputs()
  ds <- port_design()
  m <- .arm_means(ds, "seed_wt", unit = "Plot")
  expect_length(m$open, 24L)
  ref <- stats::aggregate(seed_wt ~ Site + Plot + trt, g$plants, mean)
  ref <- ref[order(ref$Site, ref$Plot), ]
  expect_equal(m$open, ref$seed_wt[ref$trt == "Open"])
  expect_equal(m$excluded, ref$seed_wt[ref$trt == "Cage"])
  expect_identical(m$unit[1:2], c("s01:p1", "s01:p2"))
})

test_that("the design path validates its arguments", {
  ds <- port_design()
  expect_error(pollin_index(ds, n_boot = 10L), "several responses")
  expect_error(pollin_index(ds, "pods"), "leave `cage` empty")
  expect_error(pollin_rii(ds, response = "rain"), "declared responses")
  expect_error(pollin_rii(ds, response = "pods", unit = "Region"), "unit col")
  three <- golden_inputs()$plants
  three$trt[three$trt == "Cage" & three$Plot == "p2"] <- "Hand"
  pd <- pollin_data(three, "pods", "trt", unit = "Site")
  sup <- pollin_design(
    pd,
    c(excluded = "Cage", open = "Open", supplemented = "Hand"),
    "multi_site"
  )
  expect_s3_class(pollin_rii(sup, n_boot = 10L, seed = 1L), "pollin_rii")
})

test_that("the design path needs an excluded and an open arm", {
  d <- golden_inputs()$plants
  d$trt[d$trt == "Open"] <- "Hand"
  pd <- pollin_data(d, "pods", "trt", unit = "Site")
  ds <- pollin_design(
    pd,
    c(excluded = "Cage", supplemented = "Hand"),
    "multi_site",
    reference = "supplemented"
  )
  expect_error(pollin_rii(ds), "missing: open")
})

test_that("units lacking an arm are dropped with a count", {
  d <- golden_inputs()$plants
  d <- d[!(d$Site == "s01" & d$trt == "Cage"), ]
  pd <- pollin_data(d, "pods", "trt", unit = "Site")
  ds <- suppressWarnings(
    pollin_design(pd, c(excluded = "Cage", open = "Open"), "multi_site")
  )
  expect_warning(m <- .arm_means(ds), "1 unit of `Site`")
  expect_length(m$open, 11L)
  tiny <- d[d$Site %in% c("s01", "s02"), ]
  pd <- pollin_data(tiny, "pods", "trt", unit = "Site")
  ds <- suppressWarnings(
    pollin_design(pd, c(excluded = "Cage", open = "Open"), "multi_site")
  )
  expect_error(suppressWarnings(.arm_means(ds)), "at least two units")
})

test_that("compensation takes the responses of a pollin_data object", {
  g <- golden_inputs()
  resp <- c("pods", "seeds_per_pod", "seed_wt")
  pd <- pollin_data(g$plants, resp, "trt", unit = "Site")
  expect_identical(
    pollin_compensation(pd),
    pollin_compensation(g$components)
  )
  one <- pollin_data(g$plants, "pods", "trt", unit = "Site")
  expect_error(pollin_compensation(one), "at least two responses")
  expect_output(print(pollin_compensation(pd)), "negatively-correlated pair")
})

# Layer L1: pollin_design() ---------------------------------------------------

# Two blocks/sites x two whole plots x two pairs x both arms, two plants each.
l1_df <- function() {
  d <- expand.grid(
    plant = 1:2,
    trt = c("Cage", "Open"),
    Pair = c("q1", "q2"),
    Plot = c("p1", "p2"),
    Site = c("a", "b"),
    stringsAsFactors = FALSE
  )
  d$pods <- seq_len(nrow(d))
  d
}
l1_arms <- c(excluded = "Cage", open = "Open")

test_that("each layout builds the random terms of the equivalent spec", {
  d <- l1_df()
  unit <- c("Site", "Plot", "Pair")
  pd <- pollin_data(d, "pods", "trt", unit = unit)
  for (lay in c("rcbd", "split_plot", "matched_pairs", "multi_site")) {
    for (s in unit) {
      ds <- pollin_design(pd, l1_arms, layout = lay, slope = s)
      sp <- pollin_spec(d, "pods", "trt", "Cage", unit, latent_unit = s)
      expect_identical(ds$random, .random_terms(sp))
    }
  }
  expect_identical(
    pollin_design(pd, l1_arms, "multi_site")$random,
    .random_terms(pollin_spec(d, "pods", "trt", "Cage", unit))
  )
})

test_that("layouts set their default slope and intercept-only terms", {
  pd <- pollin_data(l1_df(), "pods", "trt", unit = c("Site", "Plot"))
  expect_identical(
    pollin_design(pd, l1_arms, "multi_site")$random,
    c("(trt | Site)", "(1 | Site:Plot)")
  )
  expect_identical(
    pollin_design(pd, l1_arms, "multi_site", slope = FALSE)$random,
    c("(1 | Site)", "(1 | Site:Plot)")
  )
  expect_identical(
    pollin_design(pd, l1_arms, "split_plot")$random,
    c("(1 | Site)", "(1 | Site:Plot)")
  )
  expect_null(pollin_design(pd, l1_arms, "rcbd")$slope)
})

test_that("a crd needs the treatment constant within each unit", {
  d <- l1_df()
  d$Unit <- paste(d$Site, d$Plot, d$Pair, d$trt)
  pd <- pollin_data(d, "pods", "trt", unit = "Unit")
  ds <- pollin_design(pd, l1_arms, "crd")
  expect_identical(ds$random, "(1 | Unit)")
  pd_bad <- pollin_data(d, "pods", "trt", unit = "Site")
  expect_error(pollin_design(pd_bad, l1_arms, "crd"), "2 innermost unit")
})

test_that("a unit lacking an arm is counted in a warning", {
  d <- l1_df()
  d <- d[!(d$Site == "b" & d$trt == "Open"), ]
  pd <- pollin_data(d, "pods", "trt", unit = c("Site", "Plot"))
  expect_warning(
    pollin_design(pd, l1_arms, "multi_site"),
    "1 unit\\(s\\) of `Site` lack"
  )
  expect_warning(
    pollin_design(pd, l1_arms, "split_plot"),
    "2 unit\\(s\\) of `Site:Plot` lack"
  )
  expect_silent(pollin_design(
    pollin_data(l1_df(), "pods", "trt", unit = "Site"),
    l1_arms,
    "rcbd"
  ))
})

test_that("pollin_design refuses a mis-declared design", {
  pd <- pollin_data(l1_df(), "pods", "trt", unit = c("Site", "Plot"))
  expect_error(pollin_design(l1_df(), l1_arms, "rcbd"), "pollin_data")
  expect_error(pollin_design(pd, l1_arms), "layout")
  expect_error(pollin_design(pd, l1_arms, "latin"), "one of")
  expect_error(pollin_design(pd, c("Cage", "Open"), "rcbd"), "named")
  expect_error(
    pollin_design(pd, c(caged = "Cage", open = "Open"), "rcbd"),
    "Unknown arm role"
  )
  expect_error(
    pollin_design(pd, c(excluded = "Cage", open = "Cage"), "rcbd"),
    "once"
  )
  expect_error(pollin_design(pd, c(open = "Open"), "rcbd"), "two arms")
  expect_error(
    pollin_design(pd, c(excluded = "Cage", open = "Opn"), "rcbd"),
    "not found.*Opn"
  )
  expect_error(
    pollin_design(pd, c(excluded = "Cage", supplemented = "Open"), "rcbd"),
    "reference arm 'open'"
  )
  expect_error(
    pollin_design(pd, l1_arms, "rcbd", reference = "excluded"),
    "reference"
  )
  expect_error(pollin_design(pd, l1_arms, "rcbd", slope = "Pair"), "slope")
  pd1 <- pollin_data(l1_df(), "pods", "trt", unit = "Site")
  expect_error(pollin_design(pd1, l1_arms, "split_plot"), "whole")
  d3 <- l1_df()
  d3$trt[1] <- "Hand"
  pd3 <- pollin_data(d3, "pods", "trt", unit = "Site")
  expect_error(pollin_design(pd3, l1_arms, "rcbd"), "no arm: Hand")
})

test_that("supplemented reference and printing work", {
  d <- l1_df()
  d$trt[d$trt == "Open" & d$plant == 2] <- "Hand"
  pd <- pollin_data(d, "pods", "trt", unit = "Site")
  ds <- pollin_design(
    pd,
    c(excluded = "Cage", open = "Open", supplemented = "Hand"),
    "rcbd",
    reference = "supplemented"
  )
  expect_identical(ds$reference, "supplemented")
  out <- capture.output(print(ds))
  expect_match(out[1], "<pollin_design>")
  expect_true(any(grepl("supplemented = Hand", out)))
})

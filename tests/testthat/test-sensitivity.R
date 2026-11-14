# Layer L4: pollin_sensitivity() ----------------------------------------------
#
# Oracle: `EValue::evalues.RR()`, which takes the ratio and its interval and
# returns the E-values of the estimate and of the limit closer to the null.

ev_effects <- function(scale = "ratio", basis = "loss", effect = "total") {
  est <- fd_estimand(effect)
  est <- pollin_estimand(
    est$design,
    "dependence",
    scale,
    basis,
    effect = effect,
    mediators = if (effect == "direct") "veg"
  )
  pollin_effects(pollin_fit(est))
}

ev_oracle <- function(rr, lo, hi) {
  r <- EValue::evalues.RR(rr, lo, hi)
  c(point = r["E-values", "point"], ci = min(r["E-values", ], na.rm = TRUE))
}

test_that("E-values match EValue on every ratio-type scale", {
  skip_on_cran()
  skip_if_not_installed("EValue")
  for (sc in list(c("ratio", "loss"), c("ratio", "gain"), c("lnrr", "loss"),
    c("rii", "loss"))) {
    s <- pollin_sensitivity(ev_effects(sc[1L], sc[2L]))$table
    o <- ev_oracle(s$rr, s$rr.low, s$rr.high)
    expect_equal(s$evalue, unname(o["point"]), tolerance = 1e-10,
      info = sc[1L])
    expect_equal(s$evalue.ci, unname(o["ci"]), tolerance = 1e-10,
      info = sc[1L])
  }
  # The ratio is recovered from the scale: every scale gives one ratio
  r1 <- pollin_sensitivity(ev_effects("ratio", "loss"))$table$rr
  r2 <- pollin_sensitivity(ev_effects("lnrr", "loss"))$table$rr
  expect_equal(r1, r2, tolerance = 1e-10)
})

test_that("an interval containing the null has an E-value of 1", {
  skip_if_not_installed("EValue")
  expect_identical(.evalue(1), 1)
  expect_equal(.evalue(0.8), 1.25 + sqrt(1.25 * 0.25))
  o <- EValue::evalues.RR(0.8, 0.6, 1.1)
  expect_equal(.evalue(0.8), o["E-values", "point"])
})

test_that("the difference scale, a direct effect and continuous data are refused", {
  skip_on_cran()
  expect_error(
    pollin_sensitivity(ev_effects("difference", "loss")),
    "difference scale carries no ratio"
  )
  expect_error(
    pollin_sensitivity(ev_effects(effect = "direct")),
    "direct effect is refused"
  )
  e <- ev_effects()
  e$table$distribution <- "Gamma"
  expect_error(pollin_sensitivity(e), "fitted as Gamma is neither")
  expect_error(pollin_sensitivity(1), "`pollin_effects` or `pollin_moderation`")
})

test_that("a moderation term gets its link-scale bias to the null", {
  skip_on_cran()
  m <- pollin_moderation(fd_moderation_design(2L, conf = 0.5), "y")
  s <- pollin_sensitivity(m, conf_level = 0.9)
  expect_identical(s$type, "moderation")
  expect_identical(s$table$bias.estimate, m$table$estimate)
  z <- stats::qnorm(0.95)
  lo <- m$table$estimate - z * m$table$std.error
  hi <- m$table$estimate + z * m$table$std.error
  expect_equal(s$table$bias.ci, ifelse(lo > 0, lo, ifelse(hi < 0, hi, 0)))
  expect_identical(s$negative_control, m$negative_control)
  expect_output(print(s), "negative control")
})

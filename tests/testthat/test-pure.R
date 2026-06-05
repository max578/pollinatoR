# Tests for the self-contained (model-free) core ------------------------------

test_that("within-between decomposition splits correctly", {
  wb <- pollin_within_between(c(1, 3, 10, 12), c("a", "a", "b", "b"))
  expect_equal(wb$between, c(2, 2, 11, 11))
  expect_equal(wb$within, c(-1, 1, -1, 1))
  expect_equal(wb$within + wb$between, c(1, 3, 10, 12))
})

test_that("pure coupling: a constant exclusion arm makes the slope mechanical", {
  open <- c(10, 20, 30, 40, 50, 60)
  cage <- rep(10, 6)
  g <- pollin_guard_spurious(open, cage, n_perm = 50, seed = 1)
  expect_s3_class(g, "pollin_spurious")
  expect_gt(g$observed_slope, 0) # RII rises with open by construction
  expect_equal(g$perm_p, 1) # permuting a constant changes nothing
})

test_that("spurious guard accepts a decoupled axis", {
  set.seed(3)
  g <- pollin_guard_spurious(
    rpois(30, 20) + 1,
    rpois(30, 18) + 1,
    decouple = rnorm(30),
    n_perm = 100,
    seed = 3
  )
  expect_false(is.null(g$decoupled))
  expect_true(g$perm_p >= 0 && g$perm_p <= 1)
})

test_that("RII matches its definition and the bootstrap CI brackets it", {
  r <- pollin_rii(rep(8, 20), rep(6, 20), n_boot = 200, seed = 1)
  expect_equal(r$rii, (8 - 6) / (8 + 6))
  expect_true(r$conf.low <= r$rii && r$rii <= r$conf.high)
})

test_that("pollin_index returns the three indices", {
  idx <- pollin_index(
    rpois(40, 9) + 1,
    rpois(40, 6) + 1,
    n_boot = 200,
    seed = 1
  )
  expect_equal(nrow(idx), 3L)
  expect_true(all(c("RII", "lnRR", "dependence_ratio") %in% idx$index))
})

test_that("fdr wrapper matches BH", {
  p <- c(0.01, 0.02, 0.5, 0.9)
  expect_equal(pollin_fdr(p), stats::p.adjust(p, "BH"))
})

test_that("compensation flags a negative component correlation", {
  m <- data.frame(pods = c(10, 8, 6, 4), seeds = c(20, 24, 28, 32))
  cp <- pollin_compensation(m)
  expect_s3_class(cp, "pollin_compensation")
  expect_true(cp$pairs$compensation[1])
})

test_that("spec validates columns and control level", {
  df <- data.frame(
    y = rpois(20, 5),
    trt = rep(c("Cage", "Open"), 10),
    Site = rep(letters[1:4], each = 5)
  )
  sp <- pollin_spec(df, "y", "trt", "Cage", "Site")
  expect_s3_class(sp, "pollin_spec")
  expect_error(pollin_spec(df, "nope", "trt", "Cage", "Site"))
  expect_error(pollin_spec(df, "y", "trt", "Nonexist", "Site"))
})

test_that("validate reports levels and user identities", {
  df <- data.frame(
    y = c(1, 2, 3, 4),
    trt = c("Cage", "Open", "Cage", "Open"),
    Site = c("a", "a", "b", "b"),
    tot = c(3, 3, 7, 7),
    p = c(1, 2, 3, 4),
    q = c(2, 1, 4, 3)
  )
  sp <- pollin_spec(df, "y", "trt", "Cage", "Site")
  v <- pollin_validate(
    sp,
    expect_levels = list(trt = c("Cage", "Open")),
    identities = list(tot_eq = function(d) d$tot == d$p + d$q)
  )
  expect_s3_class(v, "pollin_validate")
  expect_true(any(v$table$check == "identity:tot_eq"))
  expect_equal(v$table$status[v$table$check == "identity:tot_eq"], "PASS")
})

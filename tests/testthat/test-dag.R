# Layer L2: pollin_dag() and the graph route of pollin_estimand() -------------
#
# Criterion C2. For each graph the adjustment sets and classes below are
# derived by hand (the oracle); the package must return them, and its sets
# must equal those of `dagitty::adjustmentSets()` called directly.

dag_design <- function(cols, adjust = cols, exposures = NULL) {
  set.seed(3)
  d <- data.frame(
    pods = rpois(48, 30),
    trt = rep(c("Cage", "Open"), 24),
    Site = rep(letters[1:4], each = 12)
  )
  for (v in cols) {
    d[[v]] <- rnorm(48)
  }
  pd <- pollin_data(
    d,
    "pods",
    "trt",
    unit = "Site",
    adjust = adjust,
    exposures = exposures
  )
  pollin_design(pd, c(excluded = "Cage", open = "Open"), "multi_site")
}

class_of <- function(g, target = "trt") {
  cl <- g$classes[g$classes$target == target, ]
  stats::setNames(cl$class, cl$node)[order(cl$node)]
}

same_as_dagitty <- function(g, exposure, effect, sets) {
  ref <- dagitty::adjustmentSets(g$graph, exposure, g$outcome, effect = effect)
  ref <- lapply(unclass(ref), function(z) sort(as.character(unlist(z))))
  expect_identical(unname(ref), sets)
}

test_that("confounder: c causes the arm and the yield", {
  skip_if_not_installed("dagitty")
  # Back-door path trt <- c -> yield is open; {c} blocks it and is minimal.
  # With no mediator the direct-effect set is the same. c is a confounder.
  g <- pollin_dag(
    dag_design("c"),
    "dag { c -> trt; c -> yield; trt -> yield }",
    "yield"
  )
  expect_identical(g$sets$total, list("c"))
  expect_identical(g$sets$direct, list("c"))
  expect_identical(class_of(g), c(c = "confounder"))
  same_as_dagitty(g, "trt", "total", g$sets$total)
  same_as_dagitty(g, "trt", "direct", g$sets$direct)
})

test_that("mediator: m carries part of the effect; p, d and n are bystanders", {
  skip_if_not_installed("dagitty")
  # No back-door path: the empty set is valid for the total effect. The direct
  # effect must block trt -> m -> yield, so {m}. m: descendant of trt and
  # ancestor of yield (mediator). d: descendant, not on a path to yield. p:
  # adding it keeps {} valid, independent of trt, associated with yield given
  # trt (precision). n: isolated (neutral).
  g <- pollin_dag(
    dag_design(c("m", "p", "d", "n")),
    "dag { trt -> m -> yield; trt -> yield; p -> yield; trt -> d; n }",
    "yield"
  )
  expect_identical(g$sets$total, list(character(0)))
  expect_identical(g$sets$direct, list("m"))
  expect_identical(
    class_of(g),
    c(d = "descendant", m = "mediator", n = "neutral", p = "precision")
  )
  same_as_dagitty(g, "trt", "total", g$sets$total)
  same_as_dagitty(g, "trt", "direct", g$sets$direct)
})

test_that("M-bias: conditioning on the collider m opens a path", {
  skip_if_not_installed("dagitty")
  # trt <- a -> m <- b -> yield is blocked at the collider m, so {} is valid and
  # minimal; {m} opens it (collider). a: associated with trt, and with trt's
  # out-edges cut it is independent of yield (instrument). b: independent of
  # trt, associated with yield given trt (precision).
  g <- pollin_dag(
    dag_design(c("a", "b", "m")),
    "dag { a -> trt; a -> m; b -> m; b -> yield; trt -> yield }",
    "yield"
  )
  expect_identical(g$sets$total, list(character(0)))
  expect_identical(
    class_of(g),
    c(a = "instrument", b = "precision", m = "collider")
  )
  same_as_dagitty(g, "trt", "total", g$sets$total)
})

test_that("mediator-outcome confounder: latent u leaves no direct-effect set", {
  skip_if_not_installed("dagitty")
  # Holding m fixed opens trt -> m <- u -> yield. With u latent nothing closes
  # it: no direct-effect set. With u measured, {m, u} is the set, and u is
  # independent of trt and associated with yield (precision).
  graph <- "dag { trt -> m -> yield; trt -> yield; u -> m; u -> yield }"
  g_lat <- pollin_dag(
    dag_design("m"),
    sub("u -> m", "u [latent]; u -> m", graph, fixed = TRUE),
    "yield"
  )
  expect_identical(g_lat$sets$total, list(character(0)))
  expect_identical(g_lat$sets$direct, list())
  expect_identical(class_of(g_lat), c(m = "mediator"))
  same_as_dagitty(g_lat, "trt", "direct", list())
  g_obs <- pollin_dag(dag_design(c("m", "u")), graph, "yield")
  expect_identical(g_obs$sets$direct, list(c("m", "u")))
  expect_identical(class_of(g_obs), c(m = "mediator", u = "precision"))
  same_as_dagitty(g_obs, "trt", "direct", g_obs$sets$direct)
})

test_that("instrument: z moves the arm only", {
  skip_if_not_installed("dagitty")
  # {c} blocks trt <- c -> yield. z is in no minimal set, {c, z} stays valid,
  # z is associated with trt and, with trt's out-edges cut, independent of
  # yield given c (instrument).
  g <- pollin_dag(
    dag_design(c("c", "z")),
    "dag { z -> trt; c -> trt; c -> yield; trt -> yield }",
    "yield"
  )
  expect_identical(g$sets$total, list("c"))
  expect_identical(class_of(g), c(c = "confounder", z = "instrument"))
  same_as_dagitty(g, "trt", "total", g$sets$total)
})

test_that("field graph: branch H refuses a direct effect, P is precision", {
  skip_if_not_installed("dagitty")
  # Branch H (veg weighed at harvest): trt -> veg -> yield, latent fertility
  # causes veg and yield, forest -> fertility and forest -> bees -> yield. The
  # exposure set for bees is {forest}. Branch P (veg before flowering): no
  # trt -> veg edge, veg is pre-treatment and adjusting is precision.
  ds <- dag_design(c("veg", "rain", "forest", "bees"),
    adjust = c("veg", "rain"), exposures = "bees")
  h <- paste(
    "dag { fertility [latent]; trt -> yield; trt -> veg; veg -> yield;",
    "rain -> veg; rain -> yield; fertility -> veg; fertility -> yield;",
    "forest -> fertility; forest -> bees; bees -> yield }"
  )
  g_h <- pollin_dag(ds, h, "yield")
  expect_identical(g_h$sets$direct, list())
  expect_identical(g_h$exposure_sets$bees, list("forest"))
  same_as_dagitty(g_h, "bees", "total", g_h$exposure_sets$bees)
  expect_error(
    pollin_estimand(ds, "dependence", "ratio", "loss",
      effect = "direct", mediators = "veg", dag = g_h),
    "no direct effect.*latent fertility"
  )
  est <- pollin_estimand(ds, "dependence", "ratio", "loss", dag = g_h)
  expect_identical(est$effect, "total")
  g_p <- pollin_dag(ds, sub("trt -> veg; ", "", h, fixed = TRUE), "yield")
  expect_identical(class_of(g_p)[["veg"]], "precision")
})

test_that("pollin_dag checks the nodes and the edges", {
  skip_if_not_installed("dagitty")
  ds <- dag_design("veg")
  expect_error(
    pollin_dag(ds, "dag { trt -> yield }", "yield"),
    "missing from the graph: veg"
  )
  expect_error(
    pollin_dag(ds, "dag { trt -> yield; veg -> yield; soil -> veg }", "yield"),
    "neither columns of the data nor latent: soil"
  )
  expect_error(
    pollin_dag(ds, "dag { trt -> yield; veg <-> yield }", "yield"),
    "latent node"
  )
  expect_error(pollin_dag(ds, "dag { trt -> veg }", "trt"), "missing|treatment")
  expect_error(pollin_dag(ds, 3, "yield"), "dagitty")
  g <- pollin_dag(ds, "dag { trt -> veg -> yield; trt -> yield }", "yield")
  expect_s3_class(g, "pollin_dag")
  out <- capture.output(print(g))
  expect_match(out[1], "<pollin_dag>  trt -> yield")
  expect_true(any(grepl("direct effect : \\{veg\\}", out)))
})

test_that("pollin_dag_check tests the implied independencies", {
  skip_if_not_installed("dagitty")
  ds <- dag_design(c("rain", "forest"))
  g <- pollin_dag(
    ds,
    "dag { trt -> yield; rain -> yield; forest -> yield }",
    "yield"
  )
  set.seed(5)
  site <- data.frame(rain = rnorm(30), forest = rnorm(30))
  ck <- pollin_dag_check(g, site)
  expect_identical(
    names(ck),
    c("test", "estimate", "conf.low", "conf.high", "p.value", "n")
  )
  expect_identical(ck$test, "forest _||_ rain")
  expect_identical(ck$n, 30L)
  expect_true(ck$conf.low < ck$estimate && ck$estimate < ck$conf.high)
  site$forest <- as.character(site$forest)
  expect_error(pollin_dag_check(g, site), "numeric: forest")
})

test_that("a direct effect needs mediators, the adjust set and a graph set", {
  skip_if_not_installed("dagitty")
  ds <- dag_design(c("veg", "rain"))
  g <- pollin_dag(
    ds,
    "dag { trt -> veg -> yield; trt -> yield; rain -> veg; rain -> yield }",
    "yield"
  )
  est <- pollin_estimand(ds, "dependence", "ratio", "loss",
    effect = "direct", mediators = "veg", dag = g)
  expect_identical(est$adjustment, c("rain", "veg"))
  a <- pollin_assumptions(est)
  expect_match(
    a$statement[a$assumption == "effect"],
    "Direct effect holding veg"
  )
  expect_match(a$statement[a$assumption == "adjustment"], "\\{rain, veg\\}")
  expect_error(
    pollin_estimand(ds, "dependence", "ratio", "loss",
      effect = "direct", mediators = "rain", dag = g),
    "Not on a path.*rain"
  )
  ds_v <- dag_design(c("veg", "rain"), adjust = "veg")
  g_v <- pollin_dag(
    ds_v,
    "dag { trt -> veg -> yield; trt -> yield; rain -> veg; rain -> yield }",
    "yield"
  )
  expect_error(
    pollin_estimand(ds_v, "dependence", "ratio", "loss",
      effect = "direct", mediators = "veg", dag = g_v),
    "declare one of: \\{rain, veg\\}"
  )
  expect_error(
    pollin_estimand(ds_v, "dependence", "ratio", "loss", dag = g),
    "different design"
  )
  g_c <- pollin_dag(
    dag_design("c"),
    "dag { c -> trt; c -> yield; trt -> yield }",
    "yield"
  )
  expect_error(
    pollin_estimand(dag_design("c"), "dependence", "ratio", "loss", dag = g_c),
    "common cause"
  )
})

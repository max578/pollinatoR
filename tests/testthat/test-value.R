# Layer L6: pollin_market() and pollin_value() ---------------------------------

va_design <- function(reference = "open") {
  d <- data.frame(
    pods = rep(c(20, 30, 36), 8L),
    trt = rep(c("Cage", "Open", "Hand"), 8L),
    Block = rep(sprintf("b%d", 1:4), each = 6L),
    stringsAsFactors = FALSE
  )
  pd <- pollin_data(d, "pods", "trt", unit = "Block")
  pollin_design(
    pd,
    c(excluded = "Cage", open = "Open", supplemented = "Hand"),
    "rcbd",
    reference = reference
  )
}

va_estimand <- function(
  family = "dependence",
  scale = "ratio",
  basis = "loss",
  reference = "open"
) {
  pollin_estimand(va_design(reference), family, scale, basis)
}

# Draws built from known arm means, one column per level
va_draws <- function(
  estimand = va_estimand(),
  high = cbind(overall = c(30, 32, 28, 31)),
  low = cbind(overall = c(20, 24, 21, 18)),
  population = NULL,
  by = NULL
) {
  .new_pollin_draws(
    high = high,
    low = low,
    estimand = estimand,
    response = "pods",
    by = by,
    method = "posterior",
    engine = "test",
    n = nrow(high),
    seed = NULL,
    population = population
  )
}

va_market <- function(...) {
  pollin_market(1000, 500, "t", "AUD", 2024, "Synthetic inputs", ...)
}

test_that("pollin_market records and prints its inputs", {
  m <- va_market(price_draws = c(450, 500, 550, 600))
  expect_s3_class(m, "pollin_market")
  expect_identical(m$year, 2024L)
  expect_identical(m$unit, "t")
  expect_null(m$quantity_draws)
  out <- capture.output(print(m))
  expect_true(any(grepl("Synthetic inputs", out)))
  expect_true(any(grepl("AUD per t", out)))
  expect_true(any(grepl("4 draws", out)))
})

test_that("pollin_market refuses non-positive or undeclared inputs", {
  mk <- function(...) {
    args <- list(
      quantity = 1000,
      price = 500,
      unit = "t",
      currency = "AUD",
      year = 2024,
      source = "s"
    )
    do.call(pollin_market, utils::modifyList(args, list(...)))
  }
  expect_error(mk(quantity = 0), "`quantity` must be one positive")
  expect_error(mk(quantity = -1), "`quantity` must be one positive")
  expect_error(mk(price = Inf), "`price` must be one positive")
  expect_error(mk(price = c(1, 2)), "`price` must be one positive")
  expect_error(mk(unit = ""), "`unit` must be one non-empty string")
  expect_error(mk(currency = 1), "`currency` must be one non-empty string")
  expect_error(mk(source = NA_character_), "`source`")
  expect_error(mk(year = 2024.5), "`year` must be one whole number")
  expect_error(mk(year = c(2023, 2024)), "`year`")
  expect_error(mk(price_draws = c(1, 0, 2)), "`price_draws`")
  expect_error(mk(quantity_draws = c(1, -2)), "`quantity_draws`")
  expect_error(mk(quantity_draws = c(1, NA)), "`quantity_draws`")
})

test_that("dependence value draws equal D x P x Q exactly", {
  dr <- va_draws()
  v <- pollin_value(dr, va_market(), channel = "dependence")
  d <- 1 - dr$low[, 1L] / dr$high[, 1L]
  expect_identical(v$draws, d * 500 * 1000)
  expect_s3_class(v, "pollin_value")
  expect_identical(v$level, "overall")
  expect_identical(
    names(v),
    c(
      "draws",
      "summary",
      "channel",
      "estimand",
      "level",
      "market",
      "assumptions",
      "conf_level",
      "attribution"
    )
  )
  s <- summary(v)
  expect_identical(
    names(s),
    c(
      "channel",
      "level",
      "median",
      "conf.low",
      "conf.high",
      "currency",
      "year",
      "unit",
      "basis",
      "population",
      "method",
      "n"
    )
  )
  expect_identical(s$median, stats::median(v$draws))
  expect_identical(
    c(s$conf.low, s$conf.high),
    stats::quantile(v$draws, c(0.025, 0.975), names = FALSE)
  )
  expect_identical(s$method, "posterior")
  expect_identical(s$basis, "loss")
  expect_identical(s$population, "trial")
  expect_identical(s$n, 4L)
})

test_that("deficit value draws equal (supp/open - 1) x P x Q exactly", {
  e <- va_estimand("deficit", basis = "gain")
  dr <- va_draws(
    e,
    high = cbind(overall = c(36, 38, 35)),
    low = cbind(overall = c(30, 31, 29))
  )
  v <- pollin_value(dr, va_market(), channel = "deficit")
  expect_identical(v$draws, (dr$high[, 1L] / dr$low[, 1L] - 1) * 500 * 1000)
  expect_identical(v$summary$basis, "gain")
  out <- capture.output(print(v))
  expect_true(any(grepl("closing the gap to hand pollination", out)))
})

test_that("market draws enter per draw, paired by index", {
  dr <- va_draws()
  p <- c(400, 500, 600, 700)
  q <- c(900, 1000, 1100, 1200)
  v <- pollin_value(
    dr,
    va_market(price_draws = p, quantity_draws = q),
    channel = "dependence"
  )
  expect_identical(v$draws, (1 - dr$low[, 1L] / dr$high[, 1L]) * p * q)
  v2 <- pollin_value(dr, va_market(price_draws = p), channel = "dependence")
  expect_identical(v2$draws, (1 - dr$low[, 1L] / dr$high[, 1L]) * p * 1000)
})

test_that("a gain basis is refused by the dependence channel", {
  dr <- va_draws(va_estimand(basis = "gain"))
  expect_error(
    pollin_value(dr, va_market(), channel = "dependence"),
    "overstates the loss by the factor 1/\\(1 - D\\)"
  )
  expect_error(
    pollin_value(dr, va_market(), channel = "dependence"),
    "basis = \"loss\"",
    fixed = TRUE
  )
})

test_that("the wrong family is refused, with the call for the right one", {
  dr <- va_draws(va_estimand("deficit", basis = "gain"))
  expect_error(
    pollin_value(dr, va_market(), channel = "dependence"),
    "pollin_estimand(design, \"dependence\", \"ratio\", \"loss\")",
    fixed = TRUE
  )
  dr2 <- va_draws(va_estimand())
  expect_error(
    pollin_value(dr2, va_market(), channel = "deficit"),
    "pollin_estimand(design, \"deficit\", \"ratio\", \"gain\")",
    fixed = TRUE
  )
})

test_that("a cage estimand is refused as not a pollination quantity", {
  d <- data.frame(
    pods = rep(c(20, 30, 29), 8L),
    trt = rep(c("Cage", "Open", "Sham"), 8L),
    Block = rep(sprintf("b%d", 1:4), each = 6L),
    stringsAsFactors = FALSE
  )
  pd <- pollin_data(d, "pods", "trt", unit = "Block")
  ds <- pollin_design(
    pd,
    c(excluded = "Cage", open = "Open", sham = "Sham"),
    "rcbd"
  )
  dr <- va_draws(pollin_estimand(ds, "cage", "ratio"))
  for (ch in c("dependence", "deficit")) {
    expect_error(
      pollin_value(dr, va_market(), channel = ch),
      "not a pollination quantity"
    )
  }
  # The sham-referenced dependence estimand of the same design is valued
  ds_s <- pollin_design(
    pd,
    c(excluded = "Cage", open = "Open", sham = "Sham"),
    "rcbd",
    reference = "sham"
  )
  v <- pollin_value(
    va_draws(pollin_estimand(ds_s, "dependence", "ratio", "loss")),
    va_market(),
    channel = "dependence"
  )
  expect_s3_class(v, "pollin_value")
})

test_that("a direct effect is refused; the total effect is valued", {
  e <- va_estimand()
  expect_identical(e$effect, "total")
  expect_s3_class(
    pollin_value(va_draws(e), va_market(), channel = "dependence"),
    "pollin_value"
  )
  # A direct estimand cannot be fitted yet, so its label is set by hand
  e$effect <- "direct"
  e$mediators <- "veg_weight"
  expect_error(
    pollin_value(va_draws(e), va_market(), channel = "dependence"),
    "holds `veg_weight` fixed.*value the total effect"
  )
})

test_that("a scale other than the ratio is refused", {
  for (s in c("difference", "rii", "lnrr")) {
    dr <- va_draws(va_estimand(scale = s))
    expect_error(
      pollin_value(dr, va_market(), channel = "dependence"),
      paste0("on the `", s, "` scale")
    )
  }
  dr <- va_draws(va_estimand("deficit", scale = "lnrr", basis = "gain"))
  expect_error(
    pollin_value(dr, va_market(), channel = "deficit"),
    "on the `lnrr` scale"
  )
})

test_that("a loss basis is refused by the deficit channel", {
  dr <- va_draws(va_estimand("deficit", basis = "loss"))
  expect_error(
    pollin_value(dr, va_market(), channel = "deficit"),
    "Refit the estimand with basis = \"gain\"",
    fixed = TRUE
  )
})

test_that("a descriptive index or a number is refused as draws", {
  ri <- pollin_rii(c(30, 32, 28), c(20, 22, 19), n_boot = 20L, seed = 1L)
  expect_error(
    pollin_value(ri, va_market(), channel = "dependence"),
    "only draws of a fitted estimand|Only draws of a fitted estimand"
  )
  expect_error(
    pollin_value(0.3, va_market(), channel = "dependence"),
    "has no estimand and cannot be valued"
  )
})

test_that("channel has no default and the market must be declared", {
  dr <- va_draws()
  expect_error(pollin_value(dr, va_market()), "`channel` must be one of")
  expect_error(
    pollin_value(dr, va_market(), channel = "gain"),
    "`channel` must be one of"
  )
  expect_error(
    pollin_value(dr, list(quantity = 1), channel = "dependence"),
    "`pollin_market`"
  )
})

test_that("an ambiguous level is refused and a named level is valued", {
  hi <- cbind(north = c(30, 32), south = c(28, 31))
  lo <- cbind(north = c(20, 24), south = c(21, 18))
  dr <- va_draws(high = hi, low = lo, by = "Region")
  expect_error(
    pollin_value(dr, va_market(), channel = "dependence"),
    "several columns \\(north, south\\)"
  )
  expect_error(
    pollin_value(dr, va_market(), channel = "dependence", level = "east"),
    "`level` must be one of"
  )
  v <- pollin_value(dr, va_market(), channel = "dependence", level = "south")
  expect_identical(v$draws, (1 - lo[, "south"] / hi[, "south"]) * 5e5)
  expect_identical(v$summary$population, "trial units with Region = south")
})

test_that("the population column is the default level", {
  pop <- pollin_population("Region", c(north = 3, south = 1), label = "two")
  hi <- cbind(north = c(30, 32), south = c(28, 31))
  lo <- cbind(north = c(20, 24), south = c(21, 18))
  dr <- va_draws(high = hi, low = lo, population = pop, by = "Region")
  v <- pollin_value(dr, va_market(), channel = "dependence")
  expect_identical(v$level, "population")
  ph <- 0.75 * hi[, "north"] + 0.25 * hi[, "south"]
  pl <- 0.75 * lo[, "north"] + 0.25 * lo[, "south"]
  expect_equal(v$draws, (1 - pl / ph) * 5e5)
  expect_match(v$summary$population, "two: 2 strata of Region")
})

test_that("a market draw vector of the wrong length is refused", {
  dr <- va_draws()
  expect_error(
    pollin_value(
      dr,
      va_market(price_draws = c(1, 2, 3)),
      channel = "dependence"
    ),
    "`price_draws` has 3 values and the estimand has 4 draws"
  )
  expect_error(
    pollin_value(
      dr,
      va_market(quantity_draws = rep(1000, 5L)),
      channel = "dependence"
    ),
    "`quantity_draws` has 5 values"
  )
})

test_that("a supplemented reference is labelled as potential value", {
  dr <- va_draws(va_estimand(reference = "supplemented"))
  v <- pollin_value(dr, va_market(), channel = "dependence")
  out <- capture.output(print(v))
  expect_true(any(grepl("value at full pollination, not current", out)))
  expect_true(any(grepl("full\\s+pollination", v$assumptions)))
})

test_that("the print names basis, population, year, currency and source", {
  v <- pollin_value(
    va_draws(),
    va_market(),
    channel = "dependence",
    conf_level = 0.9
  )
  out <- capture.output(print(v))
  for (pat in c(
    "basis",
    "population : trial",
    "price year : 2024",
    "currency   : AUD",
    "Synthetic inputs",
    "90% interval",
    "No price response",
    "Total-loss counterfactual",
    "independent of the market draws",
    "units in the trial"
  )) {
    expect_true(any(grepl(pat, out, fixed = TRUE)), info = pat)
  }
})

# Layer L5: attribution and pollin_guard_doublecount() -------------------------

dc_value <- function(
  agent = "all",
  product = NULL,
  crop = "canola",
  currency = "AUD",
  year = 2024,
  high = cbind(overall = c(30, 32, 28, 31)),
  low = cbind(overall = c(20, 24, 21, 18))
) {
  m <- pollin_market(1000, 500, "t", currency, year, "Synthetic", crop = crop)
  pollin_value(
    va_draws(high = high, low = low),
    m,
    channel = "dependence",
    agent = agent,
    product = product
  )
}

test_that("a value carries its attribution", {
  v <- dc_value()
  expect_identical(
    v$attribution,
    list(agent = "all", product = "pods", crop = "canola")
  )
  v2 <- dc_value(agent = "wild insects", product = "seed", crop = NULL)
  expect_identical(v2$attribution$product, "seed")
  expect_null(v2$attribution$crop)
  out <- capture.output(print(v2))
  expect_true(any(grepl("agent wild insects, product seed, crop not", out)))
  expect_error(dc_value(agent = ""), "`agent`")
  expect_error(dc_value(product = 1), "`product`")
  expect_error(
    pollin_market(1, 1, "t", "AUD", 2024, "s", crop = 2),
    "`crop`"
  )
})

test_that("(a) a part plus its total is refused; parts of two crops pass", {
  expect_error(
    pollin_guard_doublecount(dc_value(), dc_value(agent = "honey bees")),
    "The total already contains the part"
  )
  tot <- pollin_guard_doublecount(
    dc_value(),
    dc_value(agent = "honey bees", crop = "almond")
  )
  expect_s3_class(tot, "pollin_value_total")
})

test_that("(b) the same agent twice is refused; two crops pass", {
  expect_error(
    pollin_guard_doublecount(a = dc_value(), b = dc_value()),
    "`a` and `b` are both for 'pods'.*counts one contribution twice"
  )
  expect_s3_class(
    pollin_guard_doublecount(dc_value(), dc_value(crop = "faba bean")),
    "pollin_value_total"
  )
})

test_that("(c) two products of one crop and agent are refused", {
  expect_error(
    pollin_guard_doublecount(
      dc_value(product = "seed"),
      dc_value(product = "oil")
    ),
    "two products \\('seed' and 'oil'\\)"
  )
  expect_s3_class(
    pollin_guard_doublecount(
      dc_value(product = "seed"),
      dc_value(product = "oil", crop = "sunflower")
    ),
    "pollin_value_total"
  )
})

test_that("(d) two guilds need nested = FALSE, which is recorded", {
  m <- dc_value(agent = "managed honey bees")
  w <- dc_value(agent = "wild insects")
  expect_error(pollin_guard_doublecount(m, w), "nested = FALSE")
  expect_error(pollin_guard_doublecount(m, w, nested = TRUE), "two guilds")
  tot <- pollin_guard_doublecount(m, w, nested = FALSE)
  expect_s3_class(tot, "pollin_value_total")
  expect_match(tot$assertion, "separate, non-nested exclusion contrasts")
  out <- capture.output(print(tot))
  expect_true(any(grepl("asserted by the caller", out)))
  expect_true(any(grepl("managed honey bees", out)))
  expect_error(pollin_guard_doublecount(m, w, nested = NA), "`nested`")
})

test_that("a market without a crop matches any crop", {
  expect_error(
    pollin_guard_doublecount(dc_value(crop = NULL), dc_value(crop = "x")),
    "a market without `crop` matches any crop"
  )
})

test_that("agents and products compare ignoring case and spaces", {
  expect_error(
    pollin_guard_doublecount(dc_value(agent = "ALL "), dc_value()),
    "counts one contribution twice"
  )
})

test_that("the total is the per-draw sum of the components", {
  a <- dc_value(crop = "canola")
  b <- dc_value(
    crop = "faba bean",
    high = cbind(overall = c(10, 12, 11, 9)),
    low = cbind(overall = c(8, 9, 10, 7))
  )
  tot <- pollin_guard_doublecount(canola = a, faba = b, conf_level = 0.9)
  expect_identical(tot$draws, a$draws + b$draws)
  expect_identical(names(tot$components), c("canola", "faba"))
  expect_identical(
    c(tot$summary$conf.low, tot$summary$conf.high),
    stats::quantile(tot$draws, c(0.05, 0.95), names = FALSE)
  )
  expect_identical(summary(tot), tot$summary)
  expect_identical(tot$summary$components, 2L)
  expect_null(tot$assertion)
  out <- capture.output(print(tot))
  expect_true(any(grepl("canola:", out, fixed = TRUE)))
  expect_true(any(grepl("90% interval", out, fixed = TRUE)))
})

test_that("the guard refuses mixed currencies, years and draw counts", {
  expect_error(
    pollin_guard_doublecount(dc_value(), dc_value(crop = "x", currency = "NZD")),
    "differ in currency"
  )
  expect_error(
    pollin_guard_doublecount(dc_value(), dc_value(crop = "x", year = 2020)),
    "differ in price year"
  )
  short <- dc_value(
    crop = "x",
    high = cbind(overall = c(30, 32)),
    low = cbind(overall = c(20, 24))
  )
  expect_error(
    pollin_guard_doublecount(dc_value(), short),
    "different numbers of draws"
  )
  expect_error(pollin_guard_doublecount(dc_value()), "two or more")
  expect_error(pollin_guard_doublecount(dc_value(), 3), "two or more")
})

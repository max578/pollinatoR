# Economic value (layer L6) ---------------------------------------------------

#' Value a dependence or deficit estimand
#'
#' Multiplies draws of a ratio estimand by price and current production,
#' draw by draw, and summarises the value draws last: the median and a
#' percentile interval. The ratio is formed from the arm-mean draws of
#' [pollin_draws()] for each draw, so the uncertainty of the estimand (and of
#' the market inputs, when given as draws) carries through to the value.
#'
#' The quantity in [pollin_market()] is current production, the crop as it
#' is pollinated now, so each channel needs the ratio taken relative to it:
#'
#' | Channel | Estimand | Value per draw |
#' |---|---|---|
#' | `"dependence"` | `family = "dependence"`, `scale = "ratio"`, `basis = "loss"` | `D x P x Q`, with `D = 1 - low/high`: the value lost if pollinators were wholly excluded |
#' | `"deficit"` | `family = "deficit"`, `scale = "ratio"`, `basis = "gain"` | `(E[Y_supp]/E[Y_open] - 1) x P x Q`: the value of closing the gap to hand pollination on current production |
#'
#' Any other estimand is refused, with the call that declares the right one.
#' In particular a gain `g = high/low - 1` is refused by the dependence
#' channel: it relates to the loss share by `g = D / (1 - D)`, so valuing a
#' gain overstates the loss by the factor `1 / (1 - D)`.
#'
#' A dependence estimand whose reference arm is `"supplemented"` (potential
#' dependence) is valued and labelled as the value at full pollination, which
#' is not the contribution of the pollination the crop receives now.
#'
#' Draws of quantity or price in the market are paired with the estimand
#' draws by index; their number must equal the number of estimand draws.
#' Only a [pollin_draws()] object can be valued. A descriptive index from
#' [pollin_index()] or [pollin_rii()], or a plain number, carries no estimand.
#'
#' Each value carries an attribution: the pollinators the exclusion removed
#' (`agent`), the product valued (`product`) and the crop (from
#' [pollin_market()]). [pollin_guard_doublecount()] reads it before values are
#' added. For the deficit channel, `agent` names the pollinators whose
#' shortfall hand pollination closes.
#'
#' @param draws A [pollin_draws()] object.
#' @param market A [pollin_market()] object.
#' @param channel `"dependence"` or `"deficit"`. There is no default.
#' @param level The column of the draws to value: a level of the fit's `by`,
#'   `"overall"` or `"population"`. The default is `"population"` when the
#'   draws hold it, else the only column; with several columns and no
#'   population, `level` must be given.
#' @param conf_level Level of the percentile interval.
#' @param agent The pollinators the estimand's exclusion removes, one string:
#'   `"all"` (the default) for a contrast that excludes every pollinator, or
#'   the label of a guild, for example `"managed honey bees"` or
#'   `"wild insects"`.
#' @param product The product valued, one string (for example `"seed"`,
#'   `"oil"` or `"fruit"`). The default is the name of the response.
#'
#' @return An object of class `pollin_value`: a list holding
#'   * `draws`: the numeric vector of value draws, in `currency`;
#'   * `summary`: a one-row data frame with `channel`, `level`, `median`,
#'     `conf.low`, `conf.high`, `currency`, `year`, `unit`, `basis`,
#'     `population`, `method` and `n`;
#'   * `channel`, `estimand`, `level`, `market`;
#'   * `assumptions`: the statements the value rests on;
#'   * `conf_level`: the level of the interval in `summary`;
#'   * `attribution`: a list of `agent`, `product` and `crop` (`NULL` when
#'     the market names none).
#' @seealso [pollin_draws()], [pollin_market()].
#' @export
#' @examples
#' \donttest{
#' set.seed(1)
#' d <- data.frame(
#'   trt = rep(rep(c("Cage", "Open"), each = 10L), 8L),
#'   Site = rep(paste0("s", seq_len(8L)), each = 20L)
#' )
#' d$seed <- rgamma(nrow(d), 20, 20 / exp(1 + 0.3 * (d$trt == "Open")))
#' pd <- pollin_data(d, "seed", "trt", unit = "Site")
#' ds <- pollin_design(pd, c(excluded = "Cage", open = "Open"), "rcbd")
#' est <- pollin_estimand(ds, "dependence", "ratio", "loss")
#' dr <- pollin_draws(pollin_fit(est), n = 2000L, seed = 1L)
#' mk <- pollin_market(12000, 650, "t", "AUD", 2024, "Illustrative inputs")
#' pollin_value(dr, mk, channel = "dependence")
#' }
pollin_value <- function(
  draws,
  market,
  channel,
  level = NULL,
  conf_level = 0.95,
  agent = "all",
  product = NULL
) {
  .value_check_draws(draws)
  .check_class(market, "pollin_market", arg = "market")
  .check_choice(
    if (!missing(channel)) channel,
    "channel",
    c("dependence", "deficit")
  )
  .check_conf_level(conf_level)
  .market_check_string(agent, "agent")
  if (is.null(product)) {
    product <- draws$response
  }
  .market_check_string(product, "product")
  e <- draws$estimand
  .value_check_estimand(e, channel)
  level <- .value_level(level, colnames(draws$estimate))

  # Value per draw -------------------------------------------------------------
  high <- draws$high[, level]
  low <- draws$low[, level]
  ratio <- if (channel == "dependence") 1 - low / high else high / low - 1
  n <- length(ratio)
  value <- ratio *
    .market_input(market, "price", n) *
    .market_input(market, "quantity", n)

  # Summary --------------------------------------------------------------------
  a <- (1 - conf_level) / 2
  q <- stats::quantile(value, c(0.5, a, 1 - a), names = FALSE)
  summary <- data.frame(
    channel = channel,
    level = level,
    median = q[1L],
    conf.low = q[2L],
    conf.high = q[3L],
    currency = market$currency,
    year = market$year,
    unit = market$unit,
    basis = e$basis,
    population = .value_population(draws, level),
    method = draws$method,
    n = n,
    stringsAsFactors = FALSE
  )
  structure(
    list(
      draws = value,
      summary = summary,
      channel = channel,
      estimand = e,
      level = level,
      market = market,
      assumptions = .value_assumptions(draws, channel, level),
      conf_level = conf_level,
      attribution = list(agent = agent, product = product, crop = market$crop)
    ),
    class = "pollin_value"
  )
}

#' @export
print.pollin_value <- function(x, ...) {
  s <- x$summary
  e <- x$estimand
  cat("<pollin_value>\n")
  cat("  channel    :", s$channel, "-", .value_label(x$channel, e), "\n")
  cat(
    "  value      :",
    .value_format(s$median),
    s$currency,
    sprintf(
      "(%s%% interval %s to %s)",
      format(100 * x$conf_level),
      .value_format(s$conf.low),
      .value_format(s$conf.high)
    ),
    "\n"
  )
  cat(
    "  basis      : ",
    s$basis,
    " (",
    .scale_definition(
      e$scale,
      e$basis,
      e$contrast[["high"]],
      e$contrast[["low"]]
    ),
    ")\n",
    sep = ""
  )
  cat("  level      :", s$level, "\n")
  cat("  population :", s$population, "\n")
  cat("  price year :", s$year, "\n")
  cat("  currency   :", s$currency, "\n")
  cat(
    "  market     :",
    format(x$market$quantity, big.mark = ","),
    x$market$unit,
    .market_draws_note(x$market$quantity_draws),
    "at",
    format(x$market$price, big.mark = ","),
    s$currency,
    "per",
    x$market$unit,
    .market_draws_note(x$market$price_draws),
    "\n"
  )
  cat("  source     :", x$market$source, "\n")
  cat("  attribution:", .attribution_label(x$attribution), "\n")
  cat(
    "  method     :",
    .draws_method_label(s$method),
    paste0("(", s$n, " draws)"),
    "\n"
  )
  cat("  assumptions:\n")
  cat(paste0("    - ", x$assumptions, "\n"), sep = "")
  invisible(x)
}

#' Summarise a value
#'
#' @param object A [pollin_value()] object.
#' @param ... Not used.
#'
#' @return The one-row summary data frame of the value: `channel`, `level`,
#'   `median`, `conf.low`, `conf.high`, `currency`, `year`, `unit`, `basis`,
#'   `population`, `method` and `n`.
#' @export
summary.pollin_value <- function(object, ...) {
  object$summary
}

# Internal ---------------------------------------------------------------------

#' Stop unless `draws` is a `pollin_draws` object, explaining why a
#' descriptive index or a number cannot be valued
#' @noRd
#' @keywords internal
.value_check_draws <- function(draws) {
  if (!inherits(draws, "pollin_draws")) {
    stop(
      "`draws` must be a `pollin_draws` object. Only draws of a fitted ",
      "estimand carry the family, scale and basis a value rests on; a ",
      "descriptive index (from `pollin_index()` or `pollin_rii()`) or a plain ",
      "number has no estimand and cannot be valued. Declare the estimand with ",
      "`pollin_estimand()`, fit it and draw it with `pollin_draws()`.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Stop unless the estimand's family, scale and basis match the channel
#' @noRd
#' @keywords internal
.value_check_estimand <- function(estimand, channel) {
  want <- if (channel == "dependence") {
    c(family = "dependence", scale = "ratio", basis = "loss")
  } else {
    c(family = "deficit", scale = "ratio", basis = "gain")
  }
  call <- sprintf(
    "pollin_estimand(design, \"%s\", \"ratio\", \"%s\")",
    want[["family"]],
    want[["basis"]]
  )
  what <- if (channel == "dependence") {
    paste0(
      "the share of current production lost without pollinators, ",
      "D = 1 - low/high"
    )
  } else {
    paste0(
      "the gain of hand pollination over current production, ",
      "E[Y_supp]/E[Y_open] - 1"
    )
  }
  if (estimand$family != want[["family"]] || estimand$scale != "ratio") {
    stop(
      "The `",
      channel,
      "` channel needs ",
      what,
      "; the draws are of a `",
      estimand$family,
      "` estimand on the `",
      estimand$scale,
      "` scale. Declare it with `",
      call,
      "`.",
      call. = FALSE
    )
  }
  if (estimand$basis != want[["basis"]]) {
    if (channel == "dependence") {
      stop(
        "The `dependence` channel needs the share of current production ",
        "lost without pollinators, D = 1 - low/high. A gain ",
        "g = high/low - 1 relates to it as g = D / (1 - D), so valuing a ",
        "gain overstates the loss by the factor 1/(1 - D). Refit the ",
        "estimand with basis = \"loss\": `",
        call,
        "`.",
        call. = FALSE
      )
    }
    stop(
      "The `deficit` channel needs ",
      what,
      ", a gain relative to the open arm; the estimand has basis = \"",
      estimand$basis,
      "\". Refit the estimand with basis = \"gain\": `",
      call,
      "`.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' The column of the draws to value
#' @noRd
#' @keywords internal
.value_level <- function(level, columns) {
  if (is.null(level)) {
    if ("population" %in% columns) {
      return("population")
    }
    if (length(columns) == 1L) {
      return(columns)
    }
    stop(
      "The draws hold several columns (",
      paste(columns, collapse = ", "),
      "); name the one to value with `level =`.",
      call. = FALSE
    )
  }
  if (!is.character(level) || length(level) != 1L || !level %in% columns) {
    stop(
      "`level` must be one of the columns of the draws: ",
      paste(columns, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  level
}

#' The population a column of the draws refers to, in words
#' @noRd
#' @keywords internal
.value_population <- function(draws, level) {
  if (level == "population") {
    return(.population_label(draws$population))
  }
  if (level == "overall") {
    return("trial")
  }
  paste0("trial units with ", draws$by, " = ", level)
}

#' What the value measures, in words
#' @noRd
#' @keywords internal
.value_label <- function(channel, estimand) {
  if (channel == "deficit") {
    return("value of closing the gap to hand pollination")
  }
  if (estimand$design$reference == "supplemented") {
    return("value at full pollination, not current contribution")
  }
  "value lost under total exclusion of pollinators"
}

#' The assumptions a value rests on
#' @noRd
#' @keywords internal
.value_assumptions <- function(draws, channel, level) {
  e <- draws$estimand
  counterfactual <- if (channel == "dependence") {
    paste0(
      "Total-loss counterfactual: the value is the production lost if ",
      "pollinators were wholly excluded, with no substitution by other ",
      "pollinators or inputs."
    )
  } else {
    paste0(
      "Gap counterfactual: the value is the production gained if the whole ",
      "crop received hand-supplemented pollination, at no change in costs."
    )
  }
  if (channel == "dependence" && e$design$reference == "supplemented") {
    counterfactual <- paste(
      counterfactual,
      "The reference arm is supplemented, so this is the value at full",
      "pollination, not the contribution of current pollination."
    )
  }
  population <- if (level %in% c("population", "overall")) {
    e$assumptions$statement[e$assumptions$assumption == "population"]
  } else {
    paste0(
      "The value refers to the trial units with ",
      draws$by,
      " = ",
      level,
      ", and to current production of that stratum."
    )
  }
  c(
    price = paste0(
      "No price response: the price is held fixed whatever the change in ",
      "production."
    ),
    counterfactual = counterfactual,
    independence = paste0(
      "The estimand draws are independent of the market draws."
    ),
    population = population
  )
}

#' A value for printing
#' @noRd
#' @keywords internal
.value_format <- function(x) {
  format(signif(x, 4L), big.mark = ",", scientific = FALSE)
}

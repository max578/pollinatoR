# Market inputs (layer L6) ----------------------------------------------------

#' Declare the market inputs a value rests on
#'
#' The production and price that [pollin_value()] multiplies a ratio estimand
#' by. `quantity` is current production, the output of the crop as it is
#' pollinated now (for example area times yield), in `unit`; `price` is the
#' price per `unit` in `currency`, in the prices of `year`. The value is
#' `ratio x price x quantity`.
#'
#' Each input is a scalar, and may also be given as draws: a vector of
#' positive values expressing its uncertainty. When draws are given they are
#' used in place of the scalar, paired by index with the draws of the
#' estimand, so their number must equal the number of estimand draws. The
#' scalar is then the point value reported beside them. The package
#' propagates the draws given; it does not model or forecast prices.
#'
#' @param quantity Current production, one positive number in `unit`.
#' @param price Price per `unit`, one positive number in `currency`.
#' @param unit Unit of `quantity`, one string (for example `"t"`).
#' @param currency Currency of `price`, one string (for example `"AUD"`).
#' @param year Price year, one whole number.
#' @param source Where the inputs come from, free text printed with every
#'   value.
#' @param quantity_draws,price_draws Optional numeric vectors of positive
#'   draws of `quantity` and `price`.
#' @param crop Optional name of the crop, one string. It enters the
#'   attribution of each value, which [pollin_guard_doublecount()] reads; a
#'   value without a crop is taken to be of the same crop as any other.
#'
#' @return An object of class `pollin_market`: a list holding `quantity`,
#'   `price`, `unit`, `currency`, `year`, `source`, `quantity_draws` and
#'   `price_draws` (`NULL` when not given), and `crop`.
#' @seealso [pollin_value()].
#' @export
#' @examples
#' pollin_market(
#'   quantity = 12000,
#'   price = 650,
#'   unit = "t",
#'   currency = "AUD",
#'   year = 2024,
#'   source = "Illustrative inputs",
#'   crop = "canola"
#' )
#'
#' # Price uncertainty as draws
#' set.seed(1)
#' pollin_market(
#'   12000,
#'   650,
#'   "t",
#'   "AUD",
#'   2024,
#'   "Illustrative inputs",
#'   price_draws = rlnorm(4000, log(650), 0.1)
#' )
pollin_market <- function(
  quantity,
  price,
  unit,
  currency,
  year,
  source,
  quantity_draws = NULL,
  price_draws = NULL,
  crop = NULL
) {
  .market_check_positive(quantity, "quantity")
  .market_check_positive(price, "price")
  .market_check_string(unit, "unit")
  .market_check_string(currency, "currency")
  .market_check_string(source, "source")
  if (
    !is.numeric(year) ||
      length(year) != 1L ||
      !is.finite(year) ||
      year != round(year)
  ) {
    stop("`year` must be one whole number.", call. = FALSE)
  }
  .market_check_draws(quantity_draws, "quantity_draws")
  .market_check_draws(price_draws, "price_draws")
  if (!is.null(crop)) {
    .market_check_string(crop, "crop")
  }
  structure(
    list(
      quantity = quantity,
      price = price,
      unit = unit,
      currency = currency,
      year = as.integer(year),
      source = source,
      quantity_draws = quantity_draws,
      price_draws = price_draws,
      crop = crop
    ),
    class = "pollin_market"
  )
}

#' @export
print.pollin_market <- function(x, ...) {
  cat("<pollin_market>\n")
  if (!is.null(x$crop)) {
    cat("  crop     :", x$crop, "\n")
  }
  cat(
    "  quantity :",
    format(x$quantity, big.mark = ","),
    x$unit,
    .market_draws_note(x$quantity_draws),
    "\n"
  )
  cat(
    "  price    :",
    format(x$price, big.mark = ","),
    x$currency,
    "per",
    x$unit,
    .market_draws_note(x$price_draws),
    "\n"
  )
  cat("  year     :", x$year, "\n")
  cat("  source   :", x$source, "\n")
  invisible(x)
}

# Internal ---------------------------------------------------------------------

#' The draws of a market input paired with `n` estimand draws: the draw
#' vector, or the scalar when no draws were given. Stops on any other length.
#' @noRd
#' @keywords internal
.market_input <- function(market, input, n) {
  d <- market[[paste0(input, "_draws")]]
  if (is.null(d)) {
    return(market[[input]])
  }
  if (length(d) != n) {
    stop(
      "`",
      input,
      "_draws` has ",
      length(d),
      " values and the estimand has ",
      n,
      " draws. Market draws are paired with the estimand draws by index, ",
      "so give one per estimand draw, or none to use the scalar.",
      call. = FALSE
    )
  }
  d
}

#' Text noting that a market input is given as draws, for printing; NULL
#' for a scalar
#' @noRd
#' @keywords internal
.market_draws_note <- function(d) {
  if (is.null(d)) {
    return(NULL)
  }
  paste0("(", length(d), " draws)")
}

#' Stop unless `x` is one positive finite number
#' @noRd
#' @keywords internal
.market_check_positive <- function(x, arg) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0) {
    stop("`", arg, "` must be one positive finite number.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Stop unless `x` is one non-empty string
#' @noRd
#' @keywords internal
.market_check_string <- function(x, arg) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(trimws(x))) {
    stop("`", arg, "` must be one non-empty string.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Stop unless `x` is NULL or a vector of positive finite numbers
#' @noRd
#' @keywords internal
.market_check_draws <- function(x, arg) {
  if (is.null(x)) {
    return(invisible(TRUE))
  }
  if (
    !is.numeric(x) ||
      length(x) < 2L ||
      any(!is.finite(x)) ||
      any(x <= 0)
  ) {
    stop(
      "`",
      arg,
      "` must be NULL or a vector of at least two positive finite draws.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

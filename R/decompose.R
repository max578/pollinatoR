# Whole-yield estimand and its component decomposition ---------------------------

#' Decompose the effect on whole yield into its yield components
#'
#' Seed yield per plant is the product of its components (for canola, pods per
#' plant, seeds per pod and seed weight), so the effect of pollinator access on
#' whole yield is not the effect on any one component. On the log scale the
#' product becomes a sum: if \eqn{R_Y} is the ratio of the G-computation arm
#' means of yield (pollinated arm over excluded arm) and \eqn{R_j} that of
#' component \eqn{j},
#' \deqn{\log R_Y = \sum_j \log R_j + I,}
#' where the interaction term \eqn{I} is what the product of the component
#' ratios does not carry: the covariance of the components within an arm, and
#' any part of yield that is not the product of the listed components. \eqn{I}
#' is near zero when the components vary independently within each arm.
#'
#' Each term is reported with its share of \eqn{\log R_Y} (the shares, the
#' interaction included, sum to one) and a percentile interval from a
#' nonparametric bootstrap that resamples whole units (blocks, fields or sites)
#' with replacement and refits the yield model and every component model on
#' the same resampled units, so the intervals carry the correlation between the
#' components measured on the same plants. A component is flagged as
#' *compensating* when its interval excludes zero on the side opposite to the
#' effect on yield: it moves against yield, and a dependence taken from that
#' component alone would have the wrong sign.
#'
#' Every model is refitted on the rows where yield and all components are
#' observed, so the terms refer to the same plants; the estimates can then
#' differ slightly from those of [pollin_effects()] when a component has rows
#' that yield lacks. Shares are not reported when the interval of
#' \eqn{\log R_Y} contains zero, because a share of an effect that may be zero
#' has no stable reading.
#'
#' @param x A [pollin_fit()] object fitted by glmmTMB, without `by` and without
#'   a binomial denominator, whose responses include `yield` and `components`.
#'   Every response needs a log link.
#' @param yield The name of the whole-yield response.
#' @param components The names of at least two component responses whose
#'   product is yield, up to a constant factor.
#' @param n_boot Number of bootstrap replicates.
#' @param seed Optional seed for the resampling. The session's random number
#'   state is restored afterwards.
#' @param conf_level Level of the percentile intervals.
#'
#' @return An object of class `pollin_decomposition`, a list holding
#'   * `table`: a data frame with one row per component, then `interaction`
#'     and `yield`, and the columns `term`, `log_ratio`, `conf.low`,
#'     `conf.high` (the log ratio, pollinated over excluded, and its interval),
#'     `share`, `share.low`, `share.high` and `compensating`;
#'   * `draws`: the matrix of bootstrap log ratios, one row per successful
#'     replicate;
#'   * `product_check`: the coefficient of variation over plants of yield
#'     divided by the product of the components (zero when yield is an exact
#'     product);
#'   * `effect`, `contrast`, `n_boot`, `n_failed`, `units`, `seed` and
#'     `conf_level`.
#' @seealso [pollin_compensation()] for the correlations among components;
#'   [pollin_effects()] for the effect on each response.
#' @export
#' @examples
#' \donttest{
#' set.seed(2)
#' d <- data.frame(
#'   trt = rep(rep(c("Cage", "Open"), each = 6L), 10L),
#'   Site = rep(paste0("s", seq_len(10L)), each = 12L)
#' )
#' open <- d$trt == "Open"
#' d$pods <- rnbinom(nrow(d), mu = exp(4 + 0.10 * open), size = 20)
#' d$seeds <- rgamma(nrow(d), shape = 20, rate = 20 / exp(2.5 + 0.30 * open))
#' d$weight <- rgamma(nrow(d), shape = 30, rate = 30 / exp(1 - 0.10 * open))
#' d$yield <- d$pods * d$seeds * d$weight / 1000
#' pd <- pollin_data(d, c("yield", "pods", "seeds", "weight"), "trt",
#'   unit = "Site")
#' ds <- pollin_design(pd, c(excluded = "Cage", open = "Open"), "rcbd")
#' est <- pollin_estimand(ds, "dependence", "ratio", "loss")
#' fam <- list(yield = stats::Gamma(link = "log"), pods = glmmTMB::nbinom2(),
#'   seeds = stats::Gamma(link = "log"), weight = stats::Gamma(link = "log"))
#' fit <- pollin_fit(est, family = fam)
#' pollin_decompose(fit, "yield", c("pods", "seeds", "weight"),
#'   n_boot = 50L, seed = 1L)
#' }
pollin_decompose <- function(
  x,
  yield,
  components,
  n_boot = 999L,
  seed = NULL,
  conf_level = 0.95
) {
  .check_class(x, "pollin_fit")
  .decomp_check_fit(x)
  .decomp_check_names(yield, components, names(x$fits))
  .draws_check_n(n_boot)
  .check_seed(seed)
  .check_conf_level(conf_level)

  responses <- c(components, yield)
  for (r in responses) {
    if (x$families[[r]]$link != "log") {
      stop(
        "Response `",
        r,
        "` is fitted with a ",
        x$families[[r]]$link,
        " link; the decomposition needs a log link for every response.",
        call. = FALSE
      )
    }
  }
  d <- x$data[[yield]]
  d <- d[stats::complete.cases(d[responses]), , drop = FALSE]
  unit <- x$estimand$design$data$roles$unit
  groups <- .decomp_groups(x$rhs)

  product_cv <- .decomp_product_cv(d, yield, components)
  if (product_cv > 1e-6) {
    warning(
      "`",
      yield,
      "` is not a constant multiple of the product of the components ",
      "(coefficient of variation of the quotient ",
      signif(product_cv, 3),
      "); the interaction term also carries the part of yield the ",
      "components do not explain.",
      call. = FALSE
    )
  }

  point <- .decomp_log_ratios(x, d, responses)
  if (anyNA(point)) {
    stop("A model failed to fit on the common rows.", call. = FALSE)
  }

  # Unit bootstrap ---------------------------------------------------------------
  ids <- unique(as.character(d[[unit]]))
  idx <- .decomp_resample(length(ids), n_boot, seed)
  rows <- split(seq_len(nrow(d)), factor(as.character(d[[unit]]), ids))
  boot <- t(vapply(
    seq_len(n_boot),
    function(b) {
      d_b <- .decomp_resampled_data(d, rows, idx[b, ], groups)
      .decomp_log_ratios(x, d_b, responses)
    },
    numeric(length(responses))
  ))
  colnames(boot) <- responses
  ok <- stats::complete.cases(boot)
  n_failed <- sum(!ok)
  if (sum(ok) < 2L) {
    stop("Fewer than two bootstrap replicates fitted.", call. = FALSE)
  }
  if (n_failed > 0L) {
    warning(
      .n_noun(n_failed, "bootstrap replicate"),
      " of ",
      n_boot,
      " failed to fit and ",
      if (n_failed == 1L) "was" else "were",
      " dropped.",
      call. = FALSE
    )
  }
  boot <- boot[ok, , drop = FALSE]
  terms <- .decomp_terms(boot, components, yield)

  structure(
    list(
      table = .decomp_table(
        .decomp_terms(t(point), components, yield),
        terms,
        components,
        conf_level
      ),
      draws = terms,
      product_check = product_cv,
      effect = x$estimand$effect,
      contrast = x$estimand$contrast,
      n_boot = as.integer(n_boot),
      n_failed = n_failed,
      units = length(ids),
      seed = seed,
      conf_level = conf_level
    ),
    class = "pollin_decomposition"
  )
}

#' @export
print.pollin_decomposition <- function(x, ...) {
  cat("<pollin_decomposition>\n")
  cat(
    "  log ratio :",
    x$contrast[["high"]],
    "over",
    x$contrast[["low"]],
    "(G-computation arm means)\n"
  )
  cat("  effect    :", x$effect, "\n")
  cat(
    "  bootstrap :",
    x$n_boot - x$n_failed,
    "of",
    x$n_boot,
    "replicates, resampling",
    x$units,
    "units\n"
  )
  tab <- x$table
  num <- vapply(tab, is.numeric, logical(1))
  tab[num] <- lapply(tab[num], signif, digits = 3)
  print(tab, row.names = FALSE)
  comp <- tab$term[tab$compensating %in% TRUE]
  if (length(comp)) {
    cat("  compensating:", paste(comp, collapse = ", "), "\n")
  }
  if (all(is.na(tab$share))) {
    cat("  shares not reported: the yield interval contains zero\n")
  }
  invisible(x)
}

# Internal ---------------------------------------------------------------------

#' Stop unless the fit is a single-population glmmTMB fit without a
#' binomial denominator
#' @noRd
#' @keywords internal
.decomp_check_fit <- function(x) {
  if (x$engine != "glmmTMB") {
    stop("The decomposition refits models with glmmTMB; refit `x` with ",
      "`engine = \"glmmTMB\"`.",
      call. = FALSE
    )
  }
  if (!is.null(x$by)) {
    stop("The decomposition needs a fit without `by`.", call. = FALSE)
  }
  if (!is.null(x$denominator)) {
    stop(
      "The decomposition needs responses on their own scale, not ",
      "successes out of trials.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Stop unless `yield` and two or more distinct `components` are fitted
#' responses
#' @noRd
#' @keywords internal
.decomp_check_names <- function(yield, components, fitted) {
  if (!is.character(yield) || length(yield) != 1L || !yield %in% fitted) {
    stop(
      "`yield` must name one of the fitted responses: ",
      paste(fitted, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  if (
    !is.character(components) ||
      length(components) < 2L ||
      anyDuplicated(components) ||
      yield %in% components ||
      !all(components %in% fitted)
  ) {
    stop(
      "`components` must name at least two other fitted responses: ",
      paste(setdiff(fitted, yield), collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' The grouping variables of the random terms of `rhs`
#' @noRd
#' @keywords internal
.decomp_groups <- function(rhs) {
  bars <- rhs[grepl("|", rhs, fixed = TRUE)]
  unique(unlist(lapply(bars, function(b) {
    all.vars(str2lang(sub("^.*\\|", "", gsub("[()]", "", b))))
  })))
}

#' Coefficient of variation of yield over the product of the components
#' @noRd
#' @keywords internal
.decomp_product_cv <- function(d, yield, components) {
  q <- d[[yield]] / Reduce(`*`, d[components])
  if (any(!is.finite(q))) {
    return(Inf)
  }
  stats::sd(q) / mean(q)
}

#' Log ratios of the G-computation arm means of each response on `d`, NA for
#' a model that fails or warns
#' @noRd
#' @keywords internal
.decomp_log_ratios <- function(x, d, responses) {
  e <- x$estimand
  vapply(
    responses,
    function(r) {
      m <- tryCatch(
        do.call(
          .fit_glmm,
          c(list(r, x$rhs, d, x$families[[r]]), x$args)
        ),
        warning = function(w) NULL,
        error = function(err) NULL
      )
      if (is.null(m)) {
        return(NA_real_)
      }
      a <- .arm_means(
        m,
        d,
        x$treatment,
        e$contrast[["high"]],
        e$contrast[["low"]],
        NULL
      )
      unname(log(a$high / a$low))
    },
    numeric(1)
  )
}

#' Matrix of unit indices, one row per replicate, drawn with replacement
#' from `seed` when given, restoring the session's random number state
#' @noRd
#' @keywords internal
.decomp_resample <- function(n_units, n_boot, seed) {
  if (!is.null(seed)) {
    if (!exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      stats::runif(1L)
    }
    old <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", old, envir = globalenv()))
    set.seed(seed)
  }
  matrix(
    sample.int(n_units, n_units * n_boot, replace = TRUE),
    nrow = n_boot,
    byrow = TRUE
  )
}

#' The rows of the drawn units, each draw a new unit: grouping columns get
#' the draw's position appended so a unit drawn twice is two units
#' @noRd
#' @keywords internal
.decomp_resampled_data <- function(d, rows, draw, groups) {
  parts <- lapply(seq_along(draw), function(k) {
    p <- d[rows[[draw[k]]], , drop = FALSE]
    for (g in groups) {
      p[[g]] <- paste0(as.character(p[[g]]), "#", k)
    }
    p
  })
  out <- do.call(rbind, parts)
  for (g in groups) {
    out[[g]] <- factor(out[[g]])
  }
  out
}

#' Terms of the decomposition from log ratios (one row per replicate):
#' the components, the interaction and yield
#' @noRd
#' @keywords internal
.decomp_terms <- function(lr, components, yield) {
  comp <- lr[, components, drop = FALSE]
  out <- cbind(
    comp,
    interaction = lr[, yield] - rowSums(comp),
    yield = lr[, yield]
  )
  rownames(out) <- NULL
  out
}

#' The decomposition table: point log ratios and shares with percentile
#' intervals, and the compensation flag
#' @noRd
#' @keywords internal
.decomp_table <- function(point, draws, components, conf_level) {
  a <- (1 - conf_level) / 2
  ci <- apply(draws, 2L, stats::quantile, probs = c(a, 1 - a), names = FALSE)
  y_ci <- ci[, "yield"]
  share_ok <- y_ci[1L] > 0 || y_ci[2L] < 0
  shares <- draws / draws[, "yield"]
  sci <- apply(shares, 2L, stats::quantile, probs = c(a, 1 - a), names = FALSE)
  sgn <- sign(point[1L, "yield"])
  is_comp <- colnames(draws) %in% components
  opposes <- (sgn > 0 & ci[2L, ] < 0) | (sgn < 0 & ci[1L, ] > 0)
  data.frame(
    term = colnames(draws),
    log_ratio = point[1L, ],
    conf.low = ci[1L, ],
    conf.high = ci[2L, ],
    share = if (share_ok) point[1L, ] / point[1L, "yield"] else NA_real_,
    share.low = if (share_ok) sci[1L, ] else NA_real_,
    share.high = if (share_ok) sci[2L, ] else NA_real_,
    compensating = ifelse(is_comp, opposes, NA),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

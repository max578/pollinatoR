# Joint uncertainty of the fixed and random effects (glmmTMB engine) ----------

#' Estimates and joint covariance of the conditional-model fixed effects and
#' random effects of one glmmTMB fit. When the full joint precision cannot be
#' inverted, as with a variance at or near zero, the covariance holds the
#' variance parameters at their estimates (the inverse of the `(beta, b)`
#' block of the precision, Henderson's form) and a warning says so. Stops when
#' that block cannot be inverted either, or the model has a zero-inflation
#' part.
#' @noRd
#' @keywords internal
.joint_cov <- function(model, response) {
  zi <- model$modelInfo$allForm$ziformula
  if (!identical(deparse1(zi), "~0")) {
    stop(
      "The joint interval for '",
      response,
      "' does not cover a zero-inflated fit; use `interval = \"bootstrap\"`.",
      call. = FALSE
    )
  }
  env <- model$obj$env
  par <- env$last.par.best
  old <- env$last.par
  old_best <- par
  on.exit({
    env$last.par <- old
    env$last.par.best <- old_best
  })
  sdr <- TMB::sdreport(
    model$obj,
    par.fixed = model$fit$par,
    getJointPrecision = TRUE
  )
  prec <- sdr$jointPrecision
  if (is.null(prec) || nrow(prec) != length(par)) {
    stop(
      "The fit of '",
      response,
      "' gives no joint precision of its parameters; use ",
      "`interval = \"bootstrap\"`.",
      call. = FALSE
    )
  }
  idx <- which(names(par) %in% c("beta", "b"))
  v <- .joint_solve(prec)
  if (is.null(v)) {
    v <- .joint_solve(prec[idx, idx])
    if (is.null(v)) {
      stop(
        "The joint precision of '",
        response,
        "' is singular; use `interval = \"bootstrap\"`.",
        call. = FALSE
      )
    }
    warning(
      "The joint precision of '",
      response,
      "' is singular (a variance may be at zero); its interval holds the ",
      "variance parameters at their estimates.",
      call. = FALSE
    )
  } else {
    v <- v[idx, idx, drop = FALSE]
  }
  names_out <- names(par)[idx]
  beta <- glmmTMB::fixef(model)$cond
  if (sum(names_out == "beta") != length(beta)) {
    stop(
      "The fit of '",
      response,
      "' has fixed effects held at set values; use `interval = \"bootstrap\"`.",
      call. = FALSE
    )
  }
  list(
    par = unname(par[idx]),
    v = v,
    n_beta = length(beta)
  )
}

#' The inverse of a precision matrix, or NULL when it is singular or the
#' inverse has a non-positive variance
#' @noRd
#' @keywords internal
.joint_solve <- function(prec) {
  v <- tryCatch(as.matrix(Matrix::solve(prec)), error = function(err) NULL)
  if (is.null(v) || any(!is.finite(v)) || any(diag(v) <= 0)) {
    return(NULL)
  }
  v
}

#' Fixed and random design matrices of the model data with every row set to
#' the `high` arm and to the `low` arm, from glmmTMB's own model structure,
#' and the offset of each row. The model data and the two arms are stacked in
#' one structure; `x` and `z` have the columns of the fit's fixed and random
#' effects.
#' @noRd
#' @keywords internal
.joint_arm_designs <- function(fit, response) {
  d <- fit$data[[response]]
  trt <- fit$treatment
  e <- fit$estimand
  arms <- c(high = e$contrast[["high"]], low = e$contrast[["low"]])
  nd <- do.call(
    rbind,
    c(
      list(d),
      lapply(arms, function(level) {
        a <- d
        a[[trt]] <- factor(rep(level, nrow(d)), levels = levels(d[[trt]]))
        a
      })
    )
  )
  s <- do.call(
    .fit_glmm,
    c(
      list(
        .response_lhs(response, fit$denominator),
        fit$rhs,
        nd,
        fit$families[[response]]
      ),
      fit$args,
      list(doFit = FALSE)
    )
  )
  n <- nrow(d)
  rows <- list(high = n + seq_len(n), low = 2L * n + seq_len(n))
  lapply(rows, function(i) {
    list(
      x = s$data.tmb$X[i, , drop = FALSE],
      z = s$data.tmb$Z[i, , drop = FALSE],
      offset = s$data.tmb$offset[i]
    )
  })
}

#' Row-averaged arm means and their gradient with respect to `(beta, b)`, per
#' group. `group` holds one label per model row.
#' @noRd
#' @keywords internal
.joint_arm_gradient <- function(model, a, par, n_beta, group) {
  beta <- par[seq_len(n_beta)]
  b <- par[-seq_len(n_beta)]
  eta <- as.vector(a$x %*% beta + a$z %*% b) + a$offset
  fam <- stats::family(model)
  mu <- fam$linkinv(eta)
  w <- fam$mu.eta(eta)
  g <- factor(group)
  avg <- Matrix::sparseMatrix(
    i = as.integer(g),
    j = seq_along(eta),
    x = 1 / as.vector(table(g))[as.integer(g)],
    dims = c(nlevels(g), length(eta))
  )
  grad <- cbind(
    as.matrix(avg %*% (w * a$x)),
    as.matrix(avg %*% Matrix::Diagonal(x = w) %*% a$z)
  )
  rownames(grad) <- levels(g)
  mean <- as.vector(avg %*% mu)
  names(mean) <- levels(g)
  list(mean = mean, grad = grad)
}

#' Joint draws of the row-averaged arm means, one row per draw and one column
#' per group
#' @noRd
#' @keywords internal
.joint_arm_draws <- function(fit, response, group, n, seed) {
  e <- fit$estimand
  model <- fit$fits[[response]]
  jc <- .joint_cov(model, response)
  r <- tryCatch(chol(jc$v), error = function(err) NULL)
  if (is.null(r)) {
    stop(
      "The joint covariance of '",
      response,
      "' is not positive definite; use `method = \"bootstrap\"`.",
      call. = FALSE
    )
  }
  k <- length(jc$par)
  z <- matrix(.draws_rnorm(n * k, seed), nrow = n)
  draws <- z %*% r + rep(jc$par, each = n)
  g <- factor(group)
  designs <- .joint_arm_designs(fit, response)
  arm <- function(a) {
    eta <- as.matrix(
      a$x %*% t(draws[, seq_len(jc$n_beta), drop = FALSE]) +
        a$z %*% t(draws[, -seq_len(jc$n_beta), drop = FALSE])
    ) +
      a$offset
    mu <- stats::family(model)$linkinv(eta)
    out <- t(rowsum(mu, g) / as.vector(table(g)))
    colnames(out) <- levels(g)
    out
  }
  list(high = arm(designs$high), low = arm(designs$low))
}

#' Delta-method effects for one response from the joint covariance, overall
#' or per level of `by`, in the layout of `.effects_delta()`
#' @noRd
#' @keywords internal
.effects_joint <- function(fit, response, estimand, by, conf_level) {
  model <- fit$fits[[response]]
  d <- fit$data[[response]]
  jc <- .joint_cov(model, response)
  group <- if (is.null(by)) rep("all", nrow(d)) else as.character(d[[by]])
  designs <- .joint_arm_designs(fit, response)
  h <- .joint_arm_gradient(model, designs$high, jc$par, jc$n_beta, group)
  l <- .joint_arm_gradient(model, designs$low, jc$par, jc$n_beta, group)
  ratio <- estimand$scale != "difference"
  if (ratio) {
    est <- log(h$mean) - log(l$mean)
    grad <- h$grad / h$mean - l$grad / l$mean
  } else {
    est <- h$mean - l$mean
    grad <- h$grad - l$grad
  }
  se <- sqrt(rowSums((grad %*% jc$v) * grad))
  q <- stats::qnorm(1 - (1 - conf_level) / 2)
  out <- data.frame(
    estimate = unname(est),
    std.error = unname(se),
    conf.low = unname(est - q * se),
    conf.high = unname(est + q * se),
    p.value = unname(2 * stats::pnorm(-abs(est / se)))
  )
  if (ratio) {
    s <- estimand$scale
    b <- estimand$basis
    r <- exp(out$estimate)
    out$std.error <- out$std.error * .from_ratio_slope(r, s, b)
    out$estimate <- .from_ratio(r, s, b)
    out$conf.low <- .from_ratio(exp(est - q * se), s, b)
    out$conf.high <- .from_ratio(exp(est + q * se), s, b)
  }
  if (!is.null(by)) {
    lv <- data.frame(names(h$mean), stringsAsFactors = FALSE)
    names(lv) <- by
    out <- cbind(lv, out)
  }
  out
}

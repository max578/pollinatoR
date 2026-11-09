# Bayesian fits through flexyBayes (layer L3) ---------------------------------

#' glmmTMB family names and the flexyBayes family each maps to on INLA, with
#' the link both must use
#' @noRd
#' @keywords internal
.flexy_families <- data.frame(
  glmmTMB = c("nbinom2", "poisson", "binomial", "Gamma", "beta"),
  flexyBayes = c("negative_binomial", "poisson", "binomial", "gamma", "beta"),
  link = c("log", "log", "logit", "log", "logit"),
  stringsAsFactors = FALSE
)

# Availability ----------------------------------------------------------------

#' Stop unless flexyBayes and INLA are installed and flexyBayes offers
#' extension API version 1
#' @noRd
#' @keywords internal
.flexy_check_available <- function() {
  for (pkg in c("flexyBayes", "INLA")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(
        "`engine = \"flexyBayes\"` needs the ",
        pkg,
        " package, which is not installed.",
        if (pkg == "INLA") {
          " Install it from https://inla.r-inla-download.org/R/stable."
        },
        call. = FALSE
      )
    }
  }
  api <- flexyBayes::fb_extension_api()
  if (!identical(api, 1L)) {
    stop(
      "pollinatoR reads flexyBayes fits through extension API version 1; ",
      "the installed flexyBayes offers version ",
      format(api),
      ".",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

# Arguments ---------------------------------------------------------------------

#' The flexyBayes family and link for a glmmTMB family object; any family or
#' link without an INLA counterpart is refused by name
#' @noRd
#' @keywords internal
.flexy_family <- function(family, response) {
  i <- match(family$family, .flexy_families$glmmTMB)
  if (is.na(i)) {
    stop(
      "Family '",
      family$family,
      "' of '",
      response,
      "' has no flexyBayes counterpart on INLA. The flexyBayes engine fits ",
      paste(.flexy_families$glmmTMB, collapse = ", "),
      "; choose one with `family =`, or use `engine = \"glmmTMB\"`.",
      call. = FALSE
    )
  }
  if (!identical(family$link, .flexy_families$link[i])) {
    stop(
      "Family '",
      family$family,
      "' of '",
      response,
      "' has link '",
      family$link,
      "'; the flexyBayes engine fits it with link '",
      .flexy_families$link[i],
      "' only.",
      call. = FALSE
    )
  }
  list(family = .flexy_families$flexyBayes[i], link = .flexy_families$link[i])
}

#' Stop when `...` sets what the engine route sets itself, or a backend other
#' than INLA
#' @noRd
#' @keywords internal
.flexy_check_dots <- function(dots) {
  own <- intersect(
    names(dots),
    c("fixed", "random", "data", "family", "link", "syntax")
  )
  if (length(own)) {
    stop(
      "`",
      own[1L],
      "` is set by `pollin_fit()` from the estimand and cannot be passed ",
      "through `...`.",
      call. = FALSE
    )
  }
  backend <- dots$backend
  if (!is.null(backend) && !identical(backend, "inla")) {
    if (identical(backend, "brms")) {
      stop(
        "The flexyBayes engine fits on INLA only: posterior draws of a brms ",
        "fit are not in the flexyBayes extension API version 1.",
        call. = FALSE
      )
    }
    stop(
      "The flexyBayes engine takes `backend = \"inla\"` only.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

# Random structure ------------------------------------------------------------

#' The random terms of the flexyBayes route, each a column of the data: one
#' grouping column per unit level, outer to inner (a combined factor for a
#' nested level), and for the slope level a site-by-arm cell column; a
#' `pooled` column adds its own intercept and cell columns. Returns
#' the data with the columns added, the term names, and the cell terms with
#' their unit key.
#' @noRd
#' @keywords internal
.flexy_random <- function(d, unit, treatment, slope, pooled = NULL) {
  terms <- character(0)
  cells <- list()
  add_column <- function(d, key) {
    if (length(key) == 1L) {
      return(list(data = d, name = key))
    }
    name <- paste0("pollin_", paste(key, collapse = "_"))
    if (name %in% names(d)) {
      stop(
        "The data already hold a column '",
        name,
        "', which the flexyBayes engine uses for a random term.",
        call. = FALSE
      )
    }
    d[[name]] <- interaction(d[key], drop = TRUE, sep = ":")
    list(data = d, name = name)
  }
  for (i in seq_along(unit)) {
    key <- unit[seq_len(i)]
    g <- add_column(d, key)
    d <- g$data
    terms <- c(terms, g$name)
    if (identical(unit[i], slope)) {
      g <- add_column(d, c(key, treatment))
      d <- g$data
      terms <- c(terms, g$name)
      cells[[length(cells) + 1L]] <- list(term = g$name, key = key)
    }
  }
  if (!is.null(pooled)) {
    g <- add_column(d, c(pooled, treatment))
    d <- g$data
    terms <- c(terms, pooled, g$name)
    cells[[length(cells) + 1L]] <- list(term = g$name, key = pooled)
  }
  list(data = d, terms = terms, cells = cells)
}

# Fitting ---------------------------------------------------------------------

#' Fit one response with flexyBayes on INLA
#' @noRd
#' @keywords internal
.flexy_fit_one <- function(response, fixed, random, d, family, dots) {
  fam <- .flexy_family(family, response)
  args <- utils::modifyList(list(verbose = FALSE), dots)
  args <- c(
    list(
      fixed = stats::reformulate(fixed, response = response),
      random = stats::reformulate(random),
      data = d,
      family = fam$family,
      link = fam$link,
      backend = "inla",
      aggregate = FALSE
    ),
    args[setdiff(names(args), c("backend", "aggregate"))]
  )
  do.call(flexyBayes::flexybayes, args)
}

#' Stop when the convergence verdict of a flexyBayes fit is FALSE, with its
#' reasons; warn when it is NA
#' @noRd
#' @keywords internal
.flexy_check_convergence <- function(fit, response) {
  conv <- flexyBayes::fb_convergence(fit)
  reasons <- paste(conv$reasons, collapse = "; ")
  if (isFALSE(conv$ok)) {
    stop(
      "The flexyBayes fit of '",
      response,
      "' did not converge: ",
      reasons,
      call. = FALSE
    )
  }
  if (is.na(conv$ok)) {
    warning(
      "flexyBayes gives no convergence verdict for '",
      response,
      "'",
      if (nzchar(reasons)) paste0(": ", reasons),
      ".",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Build a `pollin_fit` on the flexyBayes engine: the fixed part of the
#' glmmTMB route, and its random part with the treatment slope entered as a
#' unit-by-arm cell term
#' @noRd
#' @keywords internal
.flexy_new_fit <- function(estimand, inputs, family, ..., by = NULL) {
  .flexy_check_available()
  dots <- list(...)
  .flexy_check_dots(dots)
  .check_fit_by(by, inputs$data, inputs$treatment)
  fixed <- if (is.null(by)) {
    inputs$treatment
  } else {
    paste(inputs$treatment, "*", by)
  }
  fits <- list()
  data <- list()
  families <- list()
  random <- NULL
  for (r in inputs$responses) {
    d <- inputs$data
    keep <- stats::complete.cases(d[c(r, inputs$treatment, inputs$unit, by)])
    if (!all(keep)) {
      d <- d[keep, , drop = FALSE]
    }
    if (!is.null(by)) {
      d[[by]] <- factor(d[[by]])
    }
    fam <- .family_for(family, r, d[[r]])
    .flexy_family(fam, r)
    random <- .flexy_random(
      d,
      inputs$unit,
      inputs$treatment,
      estimand$design$slope,
      inputs$pooled
    )
    fit <- tryCatch(
      .flexy_fit_one(r, fixed, random$terms, random$data, fam, dots),
      error = function(e) e
    )
    if (inherits(fit, "error")) {
      warning(
        "Outcome '",
        r,
        "' failed to fit: ",
        conditionMessage(fit),
        call. = FALSE
      )
      next
    }
    .flexy_check_convergence(fit, r)
    fits[[r]] <- fit
    data[[r]] <- random$data
    families[[r]] <- fam
  }
  if (length(fits) == 0L) {
    stop("No outcome could be fitted.", call. = FALSE)
  }
  structure(
    list(
      estimand = estimand,
      engine = "flexyBayes",
      fits = fits,
      data = data,
      rhs = c(fixed, random$terms),
      families = families,
      treatment = inputs$treatment,
      by = by,
      slope_key = inputs$slope_key,
      pooled = inputs$pooled,
      args = dots,
      flexy = list(fixed = fixed, random = random$terms, cells = random$cells)
    ),
    class = "pollin_fit"
  )
}

# Posterior draws ---------------------------------------------------------------

#' Joint posterior samples of a flexyBayes fit as a matrix, one column per
#' draw and one row per latent variable. With `seed`, R's generator (which
#' picks each draw's hyperparameter configuration) and INLA's are both set,
#' and the session's random number state is restored afterwards.
#' @noRd
#' @keywords internal
.flexy_latent <- function(model, n, seed) {
  args <- list(model, n, num.threads = "1:1")
  if (!is.null(seed)) {
    if (!exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      stats::runif(1L)
    }
    old <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", old, envir = globalenv()))
    set.seed(seed)
    args$seed <- as.integer(seed)
  }
  s <- do.call(flexyBayes::fb_inla_posterior_sample, args)
  latent <- do.call(cbind, lapply(s, `[[`, "latent"))
  rownames(latent) <- rownames(s[[1L]]$latent)
  latent
}

#' Row index into the latent matrix of each data row's level of a random
#' term, or NA for a level the fit has no effect for
#' @noRd
#' @keywords internal
.flexy_term_rows <- function(term, values, levels, latent_names) {
  match(paste0(term, ":", match(values, levels)), latent_names)
}

#' Contribution of the random terms to the linear predictor of every row and
#' draw. Cell terms take the unit-by-arm cell of `arm`, or of each row's own
#' arm when `arm` is NULL; a cell the fit has no effect for adds zero.
#' @noRd
#' @keywords internal
.flexy_random_part <- function(fit, d, latent, levels, treatment, arm) {
  cell_terms <- vapply(fit$flexy$cells, `[[`, "", "term")
  out <- matrix(0, nrow(d), ncol(latent))
  for (v in fit$flexy$random) {
    i <- match(v, cell_terms)
    values <- if (is.na(i) || is.null(arm)) {
      as.character(d[[v]])
    } else {
      key <- fit$flexy$cells[[i]]$key
      do.call(paste, c(lapply(d[key], as.character), list(arm, sep = ":")))
    }
    rows <- .flexy_term_rows(v, values, levels[[v]], rownames(latent))
    if (anyNA(rows) && (is.na(i) || is.null(arm))) {
      stop(
        "The posterior draws lack levels of the random term '",
        v,
        "' that the model data hold.",
        call. = FALSE
      )
    }
    u <- latent[rows, , drop = FALSE]
    u[is.na(rows), ] <- 0
    out <- out + u
  }
  out
}

#' Posterior draws of the G-computation arm means of one response of a
#' flexyBayes fit, overall or per level of `by`: every row of the model data
#' is predicted under each arm with its own random effects, as
#' `predict(re.form = NULL)` does on the glmmTMB engine, on the response
#' scale, and averaged. The linear predictor rebuilt for the rows as fitted
#' must equal the fit's own, draw by draw.
#' @noRd
#' @keywords internal
.flexy_arm_draws <- function(fit, response, n, seed, by) {
  e <- fit$estimand
  model <- fit$fits[[response]]
  d <- fit$data[[response]]
  trt <- fit$treatment
  latent <- .flexy_latent(model, n, seed)
  terms <- flexyBayes::fb_model_terms(model)$random_terms
  levels <- stats::setNames(
    lapply(terms, `[[`, "var_levels"),
    vapply(terms, `[[`, "", "var")
  )
  fixed <- stats::reformulate(fit$flexy$fixed)
  design <- function(data) {
    x <- flexyBayes::fb_fixed_design(fixed, data)$X
    rows <- match(paste0(colnames(x), ":1"), rownames(latent))
    if (anyNA(rows)) {
      stop(
        "The posterior draws lack the fixed effects ",
        paste(colnames(x)[is.na(rows)], collapse = ", "),
        ".",
        call. = FALSE
      )
    }
    x %*% latent[rows, , drop = FALSE]
  }

  # The fitted rows: rebuilt linear predictor against the fit's own ----------
  eta <- design(d) + .flexy_random_part(fit, d, latent, levels, trt, NULL)
  pred <- latent[paste0("Predictor:", seq_len(nrow(d))), , drop = FALSE]
  if (max(abs(eta - pred)) > 1e-8 * max(1, abs(pred))) {
    stop(
      "The posterior draws of '",
      response,
      "' do not rebuild the fit's linear predictor.",
      call. = FALSE
    )
  }

  # Each arm, every row ----------------------------------------------------
  linkinv <- fit$families[[response]]$linkinv
  g <- if (is.null(by)) rep("all", nrow(d)) else as.character(d[[by]])
  arm <- function(level) {
    nd <- d
    nd[[trt]] <- factor(rep(level, nrow(d)), levels = levels(d[[trt]]))
    eta <- design(nd) + .flexy_random_part(fit, nd, latent, levels, trt, level)
    t(rowsum(linkinv(eta), g) / as.vector(table(g)))
  }
  list(high = arm(e$contrast[["high"]]), low = arm(e$contrast[["low"]]))
}

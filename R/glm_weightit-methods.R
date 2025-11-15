#' Methods for `glm_weightit()` objects
#' @name glm_weightit-methods
#'
#' @description
#' This page documents methods for objects returned by
#' [glm_weightit()], `lm_weightit()`, `ordinal_weightit()`,
#' `multinom_weightit()`, and `coxph_weightit()`. `predict()` methods are
#' described at [predict.glm_weightit()] and `anova()` methods are described at
#' [anova.glm_weightit()].
#'
#' @inheritParams stats::vcov
#' @inheritParams stats::confint
#' @inheritParams stats::print.lm
#' @param object,x an output from one of the above modeling functions.
#' @param ci `logical`; whether to display Wald confidence intervals for
#'   estimated coefficients. Default is `FALSE`. (Note: this argument can also
#'   be supplied as `conf.int`.)
#' @param level when `ci = TRUE`, the desired confidence level.
#' @param transform the function used to transform the coefficients, e.g., `exp`
#'   (which can also be supplied as a string, e.g., `"exp"`); passed to
#'   [match.fun()] before being used on the coefficients. When `ci = TRUE`, this
#'   is also applied to the confidence interval bounds. If specified, the
#'   standard error will be omitted from the output. Default is no
#'   transformation.
#' @param thresholds `logical`; whether to include thresholds in the `summary()`
#'   output for `ordinal_weightit` objects. Default is `TRUE`.
#' @param vcov either a string indicating the method used to compute the
#'   variance of the estimated parameters for `object`, a function used to
#'   extract the variance, or the variance matrix itself. Default is to use the
#'   variance matrix already present in `object`. If a string or function,
#'   arguments passed to `...` are supplied to the method or function. (Note:
#'   for `vcov()`, can also be supplied as `type`.)
#' @param complete `logical`; whether the full variance-covariance matrix should
#'   be returned also in case of an over-determined system where some
#'   coefficients are undefined and `coef(.)` contains `NA`s correspondingly.
#'   When `complete = TRUE`, `vcov()` is compatible with `coef()` also in this
#'   singular case.
#' @param formula. changes to the model formula, passed to the `new` argument of
#'   [update.formula()].
#' @param asympt `logical`; for `estfun()`, whether to use the asymptotic
#'   empirical estimating functions that account for estimation of the weights
#'   (when `Mparts` is available). Default is `TRUE`. Set to `FALSE` to ignore
#'   estimation of the weights. Ignored when `Mparts` is not available or no
#'   argument was supplied to `weightit` in the fitting function.
#' @param \dots for `vcov()` or `summary()` or `confint()` with `vcov` supplied,
#'   other arguments used to compute the variance matrix depending on the method
#'   supplied to `vcov`, e.g., `cluster`, `R`, or `fwb.args`. For `update()`,
#'   additional arguments to the call or arguments with changed values. See
#'   [glm_weightit()] for details.
#' @param evaluate `logical`; whether to evaluate the call (`TRUE`, the default) or just
#'   return it.
#'
#' @returns
#' `summary()` returns a `summary.glm_weightit()` object, which has its
#' own `print()` method. For `coxph_weightit()` objects, the `print()` and
#' `summary()` methods are more like those for `glm` objects than for `coxph`
#' objects.
#'
#' Otherwise, all methods return the same type of object as their generics.
#'
#' @details
#' `vcov()` by default extracts the parameter covariance matrix already
#' computed by the fitting function, and `summary()` and `confint()` uses this
#' covariance matrix to compute standard errors and Wald confidence intervals
#' (internally calling [confint.lm()]), respectively. Supplying arguments to
#' `vcov` or `...` will compute a new covariance matrix. If `cluster` was
#' supplied to the original fitting function, it will be incorporated into any
#' newly computed covariance matrix unless `cluster = NULL` is specified in
#' `vcov()`, `summary()`, or `confint()`. For other arguments (e.g., `R` and
#' `fwb.args`), the defaults are those used by [glm_weightit()]. Note that for
#' `vcov = "BS"` and `vcov = "FWB"` (and `vcov = "const"` for
#' `multinom_weightit` or `ordinal_weightit` objects), the environment for the
#' fitting function is used, so any changes to that environment may affect
#' calculation. It is always safer to simply recompute the fitted object with a
#' new covariance matrix than to modify it with the `vcov` argument, but it can
#' be quicker to just request a new covariance matrix when refitting the model
#' is slow.
#'
#' `update()` updates a fitted model object with new arguments, e.g., a new
#' model formula, dataset, or variance matrix. When only arguments that control
#' the computation of the variance are supplied, only the variance will be
#' recalculated (i.e., the parameters will not be re-estimated). When `data` is
#' supplied, `weightit` is not supplied, and a `weightit` object was originally
#' passed to the model fitting function, the `weightit` object will be re-fit
#' with the new dataset before the model is refit using the new weights and new
#' data. That is, calling `update(obj, data = d)` is equivalent to calling
#' `update(obj, data = d, weightit = update(obj$weightit, data = d))` when a
#' `weightit` object was supplied to the model fitting function. Similarly,
#' supplying `s.weights` or `weights` passes the argument through to
#' `weightit()` to be refit. When `s.weights` or `weights` are supplied and no
#' `weightit` object is present, a fake one containing just the supplied weights
#' will be created.
#'
#' `estfun()` extracts the empirical estimating functions for the fitted model, optionally accounting for the estimation of the weights (if available). This, along with `bread()`, is used by [sandwich::sandwich()] to compute the robust covariance matrix of the estimated coefficients. See [glm_weightit()] and `vcov()` above for more details.
#'
#' @seealso
#' [glm_weightit()] for the page documenting `glm_weightit()`,
#' `lm_weightit()`, `ordinal_weightit()`, `multinom_weightit()`, and
#' `coxph_weightit()`. [summary.glm()], [vcov()], [confint()] for the relevant
#' methods pages. [predict.glm_weightit()] for computing predictions from the
#' models. [anova.glm_weightit()] for comparing models using a Wald test.
#'
#' [sandwich::estfun()] and [sandwich::bread()] for the `estfun()` and `bread()` generics.
#'
#' @examples
#' ## See examples at ?glm_weightit

#' @exportS3Method summary glm_weightit
#' @rdname glm_weightit-methods
summary.glm_weightit <- function(object,
                                 ci = FALSE,
                                 level = .95,
                                 transform = NULL,
                                 vcov = NULL,
                                 ...) {

  if ("conf.int" %in% ...names()) {
    conf.int <- ...elt(match("conf.int", ...names()))
    chk::chk_flag(conf.int)
    ci <- conf.int
  }
  else {
    chk::chk_flag(ci)
  }

  df.r <- object$df.residual

  fam <- object$family
  dispersion <- {
    if (is.null(fam)) NaN
    else if (!is.null(fam$dispersion) && !is.na(fam$dispersion)) fam$dispersion
    else if (fam$family %in% c("poisson", "binomial", "multinomial")) 1
    else if (df.r > 0) {
      if (any(object$weights == 0)) {
        .wrn("observations with zero weight not used for calculating dispersion")
      }
      sum((object$weights * object$residuals^2)[object$weights > 0]) / df.r
    }
    else {
      NaN
    }
  }

  if (is_null(transform)) transform <- base::identity
  else transform <- match.fun(transform)

  coef.p <- coef(object)
  aliased <- is.na(coef.p)
  p <- sum(!aliased)

  if (p > 0L) {
    p1 <- seq_len(p)
    Qr <- object$qr
    coef.p <- object$coefficients[!aliased]

    covmat.unscaled <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
    dimnames(covmat.unscaled) <- list(names(coef.p), names(coef.p))
    df.f <- NCOL(Qr$qr)

    covmat <- .vcov_glm_weightit.internal(object, vcov. = vcov, ...)

    object <- .set_vcov(object, covmat)

    if (is_null(covmat)) {
      ci <- FALSE
      s.err <- rep_with(NA_real_, coef.p)
      pvalue <- tvalue <- s.err
    }
    else {
      s.err <- sqrt(diag(covmat)[!aliased])
      tvalue <- coef.p / s.err
      pvalue <- 2 * pnorm(-abs(tvalue))
    }

    dn <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
    dimnames(coef.table) <- list(names(coef.p), dn)
  }
  else {
    coef.table <- matrix(NA_real_, nrow = 0L, ncol = 4L)
    dimnames(coef.table) <- list(NULL, c("Estimate", "Std. Error",
                                         "z value", "Pr(>|z|)"))
    covmat <- sq_matrix(n = 0L)
  }

  transformed_coefs <- transform(coef.table[, "Estimate"])
  if (!is.numeric(transformed_coefs) || length(transformed_coefs) != length(coef.table[, "Estimate"])) {
    .err("`transform` must return a numeric vector")
  }

  identity_transform <- all(transformed_coefs == coef.table[, "Estimate"])
  if (!identity_transform) {
    coef.table[, "Estimate"] <- transformed_coefs
    coef.table <- coef.table[, -2L, drop = FALSE]
  }

  if (ci) {
    conf <- confint(object, parm = which(!aliased), level = level)
    if (!identity_transform) conf[] <- transform(conf)
    coef.table <- cbind(coef.table, conf)
  }

  keep <- match(c("call", "terms", "family", "deviance", "aic", "contrasts", "df.residual",
                  "null.deviance", "df.null", "iter", "na.action", "vcov_type"),
                names(object), 0L)

  ans <- c(object[keep],
           list(deviance.resid = residuals(object, type = "deviance"),
                coefficients = coef.table,
                aliased = aliased,
                dispersion = dispersion,
                df = c(object$rank, df.r, df.f),
                cov.unscaled = covmat.unscaled,
                cov.scaled = covmat,
                cluster = .attr(object, "cluster"),
                transformed = !identity_transform))

  class(ans) <- c("summary.glm_weightit", "summary.glm")

  ans
}

#' @exportS3Method summary multinom_weightit
#' @rdname glm_weightit-methods
summary.multinom_weightit <- function(object,
                                      ci = FALSE,
                                      level = .95,
                                      transform = NULL,
                                      vcov = NULL,
                                      ...) {
  if ("conf.int" %in% ...names()) {
    conf.int <- ...elt(match("conf.int", ...names()))
    chk::chk_flag(conf.int)
    ci <- conf.int
  }
  else {
    chk::chk_flag(ci)
  }

  df.r <- object$df.residual

  dispersion <-  1

  if (is_null(transform)) transform <- base::identity
  else transform <- match.fun(transform)

  coef.p <- coef(object)
  aliased <- is.na(coef.p)
  p <- sum(!aliased)
  if (p > 0L) {
    Qr <- object$qr
    coef.p <- object$coefficients[!aliased]

    covmat.unscaled <- NULL
    df.f <- NCOL(Qr$qr)

    covmat <- .vcov_glm_weightit.internal(object, vcov. = vcov, ...)

    object <- .set_vcov(object, covmat)

    if (is_null(covmat)) {
      ci <- FALSE
      s.err <- rep_with(NA_real_, coef.p)
      pvalue <- tvalue <- s.err
    }
    else {
      s.err <- sqrt(diag(covmat)[!aliased])
      tvalue <- coef.p / s.err
      pvalue <- 2 * pnorm(-abs(tvalue))
    }

    dn <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
    dimnames(coef.table) <- list(names(coef.p), dn)
  }
  else {
    coef.table <- matrix(NA_real_, nrow = 0L, ncol = 4L)
    dimnames(coef.table) <- list(NULL, c("Estimate", "Std. Error",
                                         "z value", "Pr(>|z|)"))
    covmat <- sq_matrix(n = 0L)
  }

  transformed_coefs <- transform(coef.table[, "Estimate"])
  if (!is.numeric(transformed_coefs) || length(transformed_coefs) != length(coef.table[, "Estimate"])) {
    .err("`transform` must return a numeric vector")
  }

  identity_transform <- all(transformed_coefs == coef.table[, "Estimate"])
  if (!identity_transform) {
    coef.table[, "Estimate"] <- transformed_coefs
    coef.table <- coef.table[, -2L, drop = FALSE]
  }

  if (ci) {
    conf <- confint(object, parm = which(!aliased), level = level)
    if (!identity_transform) conf[] <- transform(conf)
    coef.table <- cbind(coef.table, conf)
  }

  keep <- match(c("call", "terms", "family", "deviance", "aic", "contrasts", "df.residual",
                  "null.deviance", "df.null", "iter", "na.action", "vcov_type"),
                names(object), 0L)

  ans <- c(object[keep],
           list(deviance.resid = residuals(object, type = "deviance"),
                coefficients = coef.table,
                aliased = aliased,
                dispersion = dispersion,
                df = c(object$rank, df.r, df.f),
                cov.unscaled = covmat.unscaled,
                cov.scaled = covmat,
                cluster = .attr(object, "cluster"),
                transformed = !identity_transform))

  class(ans) <- c("summary.glm_weightit", "summary.glm")

  ans
}

#' @exportS3Method summary ordinal_weightit
#' @rdname glm_weightit-methods
summary.ordinal_weightit <- function(object,
                                     ci = FALSE,
                                     level = .95,
                                     transform = NULL,
                                     thresholds = TRUE,
                                     vcov = NULL,
                                     ...) {

  chk::chk_flag(thresholds)

  out <- summary.multinom_weightit(object, ci = ci, level = level,
                                   transform = transform, vcov = vcov, ...)

  nthreshold <- ncol(object$fitted.values) - 1L

  thresholds_ind <- seq(nrow(out$coefficients) - nthreshold + 1L, nrow(out$coefficients))

  if (thresholds) {
    attr(out, "thresholds") <- rownames(out$coefficients)[thresholds_ind]
  }
  else {
    out$coefficients <- out$coefficients[-thresholds_ind, , drop = FALSE]
  }

  out
}

#' @exportS3Method summary coxph_weightit
#' @rdname glm_weightit-methods
summary.coxph_weightit <- function(object,
                                   ci = FALSE,
                                   level = .95,
                                   transform = NULL,
                                   vcov = NULL,
                                   ...) {
  if ("conf.int" %in% ...names()) {
    conf.int <- ...elt(match("conf.int", ...names()))
    chk::chk_flag(conf.int)
    ci <- conf.int
  }
  else {
    chk::chk_flag(ci)
  }

  if (is.null(transform)) transform <- base::identity
  else transform <- match.fun(transform)

  coef.p <- coef(object)
  aliased <- is.na(coef.p)
  p <- sum(!aliased)
  if (p > 0L) {
    coef.p <- object$coefficients[!aliased]

    covmat.unscaled <- NULL
    df.f <- NULL

    covmat <- .vcov_glm_weightit.internal(object, vcov. = vcov, ...)

    object <- .set_vcov(object, covmat)

    if (is_null(covmat)) {
      ci <- FALSE
      s.err <- rep_with(NA_real_, coef.p)
      pvalue <- tvalue <- s.err
    }
    else {
      s.err <- sqrt(diag(covmat))
      tvalue <- coef.p / s.err
      pvalue <- 2 * pnorm(-abs(tvalue))
    }
    dn <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
    dimnames(coef.table) <- list(names(coef.p), dn)
  }
  else {
    coef.table <- matrix(NA_real_, nrow = 0L, ncol = 4L)
    colnames(coef.table) <- c("Estimate", "Std. Error",
                              "z value", "Pr(>|z|)")

    covmat <- sq_matrix(n = 0L)
  }

  transformed_coefs <- transform(coef.table[, "Estimate"])
  if (!is.numeric(transformed_coefs) || length(transformed_coefs) != length(coef.table[, "Estimate"])) {
    .err("`transform` must return a numeric vector")
  }

  identity_transform <- all(transformed_coefs == coef.table[, "Estimate"])
  if (!identity_transform) {
    coef.table[, "Estimate"] <- transformed_coefs
    coef.table <- coef.table[, -2L, drop = FALSE]
  }

  if (ci) {
    conf <- confint(object, parm = which(!aliased), level = level)
    if (!identity_transform) conf[] <- transform(conf)
    coef.table <- cbind(coef.table, conf)
  }

  keep <- match(c("call", "terms", "family", "deviance", "aic", "contrasts", "df.residual",
                  "null.deviance", "df.null", "iter", "na.action", "vcov_type"),
                names(object), 0L)

  ans <- c(object[keep],
           list(coefficients = coef.table,
                aliased = aliased,
                cov.scaled = covmat,
                cluster = .attr(object, "cluster"),
                transformed = !identity_transform))

  class(ans) <- c("summary.glm_weightit", "summary.glm")

  ans
}

#' @exportS3Method print summary.glm_weightit
print.summary.glm_weightit <- function(x, digits = max(3L, getOption("digits") - 3L),
                                       signif.stars = getOption("show.signif.stars"),
                                       call = TRUE,
                                       ...) {
  chk::chk_count(digits)
  chk::chk_flag(call)

  if (call) {
    cat0("\n", .ul("Call:"), "\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
         "\n")
  }

  if (is_null(x$coefficients)) {
    cat("\nNo Coefficients\n\n")
    return(invisible(x))
  }

  cat0("\n", .ul(sprintf("Coefficients%s:",
                         if (x$transformed) " (transformed)" else "")),
       "\n")

  coefs <- x$coefficients

  if (is_not_null(.attr(x, "thresholds"))) {
    coefs <- coefs[-match(.attr(x, "thresholds"), rownames(x$coefficients)), , drop = FALSE]
  }

  aliased <- x$aliased
  if (is_not_null(aliased) && any(aliased)) {
    if (is_not_null(.attr(x, "thresholds"))) {
      aliased <- aliased[-match(.attr(x, "thresholds"), names(aliased))]
    }

    cn <- names(aliased)
    coefs <- matrix(NA_real_, nrow = length(aliased), ncol = ncol(coefs),
                    dimnames = list(cn, colnames(coefs)))
    coefs[!aliased, ] <- x$coefficients
  }

  .printCoefmat_glm_weightit(coefs, digits = digits,
                             has.Pvalue = TRUE,
                             signif.stars = signif.stars,
                             ...)

  cat(.it(sprintf("Standard error: %s\n",
                  .vcov_to_phrase(x$vcov_type,
                                  is_not_null(x[["cluster"]])))))

  if (is_not_null(.attr(x, "thresholds"))) {
    thresholds <- x$coefficients[.attr(x, "thresholds"), , drop = FALSE]

    cat0("\n", .ul(sprintf("Thresholds%s:", if (x$transformed) " (transformed)" else "")),
         "\n")

    .printCoefmat_glm_weightit(thresholds, digits = digits,
                               has.Pvalue = TRUE,
                               signif.stars = signif.stars,
                               ...)
  }

  invisible(x)
}

# print() methods
#' @exportS3Method print glm_weightit
#' @rdname glm_weightit-methods
print.glm_weightit <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {

  cat0("\n", .ul("Call:"), "\n",
       paste(deparse(x$call), collapse = "\n"),
       "\n")

  if (is_null(x$coefficients)) {
    cat("\nNo Coefficients\n\n")
    return(invisible(x))
  }

  co <- x$contrasts
  cat0("\n", .ul(sprintf("Coefficients%s:",
                         if (is.character(co))
                           paste("  [contrasts: ", apply(cbind(names(co), co), 1L, paste, collapse = "="), "]")
                         else "")),
       "\n")

  format(x$coefficients, digits = digits) |>
    print.default(print.gap = 2L, quote = FALSE)

  cat(.it(sprintf("Standard error: %s\n",
                  .vcov_to_phrase(x$vcov_type,
                                  is_not_null(.attr(x, "cluster"))))))

  invisible(x)
}

#' @exportS3Method print multinom_weightit
print.multinom_weightit <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  print.glm_weightit(x, digits = digits, ...)
}

#' @exportS3Method print ordinal_weightit
print.ordinal_weightit <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  print.glm_weightit(x, digits = digits, ...)
}

#' @exportS3Method print coxph_weightit
print.coxph_weightit <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  print.glm_weightit(x, digits = digits, ...)
}

# vcov() methods
#' @exportS3Method stats::vcov glm_weightit
#' @rdname glm_weightit-methods
vcov.glm_weightit <- function(object, complete = TRUE, vcov = NULL, ...) {
  chk::chk_flag(complete)

  .vcov_glm_weightit.internal(object, vcov. = vcov, ...) |>
    .modify_vcov_info() |>
    .vcov.aliased(aliased = is.na(object$coefficients),
                  complete = complete)
}

#' @exportS3Method stats::vcov multinom_weightit
vcov.multinom_weightit <- function(object, complete = TRUE, vcov = NULL, ...) {
  vcov.glm_weightit(object = object, complete = complete, vcov = vcov, ...)
}

#' @exportS3Method stats::vcov ordinal_weightit
vcov.ordinal_weightit <- function(object, complete = TRUE, vcov = NULL, ...) {
  vcov.glm_weightit(object = object, complete = complete, vcov = vcov, ...)
}

#' @exportS3Method stats::vcov coxph_weightit
vcov.coxph_weightit <- function(object, complete = TRUE, vcov = NULL, ...) {
  vcov.glm_weightit(object = object, complete = complete, vcov = vcov, ...)
}

# confint() methods
#' @exportS3Method stats::confint glm_weightit
confint.glm_weightit <- function(object, parm, level = 0.95, vcov = NULL, ...) {
  chk::chk_number(level)
  chk::chk_gt(level, .5)
  chk::chk_lt(level, 1)

  object$df.residual <- Inf

  covmat <- .vcov_glm_weightit.internal(object, vcov. = vcov, ...)

  object |>
    .set_vcov(covmat) |>
    confint.lm(parm = parm, level = level, ...)
}

#' @exportS3Method stats::confint multinom_weightit
confint.multinom_weightit <- function(object, parm, level = 0.95, vcov = NULL, ...) {
  confint.glm_weightit(object = object, parm = parm, level = level, vcov = vcov, ...)
}

#' @exportS3Method stats::confint ordinal_weightit
confint.ordinal_weightit <- function(object, parm, level = 0.95, vcov = NULL, ...) {
  confint.glm_weightit(object = object, parm = parm, level = level, vcov = vcov, ...)
}

#' @exportS3Method stats::confint coxph_weightit
confint.coxph_weightit <- function(object, parm, level = 0.95, vcov = NULL, ...) {
  confint.glm_weightit(object = object, parm = parm, level = level, vcov = vcov, ...)
}

#' @exportS3Method stats::model.matrix multinom_weightit
model.matrix.multinom_weightit <- function(object, ...) {
  class(object) <- "lm"
  model.matrix(object, ...)
}

#' @exportS3Method stats::model.matrix ordinal_weightit
model.matrix.ordinal_weightit <- function(object, ...) {
  x <- model.matrix.multinom_weightit(object, ...)

  x[, colnames(x) != "(Intercept)", drop = FALSE]
}

#' @exportS3Method stats::nobs ordinal_weightit
nobs.ordinal_weightit <- function(object, ...) {
  if (is_not_null(object[["weights"]])) {
    sum(object[["weights"]] != 0)
  }
  else {
    length(object[["residuals"]])
  }
}

#' @exportS3Method stats::nobs multinom_weightit
nobs.multinom_weightit <- function(object, ...) {
  nobs.ordinal_weightit(object, ...)
}

#' @importFrom sandwich estfun
#' @exportS3Method sandwich::estfun glm_weightit
#' @rdname glm_weightit-methods
estfun.glm_weightit <- function(x, asympt = TRUE, ...) {

  # Check missing
  if (is_not_null(x[["na.action"]])) {
    .err("missing values are not allowed in the model variables")
  }

  bout <- x[["coefficients"]]
  aliased <- is.na(bout)

  Xout <- x[["x"]] %or% model.matrix(x)
  Y <- x[["y"]] %or% model.response(model.frame(x))

  if (is_not_null(x[["weightit"]])) {
    W <- x[["weightit"]][["weights"]]
    SW <- x[["weightit"]][["s.weights"]]
  }
  else {
    W <- SW <- NULL
  }

  if (is_null(W)) {
    W <- rep_with(1, Y)
  }

  if (is_null(SW)) {
    SW <- rep_with(1, Y)
  }

  offset <- x[["offset"]] %or% rep_with(0, Y)

  if (any(aliased)) {
    if (is_not_null(.attr(x[["qr"]][["qr"]], "aliased"))) {
      Xout <- Xout[, !.attr(x[["qr"]][["qr"]], "aliased"), drop = FALSE]
    }
    else {
      Xout <- make_full_rank(Xout, with.intercept = FALSE)
    }

    bout <- bout[!aliased]
  }

  Mparts <- .attr(x[["weightit"]], "Mparts")
  Mparts.list <- .attr(x[["weightit"]], "Mparts.list")

  if (is_not_null(Mparts) || is_not_null(Mparts.list)) {
    chk::chk_flag(asympt)
  }
  else {
    asympt <- FALSE
  }

  if (!asympt) {
    Mparts <- Mparts.list <- NULL
  }

  psi_out <- function(Bout, w, Y, Xout, SW, offset) {
    x$psi(Bout, Xout, Y, w * SW, offset = offset)
  }

  psi_b <- x[["gradient"]] %or% psi_out(bout, W, Y, Xout, SW, offset)

  if (is_not_null(Mparts)) {
    # Mparts from weightit()
    psi_treat <- Mparts[["psi_treat"]]
    wfun <- Mparts[["wfun"]]
    Xtreat <- Mparts[["Xtreat"]]
    A <- Mparts[["A"]]
    btreat <- Mparts[["btreat"]]
    hess_treat <- Mparts[["hess_treat"]]
    dw_dBtreat <- Mparts[["dw_dBtreat"]]

    H_treat <- {
      if (is_not_null(hess_treat)) {
        hess_treat(btreat, Xtreat, A, SW)
      }
      else {
        .gradient(function(Btreat) {
          colSums(psi_treat(Btreat, Xtreat, A, SW))
        }, .x = btreat)
      }
    }

    H_out_treat <- {
      if (is_not_null(dw_dBtreat)) {
        crossprod(psi_out(bout, 1, Y, Xout, SW, offset),
                  dw_dBtreat(btreat, Xtreat, A, SW))
      }
      else {
        .gradient(function(Btreat) {
          w <- wfun(Btreat, Xtreat, A)
          colSums(psi_out(bout, w, Y, Xout, SW, offset))
        }, .x = btreat)
      }
    }

    #Using formula from Wooldridge (2010) p. 419
    psi_b <- psi_b + psi_treat(btreat, Xtreat, A, SW) %*%
      .solve_hessian(-H_treat, t(H_out_treat), model = "weights")
  }
  else if (is_not_null(Mparts.list)) {
    # Mparts.list from weightitMSM() or weightit()
    psi_treat.list <- grab(Mparts.list, "psi_treat")
    wfun.list <- grab(Mparts.list, "wfun")
    Xtreat.list <- grab(Mparts.list, "Xtreat")
    A.list <- grab(Mparts.list, "A")
    btreat.list <- grab(Mparts.list, "btreat")
    hess_treat.list <- grab(Mparts.list, "hess_treat")
    dw_dBtreat.list <- grab(Mparts.list, "dw_dBtreat")

    psi_treat <- function(Btreat.list, Xtreat.list, A.list, SW) {
      do.call("cbind", lapply(seq_along(Btreat.list), function(i) {
        psi_treat.list[[i]](Btreat.list[[i]], Xtreat.list[[i]], A.list[[i]], SW)
      }))
    }

    wfun <- function(Btreat.list, Xtreat.list, A.list) {
      Reduce("*", lapply(seq_along(Btreat.list), function(i) {
        wfun.list[[i]](Btreat.list[[i]], Xtreat.list[[i]], A.list[[i]])
      }), init = 1)
    }

    H_treat <- {
      if (all(lengths(hess_treat.list) > 0L)) {
        .block_diag(lapply(seq_along(hess_treat.list), function(i) {
          hess_treat.list[[i]](btreat.list[[i]], Xtreat.list[[i]], A.list[[i]], SW)
        }))
      }
      else {
        .gradient(function(Btreat) {
          Btreat.list <- .vec2list(Btreat, lengths(btreat.list))
          colSums(psi_treat(Btreat.list, Xtreat.list, A.list, SW))
        }, .x = unlist(btreat.list))
      }
    }

    if (all(lengths(dw_dBtreat.list) > 0L)) {
      w.list <- c(lapply(seq_along(btreat.list), function(i) {
        wfun.list[[i]](btreat.list[[i]], Xtreat.list[[i]], A.list[[i]])
      }), list(rep_with(1, A.list[[1L]])))

      dw_dBtreat <- do.call("cbind", lapply(seq_along(btreat.list), function(i) {
        dw_dBtreat.list[[i]](btreat.list[[i]], Xtreat.list[[i]], A.list[[i]], SW) *
          Reduce("*", w.list[-i])
      }))

      H_out_treat <- crossprod(psi_out(bout, 1, Y, Xout, SW, offset), dw_dBtreat)
    }
    else {
      H_out_treat <- .gradient(function(Btreat) {
        Btreat.list <- .vec2list(Btreat, lengths(btreat.list))
        w <- wfun(Btreat.list, Xtreat.list, A.list)
        colSums(psi_out(bout, w, Y, Xout, SW, offset))
      }, .x = unlist(btreat.list))
    }

    #Using formula from Wooldridge (2010) p. 419
    psi_b <- psi_b + psi_treat(btreat.list, Xtreat.list, A.list, SW) %*%
      .solve_hessian(-H_treat, t(H_out_treat), model = "weights")
  }

  colnames(psi_b) <- names(aliased)[!aliased]

  psi_b
}

#' @exportS3Method sandwich::estfun multinom_weightit
estfun.multinom_weightit <- function(x, asympt = TRUE, ...) {
  estfun.glm_weightit(x, asympt = asympt, ...)
}

#' @exportS3Method sandwich::estfun ordinal_weightit
estfun.ordinal_weightit <- function(x, asympt = TRUE, ...) {
  estfun.glm_weightit(x, asympt = asympt, ...)
}

#' @importFrom sandwich bread
#' @exportS3Method sandwich::bread glm_weightit
bread.glm_weightit <- function(x, ...) {

  bout <- x[["coefficients"]]
  aliased <- is.na(bout)

  if (is_not_null(x[["hessian"]])) {
    H <- x$hessian
  }
  else if (inherits(x, "ordinal_weightit")) {
    H <- .get_hess_ordinal(x)
  }
  else if (inherits(x, "multinom_weightit")) {
    H <- .get_hess_multinom(x)
  }
  else if (inherits(x, "glm")) {
    H <- .get_hess_glm(x)
  }
  else {
    Xout <- x[["x"]] %or% model.matrix(x)
    Y <- x[["y"]] %or% model.response(model.frame(x))

    if (is_not_null(x[["weightit"]])) {
      W <- x[["weightit"]][["weights"]]
      SW <- x[["weightit"]][["s.weights"]]
    }
    else {
      W <- SW <- NULL
    }

    if (is_null(W)) {
      W <- rep_with(1, Y)
    }

    if (is_null(SW)) {
      SW <- rep_with(1, Y)
    }

    offset <- x[["offset"]] %or% rep_with(0, Y)

    if (any(aliased)) {
      if (is_not_null(.attr(x[["qr"]][["qr"]], "aliased"))) {
        Xout <- Xout[, !.attr(x[["qr"]][["qr"]], "aliased"), drop = FALSE]
      }
      else {
        Xout <- make_full_rank(Xout, with.intercept = FALSE)
      }
      bout <- bout[!aliased]
    }

    psi_out <- function(Bout, w, Y, Xout, SW, offset) {
      x$psi(Bout, Xout, Y, w * SW, offset = offset)
    }

    H <- .gradient(function(Bout) {
      colSums(psi_out(Bout, W, Y, Xout, SW, offset))
    }, .x = bout)
  }

  A1 <- -nobs(x) * .solve_hessian(H)

  colnames(A1) <- rownames(A1) <- names(aliased)[!aliased]

  A1
}

#' @exportS3Method sandwich::bread multinom_weightit
bread.multinom_weightit <- function(x, ...) {
  bread.glm_weightit(x, ...)
}

#' @exportS3Method sandwich::bread ordinal_weightit
bread.ordinal_weightit <- function(x, ...) {
  bread.glm_weightit(x, ...)
}

#' @importFrom generics tidy
#' @exportS3Method generics::tidy glm_weightit
tidy.glm_weightit <- function(x, conf.int = FALSE, conf.level = 0.95, exponentiate = FALSE, ...) {
  s <- summary(x, ci = conf.int, level = conf.level,
               transform = if (exponentiate) exp else NULL,
               ...)

  ret <- cbind(rownames(s$coefficients),
               as.data.frame(s$coefficients)) |>
    setNames(c("term", "estimate", "std.error", "statistic",
               "p.value"))

  class(ret) <- c("tbl_df", "tbl", "data.frame")

  ret
}

#' @exportS3Method generics::tidy multinom_weightit
tidy.multinom_weightit <- function(x, conf.int = FALSE, conf.level = 0.95, exponentiate = FALSE, ...) {
  tidy.glm_weightit(x = x, conf.int = conf.int, conf.level = conf.level, exponentiate = exponentiate, ...)
}

#' @exportS3Method generics::tidy ordinal_weightit
tidy.ordinal_weightit <- function(x, conf.int = FALSE, conf.level = 0.95, exponentiate = FALSE, ...) {
  tidy.glm_weightit(x = x, conf.int = conf.int, conf.level = conf.level, exponentiate = exponentiate, ...)
}

#' @exportS3Method generics::tidy coxph_weightit
tidy.coxph_weightit <- function(x, conf.int = FALSE, conf.level = 0.95, exponentiate = FALSE, ...) {
  tidy.glm_weightit(x = x, conf.int = conf.int, conf.level = conf.level, exponentiate = exponentiate, ...)
}

#' @importFrom generics glance
#' @exportS3Method generics::glance glm_weightit
glance.glm_weightit <- function(x, ...) {
  ret <- data.frame(nobs = nobs(x))

  class(ret) <- c("tbl_df", "tbl", "data.frame")
  ret
}

#' @exportS3Method generics::glance multinom_weightit
glance.multinom_weightit <- function(x, ...) {
  glance.glm_weightit(x, ...)
}

#' @exportS3Method generics::glance ordinal_weightit
glance.ordinal_weightit <- function(x, ...) {
  glance.glm_weightit(x, ...)
}

#' @exportS3Method generics::glance coxph_weightit
glance.coxph_weightit <- function(x, ...) {
  glance.glm_weightit(x, ...)
}

#' @exportS3Method stats::update glm_weightit
#' @rdname glm_weightit-methods
update.glm_weightit <- function(object, formula. = NULL, ..., evaluate = TRUE) {
  obj_call <- getCall(object)

  if (is_null(obj_call)) {
    .err("need an object with `call` component")
  }

  chk::chk_flag(evaluate)

  refit <- TRUE

  update_call <- match.call(expand.dots = FALSE)
  extras <- update_call[["..."]]

  if (is_not_null(formula.)) {
    extras$formula <- update(formula(object), formula.)
  }

  if (all(c("weights", "s.weights") %in% names(extras))) {
    .err("`weights` and `s.weights` can not both be supplied to `update()`")
  }

  if (utils::hasName(extras, "weights")) {
    names(extras)[names(extras) == "weights"] <- "s.weights"
  }

  vcov_args <- c("vcov", "R", "fwb.args", "cluster")

  weightit_args <- c("data", "s.weights")
  orig_weightit <- NULL

  if (is_not_null(extras)) {
    if (evaluate && all(names(extras) %in% vcov_args)) {
      #Just re-estimate vcov, don't change anything else
      update_call[[1L]] <- .vcov_glm_weightit.internal

      vn_in_extras <- which(names(extras) %in% vcov_args)
      vn_in_call <- which(names(obj_call) %in% vcov_args & names(obj_call) %nin% names(extras))

      update_call <- as.call(c(as.list(update_call[c(1L, match("object", names(update_call)))]),
                               as.list(extras[vn_in_extras]),
                               as.list(obj_call[vn_in_call])))

      if (any(names(extras) == "vcov")) {
        names(update_call)[names(update_call) == "vcov"] <- "vcov."
      }
      else {
        update_call[["vcov."]] <- object[["vcov_type"]]
      }

      V <- eval.parent(update_call)

      object <- .set_vcov(object, V)

      refit <- FALSE
    }
    else if (any(names(extras) %in% weightit_args)) {
      if (utils::hasName(extras, "weightit")) {
        .err(sprintf("when `weightit` is supplied, %s cannot be supplied",
                     word_list(intersect(weightit_args, names(extras)), and.or = "and",
                               is.are = TRUE, quotes = "`")))
      }

      if (!evaluate) {
        .err(sprintf("`evaluate` must be `TRUE` when %s supplied",
                     word_list(intersect(weightit_args, names(extras)), and.or = "or",
                               is.are = TRUE, quotes = "`")))
      }

      weightit_args <- intersect(weightit_args, names(extras))

      if (is_not_null(object[["weightit"]]) && !inherits(object[["weightit"]], "fake_weightit")) {
        #Refit weightit object with new args
        wucall <- update_call
        wucall[[1L]] <- quote(stats::update)
        wucall[["object"]] <- object[["weightit"]]
        wucall[weightit_args] <- extras[weightit_args]
        wucall <- wucall[c(1L, match(c("object", weightit_args), names(wucall)))]

        extras[["weightit"]] <- eval(wucall, envir = object[["weightit"]]$env)

        wucall[[1L]] <- quote(update)
        wucall[["object"]] <- obj_call[["weightit"]]

        orig_weightit <- wucall
      }
      else if (utils::hasName(extras, "s.weights")) {
        #Construct a fake weightit object with just the new components
        data <- {
          if (utils::hasName(extras, "data")) eval.parent(extras[["data"]])
          else object[["data"]]
        }

        s.weights <- eval.parent(extras[["s.weights"]]) |>
          .process.s.weights(data)

        if (is_not_null(s.weights)) {
          extras[["weightit"]] <- list(s.weights = s.weights,
                                       weights = rep_with(1, s.weights),
                                       method = NULL)

          class(extras[["weightit"]]) <- c("weightit", "fake_weightit")

          orig_weightit <- str2lang("(fake_weightit)")
        }
        else if (is_not_null(object[["weightit"]])) {
          obj_call[["weightit"]] <- NULL
          object[["weightit"]] <- NULL
        }
      }

      if (utils::hasName(extras, "s.weights")) {
        extras[["s.weights"]] <- NULL
      }
    }

    existing <- names(extras) %in% names(obj_call)

    for (a in names(extras)[existing]) {
      obj_call[[a]] <- extras[[a]]
    }

    if (!all(existing)) {
      obj_call <- c(as.list(obj_call), extras[!existing])
      obj_call <- as.call(obj_call)
    }
  }

  if (!evaluate) {
    return(obj_call)
  }

  if (!refit) {
    object$call <- obj_call
    return(object)
  }

  if (is_not_null(object[["weightit"]]) &&
      inherits(object[["weightit"]], "fake_weightit") &&
      utils::hasName(obj_call, "weightit")) {
    obj_call[["weightit"]] <- object[["weightit"]]
    orig_weightit <- as.name("<fake_weightit>")
  }

  object <- eval.parent(obj_call)

  if (is_not_null(orig_weightit)) {
    object[["call"]][["weightit"]] <- orig_weightit
  }

  object
}

#Note: the following need to be fun.foo <- fun.bar for environment reasons
#' @exportS3Method stats::update multinom_weightit
update.multinom_weightit <- update.glm_weightit

#' @exportS3Method stats::update ordinal_weightit
update.ordinal_weightit <- update.glm_weightit

#' @exportS3Method stats::update coxph_weightit
update.coxph_weightit <- update.glm_weightit

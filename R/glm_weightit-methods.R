#' Methods for `glm_weightit()` objects
#' @name glm_weightit-methods
#'
#' @description
#' This page documents methods for objects returned by [glm_weightit()], `lm_weightit()`, `ordinal_weightit()`, `multinom_weightit()`, and `coxph_weightit()`. `predict()` methods are described at [predict.glm_weightit()] and `anova()` methods are described at [anova.glm_weightit()].
#'
#' @inheritParams stats::vcov
#' @inheritParams stats::confint
#' @inheritParams stats::print.lm
#' @param object,x an output from one of the above modeling functions.
#' @param ci `logical`; whether to display Wald confidence intervals for estimated coefficients. Default is `FALSE`. (Note: this argument can also be supplied as `conf.int`.)
#' @param level when `ci = TRUE`, the desired confidence level.
#' @param transform the function used to transform the coefficients, e.g., `exp` (which can also be supplied as a string, e.g., `"exp"`); passed to [match.fun()] before being used on the coefficients. When `ci = TRUE`, this is also applied to the confidence interval bounds. If specified, the standard error will be omitted from the output. Default is no transformation.
#' @param thresholds `logical`; whether to include thresholds in the `summary()` output for `ordinal_weightit` objects. Default is `TRUE`.
#' @param vcov either a string indicating the method used to compute the variance of the estimated parameters for `object`, a function used to extract the variance, or the variance matrix itself. Default is to use the variance matrix already present in `object`. If a string or function, arguments passed to `...` are supplied to the method or function. (Note: for `vcov()`, can also be supplied as `type`.)
#' @param complete `logical`; whether the full variance-covariance matrix should be returned also in case of an over-determined system where some coefficients are undefined and `coef(.)` contains `NA`s correspondingly. When `complete = TRUE`, `vcov()` is compatible with `coef()` also in this singular case.
#' @param formula. changes to the model formula, passed to the `new` argument of [update.formula()].
#' @param \dots for `vcov()` or `summary()` or `confint()` with `vcov` supplied, other arguments used to compute the variance matrix depending on the method supplied to `vcov`, e.g., `cluster`, `R`, or `fwb.args`. For `update()`, additional arguments to the call or arguments with changed values. See [glm_weightit()] for details.
#' @param evaluate whether to evaluate the call (`TRUE`, the default) or just return it.
#'
#' @returns
#' `summary()` returns a `summary.glm_weightit()` object, which has its own `print()` method. For `coxph_weightit()` objects, the `print()` and `summary()` methods are more like those for `glm` objects than for `coxph` objects.
#'
#' Otherwise, all methods return the same type of object as their generics.
#'
#' @details
#' `vcov()` by default extracts the parameter covariance matrix already computed by the fitting function, and `summary()` and `confint()` uses this covariance matrix to compute standard errors and Wald confidence intervals (internally calling [confint.lm()]), respectively. Supplying arguments to `vcov` or `...` will compute a new covariance matrix. If `cluster` was supplied to the original fitting function, it will be incorporated into any newly computed covariance matrix unless `cluster = NULL` is specified in `vcov()`, `summary()`, or `confint()`. For other arguments (e.g., `R` and `fwb.args`), the defaults are those used by [glm_weightit()]. Note that for `vcov = "BS"` and `vcov = "FWB"` (and `vcov = "const"` for `multinom_weightit` or `ordinal_weightit` objects), the environment for the fitting function is used, so any changes to that environment may affect calculation. It is always safer to simply recompute the fitted object with a new covariance matrix than to modify it with the `vcov` argument, but it can be quicker to just request a new covariance matrix when refitting the model is slow.
#'
#' `update()` updates a fitted model object with new arguments, e.g., a new model formula, dataset, or variance matrix. When only arguments that control the computation of the variance are supplied, only the variance will be recalculated (i.e., the parameters will not be re-estimated). When `data` is supplied, `weightit` is not supplied, and a `weightit` object was originally passed to the model fitting function, the `weightit` object will be re-fit with the new dataset before the model is refit using the new weights and new data. That is, calling `update(obj, data = d)` is equivalent to calling `update(obj, data = d, weightit = update(obj$weightit, data = d))` when a `weightit` object was supplied to the model fitting function.
#'
#' The `estfun()` method for `multinom_weightit` and `ordinal_weightit` objects (which is used by function in the \pkg{sandwich} package to compute coefficient covariance matrices) simply extracts the `gradient` component of the object. For `glm_weightit` and `coxph_weightit` objects, the `glm` and `coxph` methods are dispatched instead.
#'
#' @seealso
#' [glm_weightit()] for the page documenting `glm_weightit()`, `lm_weightit()`, `ordinal_weightit()`, `multinom_weightit()`, and `coxph_weightit()`. [summary.glm()], [vcov()], [confint()] for the relevant methods pages. [predict.glm_weightit()] for computing predictions from the models. [anova.glm_weightit()] for comparing models using a Wald test.
#'
#' @examples
#' ## See more examples at ?glm_weightit

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
      sum((object$weights * object$residuals^2)[object$weights > 0])/df.r
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

  if (p > 0) {
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
      s.err <- rep.int(NA_real_, length(coef.p))
      pvalue <- tvalue <- s.err
    }
    else {
      s.err <- sqrt(diag(covmat)[!aliased])
      tvalue <- coef.p/s.err
      pvalue <- 2 * pnorm(-abs(tvalue))
    }

    dn <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
    dimnames(coef.table) <- list(names(coef.p), dn)
  }
  else {
    coef.table <- matrix(NA_real_, 0L, 4L)
    dimnames(coef.table) <- list(NULL, c("Estimate", "Std. Error",
                                         "z value", "Pr(>|z|)"))
    covmat <- matrix(NA_real_, 0L, 0L)
  }

  transformed_coefs <- transform(coef.table[,"Estimate"])
  if (!is.numeric(transformed_coefs) || length(transformed_coefs) != length(coef.table[,"Estimate"])) {
    .err("`transform` must return a numeric vector")
  }

  identity_transform <- all(transformed_coefs == coef.table[,"Estimate"])
  if (!identity_transform) {
    coef.table[,"Estimate"] <- transformed_coefs
    coef.table <- coef.table[, -2, drop = FALSE]
  }

  if (ci) {
    conf <- confint(object, parm = which(!aliased), level = level)
    if (!identity_transform) conf[] <- transform(conf)
    coef.table <- cbind(coef.table, conf)
  }

  keep <- match(c("call", "terms", "family", "deviance", "aic", "contrasts", "df.residual",
                  "null.deviance", "df.null", "iter", "na.action", "vcov_type"), names(object), 0L)
  ans <- c(object[keep], list(deviance.resid = residuals(object, type = "deviance"),
                              coefficients = coef.table, aliased = aliased,
                              dispersion = dispersion, df = c(object$rank, df.r, df.f),
                              cov.unscaled = covmat.unscaled, cov.scaled = covmat,
                              cluster = attr(object, "cluster"),
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
  if (p > 0) {
    Qr <- object$qr
    coef.p <- object$coefficients[!aliased]

    covmat.unscaled <- NULL
    df.f <- NCOL(Qr$qr)

    covmat <- .vcov_glm_weightit.internal(object, vcov. = vcov, ...)

    object <- .set_vcov(object, covmat)

    if (is_null(covmat)) {
      ci <- FALSE
      s.err <- rep.int(NA_real_, length(coef.p))
      pvalue <- tvalue <- s.err
    }
    else {
      s.err <- sqrt(diag(covmat)[!aliased])
      tvalue <- coef.p/s.err
      pvalue <- 2 * pnorm(-abs(tvalue))
    }

    dn <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
    dimnames(coef.table) <- list(names(coef.p), dn)
  }
  else {
    coef.table <- matrix(NA_real_, 0L, 4L)
    dimnames(coef.table) <- list(NULL, c("Estimate", "Std. Error",
                                         "z value", "Pr(>|z|)"))
    covmat <- matrix(NA_real_, 0L, 0L)
  }

  transformed_coefs <- transform(coef.table[,"Estimate"])
  if (!is.numeric(transformed_coefs) || length(transformed_coefs) != length(coef.table[,"Estimate"])) {
    .err("`transform` must return a numeric vector")
  }

  identity_transform <- all(transformed_coefs == coef.table[,"Estimate"])
  if (!identity_transform) {
    coef.table[,"Estimate"] <- transformed_coefs
    coef.table <- coef.table[, -2, drop = FALSE]
  }

  if (ci) {
    conf <- confint(object, parm = which(!aliased), level = level)
    if (!identity_transform) conf[] <- transform(conf)
    coef.table <- cbind(coef.table, conf)
  }

  keep <- match(c("call", "terms", "family", "deviance", "aic", "contrasts", "df.residual",
                  "null.deviance", "df.null", "iter", "na.action", "vcov_type"), names(object), 0L)
  ans <- c(object[keep], list(deviance.resid = residuals(object, type = "deviance"),
                              coefficients = coef.table, aliased = aliased,
                              dispersion = dispersion, df = c(object$rank, df.r, df.f),
                              cov.unscaled = covmat.unscaled, cov.scaled = covmat,
                              cluster = attr(object, "cluster"),
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

  thresholds_ind <- seq(nrow(out$coefficients) - nthreshold + 1, nrow(out$coefficients))

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

  df.r <- object$df.residual

  dispersion <- NaN

  if (is.null(transform)) transform <- base::identity
  else transform <- match.fun(transform)

  coef.p <- coef(object)
  aliased <- is.na(coef.p)
  p <- sum(!aliased)
  if (p > 0) {
    coef.p <- object$coefficients[!aliased]

    covmat.unscaled <- NULL
    df.f <- NULL

    covmat <- .vcov_glm_weightit.internal(object, vcov. = vcov, ...)

    object <- .set_vcov(object, covmat)

    if (is_null(covmat)) {
      ci <- FALSE
      s.err <- rep.int(NA_real_, length(coef.p))
      pvalue <- tvalue <- s.err
    }
    else {
      s.err <- sqrt(diag(covmat)[!aliased])
      tvalue <- coef.p/s.err
      pvalue <- 2 * pnorm(-abs(tvalue))
    }
    dn <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
    dimnames(coef.table) <- list(names(coef.p), dn)
  }
  else {
    coef.table <- matrix(NA_real_, 0L, 4L)
    dimnames(coef.table) <- list(NULL, c("Estimate", "Std. Error",
                                         "z value", "Pr(>|z|)"))
    covmat <- matrix(NA_real_, 0L, 0L)
  }

  transformed_coefs <- transform(coef.table[,"Estimate"])
  if (!is.numeric(transformed_coefs) || length(transformed_coefs) != length(coef.table[,"Estimate"])) {
    .err("`transform` must return a numeric vector")
  }

  identity_transform <- all(transformed_coefs == coef.table[,"Estimate"])
  if (!identity_transform) {
    coef.table[,"Estimate"] <- transformed_coefs
    coef.table <- coef.table[, -2, drop = FALSE]
  }

  if (ci) {
    conf <- confint(object, parm = which(!aliased), level = level)
    if (!identity_transform) conf[] <- transform(conf)
    coef.table <- cbind(coef.table, conf)
  }

  keep <- match(c("call", "terms", "family", "deviance", "aic", "contrasts", "df.residual",
                  "null.deviance", "df.null", "iter", "na.action", "vcov_type"), names(object), 0L)
  ans <- c(object[keep], list(deviance.resid = residuals(object, type = "deviance"),
                              coefficients = coef.table, aliased = aliased,
                              dispersion = dispersion, df = c(object$rank, df.r, df.f),
                              cov.unscaled = covmat.unscaled, cov.scaled = covmat,
                              cluster = attr(object, "cluster"),
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
    cat0("\n", underline("Call:"), "\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
         "\n")
  }

  if (is_null(x$coefficients)) {
    cat("\nNo Coefficients\n\n")
    return(invisible(x))
  }

  cat0("\n", underline(sprintf("Coefficients%s:",
                               if (x$transformed) " (transformed)" else "")),
       "\n")

  coefs <- x$coefficients

  if (is_not_null(attr(x, "thresholds"))) {
    coefs <- coefs[-match(attr(x, "thresholds"), rownames(x$coefficients)),, drop = FALSE]
  }

  if (is_not_null(aliased <- x$aliased) && any(aliased)) {
    if (is_not_null(attr(x, "thresholds"))) {
      aliased <- aliased[-match(attr(x, "thresholds"), names(aliased))]
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

  cat(italic(sprintf("Standard error: %s\n",
                     .vcov_to_phrase(x$vcov_type,
                                     is_not_null(x[["cluster"]])))))

  if (is_not_null(attr(x, "thresholds"))) {
    thresholds <- x$coefficients[attr(x, "thresholds"),, drop = FALSE]

    cat0("\n", underline(sprintf("Thresholds%s:", if (x$transformed) " (transformed)" else "")),
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

  cat0("\n", underline("Call:"), "\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
       "\n")

  if (is_not_null(coef(x))) {
    cat0("\n", underline(sprintf("Coefficients%s:",
                                 if (is.character(co <- x$contrasts))
                                   paste("  [contrasts: ", apply(cbind(names(co), co), 1L, paste, collapse = "="), "]")
                                 else "")),
         "\n")

    print.default(format(x$coefficients, digits = digits),
                  print.gap = 2, quote = FALSE)
    cat(italic(sprintf("Standard error: %s\n",
                       .vcov_to_phrase(x$vcov_type,
                                       is_not_null(attr(x, "cluster"))))))
  }
  else {
    cat("No coefficients\n\n")
  }

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

  V <- .vcov_glm_weightit.internal(object, vcov. = vcov, ...)

  V <- .modify_vcov_info(V)

  .vcov.aliased(is.na(object$coefficients), V,
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

  object <- .set_vcov(object, covmat)

  confint.lm(object, parm = parm, level = level, ...)
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

  x[,colnames(x) != "(Intercept)", drop = FALSE]
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

#' @importFrom generics estfun
#' @exportS3Method generics::estfun multinom_weightit
estfun.multinom_weightit <- function(x, ...) {
  x$gradient
}

#' @exportS3Method generics::estfun ordinal_weightit
estfun.ordinal_weightit <- function(x, ...) {
  estfun.multinom_weightit(x, ...)
}

#' @importFrom generics tidy
#' @exportS3Method generics::tidy glm_weightit
tidy.glm_weightit <- function(x, conf.int = FALSE, conf.level = 0.95, exponentiate = FALSE, ...) {
  s <- summary(x, ci = conf.int, level = conf.level,
               transform = if (exponentiate) exp else NULL,
               ...)

  ret <- cbind(rownames(s$coefficients),
               as.data.frame(s$coefficients))
  names(ret) <- c("term", "estimate", "std.error", "statistic",
                  "p.value")

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
  ret <- data.frame(nobs = stats::nobs(x))

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
  call <- getCall(object)

  if (is_null(call)) {
    .err("need an object with `call` component")
  }

  chk::chk_flag(evaluate)

  refit <- TRUE

  ucall <- match.call(expand.dots = FALSE)
  extras <- ucall$...

  if (is_not_null(formula.)) {
    extras$formula <- update(formula(object), formula.)
  }

  vcov_args <- c("vcov", "R", "fwb.args", "cluster")

  weightit_args <- c("data")
  orig_weightit <- NULL

  if (is_not_null(extras)) {
    if (evaluate && all(names(extras) %in% vcov_args)) {
      ucall[[1L]] <- .vcov_glm_weightit.internal

      vn_in_extras <- which(names(extras) %in% vcov_args)
      vn_in_call <- which(names(call) %in% vcov_args & names(call) %nin% names(extras))

      ucall <- as.call(c(as.list(ucall[c(1L, match("object", names(ucall)))]),
                         as.list(extras[vn_in_extras]),
                         as.list(call[vn_in_call])))

      if ("vcov" %in% names(extras)) {
        names(ucall)[names(ucall) == "vcov"] <- "vcov."
      }
      else {
        ucall[["vcov."]] <- object[["vcov_type"]]
      }

      V <- eval.parent(ucall)

      object <- .set_vcov(object, V)

      refit <- FALSE
    }

    if (any(names(extras) %in% weightit_args)) {
      if (is_not_null(object[["weightit"]]) && "weightit" %nin% names(extras)) {
        if (!evaluate) {
          .err(sprintf("`evaluate` must be `TRUE` when %s supplied, `weightit` is not supplied, and the object has a `weightit` component",
                       word_list(intersect(weightit_args, names(extras)), and.or = "or",
                                 is.are = TRUE, quotes = "`")))
        }

        wucall <- ucall
        wucall[[1L]] <- update
        wucall[["object"]] <- object[["weightit"]]
        wucall[weightit_args] <- extras[weightit_args]
        wucall <- wucall[c(1L, match(c("object", weightit_args), names(wucall)))]

        orig_weightit <- call[["weightit"]]

        call[["weightit"]] <- eval.parent(wucall)
      }
    }

    existing <- names(extras) %in% names(call)

    for (a in names(extras)[existing]) {
      call[[a]] <- extras[[a]]
    }

    if (any(!existing)) {
      call <- c(as.list(call), extras[!existing])
      call <- as.call(call)
    }
  }

  if (!evaluate) {
    return(call)
  }

  if (!refit) {
    object$call <- call
    return(object)
  }

  object <- eval.parent(call)

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
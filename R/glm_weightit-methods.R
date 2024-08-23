#' Methods for `glm_weightit()` objects
#' @name glm_weightit-methods
#'
#' @description
#' This page documents methods for objects returned by [glm_weightit()], `lm_weightit()`, `ordinal_weightit()`, `multinom_weightit()`, and `coxph_weightit()`. `predict()` methods are described at [predict.glm_weightit()].
#'
#' @inheritParams stats::vcov
#' @inheritParams stats::confint
#' @inheritParams stats::print.lm
#' @param object,object2,x an output from one of the above modeling functions. For `anova()`, `object2` is required.
#' @param ci `logical`; whether to display Wald confidence intervals for estimated coefficients. Default is `FALSE`. (Note: this argument can also be supplied as `conf.int`.)
#' @param level when `ci = TRUE`, the desired confidence level.
#' @param transform the function used to transform the coefficients, e.g., `exp` (which can also be supplied as a string, e.g., `"exp"`); passed to [match.fun()] before being used on the coefficients. When `ci = TRUE`, this is also applied to the confidence interval bounds. If specified, the standard error will be omitted from the output. Default is no transformation.
#' @param thresholds `logical`; whether to include thresholds in the `summary()` output for `ordinal_weightit` objects. Default is `TRUE`.
#' @param complete `logical`; whether the full variance-covariance matrix should be returned also in case of an over-determined system where some coefficients are undefined and `coef(.)` contains `NA`s correspondingly. When `complete = TRUE`, `vcov()` is compatible with `coef()` also in this singular case.
#' @param \dots ignored.
#' @param test the type of test statistic used to compare models. Currently only `"Chisq"` (the chi-square statistic) is allowed.
#' @param method the kind of test used to compare models. Currently only `"Wald"` is allowed.
#' @param tolerance for the Wald test, the tolerance used to determine if models are symbolically nested.
#'
#' @returns
#' `summary()` returns a `summary.glm_weightit()` object, which has its own `print()` method. For `coxph_weightit()` objects, the `print()` and `summary()` methods are more like those for `glm` objects then for `coxph` objects.
#'
#' Otherwise, all methods return the same type of object as their generics.
#'
#' @details
#' `vcov()` (which is called by `summary()`) simply extracts the covariance matrix already computed by the fitting function. `confint()` computes Wald confidence intervals (internally calling [confint.lm()]). The `estfun()` method for `multinom_weightit` and `ordinal_weightit` objects (which is used by function in the \pkg{sandwich} package to compute coefficient covariance matrices) simply extracts the `gradient` component of the object. For `glm_weightit` and `coxph_weightit` objects, the `glm` and `coxph` methods are dispatched instead.
#'
#' `anova()` performs a Wald test to compare two fitted models. The models must be nested, but they don't have to be nested symbolically (i.e., the names of the coefficients of the smaller model do not have to be a subset of the names of the coefficients of the larger model). The larger model must be supplied to `object` and the smaller to `object2`. Both models must contain the same units, weights (if any), and outcomes. The variance-covariance matrix of the coefficients of the smaller model is not used, so it can be specified as `"none"` in the original model call. Otherwise, a warning be thrown if the covariances were computed using different methods.
#'
#' @seealso
#' [glm_weightit()] for the page documenting `glm_weightit()`, `lm_weightit()`, `ordinal_weightit()`, `multinom_weightit()`, and `coxph_weightit()`. [summary.glm()], [vcov], [confint()] for the relevant methods pages. [predict.glm_weightit()] for computing predictions from the models.
#'
#' @examples
#' ## See more examples at ?glm_weightit
#'
#' data("lalonde", package = "cobalt")
#'
#' # Model comparison for any relationship between `treat`
#' # and `re78` (not the same as testing for the ATE)
#' fit1 <- glm_weightit(
#'   re78 ~ treat * (age + educ + race + married + nodegree +
#'                     re74 + re75), data = lalonde
#' )
#'
#' fit2 <- glm_weightit(
#'   re78 ~ age + educ + race + married + nodegree +
#'     re74 + re75, data = lalonde
#' )
#'
#' anova(fit1, fit2)
#'
#' @examplesIf requireNamespace("splines", quietly = TRUE)
#' # Model comparison between spline model and linear
#' # model; note they are nested but not symbolically
#' # nested
#' fit_s <- glm_weightit(
#'   re78 ~ splines::ns(age, df = 4), data = lalonde
#' )
#'
#' fit_l <- glm_weightit(
#'   re78 ~ age, data = lalonde
#' )
#'
#' anova(fit_s, fit_l)

#' @exportS3Method summary glm_weightit
#' @rdname glm_weightit-methods
summary.glm_weightit <- function(object,
                                 ci = FALSE,
                                 level = .95,
                                 transform = NULL,
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
    else if (fam$family %in% c("poisson", "binomial")) 1
    else if (df.r > 0) {
      if (any(object$weights == 0))
        .wrn("observations with zero weight not used for calculating dispersion")
      sum((object$weights * object$residuals^2)[object$weights > 0])/df.r
    }
    else {
      NaN
    }
  }

  if (is.null(transform)) transform <- base::identity
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

    covmat <- suppressWarnings(vcov(object, complete = TRUE))
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
summary.multinom_weightit <- function(object, ci = FALSE, level = .95, transform = NULL, ...) {
  chk::chk_flag(ci)

  df.r <- object$df.residual

  dispersion <-  1

  if (is.null(transform)) transform <- base::identity
  else transform <- match.fun(transform)

  coef.p <- coef(object)
  aliased <- is.na(coef.p)
  p <- sum(!aliased)
  if (p > 0) {
    p1 <- seq_len(p)
    Qr <- object$qr
    coef.p <- object$coefficients[!aliased]

    covmat.unscaled <- NULL
    df.f <- NCOL(Qr$qr)

    covmat <- suppressWarnings(vcov(object, complete = TRUE))
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
summary.ordinal_weightit <- function(object, ci = FALSE, level = .95, transform = NULL, thresholds = TRUE, ...) {
  chk::chk_flag(thresholds)

  out <- summary.multinom_weightit(object, ci = ci, level = level, transform = transform, ...)

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
summary.coxph_weightit <- function(object, ci = FALSE, level = .95, transform = NULL, ...) {
  chk::chk_flag(ci)

  df.r <- object$df.residual

  dispersion <- NaN

  if (is.null(transform)) transform <- base::identity
  else transform <- match.fun(transform)

  coef.p <- coef(object)
  aliased <- is.na(coef.p)
  p <- sum(!aliased)
  if (p > 0) {
    p1 <- seq_len(p)
    coef.p <- object$coefficients[!aliased]

    covmat.unscaled <- NULL
    df.f <- NULL

    covmat <- suppressWarnings(vcov(object, complete = TRUE))
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
                                       ...) {
  cat0("\n", underline("Call:"), "\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
       "\n")

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
                                     is_not_null(x$cluster)))))

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
                                       is_not_null(x$cluster)))))
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

.vcov_to_phrase <- function(vcov_type, cluster = FALSE) {
  switch(vcov_type,
         "const" = "maximum likelihood",
         "asympt" = {
           if (!cluster)
             "HC0 robust (adjusted for estimation of weights)"
           else
             "HC0 cluster-robust (adjusted for estimation of weights)"
         },
         "HC0" = {
           if (!cluster)
             "HC0 robust"
           else
             "HC0 cluster-robust"
         },
         "BS" = {
           if (!cluster)
             "traditional bootstrap"
           else
             "traditional cluster bootstrap"
         },
         "FWB" = {
           if (!cluster)
             "fractional weighted bootstrap"
           else
             "fractional weighted cluster bootstrap"
         })
}

# vcov() methods
#' @exportS3Method stats::vcov glm_weightit
#' @rdname glm_weightit-methods
vcov.glm_weightit <- function(object, complete = TRUE, ...) {
  if (is_null(object[["vcov"]])) {
    if (!identical(object[["vcov_type"]], "none")) {
      .err("no variance-covariance matrix was found in the supplied object; this is likely a bug")
    }

    .wrn('`vcov` was specified as `"none"` in the original fitting call, so no variance-covariance matrix will be returned')

    return(NULL)
  }

  .vcov.aliased(is.na(object$coefficients), object[["vcov"]],
                complete = complete)
}

#' @exportS3Method stats::vcov multinom_weightit
vcov.multinom_weightit <- function(object, complete = TRUE, ...) {
  vcov.glm_weightit(object, complete = complete, ...)
}

#' @exportS3Method stats::vcov ordinal_weightit
vcov.ordinal_weightit <- function(object, complete = TRUE, ...) {
  vcov.glm_weightit(object, complete = complete, ...)
}

#' @exportS3Method stats::vcov coxph_weightit
vcov.coxph_weightit <- function(object, complete = TRUE, ...) {
  vcov.glm_weightit(object, complete = complete, ...)
}

# confint() methods
#' @exportS3Method stats::confint glm_weightit
confint.glm_weightit <- function(object, parm, level = 0.95, ...) {
  chk::chk_number(level)
  chk::chk_gt(level, .5)
  chk::chk_lt(level, 1)
  object$df.residual <- Inf
  confint.lm(object, parm = parm, level = level, ...)
}

#' @exportS3Method stats::confint multinom_weightit
confint.multinom_weightit <- function(object, parm, level = 0.95, ...) {
  confint.glm_weightit(object, parm = parm, level = level, ...)
}

#' @exportS3Method stats::confint ordinal_weightit
confint.ordinal_weightit <- function(object, parm, level = 0.95, ...) {
  confint.glm_weightit(object, parm = parm, level = level, ...)
}

#' @exportS3Method stats::confint coxph_weightit
confint.coxph_weightit <- function(object, parm, level = 0.95, ...) {
  confint.glm_weightit(object, parm = parm, level = level, ...)
}

#' @exportS3Method stats::model.matrix multinom_weightit
model.matrix.multinom_weightit <- function(object, ...) {
  class(object) <- "lm"
  model.matrix(object, ...)
}

#' @exportS3Method stats::model.matrix ordinal_weightit
model.matrix.ordinal_weightit <- function(object, ...) {
  model.matrix.multinom_weightit(object, ...)
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
               transform = if (exponentiate) exp else NULL)

  ret <- cbind(rownames(s$coefficients),
               as.data.frame(s$coefficients))
  names(ret) <- c("term", "estimate", "std.error", "statistic",
                  "p.value")

  class(ret) <- c("tbl_df", "tbl", "data.frame")
  ret
}

#' @exportS3Method generics::tidy multinom_weightit
tidy.multinom_weightit <- tidy.glm_weightit

#' @exportS3Method generics::tidy ordinal_weightit
tidy.ordinal_weightit <- tidy.glm_weightit

#' @exportS3Method generics::tidy coxph_weightit
tidy.coxph_weightit <- tidy.glm_weightit

#' @importFrom generics glance
#' @exportS3Method generics::glance glm_weightit
glance.glm_weightit <- function(x, ...) {
  ret <- data.frame(nobs = stats::nobs(x))

  class(ret) <- c("tbl_df", "tbl", "data.frame")
  ret
}

#' @exportS3Method generics::glance multinom_weightit
glance.multinom_weightit <- glance.glm_weightit

#' @exportS3Method generics::glance ordinal_weightit
glance.ordinal_weightit <- glance.glm_weightit

#' @exportS3Method generics::glance coxph_weightit
glance.coxph_weightit <- glance.glm_weightit

#' @exportS3Method stats::anova glm_weightit
#' @rdname glm_weightit-methods
anova.glm_weightit <- function(object, object2, test = "Chisq",
                               method = "Wald", tolerance = 1e-7, ...) {

  chk::chk_not_missing(object, "`object`")

  chk::chk_not_missing(object2, "`object2`")
  chk::chk_is(object2, class(object)[1])

  if (!identical(nobs(object), nobs(object2)) ||
      !identical(weights(object), weights(object2))) {
    .err("models must be fit with the same units to be compared")
  }

  if (is_null(object[["y"]]) || is_null(object2[["y"]])) {
    .err("models must be fit with `y = TRUE` to be compared")
  }

  if (!identical(object[["y"]], object2[["y"]])) {
    .err("models must be fit with the same outcomes to be compared")
  }

  chk::chk_string(test)
  test <- match_arg(test, c("Chisq"))

  chk::chk_string(method)
  method <- match_arg(method, c("Wald"))

  if (!identical(attr(object, "vcov_type"), attr(object2, "vcov_type")) &&
      !identical(attr(object2, "vcov_type"), "none")) {
    .wrn("different `vcov` types detected for each model; using the `vcov` from the larger model")
  }

  b1 <- coef(object, complete = FALSE)
  b2 <- coef(object2, complete = FALSE)

  df1 <- nobs(object) - length(b1)
  df2 <- nobs(object2) - length(b2)

  if (df1 >= df2) {
    .err("`object2` does not appear to be nested within `object`")
  }

  Z1 <- .lm.fit(x = model.matrix(object2)[, names(b2), drop = FALSE],
                y = model.matrix(object)[, names(b1), drop = FALSE])$residuals

  Z1_svd <- svd(Z1)
  keep <- Z1_svd$d > tolerance
  q <- sum(keep)

  if (q > df2 - df1) {
    .err("`object2` does not appear to be nested within `object`")
  }

  L <- t(Z1_svd$v[, keep, drop = FALSE])

  V <- vcov(object, complete = FALSE)

  value.hyp <- L %*% b1
  vcov.hyp <- L %*% V %*% t(L)

  SSH <- drop(crossprod(value.hyp, solve(vcov.hyp, value.hyp)))

  title <- "Wald test\n"
  topnote <- sprintf("Model 1: %s\nModel 2: %s",
                     deparse1(formula(object)),
                     deparse1(formula(object2)))

  rval <- make_df(c("Res.Df", "Df", test, sprintf("Pr(>%s)", test)),
                  c("1", "2"))

  rval[[1]] <- c(df1, df2)
  rval[[2]] <- c(NA_integer_, as.integer(q))
  rval[[3]][2] <- SSH
  rval[[4]][2] <- pchisq(SSH, q, lower.tail = FALSE)

  result <- structure(rval,
                      heading = c(title, topnote),
                      class = c("anova", "data.frame"))

  attr(result, "value") <- value.hyp
  attr(result, "vcov") <- vcov.hyp

  result
}

#' @exportS3Method stats::anova ordinal_weightit
anova.ordinal_weightit <- function(object, object2, test = "Chisq",
                                   method = "Wald", tolerance = 1e-7, ...) {

  chk::chk_not_missing(object, "`object`")

  chk::chk_not_missing(object2, "`object2`")
  chk::chk_is(object2, class(object)[1])

  if (!identical(nobs(object), nobs(object2)) ||
      !identical(weights(object), weights(object2))) {
    .err("models must be fit with the same units to be compared")
  }

  if (is_null(object[["y"]]) || is_null(object2[["y"]])) {
    .err("models must be fit with `y = TRUE` to be compared")
  }

  if (!identical(object[["y"]], object2[["y"]])) {
    .err("models must be fit with the same outcomes to be compared")
  }

  nthreshold <- ncol(object$fitted.values) - 1L

  chk::chk_string(test)
  test <- match_arg(test, c("Chisq"))

  chk::chk_string(method)
  method <- match_arg(method, c("Wald"))

  if (!identical(attr(object, "vcov_type"), attr(object2, "vcov_type")) &&
      !identical(attr(object2, "vcov_type"), "none")) {
    .wrn("different `vcov` types detected for each model; using the `vcov` from the larger model")
  }

  b1 <- coef(object)
  b2 <- coef(object2)

  df1 <- nobs(object) - sum(!is.na(b1))
  df2 <- nobs(object2) - sum(!is.na(b2))

  if (df1 >= df2) {
    .err("`object2` does not appear to be nested within `object`")
  }

  b1 <- b1[seq_len(length(b1) - nthreshold)]
  b2 <- b2[seq_len(length(b2) - nthreshold)]

  X1 <- model.matrix(object)
  X2 <- model.matrix(object2)

  nm <- names(na.rem(b1))[names(na.rem(b1)) %in% colnames(X1)[-1]]
  nm2 <- names(na.rem(b2))[names(na.rem(b2)) %in% colnames(X2)[-1]]

  #Drop intercept column
  Z1 <- .lm.fit(x = cbind(1, X2[, nm2, drop = FALSE]),
                y = cbind(1, X1[, nm, drop = FALSE]))$residuals[, -1, drop = FALSE]

  Z1_svd <- svd(Z1)
  keep <- Z1_svd$d > tolerance
  q <- sum(keep)

  if (q > df2 - df1) {
    .err("`object2` does not appear to be nested within `object`")
  }

  L <- t(Z1_svd$v[, keep, drop = FALSE])

  V <- vcov(object)

  value.hyp <- L %*% b1[nm]
  vcov.hyp <- L %*% V[nm, nm, drop = FALSE] %*% t(L)

  SSH <- drop(crossprod(value.hyp, solve(vcov.hyp, value.hyp)))

  title <- "Wald test\n"
  topnote <- sprintf("Model 1: %s\nModel 2: %s",
                     deparse1(formula(object)),
                     deparse1(formula(object2)))

  rval <- make_df(c("Res.Df", "Df", test, sprintf("Pr(>%s)", test)),
                  c("1", "2"))

  rval[[1]] <- c(df1, df2)
  rval[[2]] <- c(NA_integer_, as.integer(q))
  rval[[3]][2] <- SSH
  rval[[4]][2] <- pchisq(SSH, q, lower.tail = FALSE)

  result <- structure(rval,
                      heading = c(title, topnote),
                      class = c("anova", "data.frame"))

  attr(result, "value") <- value.hyp
  attr(result, "vcov") <- vcov.hyp

  result
}

#' @exportS3Method stats::anova multinom_weightit
anova.multinom_weightit <- function(object, object2, test = "Chisq",
                                    method = "Wald", tolerance = 1e-7, ...) {

  chk::chk_not_missing(object, "`object`")

  chk::chk_not_missing(object2, "`object2`")
  chk::chk_is(object2, class(object)[1])

  if (!identical(nobs(object), nobs(object2)) ||
      !identical(weights(object), weights(object2))) {
    .err("models must be fit with the same units to be compared")
  }

  if (is_null(object[["y"]]) || is_null(object2[["y"]])) {
    .err("models must be fit with `y = TRUE` to be compared")
  }

  if (!identical(object[["y"]], object2[["y"]])) {
    .err("models must be fit with the same outcomes to be compared")
  }

  chk::chk_string(test)
  test <- match_arg(test, c("Chisq"))

  chk::chk_string(method)
  method <- match_arg(method, c("Wald"))

  if (!identical(attr(object, "vcov_type"), attr(object2, "vcov_type")) &&
      !identical(attr(object2, "vcov_type"), "none")) {
    .wrn("different `vcov` types detected for each model; using the `vcov` from the larger model")
  }

  b1 <- coef(object)
  b2 <- coef(object2)

  df1 <- nobs(object) - sum(!is.na(b1))
  df2 <- nobs(object2) - sum(!is.na(b2))

  if (df1 >= df2) {
    .err("`object2` does not appear to be nested within `object`")
  }

  X1 <- model.matrix(object)
  X2 <- model.matrix(object2)

  k <- nlevels(object[["y"]]) - 1L

  Z1 <- matrix(0, nrow = nrow(X1), ncol = sum(!is.na(b1)),
               dimnames = list(NULL, names(na.rem(b1))))

  j <- 0
  for (i in 1L + seq_len(k)) {
    in_i <- which(object[["y"]] == levels(object[["y"]])[i])
    b1_i <- b1[(i - 2L) * ncol(X1) + seq_col(X1)]
    b2_i <- b2[(i - 2L) * ncol(X2) + seq_col(X2)]

    Z1[in_i, j + seq_len(sum(!is.na(b1_i)))] <- .lm.fit(x = X2[in_i, !is.na(b2_i), drop = FALSE],
                                                        y = X1[in_i, !is.na(b1_i), drop = FALSE])$residuals
    j <- j + sum(!is.na(b1_i))
  }

  Z1_svd <- svd(Z1)
  keep <- which(Z1_svd$d > tolerance)
  q <- length(keep)

  if (q > df2 - df1) {
    .err("`object2` does not appear to be nested within `object`")
  }

  L <- t(Z1_svd$v[, keep, drop = FALSE])

  V <- vcov(object, complete = FALSE)

  value.hyp <- L %*% na.rem(b1)
  vcov.hyp <- L %*% V %*% t(L)

  SSH <- drop(crossprod(value.hyp, solve(vcov.hyp, value.hyp)))

  title <- "Wald test\n"
  topnote <- sprintf("Model 1: %s\nModel 2: %s",
                     deparse1(formula(object)),
                     deparse1(formula(object2)))

  rval <- make_df(c("Res.Df", "Df", test, sprintf("Pr(>%s)", test)),
                  c("1", "2"))

  rval[[1]] <- c(df1, df2)
  rval[[2]] <- c(NA_integer_, as.integer(q))
  rval[[3]][2] <- SSH
  rval[[4]][2] <- pchisq(SSH, q, lower.tail = FALSE)

  result <- structure(rval,
                      heading = c(title, topnote),
                      class = c("anova", "data.frame"))

  attr(result, "value") <- value.hyp
  attr(result, "vcov") <- vcov.hyp

  result
}

#' @exportS3Method stats::anova coxph_weightit
anova.coxph_weightit <- anova.glm_weightit

.printCoefmat_glm_weightit <- function(x,
                                      digits = max(3L, getOption("digits") - 2L),
                                      signif.stars = TRUE,
                                      signif.legend = FALSE,
                                      dig.tst = max(1L, min(5L, digits - 1L)),
                                      cs.ind = NULL,
                                      tst.ind = NULL,
                                      p.ind = NULL,
                                      zap.ind = integer(),
                                      P.values = NULL,
                                      has.Pvalue = TRUE,
                                      eps.Pvalue = 1e-6,
                                      na.print = ".",
                                      quote = FALSE,
                                      right = TRUE,
                                      ...) {
  if (is_null(d <- dim(x)) || length(d) != 2L) {
    .err("'x' must be coefficient matrix/data frame")
  }

  nm <- colnames(x)

  chk::chk_flag(has.Pvalue)

  if (has.Pvalue) {
    if (is_null(p.ind)) {
      if (is_null(nm)) {
        .err("`has.Pvalue` set to `TRUE` but `p.ind` is `NULL` and no colnames present")
      }

      p.ind <- which(substr(nm, 1L, 3L) %in% c("Pr(", "p-v"))

      if (is_null(p.ind)) {
        .err("`has.Pvalue` set to `TRUE` but `p.ind` is `NULL` and no colnames match p-value strings")
      }
    }
    else {
      chk::chk_whole_number(p.ind)
      chk::chk_subset(p.ind, seq_col(x))
    }
  }
  else {
    if (is_not_null(p.ind)) {
      .err("`has.Pvalue` set to `FALSE` but `p.ind` is not `NULL`")
    }
  }

  if (is_null(P.values)) {
    scp <- getOption("show.coef.Pvalues")
    if (!is.logical(scp) || is.na(scp)) {
      .wrn("option \"show.coef.Pvalues\" is invalid: assuming `TRUE`")
      scp <- TRUE
    }
    P.values <- has.Pvalue && scp
  }
  else {
    chk::chk_flag(P.values)
  }

  if (P.values && !has.Pvalue) {
    .err("`P.values` is `TRUE`, but `has.Pvalue` is not")
  }

  if (is_null(cs.ind)) {
    cs.ind <- which(nm %in% c("Estimate", "Std. Error") | endsWith(nm, " %"))
  }
  else {
    chk::chk_whole_numeric(cs.ind)
    chk::chk_subset(cs.ind, seq_col(x))
  }

  cs.ind <- setdiff(cs.ind, p.ind)

  if (is_null(tst.ind)) {
    tst.ind <- which(endsWith(nm, "value"))
  }
  else {
    chk::chk_whole_numeric(tst.ind)
    chk::chk_subset(tst.ind, seq_col(x))
  }

  tst.ind <- setdiff(tst.ind, p.ind)

  if (any(tst.ind %in% cs.ind)) {
    .err("`tst.ind` must not overlap with `cs.ind`")
  }

  xm <- data.matrix(x)

  if (is_null(tst.ind)) {
    tst.ind <- setdiff(which(endsWith(nm, "value")), p.ind)
  }

  Cf <- array("", dim = d, dimnames = dimnames(xm))

  ok <- !(ina <- is.na(xm))

  for (i in zap.ind) {
    xm[, i] <- zapsmall(xm[, i], digits)
  }

  if (is_not_null(cs.ind)) {
    acs <- abs(coef.se <- xm[, cs.ind, drop = FALSE])
    if (any(ia <- is.finite(acs))) {
      digmin <- 1 + if (length(acs <- acs[ia & acs != 0]))
        floor(log10(range(acs[acs != 0], finite = TRUE)))
      else 0

      Cf[, cs.ind] <- format(round(coef.se, max(1L, digits - digmin)), digits = digits)
    }
  }

  if (is_not_null(tst.ind)) {
    Cf[, tst.ind] <- format(round(xm[, tst.ind], digits = dig.tst),
                            digits = digits)
  }

  if (is_not_null(r.ind <- setdiff(seq_col(xm), c(cs.ind, tst.ind, p.ind)))) {
    for (i in r.ind) {
      Cf[, i] <- format(xm[, i], digits = digits)
    }
  }

  ok[, tst.ind] <- FALSE
  okP <- if (has.Pvalue) ok[, -p.ind] else ok

  x1 <- Cf[okP]
  dec <- getOption("OutDec")

  if (dec != ".") {
    x1 <- chartr(dec, ".", x1)
  }

  x0 <- (xm[okP] == 0) != (as.numeric(x1) == 0)

  if (is_not_null(not.both.0 <- which(x0 & !is.na(x0)))) {
    Cf[okP][not.both.0] <- format(xm[okP][not.both.0], digits = max(1L, digits - 1L))
  }

  if (any(ina)) {
    Cf[ina] <- na.print
  }

  if (any(inan <- is.nan(xm))) {
    Cf[inan] <- "NaN"
  }

  if (P.values) {
    chk::chk_flag(signif.stars)

    if (any(okP <- ok[, p.ind])) {
      pv <- as.vector(xm[, p.ind])
      Cf[okP, p.ind] <- format.pval(pv[okP], digits = dig.tst,
                                    eps = eps.Pvalue)

      signif.stars <- signif.stars && any(pv[okP] < 0.1)

      if (signif.stars) {
        Signif <- symnum(pv, corr = FALSE, na = FALSE,
                         cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                         symbols = c("***", "**", "*", ".", " "))
        Cf <- cbind(Cf, format(Signif))
      }
    }
    else signif.stars <- FALSE
  }
  else {
    if (has.Pvalue) {
      Cf <- Cf[, -p.ind, drop = FALSE]
    }
    signif.stars <- FALSE
  }

  print.default(Cf, quote = quote, right = right, na.print = na.print, ...)

  if (signif.stars) {
    chk::chk_flag(signif.legend)

    if (signif.legend) {
      if ((w <- getOption("width")) < nchar(sleg <- attr(Signif, "legend"))) {
        sleg <- strwrap(sleg, width = w - 2, prefix = space(2))
      }

      cat0("---\nSignif. codes:  ", sleg, fill = w + 4 + max(nchar(sleg, "bytes") - nchar(sleg)))
    }
  }

  invisible(x)
}
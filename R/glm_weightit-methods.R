#' Methods for `glm_weightit()` objects
#' @name glm_weightit-methods
#'
#' @description
#' This page documents methods for objects returned by [glm_weightit()], `lm_weightit()`, `ordinal_weightit()`, `multinom_weightit()`, and `coxph_weightit()`. `predict()` methods are described at [predict.glm_weightit()].
#'
#' @inheritParams stats::vcov
#' @inheritParams stats::confint
#' @inheritParams stats::print.lm
#' @param object,x an output from one of the above modeling functions.
#' @param ci `logical`; whether to display Wald confidence intervals for estimated coefficients. Default is `FALSE`.
#' @param level when `ci = TRUE`, the desired confidence level.
#' @param transform the function used to transform the coefficients, e.g., `exp` (which can also be supplied as a string, e.g., `"exp"`); passed to [match.fun()] before being used on the coefficients. When `ci = TRUE`, this is also applied to the confidence interval bounds. If specified, the standard error will be omitted from the output. Default is no transformation.
#' @param complete `logical`; whether the full variance-covariance matrix should be returned also in case of an over-determined system where some coefficients are undefined and `coef(.)` contains `NA`s correspondingly. When `complete = TRUE`, `vcov()` is compatible with `coef()` also in this singular case.
#' @param \dots ignored.
#'
#' @returns
#' `summary()` returns a `summary.glm_weightit()` object, which has its own print method. For `coxph_weightit()` objects, the `print()` and `summary()` methods are more like those for `glm` objects then for `coxph` objects.
#'
#' Otherwise, all methods return the same type of object as their generics.
#'
#' @details
#' `vcov()` (which is called by `summary()`) simply extracts the covariance matrix already computed by the fitting function. `confint()` computes Wald confidence intervals (internally calling [confint.lm()]). The `estfun()` method for `multinom_weightit` and `ordinal_weightit` objects (which is used by function in the \pkg{sandwich} package to compute coefficient covariance matrices) simply extracts the `gradient` component of the object. For `glm_weightit` and `coxph_weightit` objects, the `glm` and `coxph` methods are dispatched instead.
#'
#' @seealso
#' [glm_weightit()] for the page documenting `glm_weightit()`, `lm_weightit()`, `ordinal_weightit()`, `multinom_weightit()`, and `coxph_weightit()`. [summary.glm()], [vcov], [confint()] for the relevant methods pages. [predict.glm_weightit()] for computing predictions from the models.
#'
#' @examples
#' ## See examples at ?glm_weightit
#'

#' @exportS3Method summary glm_weightit
#' @rdname glm_weightit-methods
summary.glm_weightit <- function(object,
                                 ci = FALSE,
                                 level = .95,
                                 transform = NULL,
                                 ...) {
  chk::chk_flag(ci)

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
summary.ordinal_weightit <- function(object, ci = FALSE, level = .95, transform = NULL, ...) {
  out <- summary.multinom_weightit(object, ci = ci, level = level, transform = transform, ...)

  nthreshold <- ncol(object$fitted.values) - 1

  attr(out, "thresholds") <- rownames(out$coefficients)[-seq(1, nrow(out$coefficients) - nthreshold)]

  out
}

#' @exportS3Method summary coxph_weightit
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
  cat("\n", underline("Call:"), "\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n", sep = "")

  if (length(x$aliased) == 0L) {
    cat("\nNo Coefficients\n\n")
    return(invisible(x))
  }

  cat("\n", underline(paste0("Coefficients", if (x$transformed) " (transformed):" else ":")),
      "\n", sep = "")

  coefs <- x$coefficients

  if (is_not_null(attr(x, "thresholds"))) {
    coefs <- coefs[-match(attr(x, "thresholds"), rownames(x$coefficients)),, drop = FALSE]
  }

  if (!is.null(aliased <- x$aliased) && any(aliased)) {
    if (is_not_null(attr(x, "thresholds"))) {
      aliased <- aliased[-match(attr(x, "thresholds"), names(aliased))]
    }

    cn <- names(aliased)
    coefs <- matrix(NA, length(aliased), ncol(coefs), dimnames = list(cn, colnames(coefs)))
    coefs[!aliased, ] <- x$coefficients
  }

  printCoefmat(coefs, digits = digits, signif.legend = FALSE,
               na.print = ".",
               cs.ind = if (x$transformed) -3L else -(3:4),
               tst.ind = if (x$transformed) 2L else 3L,
               has.Pvalue = TRUE,
               ...)

  cat(italic(sprintf("Standard error: %s\n",
                     .vcov_to_phrase(x$vcov_type,
                                     is_not_null(x$cluster)))))

  if (is_not_null(attr(x, "thresholds"))) {
    thresholds <- x$coefficients[attr(x, "thresholds"),, drop = FALSE]

    cat("\n", underline(paste0("Thresholds", if (x$transformed) " (transformed):" else ":")),
        "\n", sep = "")

    printCoefmat(thresholds, digits = digits, signif.legend = FALSE,
                 na.print = ".",
                 cs.ind = if (x$transformed) -3L else -(3:4),
                 tst.ind = if (x$transformed) 2L else 3L,
                 has.Pvalue = TRUE,
                 ...)
  }

  invisible(x)
}

# print() methods
#' @exportS3Method print glm_weightit
#' @rdname glm_weightit-methods
print.glm_weightit <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {

  cat("\n", underline("Call:"), "\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n", sep = "")

  if (is_not_null(coef(x))) {
    cat("\n", underline(paste0("Coefficients",
                               if (is.character(co <- x$contrasts))
                                 paste("  [contrasts: ", apply(cbind(names(co), co), 1L, paste, collapse = "="), "]"), ":\n")), sep = "")

    print.default(format(x$coefficients, digits = digits),
                  print.gap = 2, quote = FALSE)
    cat(italic(sprintf("Standard error: %s\n",
                       .vcov_to_phrase(x$vcov_type,
                                       is_not_null(x$cluster)))))
  }
  else cat("No coefficients\n\n")

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

    .wrn("`vcov` was specified as `\"none\"` in the original fitting call, so no variance-covariance matrix will be returned")

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
#' @rdname glm_weightit-methods
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

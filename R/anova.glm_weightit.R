#' Methods for `glm_weightit()` objects
#' @name anova.glm_weightit
#'
#' @description
#' `anova()` is used to compare nested models fit with
#' `glm_weightit()`, `mutinom_weightit()`, `ordinal_weightit()`, or
#' `coxph_weightit()` using a Wald test that incorporates uncertainty in
#' estimating the weights (if any).
#'
#' @inheritParams glm_weightit-methods
#' @param object,object2 an output from one of the above modeling functions.
#'   `object2` is required.
#' @param test the type of test statistic used to compare models. Currently only
#'   `"Chisq"` (the chi-square statistic) is allowed.
#' @param method the kind of test used to compare models. Currently only
#'   `"Wald"` is allowed.
#' @param tolerance for the Wald test, the tolerance used to determine if models
#'   are symbolically nested.
#' @param \dots other arguments passed to the function used for computing the
#'   parameter variance matrix, if supplied as a string or function, e.g.,
#'   `cluster`, `R`, or `fwb.args`.
#'
#' @returns
#' An object of class `"anova"` inheriting from class `"data.frame"`.
#'
#' @details
#' `anova()` performs a Wald test to compare two fitted models. The
#' models must be nested, but they don't have to be nested symbolically (i.e.,
#' the names of the coefficients of the smaller model do not have to be a subset
#' of the names of the coefficients of the larger model). The larger model must
#' be supplied to `object` and the smaller to `object2`. Both models must
#' contain the same units, weights (if any), and outcomes. The
#' variance-covariance matrix of the coefficients of the smaller model is not
#' used.
#'
#' @seealso
#' [glm_weightit()] for the page documenting `glm_weightit()`,
#' `lm_weightit()`, `ordinal_weightit()`, `multinom_weightit()`, and
#' `coxph_weightit()`. [anova.glm()] for model comparison of `glm` objects.
#'
#' @examples
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
#' # Using the usual maximum likelihood variance matrix
#' anova(fit1, fit2, vcov = "const")
#'
#' # Using a bootstrapped variance matrix
#' anova(fit1, fit2, vcov = "BS", R = 100)
#'
#' @examplesIf rlang::is_installed("splines")
#' # Model comparison between spline model and linear
#' # model; note they are nested but not symbolically
#' # nested
#' fit_s <- glm_weightit(re78 ~ splines::ns(age, df =4),
#'                       data = lalonde)
#'
#' fit_l <- glm_weightit(re78 ~ age,
#'                       data = lalonde)
#'
#' anova(fit_s, fit_l)

#' @exportS3Method stats::anova glm_weightit
anova.glm_weightit <- function(object, object2, test = "Chisq",
                               method = "Wald", tolerance = 1e-7, vcov = NULL, ...) {

  chk::chk_not_missing(object, "`object`")

  chk::chk_not_missing(object2, "`object2`")
  chk::chk_is(object2, class(object)[1L])

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
  test <- match_arg(test, "Chisq")

  chk::chk_string(method)
  method <- match_arg(method, "Wald")

  b1 <- coef(object, complete = FALSE)
  b2 <- coef(object2, complete = FALSE)

  df1 <- nobs(object) - length(b1)
  df2 <- nobs(object2) - length(b2)

  if (df1 >= df2) {
    .err("`object2` does not appear to be nested within `object`")
  }

  Z1 <- .lm.fit(x = model.matrix(object2)[, names(b2), drop = FALSE],
                y = model.matrix(object)[, names(b1), drop = FALSE])[["residuals"]]

  Z1_svd <- svd(Z1)
  keep <- Z1_svd[["d"]] > tolerance
  .q <- sum(keep)

  if (.q > df2 - df1) {
    .err("`object2` does not appear to be nested within `object`")
  }

  L <- t(Z1_svd[["v"]][, keep, drop = FALSE])

  V <- .process_vcov_anova(object = object,
                           object2 = object2,
                           vcov. = vcov,
                           b1 = b1, ...)

  object <- .set_vcov(object, V)

  value.hyp <- L %*% b1
  vcov.hyp <- L %*% V %*% t(L)

  SSH <- drop(crossprod(value.hyp, solve(vcov.hyp, value.hyp)))

  .title <- paste0("\n", .ul("Wald test"))
  .topnote <- sprintf("Model 1: %s\nModel 2: %s\n",
                      deparse1(formula(object)),
                      deparse1(formula(object2)))

  .varnote <- .it(sprintf("Variance: %s\n",
                          .vcov_to_phrase(object[["vcov_type"]],
                                          is_not_null(.attr(object, "cluster")))))

  result <- make_df(c("Res.Df", "Df", test, sprintf("Pr(>%s)", test)),
                    c("1", "2"))

  result[[1L]] <- c(df1, df2)
  result[[2L]] <- c(NA_integer_, as.integer(.q))
  result[[3L]][2L] <- SSH
  result[[4L]][2L] <- pchisq(SSH, .q, lower.tail = FALSE)

  attr(result, "heading") <- c(.title, .varnote, .topnote)
  attr(result, "value") <- value.hyp
  attr(result, "vcov") <- vcov.hyp

  class(result) <- c("anova", class(result))

  result
}

#' @exportS3Method stats::anova ordinal_weightit
anova.ordinal_weightit <- function(object, object2, test = "Chisq",
                                   method = "Wald", tolerance = 1e-7, vcov = NULL, ...) {

  chk::chk_not_missing(object, "`object`")

  chk::chk_not_missing(object2, "`object2`")
  chk::chk_is(object2, class(object)[1L])

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
  test <- match_arg(test, "Chisq")

  chk::chk_string(method)
  method <- match_arg(method, "Wald")

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

  nm <- names(na.rem(b1))[names(na.rem(b1)) %in% colnames(X1)[-1L]]
  nm2 <- names(na.rem(b2))[names(na.rem(b2)) %in% colnames(X2)[-1L]]

  #Drop intercept column
  Z1 <- .lm.fit(x = cbind(1, X2[, nm2, drop = FALSE]),
                y = cbind(1, X1[, nm, drop = FALSE]))$residuals[, -1L, drop = FALSE]

  Z1_svd <- svd(Z1)
  keep <- Z1_svd$d > tolerance
  .q <- sum(keep)

  if (.q > df2 - df1) {
    .err("`object2` does not appear to be nested within `object`")
  }

  L <- t(Z1_svd$v[, keep, drop = FALSE])

  b1 <- b1[nm]

  V <- .process_vcov_anova(object = object,
                           object2 = object2,
                           vcov. = vcov,
                           b1 = b1, ...)

  object <- .set_vcov(object, V)

  value.hyp <- L %*% b1
  vcov.hyp <- L %*% V %*% t(L)

  SSH <- drop(crossprod(value.hyp, solve(vcov.hyp, value.hyp)))

  .title <- paste0("\n", .ul("Wald test"))
  .topnote <- sprintf("Model 1: %s\nModel 2: %s\n",
                      deparse1(formula(object)),
                      deparse1(formula(object2)))

  .varnote <- .it(sprintf("Variance: %s\n",
                          .vcov_to_phrase(object$vcov_type,
                                          is_not_null(.attr(object, "cluster")))))

  result <- make_df(c("Res.Df", "Df", test, sprintf("Pr(>%s)", test)),
                    c("1", "2"))

  result[[1L]] <- c(df1, df2)
  result[[2L]] <- c(NA_integer_, as.integer(.q))
  result[[3L]][2L] <- SSH
  result[[4L]][2L] <- pchisq(SSH, .q, lower.tail = FALSE)

  attr(result, "heading") <- c(.title, .varnote, .topnote)
  attr(result, "value") <- value.hyp
  attr(result, "vcov") <- vcov.hyp

  class(result) <- c("anova", class(result))

  result
}

#' @exportS3Method stats::anova multinom_weightit
anova.multinom_weightit <- function(object, object2, test = "Chisq",
                                    method = "Wald", tolerance = 1e-7, vcov = NULL, ...) {

  chk::chk_not_missing(object, "`object`")

  chk::chk_not_missing(object2, "`object2`")
  chk::chk_is(object2, class(object)[1L])

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
  test <- match_arg(test, "Chisq")

  chk::chk_string(method)
  method <- match_arg(method, "Wald")

  if (!identical(.attr(object, "vcov_type"), .attr(object2, "vcov_type")) &&
      !identical(.attr(object2, "vcov_type"), "none")) {
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

  j <- 0L
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
  .q <- length(keep)

  if (.q > df2 - df1) {
    .err("`object2` does not appear to be nested within `object`")
  }

  L <- t(Z1_svd$v[, keep, drop = FALSE])

  b1 <- na.rem(b1)

  V <- .process_vcov_anova(object = object,
                           object2 = object2,
                           vcov. = vcov,
                           b1 = b1, ...)

  object <- .set_vcov(object, V)

  value.hyp <- L %*% b1
  vcov.hyp <- L %*% V %*% t(L)

  SSH <- drop(crossprod(value.hyp, solve(vcov.hyp, value.hyp)))

  .title <- paste0("\n", .ul("Wald test"))
  .topnote <- sprintf("Model 1: %s\nModel 2: %s\n",
                      deparse1(formula(object)),
                      deparse1(formula(object2)))

  .varnote <- .it(sprintf("Variance: %s\n",
                          .vcov_to_phrase(object$vcov_type,
                                          is_not_null(.attr(object, "cluster")))))

  result <- make_df(c("Res.Df", "Df", test, sprintf("Pr(>%s)", test)),
                    c("1", "2"))

  result[[1L]] <- c(df1, df2)
  result[[2L]] <- c(NA_integer_, as.integer(.q))
  result[[3L]][2L] <- SSH
  result[[4L]][2L] <- pchisq(SSH, .q, lower.tail = FALSE)

  attr(result, "heading") <- c(.title, .varnote, .topnote)
  attr(result, "value") <- value.hyp
  attr(result, "vcov") <- vcov.hyp

  class(result) <- c("anova", class(result))

  result
}

#' @exportS3Method stats::anova coxph_weightit
anova.coxph_weightit <- anova.glm_weightit

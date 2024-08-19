#' Predictions for `glm_weightit` objects
#'
#' @description
#' `predict()` generates predictions for models fit using `glm_weightit()`, `ordinal_weightit()`, `multinom_weightit()`, or `coxph_weightit()`. This page only details the `predict()` methods after using `glm_weightit()`, `ordinal_weightit()`, or `multinom_weightit()`. See [survival::predict.coxph()] for predictions when fitting Cox proportional hazards models using `coxph_weightit()`.
#'
#' @param object a `glm_weightit` object.
#' @param newdata optionally, a data frame in which to look for variables with which to predict. If omitted, the fitted values applied to the original dataset are used.
#' @param type the type of prediction desired. Allowable options include `"response"`, predictions on the scale of the original response variable (also `"probs"`); `"link"`, predictions on the scale of the linear predictor (also `"lp"`); `"class"`, the modal predicted category for ordinal and multinomial models; and `"mean"`, the expected value of the outcome for ordinal and multinomial models. See Details for more information. The default is `"response"` for all models, which differs from [stats::predict.glm()].
#' @param na.action function determining what should be done with missing values in `newdata`. The default is to predict `NA`.
#' @param values when `type = "mean"`, the numeric values each level corresponds to. Should be supplied as a named vector with outcome levels as the names. If `NULL` and the outcome levels can be converted to numeric, those will be used. See Details.
#' @param \dots further arguments passed to or from other methods.
#'
#' @returns
#' A numeric vector containing the desired predictions, except for the following circumstances when an ordinal or multinomial model was fit:
#' * when `type = "response"`, a numeric matrix with a row for each unit and a column for each level of the outcome with the predicted probability of the corresponding outcome in the cells
#' * when `type = "class"`, a factor with the model predicted class for each unit; for ordinal models, this will be an ordered factor.
#'
#' @details
#' For generalized linear models other than ordinal and multinomial models, see [stats::predict.glm()] for more information on how predictions are computed and which arguments can be specified. Note that standard errors cannot be computed for the predictions using `predict.glm_weightit()`.
#'
#' For ordinal and multinomial models, setting `type = "mean"` computes the expected value of the outcome for each unit; this corresponds to the sum of the values supplied in `values` weighted by the predicted probability of those values. If `values` is omitted, `predict()` will attempt to convert the outcome levels to numeric values, and if this cannot be done, an error will be thrown. `values` should be specified as a named vector, e.g., `values = c(one = 1, two = 2, three = 3)`, where `"one"`, `"two"`, and `"three"` are the original outcome levels and 1, 2, and 3 are the numeric values they correspond to. This method only makes sense to use if the outcome levels meaningfully correspond to numeric values.
#'
#' For ordinal models, setting `type = "link"` (also `"lp"`) computes the linear predictor without including the thresholds. This can be interpreted as the prediction of the latent variable underlying the ordinal response. This cannot be used with multinomial models.
#'
#' @seealso
#' [stats::predict.glm()] for predictions from generalized linear models. [glm_weightit()] for the fitting function. [survival::predict.coxph()] for predictions from Cox proportional hazards models.
#'
#' @examples
#' data("lalonde", package = "cobalt")
#'
#' # Logistic regression model
#' fit1 <- glm_weightit(
#'   re78 > 0 ~ treat * (age + educ + race + married +
#'                         re74 + re75),
#'   data = lalonde, family = binomial, vcov = "none")
#'
#' summary(predict(fit1))
#'
#' # G-computation using predicted probabilities
#' p0 <- predict(fit1, type = "response",
#'               newdata = transform(lalonde,
#'                                   treat = 0))
#'
#' p1 <- predict(fit1, type = "response",
#'               newdata = transform(lalonde,
#'                                   treat = 1))
#'
#' mean(p1) - mean(p0)
#'
#' # Multinomial logistic regression model
#' lalonde$re78_3 <- factor(findInterval(lalonde$re78,
#'                                       c(0, 5e3, 1e4)),
#'                          labels = c("low", "med", "high"))
#'
#' fit2 <- multinom_weightit(
#'   re78_3 ~ treat * (age + educ + race + married +
#'                       re74 + re75),
#'   data = lalonde, vcov = "none")
#'
#' # Predicted probabilities
#' head(predict(fit2))
#'
#' # Class assignment accuracy
#' mean(predict(fit2, type = "class") == lalonde$re78_3)
#'
#' # G-computation using expected value of the outcome
#' values <- c("low" = 2500,
#'             "med" = 7500,
#'             "high" = 12500)
#'
#' p0 <- predict(fit2, type = "mean", values = values,
#'               newdata = transform(lalonde,
#'                                   treat = 0))
#'
#' p1 <- predict(fit2, type = "mean", values = values,
#'               newdata = transform(lalonde,
#'                                   treat = 1))
#'
#' mean(p1) - mean(p0)
#' \donttest{
#' # Ordinal logistic regression
#' fit3 <- ordinal_weightit(
#'   re78 ~ treat * (age + educ + race + married +
#'                     re74 + re75),
#'   data = lalonde, vcov = "none")
#'
#' # G-computation using expected value of the outcome;
#' # using original outcome values
#' p0 <- predict(fit3, type = "mean",
#'               newdata = transform(lalonde,
#'                                   treat = 0))
#'
#' p1 <- predict(fit3, type = "mean",
#'               newdata = transform(lalonde,
#'                                   treat = 1))
#'
#' mean(p1) - mean(p0)
#' }

#' @exportS3Method predict glm_weightit
#' @name predict.glm_weightit
predict.glm_weightit <- function(object, newdata = NULL, type = "response",
                                 na.action = na.pass, ...) {
  chk::chk_string(type)
  type <- switch(type, "probs" = "response", "lp" = "link", type)
  type <- match_arg(type, c("response", "link"))

  stats::predict.glm(object, newdata = newdata, type = type,
                     na.action = na.action, se.fit = FALSE,
                     dispersion = NULL, terms = NULL, ...)
}

#' @exportS3Method predict ordinal_weightit
#' @rdname predict.glm_weightit
predict.ordinal_weightit <- function(object, newdata = NULL, type = "response",
                                  na.action = na.pass, values = NULL, ...) {

  chk::chk_string(type)
  type <- switch(type, "probs" = "response", "lp" = "link", type)
  type <- match_arg(type, c("response", "link", "class", "mean"))

  na.act <- object$na.action
  object$na.action <- NULL

  if (is_null(newdata)) {
    if (type == "response") {
      out <- object$fitted.values
    }
    else if (type == "link") {
      out <- object$linear.predictors
    }
    else if (type == "class") {
      out <- factor(max.col(object$fitted.values, ties.method = "first"),
                    levels = seq_col(object$fitted.values),
                    labels = colnames(object$fitted.values),
                    ordered = inherits(object, "ordinal_weightit"))
    }
    else if (type == "mean") {
      if (is_null(values)) {
        if (!can_str2num(colnames(object$fitted.values))) {
          .err("when `type = \"mean\"` and `values` is not set, the outcome levels must be able to be read as numbers")
        }
        values <- setNames(str2num(colnames(object$fitted.values)),
                           colnames(object$fitted.values))
      }
      else {
        chk::chk_numeric(values)
        chk::chk_named(values)
        if (!all(colnames(object$fitted.values) %in% names(values))) {
          .err("when `type = \"mean\"`, all outcome levels must be named in `values`")
        }
      }

      out <- drop(object$fitted.values %*% values[colnames(object$fitted.values)])
    }

    if (is_not_null(na.act))
      out <- napredict(na.act, out)

    return(out)
  }

  tt <- terms(object)

  Terms <- delete.response(tt)
  m <- model.frame(tt, newdata, na.action = na.action,
                   xlev = object$xlevels)

  if (is_not_null(cl <- attr(Terms, "dataClasses")))
    .checkMFClasses(cl, m)
  x <- model.matrix(Terms, m, contrasts.arg = object$contrasts)

  offset <- model.offset(m)
  if (is_not_null(addO <- object$call$offset)) {
    addO <- eval(addO, newdata, environment(tt))
    offset <- {
      if (is_null(offset)) addO
      else offset + addO
    }
  }
  if (is_null(offset)) offset <- rep.int(0, nrow(x))

  x <- x[,colnames(x) != "(Intercept)", drop = FALSE]

  if (type == "link") {
    return(offset + drop(x %*% object$coefficients[seq_col(x)]))
  }

  p <- object$get_p(object$coefficients, x, offset)

  if (type == "response") {
    return(p)
  }

  if (type == "class") {
    out <- factor(max.col(p, ties.method = "first"),
                  levels = seq_col(p),
                  labels = colnames(p),
                  ordered = inherits(object, "ordinal_weightit"))
  }
  else if (type == "mean") {
    if (is_null(values)) {
      if (!can_str2num(colnames(p))) {
        .err("when `type = \"mean\"` and `values` is not set, the outcome levels must be able to be read as numbers")
      }
      values <- setNames(str2num(colnames(p)), colnames(p))
    }
    else {
      chk::chk_numeric(values)
      chk::chk_named(values)
      if (!all(colnames(p) %in% names(values))) {
        .err("when `type = \"mean\"`, all outcome levels must be named in `values`")
      }
    }

    out <- drop(p %*% values[colnames(p)])
  }

  out
}

#' @exportS3Method predict multinom_weightit
#' @rdname predict.glm_weightit
predict.multinom_weightit <- function(object, newdata = NULL, type = "response",
                                    na.action = na.pass, values = NULL, ...) {

  chk::chk_string(type)
  type <- switch(type, "probs" = "response", type)
  type <- match_arg(type, c("response", "class", "mean"))

  na.act <- object$na.action
  object$na.action <- NULL

  if (is_null(newdata)) {
    if (type == "response") {
      out <- object$fitted.values
    }
    else if (type == "link") {
      out <- object$linear.predictors
    }
    else if (type == "class") {
      out <- factor(max.col(object$fitted.values, ties.method = "first"),
                    levels = seq_col(object$fitted.values),
                    labels = colnames(object$fitted.values))
    }
    else if (type == "mean") {
      if (is_null(values)) {
        if (!can_str2num(colnames(object$fitted.values))) {
          .err("when `type = \"mean\"` and `values` is not set, the outcome levels must be able to be read as numbers")
        }
        values <- setNames(str2num(colnames(object$fitted.values)),
                           colnames(object$fitted.values))
      }
      else {
        chk::chk_numeric(values)
        chk::chk_named(values)
        if (!all(colnames(object$fitted.values) %in% names(values))) {
          .err("when `type = \"mean\"`, all outcome levels must be named in `values`")
        }
      }

      out <- drop(object$fitted.values %*% values[colnames(object$fitted.values)])
    }

    if (is_not_null(na.act))
      out <- napredict(na.act, out)

    return(out)
  }

  tt <- terms(object)

  Terms <- delete.response(tt)
  m <- model.frame(tt, newdata, na.action = na.action,
                   xlev = object$xlevels)

  if (is_not_null(cl <- attr(Terms, "dataClasses")))
    .checkMFClasses(cl, m)
  x <- model.matrix(Terms, m, contrasts.arg = object$contrasts)

  offset <- model.offset(m)
  if (is_not_null(addO <- object$call$offset)) {
    addO <- eval(addO, newdata, environment(tt))
    offset <- {
      if (is_null(offset)) addO
      else offset + addO
    }
  }

  if (is_null(offset)) offset <- rep.int(0, nrow(x))

  p <- object$get_p(object$coefficients, x, offset)

  if (type == "response") {
    return(p)
  }

  if (type == "class") {
    out <- factor(max.col(p, ties.method = "first"),
                  levels = seq_col(p),
                  labels = colnames(p))
  }
  else if (type == "mean") {
    if (is_null(values)) {
      if (!can_str2num(colnames(p))) {
        .err("when `type = \"mean\"` and `values` is not set, the outcome levels must be able to be read as numbers")
      }
      values <- setNames(str2num(colnames(p)), colnames(p))
    }
    else {
      chk::chk_numeric(values)
      chk::chk_named(values)
      if (!all(colnames(p) %in% names(values))) {
        .err("when `type = \"mean\"`, all outcome levels must be named in `values`")
      }
    }

    out <- drop(p %*% values[colnames(p)])
  }

  out
}

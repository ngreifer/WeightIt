#' Calibrate Propensity Score Weights
#' @name calibrate
#'
#' @description
#' `calibrate()` performs Platt scaling to calibrate propensity scores as recommended by Gutman et al. (2022). This involves fitting a new propensity score model using logistic regression with the previously estimated propensity score as the sole predictor. Weights are computed using this new propensity score.
#'
#' @param x A `weightit` object or a vector of propensity scores. Only binary treatments are supported.
#' @param treat A vector of treatment status for each unit. Only binary treatments are supported.
#' @param s.weights A vector of sampling weights or the name of a variable in
#' `data` that contains sampling weights.
#' @param data An optional data frame containing the variable named in `s.weights` when supplied as a string.
#' @param \dots Not used.
#'
#' @returns
#' If the input is a `weightit` object, the output will be a
#' `weightit` object with the propensity scores replaced with the calibrated propensity scores and the weights replaced by weights computed from the calibrated propensity scores.
#'
#' If the input is a numeric vector of weights, the output will be a numeric
#' vector of the calibrated propensity scores.
#'
#' @seealso
#' [weightit()], [weightitMSM()]
#'
#' @references
#' Gutman, R., Karavani, E., & Shimoni, Y. (2022). Propensity score models are better when post-calibrated (arXiv:2211.01221). arXiv. \url{http://arxiv.org/abs/2211.01221}
#'
#' @examplesIf requireNamespace("gbm", quietly = TRUE)
#' library("cobalt")
#' data("lalonde", package = "cobalt")
#'
#' #Using GBM to estimate weights
#' (W <- weightit(treat ~ age + educ + married +
#'                  nodegree + re74, data = lalonde,
#'                method = "gbm", estimand = "ATT",
#'                criterion = "smd.max"))
#' summary(W)
#'
#' #Calibrating the GBM propensity scores
#' Wc <- calibrate(W)
#'
#' #Calibrating propensity scores directly
#' PSc <- calibrate(W$ps, treat = lalonde$treat)
#'

#' @export
calibrate <- function(x, ...) {
  UseMethod("calibrate")
}

#' @exportS3Method calibrate default
#' @rdname calibrate
calibrate.default <- function(x, treat, s.weights = NULL, data = NULL, ...) {
  chk::chk_not_missing(treat, "`treat`")
  if (length(unique(treat)) != 2L) {
    .err("`calibrate()` can only be used with binary treatments")
  }
  chk::chk_numeric(x)

  s.weights <- process.s.weights(s.weights, data)
  if (is_null(s.weights)) s.weights <- rep(1, length(x))

  p <- glm.fit(cbind(1, x), treat, weights = s.weights,
               family = quasibinomial())$fitted.values

  nm <- if (!is_null(names(x))) names(x) else if (!is_null(data)) rownames(data) else names(treat)
  setNames(p, nm)
}

#' @exportS3Method calibrate weightit
#' @rdname calibrate
calibrate.weightit <- function(x, ...) {
  if (is_null(x[["ps"]])) {
    .err("`calibrate()` can only be used on `weightit` objects when propensity scores have been estimated")
  }

  if (get_treat_type(x[["treat"]]) != "binary") {
    .err("`calibrate()` can only be used with binary treatments")
  }

  x$ps[] <- calibrate.default(x[["ps"]], treat = x[["treat"]],
                              s.weights = x[["s.weights"]])

  x$weights[] <- get_w_from_ps(x$ps, x[["treat"]],
                               estimand = x[["estimand"]],
                               focal = x[["focal"]])

  attr(x, "Mparts") <- NULL

  x
}
#' Calibrate Propensity Score Weights
#' @name calibrate
#'
#' @description `calibrate()` calibrates propensity scores used in weights. This
#' involves fitting a new propensity score model using logistic or isotonic
#' regression with the previously estimated propensity score as the sole
#' predictor. Weights are computed using this new propensity score.
#'
#' @param x A `weightit` object or a vector of propensity scores. Only binary
#'   treatments are supported.
#' @param treat A vector of treatment status for each unit. Only binary
#'   treatments are supported.
#' @param s.weights A vector of sampling weights or the name of a variable in
#'   `data` that contains sampling weights.
#' @param data An optional data frame containing the variable named in
#'   `s.weights` when supplied as a string.
#' @param method `character`; the method of calibration used. Allowable options
#'   include `"platt"` (default) for Platt scaling as described by Gutman et al.
#'   (2024) and `"isoreg"` for isotonic regression as described by van der Laan
#'   et al. (2024) and implemented in [isoreg()].
#' @param \dots Not used.
#'
#' @returns If the input is a `weightit` object, the output will be a `weightit`
#' object with the propensity scores replaced with the calibrated propensity
#' scores and the weights replaced by weights computed from the calibrated
#' propensity scores.
#'
#' If the input is a numeric vector of weights, the output will be a numeric
#' vector of the calibrated propensity scores.
#'
#' @seealso [weightit()], [weightitMSM()]
#'
#' @references Gutman, R., Karavani, E., & Shimoni, Y. (2024). Improving Inverse
#' Probability Weighting by Post-calibrating Its Propensity Scores.
#' *Epidemiology*, 35(4). \doi{10.1097/EDE.0000000000001733}
#'
#' van der Laan, L., Lin, Z., Carone, M., & Luedtke, A. (2024). Stabilized
#' Inverse Probability Weighting via Isotonic Calibration (arXiv:2411.06342).
#' arXiv. \url{http://arxiv.org/abs/2411.06342}
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
calibrate.default <- function(x, treat, s.weights = NULL, data = NULL, method = "platt", ...) {
  chk::chk_not_missing(treat, "`treat`")
  if (length(unique(treat)) != 2L) {
    .err("`calibrate()` can only be used with binary treatments")
  }
  chk::chk_numeric(x)

  method <- match_arg(method, c("platt", "isoreg"))

  s.weights <- .process.s.weights(s.weights, data)
  if (is_null(s.weights)) s.weights <- rep.int(1, length(x))

  if (method == "platt") {
    p <- glm.fit(cbind(1, x), treat, weights = s.weights,
                 family = quasibinomial())$fitted.values
  }
  else {
    if (!all_the_same(s.weights)) {
      .wrn("sampling weights will not be incorporated into isotonic regression. Use with caution")
    }

    i0 <- stats::isoreg(1 - x, 1 - treat)
    i1 <- stats::isoreg(x, treat)

    p0 <- 1 - i0$yf[order(i0$ord)]
    p1 <- i1$yf[order(i1$ord)]

    p0 <- squish(p0, lo = min(p0[treat == 0]), hi = Inf)
    p1 <- squish(p1, lo = min(p1[treat == 1]), hi = Inf)

    p <- p0
    p[treat == 1] <- p1[treat == 1]
  }

  nm <- {
    if (is_not_null(names(x))) names(x)
    else if (is_null(data)) names(treat)
    else rownames(data)
  }

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
                              s.weights = x[["s.weights"]], ...)

  x$weights[] <- get_w_from_ps(x$ps, x[["treat"]],
                               estimand = x[["estimand"]],
                               focal = x[["focal"]])

  attr(x, "Mparts") <- NULL

  x
}

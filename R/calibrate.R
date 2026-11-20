#' Calibrate Propensity Score Weights
#' @name calibrate
#'
#' @description
#' `calibrate()` calibrates propensity scores used in weights. This
#' involves fitting a new propensity score model using logistic or isotonic
#' regression with the previously estimated propensity score as the sole
#' predictor. Weights are computed using this new propensity score.
#'
#' @param x a `weightit` object or a vector of propensity scores. Only binary
#'   treatments are supported.
#' @param treat a vector of treatment status for each unit. Only binary
#'   treatments are supported.
#' @param s.weights a vector of sampling weights or the name of a variable in
#'   `data` that contains sampling weights.
#' @param data an optional data frame containing the variable named in
#'   `s.weights` when supplied as a string.
#' @param method `character`; the method of calibration used. Allowable options
#'   include `"platt"` (default) for Platt scaling as described by Gutman et al.
#'   (2024) and `"isoreg"` for isotonic regression as described by van der Laan
#'   et al. (2024).
#' @param treated when `x` is a vector of propensity scores, the value of
#'   `treat` that is considered the "treated" group (i.e., the group for which
#'   the propensity scores are the probability of being in). If `NULL`,
#'   `calibrate()` will attempt to figure it out on its own using the same
#'   heuristics as [get_w_from_ps()]. This only matters when `treat` has values
#'   other than 0 and 1.
#' @param \dots not used.
#'
#' @returns
#' If the input is a `weightit` object, the output will be a `weightit`
#' object with the propensity scores replaced with the calibrated propensity
#' scores and the weights replaced by weights computed from the calibrated
#' propensity scores.
#'
#' If the input is a numeric vector of propensity scores, the output will be a
#' numeric vector of the calibrated propensity scores.
#'
#' @seealso [weightit()], [weightitMSM()]
#'
#' `vignette("weighting-methods")` for the role of calibration with machine learning methods.
#'
#' @references
#' Gutman, R., Karavani, E., & Shimoni, Y. (2024). Improving Inverse
#' Probability Weighting by Post-calibrating Its Propensity Scores.
#' *Epidemiology*, 35(4). \doi{10.1097/EDE.0000000000001733}
#'
#' van der Laan, L., Lin, Z., Carone, M., & Luedtke, A. (2024). Stabilized
#' Inverse Probability Weighting via Isotonic Calibration.
#' arXiv. \url{https://arxiv.org/abs/2411.06342}
#'
#' @examplesIf rlang::is_installed("gbm")
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
calibrate.default <- function(x, treat, s.weights = NULL, data = NULL, method = "platt",
                              treated = NULL, ...) {
  arg::arg_supplied(treat)

  if (length(unique(treat)) != 2L) {
    arg::err("{.fun calibrate} can only be used with binary treatments")
  }

  #The propensity scores are the probabilities of being in the treated group,
  #which is determined as it is by `get_w_from_ps()`
  treated <- .process_focal_and_estimand(NULL, "ATE", treat, treated)[["treated"]]
  treat <- binarize(treat, one = treated)

  arg::arg_numeric(x)

  method <- arg::match_arg(method, c("platt", "isoreg"))

  s.weights <- .process.s.weights(s.weights, data) %or% rep_with(1.0, x)

  if (method == "platt") {
    p <- glm.fit(cbind(1, x), treat, weights = s.weights,
                 family = quasibinomial())$fitted.values
  }
  else {
    p0 <- 1 - .isoregw(1 - x, 1 - treat, s.weights)
    p1 <- .isoregw(x, treat, s.weights)

    p0 <- squish(p0, lo = min(p0[treat == 0]), hi = Inf)
    p1 <- squish(p1, lo = min(p1[treat == 1]), hi = Inf)

    p <- p0
    p[treat == 1] <- p1[treat == 1]
  }

  nm <- names(x) %or% rownames(data) %or% names(treat)

  setNames(p, nm)
}

#' @exportS3Method calibrate weightit
#' @rdname calibrate
calibrate.weightit <- function(x, method = "platt", ...) {
  if (is_null(x[["ps"]])) {
    arg::err("{.fun calibrate} can only be used on {.cls weightit} objects when propensity scores have been estimated")
  }

  if (!identical(get_treat_type(x[["treat"]]), "binary")) {
    arg::err("{.fun calibrate} can only be used with binary treatments")
  }

  treated <- .recorded_treated_level(x[["treat"]]) %or%
    .get_treated_level(x[["treat"]], x[["estimand"]], x[["focal"]])

  x$ps[] <- calibrate.default(x[["ps"]], treat = x[["treat"]],
                              s.weights = x[["s.weights"]],
                              method = method, treated = treated, ...)

  x$weights[] <- get_w_from_ps(x$ps, x[["treat"]],
                               estimand = x[["estimand"]],
                               focal = x[["focal"]],
                               treated = treated)

  attr(x, "calibrate") <- list(method = method)

  attr(x, "trim") <- NULL

  #The M-estimation parts must always reproduce `x$weights`, and they describe the
  #uncalibrated propensity score model. `Mparts.list` (from `by` or stabilization)
  #has to go too, or it would be used in place of `Mparts`.
  attr(x, "Mparts") <- NULL
  attr(x, "Mparts.list") <- NULL

  x
}

#Weighted isotonic regression by the pool-adjacent-violators algorithm. The blocks
#are kept on a stack of preallocated vectors (fitted value, total weight, number
#of units), so each unit is pushed once and merged at most once, and the whole
#fit takes linear time after sorting.
.isoregw <- function(x, y, w = rep.int(1, length(y))) {
  ord <- order(x)
  y <- y[ord]
  w <- w[ord]

  n <- length(y)

  val <- wt <- numeric(n)
  cnt <- integer(n)
  k <- 0L

  for (i in seq_len(n)) {
    k <- k + 1L
    val[k] <- y[i]
    wt[k] <- w[i]
    cnt[k] <- 1L

    #Pool the newest block with the one before it while they violate monotonicity
    while (k > 1L && val[k - 1L] > val[k]) {
      tw <- wt[k - 1L] + wt[k]

      val[k - 1L] <- {
        if (tw > 0) (val[k - 1L] * wt[k - 1L] + val[k] * wt[k]) / tw
        else (val[k - 1L] * cnt[k - 1L] + val[k] * cnt[k]) / (cnt[k - 1L] + cnt[k])
      }

      wt[k - 1L] <- tw
      cnt[k - 1L] <- cnt[k - 1L] + cnt[k]
      k <- k - 1L
    }
  }

  #Fitted values in `x` order, returned in the original order
  yhat <- numeric(n)
  yhat[ord] <- rep.int(val[seq_len(k)], cnt[seq_len(k)])

  yhat
}

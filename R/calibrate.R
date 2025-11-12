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
#' @param \dots not used.
#'
#' @returns
#' If the input is a `weightit` object, the output will be a `weightit`
#' object with the propensity scores replaced with the calibrated propensity
#' scores and the weights replaced by weights computed from the calibrated
#' propensity scores.
#'
#' If the input is a numeric vector of weights, the output will be a numeric
#' vector of the calibrated propensity scores.
#'
#' @seealso [weightit()], [weightitMSM()]
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
calibrate.default <- function(x, treat, s.weights = NULL, data = NULL, method = "platt", ...) {
  chk::chk_not_missing(treat, "`treat`")
  if (length(unique(treat)) != 2L) {
    .err("`calibrate()` can only be used with binary treatments")
  }
  chk::chk_numeric(x)

  chk::chk_string(method)
  method <- match_arg(method, c("platt", "isoreg"))

  s.weights <- .process.s.weights(s.weights, data) %or% rep_with(1, x)

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
    .err("`calibrate()` can only be used on `weightit` objects when propensity scores have been estimated")
  }

  if (!identical(get_treat_type(x[["treat"]]), "binary")) {
    .err("`calibrate()` can only be used with binary treatments")
  }

  x$ps[] <- calibrate.default(x[["ps"]], treat = x[["treat"]],
                              s.weights = x[["s.weights"]],
                              method = method, ...)

  x$weights[] <- get_w_from_ps(x$ps, x[["treat"]],
                               estimand = x[["estimand"]],
                               focal = x[["focal"]])

  attr(x, "calibrate") <- list(method = method)

  attr(x, "trim") <- NULL
  attr(x, "Mparts") <- NULL

  x
}

.isoregw <- function(x, y, w = rep(1, length(y))) {
  stopifnot(length(x) == length(y), length(y) == length(w))

  # Order by x
  ord <- order(x)
  x <- x[ord]
  y <- y[ord]
  w <- w[ord]

  # Start with each point as its own block
  blocks <- lapply(seq_along(y), function(i) list(
    idx = i,
    value = y[i],
    weight = w[i]
  ))

  i <- 1
  while (i < length(blocks)) {
    if (blocks[[i]]$value > blocks[[i + 1]]$value) {
      # merge blocks i and i+1
      new_weight <- blocks[[i]]$weight + blocks[[i + 1L]]$weight
      new_value <- (blocks[[i]]$value * blocks[[i]]$weight +
                      blocks[[i + 1L]]$value * blocks[[i + 1L]]$weight) / new_weight
      new_idx <- c(blocks[[i]]$idx, blocks[[i + 1]]$idx)

      blocks[[i]] <- list(idx = new_idx, value = new_value, weight = new_weight)
      blocks[[i + 1L]] <- NULL

      if (i > 1L) {
        i <- i - 1L
      }
    }
    else {
      i <- i + 1
    }
  }

  # Build fitted values
  fit <- numeric(length(y))
  for (b in blocks) {
    fit[b$idx] <- b$value
  }

  # Return fitted values in original x order
  yhat <- numeric(length(y))
  yhat[ord] <- fit

  yhat
}

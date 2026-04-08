#' Fitting (Weighted) Multinomial Regression Models
#'
#' @description
#' `multinom_weightit()` fits a multinomial logistic regression model with a
#' covariance matrix that accounts for estimation of weights, if supplied. By default, this function uses M-estimation to construct a robust covariance
#' matrix using the estimating equations for the weighting model and the outcome
#' model when available.
#'
#' @inheritParams glm_weightit
#' @param link a string corresponding to the desired link function. Currently, only
#'   `"logit"` is allowed.
#'
#' @returns
#' A `multinom_weightit` object.
#'
#' Unless `vcov = "none"`, the `vcov` component contains the covariance matrix
#' adjusted for the estimation of the weights if requested and a compatible
#' `weightit` object was supplied. The `vcov_type` component contains the type
#' of variance matrix requested. If `cluster` is supplied, it will be stored in
#' the `"cluster"` attribute of the output object, even if not used.
#'
#' The `model` component of the output object (also the `model.frame()` output)
#' will include two extra columns when `weightit` is supplied: `(weights)`
#' containing the weights used in the model (the product of the estimated
#' weights and the sampling weights, if any) and `(s.weights)` containing the
#' sampling weights, which will be 1 if `s.weights` is not supplied in the
#' original `weightit()` call.
#'
#' @details
#' `multinom_weightit()` implements multinomial logistic regression using a
#' custom function in \pkg{WeightIt} that optionally computes a coefficient variance matrix that can be adjusted to
#' account for estimation of the weights if a `weightit` or `weightitMSM` object
#' is supplied to the `weightit` argument. This implementation is less robust to
#' failures than other multinomial logistic regression solvers and should be
#' used with caution. Estimation of coefficients should align with that from
#' `mlogit::mlogit()` and `mclogit::mblogit()` but might differ from `nnet::multinom()` due to the relaxed convergence thresholds of the latter.
#'
#' When no argument is supplied to
#' `weightit` or there is no `"Mparts"` attribute in the supplied object, the
#' default variance matrix returned will be the "HC0" sandwich variance matrix,
#' which is robust to misspecification of the outcome family (including
#' heteroscedasticity). Otherwise, the default variance matrix uses M-estimation
#' to additionally adjust for estimation of the weights. When possible, this
#' often yields smaller (and more accurate) standard errors. See the individual
#' methods pages to see whether and when an `"Mparts"` attribute is included in
#' the supplied object. To request that a variance matrix be computed that
#' doesn't account for estimation of the weights even when a compatible
#' `weightit` object is supplied, set `vcov = "HC0"`, which treats the weights
#' as fixed.
#'
#' Bootstrapping can also be used to compute the coefficient variance matrix;
#' when `vcov = "BS"` or `vcov = "FWB"`, which implement the traditional
#' resampling-based and fractional weighted bootstrap, respectively, the entire
#' process of estimating the weights and fitting the outcome model is repeated
#' in bootstrap samples (if a `weightit` object is supplied). This accounts for
#' estimation of the weights and can be used with any weighting method. It is
#' important to set a seed using `set.seed()` to ensure replicability of the
#' results. The fractional weighted bootstrap is more reliable but requires the
#' weighting method to accept sampling weights (which most do, and you'll get an
#' error if it doesn't). Setting `vcov = "FWB"` and supplying `fwb.args = list(wtype = "multinom")`
#' also performs the resampling-based bootstrap but
#' with the additional features \pkg{fwb} provides (e.g., a progress bar and
#' parallelization).
#'
#' @seealso
#' * [glm_weightit()] for fitting generalized linear models that adjust for estimation of the weights.
#' * [ordinal_weightit()] for fitting ordinal regression models that adjust for estimation of the weights.
#' * [coxph_weightit()] for fitting Cox proportional hazards models that adjust for estimation of the weights.
#' * \pkgfun{mclogit}{mblogit} for fitting multinomial regression models that do not account for estimation of the weights.
#'
#' @examples
#' data("lalonde", package = "cobalt")
#'
#' # Logistic regression ATT weights
#' w.out <- weightit(treat ~ age + educ + married + re74,
#'                   data = lalonde, method = "glm",
#'                   estimand = "ATT")
#'
#' # Multinomial logistic regression outcome model
#' # that adjusts for estimation of weights
#' lalonde$re78_3 <- factor(findInterval(lalonde$re78,
#'                                       c(0, 5e3, 1e4)))
#'
#' fit <- multinom_weightit(re78_3 ~ treat,
#'                          data = lalonde,
#'                          weightit = w.out)
#'
#' summary(fit)

#' @export
multinom_weightit <- function(formula, data, link = "logit", weightit = NULL,
                              vcov = NULL, cluster, R = 500L,
                              offset, start = NULL,
                              control = list(...),
                              x = FALSE, y = TRUE,
                              contrasts = NULL, fwb.args = list(), ...) {

  vcov <- .process_vcov(vcov, weightit, R, fwb.args)

  if (missing(cluster)) {
    cluster <- NULL
  }

  model_call <- match.call()

  ###
  if (is_not_null(...get("family"))) {
    arg::err("{.arg family} cannot be used with {.fun multinom_weightit}")
  }

  internal_model_call <- .build_internal_model_call(model = "multinom",
                                                    model_call = model_call,
                                                    weightit = weightit,
                                                    vcov = vcov)

  fit <- .eval_fit(internal_model_call,
                   errors = c("missing values in object" = "missing values are not allowed in the model variables"),
                   from = FALSE)

  fit$family <- list(family = "multinomial",
                     link = "logit")
  ###

  fit$vcov <- .compute_vcov(fit, weightit, vcov, cluster, model_call, internal_model_call)

  fit <- .process_fit(fit, weightit, vcov, model_call, x, y)

  class(fit) <- "multinom_weightit"

  fit
}

# Multinomial logistic regression
.multinom_weightit.fit <- function(x, y, weights = NULL, offset = NULL, start = NULL,
                                   hess = TRUE, control = list(), ...) {
  arg::arg_atomic(y)
  y <- as.factor(y)
  arg::arg_numeric(x)
  arg::arg_matrix(x)

  if (is_null(colnames(x))) {
    colnames(x) <- paste0("x", seq_col(x))
  }

  N <- length(y)

  if (is_null(weights)) weights <- rep.int(1, N)
  else arg::arg_numeric(weights)

  if (is_null(offset)) offset <- rep.int(0, N)
  else arg::arg_numeric(offset)

  if (!all_the_same(c(length(y), nrow(x), length(weights), length(offset)))) {
    arg::err('{.arg {c("y", "x", "weights", "offset")}} must all have the same number of units')
  }

  K <- nlevels(y) - 1L

  aliased_X <- colnames(x) %nin% colnames(make_full_rank(x, with.intercept = FALSE))
  aliased_B <- rep.int(aliased_X, K)

  k0 <- K * ncol(x)

  if (is_null(start)) {
    start <- rep.int(0, k0)
  }
  else {
    arg::arg_numeric(start)
    arg::arg_length(start, k0)
  }

  nm <- unlist(lapply(levels(y)[-1L], function(i) paste(i, colnames(x), sep = "~")))
  names(start) <- nm

  x_ <- x[, !aliased_X, drop = FALSE]

  get_pp <- function(B, X, offset = NULL) {
    if (is_null(offset)) {
      offset <- 0
    }

    qq <- exp(offset + X %*% matrix(B, nrow = ncol(X)))

    pp <- cbind(1, qq) / (1 + rowSums(qq))

    colnames(pp) <- levels(y)
    rownames(pp) <- rownames(X)

    pp
  }

  #Multinomial logistic regression score
  psi <- function(B, X, y, weights, offset = NULL) {
    pp <- get_pp(B, X, offset)

    out <- do.call("cbind", lapply(levels(y)[-1L], function(i) {
      weights * ((y == i) - pp[, i]) * X
    }))

    if (is_not_null(names(B))) {
      colnames(out) <- names(B)
    }

    out
  }

  gr <- function(B, X, y, weights, offset) {
    colSums(psi(B, X, y, weights, offset))
  }

  ind_mat <- cbind(seq_along(y), as.integer(y))

  ll <- function(B, X, y, weights, offset) {
    p <- get_pp(B, X, offset)[ind_mat]

    sum(weights * log(p))
  }

  m_control <- list(fnscale = -1, #maximize likelihood; optim() minimizes by default
                    trace = 0,
                    maxit = 1e3L,
                    reltol = 1e-12)

  control <- utils::modifyList(m_control, control)

  out <- optim(par = start[!aliased_B],
               ll,
               X = x_,
               y = y,
               weights = weights,
               offset = offset,
               gr = gr,
               method = "BFGS",
               control = control)

  grad <- psi(out$par, X = x_, y = y,
              weights = weights, offset = offset)

  pp <- get_pp(out$par, x_, offset)

  res <- setNames(1 - pp[ind_mat], rownames(x))

  coefs <- rep_with(NA_real_, start)
  coefs[!aliased_B] <- out$par

  fit <- list(coefficients = coefs,
              residuals = res,
              fitted.values = pp,
              solve = out,
              psi = psi,
              f = gr,
              get_p = get_pp,
              df.residual = length(res) - sum(!is.na(coefs)),
              x = x,
              y = y,
              weights = weights,
              gradient = grad)

  if (hess) {
    hessian <- sq_matrix(NA_real_, n = sum(!aliased_B),
                         names = nm[!aliased_B])

    for (i in seq_len(K)) {
      i_ind <- (i - 1L) * ncol(x_) + seq_len(ncol(x_))
      for (j in seq_len(i)) {
        if (i == j) {
          hessian[i_ind, i_ind] <- -crossprod(x_ * ((1 - pp[, i + 1L]) * pp[, i + 1L] * weights), x_)
        }
        else {
          j_ind <- (j - 1L) * ncol(x_) + seq_len(ncol(x_))

          hessian[i_ind, j_ind] <- -crossprod(x_ * (-pp[, i + 1L] * pp[, j + 1L] * weights), x_)
          hessian[j_ind, i_ind] <- t(hessian[i_ind, j_ind])
        }
      }
    }

    fit$hessian <- hessian
  }

  fit
}

.multinom_weightit <- function(formula, data, weights, subset, start = NULL, na.action,
                               hess = TRUE, control = list(), model = TRUE,
                               x = FALSE, y = TRUE, contrasts = NULL, ...) {
  cal <- match.call()

  arg::arg_supplied(formula)
  arg::arg_formula(formula, one_sided = FALSE)
  arg::arg_flag(hess)
  arg::arg_flag(model)
  arg::arg_flag(x)
  arg::arg_flag(y)

  if (missing(data)) {
    data <- environment(formula)
  }

  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"),
             names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())

  mt <- .attr(mf, "terms")

  Y <- model.response(mf, "any")
  if (length(dim(Y)) == 1L) {
    nm <- rownames(Y)
    dim(Y) <- NULL

    if (is_not_null(nm)) {
      names(Y) <- nm
    }
  }

  X <- {
    if (is.empty.model(mt)) matrix(NA_real_, nrow = NROW(Y), ncol = 0L)
    else model.matrix(mt, data = mf, contrasts.arg = contrasts)
  }

  weights <- as.vector(model.weights(mf))

  if (is_not_null(weights)) {
    arg::arg_numeric(weights)
    arg::arg_gte(weights, 0)
  }

  offset <- as.vector(model.offset(mf))
  if (is_not_null(offset)) {
    arg::arg_numeric(offset)

    if (length(offset) != NROW(Y)) {
      arg::err("number of offsets is {length(offset)}; should equal {NROW(Y)} (number of observations)")
    }
  }

  fit <- eval(call(".multinom_weightit.fit",
                   x = X, y = Y, weights = weights,
                   offset = offset, start = start,
                   hess = hess, control = control))

  if (model) fit$model <- mf
  fit$na.action <- .attr(mf, "na.action")
  if (!x) fit$x <- NULL
  if (!y) fit$y <- NULL

  c(fit,
    list(call = cal, formula = formula, terms = mt,
         data = data, offset = offset,
         contrasts = .attr(X, "contrasts"),
         xlevels = .getXlevels(mt, mf)))
}

.get_hess_multinom <- function(fit) {
  x <- fit[["x"]] %or% model.matrix(fit)
  y <- fit[["y"]] %or% model.response(model.frame(fit))
  weights <- weights(fit)
  coefs <- coef(fit)

  y <- as.factor(y)

  if (is_null(colnames(x))) {
    colnames(x) <- paste0("x", seq_col(x))
  }

  if (is_null(weights)) {
    weights <- rep_with(1, y)
  }

  K <- nlevels(y) - 1L

  aliased_X <- colnames(x) %nin% colnames(make_full_rank(x, with.intercept = FALSE))

  x_ <- x[, !aliased_X, drop = FALSE]

  theta0 <- na.rem(coefs)

  pp <- fit$fitted.values

  hessian <- sq_matrix(NA_real_, n = length(theta0),
                       names = names(theta0))

  for (i in seq_len(K)) {
    i_ind <- (i - 1L) * ncol(x_) + seq_len(ncol(x_))
    for (j in seq_len(i)) {
      if (i == j) {
        hessian[i_ind, i_ind] <- -crossprod(x_ * ((1 - pp[, i + 1L]) * pp[, i + 1L] * weights), x_)
      }
      else {
        j_ind <- (j - 1L) * ncol(x_) + seq_len(ncol(x_))

        hessian[i_ind, j_ind] <- -crossprod(x_ * (-pp[, i + 1L] * pp[, j + 1L] * weights), x_)
        hessian[j_ind, i_ind] <- t(hessian[i_ind, j_ind])
      }
    }
  }

  hessian
}

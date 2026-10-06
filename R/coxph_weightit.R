#' Fitting (Weighted) Cox Proportional Hazards Models
#'
#' @description
#' `coxph_weightit()` fits a Cox proportional hazards model with a
#' covariance matrix that accounts for estimation of weights, if supplied, and is a wrapper for functions in the \pkg{survival} package. By default, this function uses M-estimation to construct a robust covariance
#' matrix using the estimating equations for the weighting model and the outcome
#' model when available.
#'
#' @inheritParams glm_weightit
#' @param formula an object of class [`formula`] (or one that can be coerced to
#'   that class): a symbolic description of the model to be fitted. Should include a \pkgfun2{survival}{Surv}{Surv} term as the response. See \pkgfun{survival}{coxph} for how this should be specified.
#' @param control a list of parameters for controlling the fitting process, passed to \pkgfun{survival}{coxph.control}.
#' @param br `logical`; whether to use bias reduction, i.e., to maximize the
#'   partial likelihood penalized by the Jeffreys invariant prior as described by
#'   Heinze and Schemper (2001) rather than the partial likelihood. This yields
#'   estimates with smaller asymptotic bias that are always finite, even when the
#'   maximum partial likelihood estimates are not (e.g., under a monotone
#'   likelihood). Default is `FALSE`. See Details.
#' @param \dots other arguments passed to \pkgfun{survival}{coxph.control}.
#'
#' @returns
#' A `coxph_weightit` object, which inherits from `coxph`. See \pkgfun{survival}{coxph} for details.
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
#' original `weightit()` call. When `weights` is supplied instead, only `(weights)`
#' is included, containing those weights.
#'
#' @details
#' `coxph_weightit()` is essentially a simplified version of \pkgfun{survival}{coxph} to fit weighted
#' survival models that optionally computes a coefficient variance matrix that can be adjusted to
#' account for estimation of the weights if a `weightit` or `weightitMSM` object
#' is supplied to the `weightit` argument. It differs from `coxph()` in a few ways:
#'
#' * the `cluster` argument (if used) should be specified as a one-sided formula (which can include multiple
#' clustering variables) and uses a small sample correction for cluster variance
#' estimates when specified
#' * Special formula components, such as `strata()`, `cluster()`, `pspline()`, `frailty()`, `ridge()`, and `tt()` are not allowed
#' * Only right censoring is allowed, and only two-state models are allowed (i.e., the `Surv()` component of `formula` must be of the form `Surv(time, event)`)
#' * Time-varying predictors are not allowed and there must be one observation per unit (and the `id` and `istate` arguments to `coxph()` are ignored)
#' * Weights of 0 are allowed (`coxph()` rejects them). Such units are omitted from the model fit, which is exact rather than an approximation: a unit with a weight of 0 contributes neither its own term nor anything to any risk set denominator, both of which are weighted by the same weights. Missing values in the model variables are tolerated for these units, which makes it possible to fit an outcome model after censoring, where the event time is unascertained for censored units (see [.cens()]). Missing values in units with a nonzero weight produce an error.
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
#' ## Bias reduction
#'
#' When `br = TRUE`, the coefficients maximize the partial likelihood penalized by
#' the Jeffreys invariant prior, which is equivalent to solving the bias-reducing
#' adjusted score equations of Firth (1993) and removes the first-order term in the
#' asymptotic bias of the estimates. The estimates are finite even when the
#' likelihood is monotone, i.e., when some maximum partial likelihood estimates are
#' infinite, which happens most often in small samples with heavy censoring and
#' strongly predictive covariates (Heinze & Schemper, 2001). Without weights,
#' estimation should align with that from \pkgfun{coxphf}{coxphf}. The penalized
#' likelihood is maximized by Newton-Raphson iterations using the information
#' matrix, which converge more slowly than those for the unpenalized likelihood, so
#' the maximum number of iterations is 100 rather than 20 unless `iter.max` is
#' supplied.
#'
#' As in [multinom_weightit()] and [ordinal_weightit()], the weights are scaled to
#' have a mean of 1 among the units with a nonzero weight before fitting
#' (Mukhopadhyay, 2020), which makes the estimates invariant to multiplying the
#' weights by a constant, and the reported variance matrix uses the information
#' matrix at the estimates, i.e., the adjustment is treated as fixed. M-estimation
#' and bootstrapping can be used with `br = TRUE` just as they can without it. Note
#' that under a monotone likelihood, the distribution of the estimates can be far
#' from normal, so Wald confidence intervals can perform poorly (Heinze & Schemper,
#' 2001); bootstrapping may be preferable in that case.
#'
#' @seealso
#' * \pkgfun{survival}{coxph} for fitting Cox proportional hazards models without adjusting standard errors
#' for estimation of the weights.
#' * \pkgfun{coxphf}{coxphf} for fitting bias-reduced Cox proportional hazards models that do not account for estimation of the weights.
#' * [glm_weightit()] for fitting generalized linear models that adjust for estimation of the weights.
#' * [ordinal_weightit()] and [multinom_weightit()] for fitting ordinal and multinomial regression models that adjust for estimation of the weights.
#'
#' @references
#' Firth, D. (1993). Bias reduction of maximum likelihood estimates.
#' *Biometrika*, 80(1), 27–38. \doi{10.1093/biomet/80.1.27}
#'
#' Heinze, G., & Schemper, M. (2001). A solution to the problem of monotone
#' likelihood in Cox regression. *Biometrics*, 57(1), 114–119.
#' \doi{10.1111/j.0006-341X.2001.00114.x}
#'
#' Mukhopadhyay, P. K. (2020). Firth's penalized likelihood for proportional
#' hazards regressions for complex surveys. *Survey Methodology*, 46(2), 215–241.
#'
#' @examples
#' # See `vignette("estimating-effects")` for an example

#' @export
coxph_weightit <- function(formula, data, weightit = NULL,
                           vcov = NULL, cluster, R = 500L,
                           control = list(...),
                           x = FALSE, y = TRUE,
                           fwb.args = list(), br = FALSE, weights, ...) {

  rlang::check_installed("survival")

  model_call <- match.call()

  if (!missing(weights)) {
    w_out <- .process_weights_arg(substitute(weights), weightit,
                                  data = if (!missing(data)) data,
                                  env = .formula_env(formula),
                                  model_call = model_call)

    weightit <- w_out[["weightit"]]
    model_call <- w_out[["model_call"]]
  }

  vcov <- .process_vcov(vcov, weightit, R, fwb.args,
                        m_est_supported = TRUE)

  if (missing(cluster)) {
    cluster <- NULL
  }

  ##

  arg::arg_flag(br)

  internal_model_call <- .build_internal_model_call(model = "coxph",
                                                    model_call = model_call,
                                                    weightit = weightit,
                                                    vcov = vcov,
                                                    br = br)

  fit <- .eval_fit(internal_model_call,
                   errors = c("missing values in object" = "missing values are not allowed in the model variables"),
                   from = FALSE)

  #No `gradient` component is stored, so `.compute_vcov()` and `estfun()` evaluate
  #`psi` themselves. `residuals(fit, type = "score", weighted = TRUE)` would give
  #the same values (they agree to ~1e-14), but it returns one column per
  #model.matrix column, whereas `.compute_vcov()` reduces the model matrix to the
  #estimable columns; supplying it made rank-deficient fits fail with
  #"non-conformable arguments". It is also undefined for a fit whose row-indexed
  #components were scattered back from a subset fit (see `.coxph_weightit()`).
  fit[["psi"]] <- .get_coxph_psi(fit)

  if (isTRUE(fit[["nevent"]] == 0)) {
    # With no events, the coefficients are all `NA`, so there is nothing to
    # differentiate or invert
    arg::wrn("no events were observed, so no coefficients could be estimated")
  }
  else if (inherits(fit, "coxph.null")) {
    # `coxph.fit()` returns a `coxph.null` object when the model has no
    # covariates, which has no coefficients and so no Hessian
    fit[["coefficients"]] <- setNames(numeric(0L), character(0L))
  }
  else {
    aliased <- is.na(fit[["coefficients"]])

    if (any(aliased)) {
      # `coxph.fit()` returns `NA` for the coefficients of collinear columns of
      # the model matrix. Record which those are so that everything downstream
      # drops the same columns rather than inferring them.
      fit[["aliased"]] <- aliased
    }

    # The model-based ("naive") variance returned by `coxph.fit()` in `var` is
    # the only source for the outcome model Hessian, so store the Hessian
    # before `var` is cleared below.
    fit[["hessian"]] <- .get_hess_coxph(fit)
  }

  fit[["var"]] <- NULL

  ##

  fit[["vcov"]] <- .compute_vcov(fit, weightit, vcov, cluster, model_call, internal_model_call)

  fit <- .process_fit(fit, weightit, vcov, model_call, x, y)

  class(fit) <- c("coxph_weightit", class(fit))

  fit
}

.coxph_weightit <- function(formula, data, weights, subset, na.action,
                            control = list(), model = TRUE,
                            x = FALSE, y = TRUE, contrasts = NULL, br = FALSE, ...) {

  rlang::check_installed("survival")

  #`base::missing()` because `rlang::is_missing()` is `FALSE` for an argument whose
  #default has been evaluated, which `arg::when_not_null()` below does
  control_missing <- missing(control)

  method <- "breslow"

  cal <- match.call()

  arg::arg_supplied(formula)
  arg::arg_formula(formula, one_sided = FALSE)
  arg::arg_flag(model)
  arg::arg_flag(x)
  arg::arg_flag(y)
  arg::arg_flag(br)

  if (...length() > 0L) {
    controlargs <- names(formals(survival::coxph.control))
    indx <- pmatch(...names(), controlargs, nomatch = 0L)

    if (any(indx == 0L)) {
      bad_args <- ...names()[indx == 0L]
      arg::err("argument{?s} {.arg {bad_args}} not matched")
    }
  }

  arg::when_not_null(control, arg::arg_list)

  #Bias-reduced fits converge more slowly than unpenalized ones (see
  #`.coxph_firth.fit()`), so they get more iterations unless `iter.max` is set.
  #Arguments to `coxph.control()` can be abbreviated, so any prefix counts.
  control_names <- as.character({
    if (control_missing) ...names()
    else names(control)
  })

  iter.max_set <- any(startsWith("iter.max", control_names[nzchar(control_names)]))

  if (control_missing) {
    control <- survival::coxph.control(...)
  }
  else {
    control <- do.call(survival::coxph.control, control)
  }

  if (br && !iter.max_set) {
    control$iter.max <- 100L
  }

  newform <- .removeDoubleColonSurv(formula)
  if (is_not_null(newform)) {
    formula <- newform$formula

    if (newform$newcall) {
      cal$formula <- formula
    }
  }

  ss <- "cluster"
  Terms <- {
    if (missing(data)) terms(formula, specials = ss)
    else terms(formula, specials = ss, data = data)
  }

  if (is_not_null(attr(Terms, "specials")$cluster)) {
    arg::err("{.fun cluster} cannot be used in the model formula")
  }

  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action", "offset",
               "id", "istate"),
             names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)

  special <- c("strata", "tt", "frailty", "ridge", "pspline")
  mf$formula <- {
    if (missing(data)) terms(formula, special)
    else terms(formula, special, data = data)
  }

  mf <- eval(mf, parent.frame())
  Terms <- terms(mf)

  specials <- as.list(attr(Terms, "specials"))
  if (any(lengths(specials) > 0L)) {
    arg::err('special terms ({.fun {names(specials)[lengths(specials) > 0L]}}) cannot be used with {.fun coxph_weightit}')
  }

  for (i in c("id", "istate")) {
    if (is_not_null(model.extract(mf, i))) {
      arg::err("{.arg {i}} cannot be used with {.fun coxph_weightit}")
    }
  }

  n <- nrow(mf)
  if (n == 0) {
    arg::err("no (non-missing) observations")
  }

  # Process Y
  Y <- model.response(mf)

  if (!survival::is.Surv(Y) || inherits(Y, "Surv2")) {
    arg::err("the response must be a survival ({.cls Surv}) object")
  }

  type <- attr(Y, "type")

  if (type != "right") {
    arg::err("{.fun coxph_weightit} only supports right-censoring")
  }

  data.n <- nrow(Y)

  if (control$timefix) {
    Y <- survival::aeqSurv(Y)
  }

  # Process X
  xlevels <- .getXlevels(Terms, mf)

  attr(Terms, "intercept") <- 1

  X <- model.matrix(Terms, mf, contrasts.arg = contrasts)

  Xatt <- attributes(X)
  xdrop <- Xatt$assign == 0
  X <- X[, !xdrop, drop = FALSE]

  if (!all(is.finite(X))) {
    arg::err("all predictors must be finite")
  }

  attr(X, "assign") <- Xatt$assign[!xdrop]
  attr(X, "contrasts") <- Xatt$contrasts

  # Process weights
  weights <- as.vector(model.weights(mf))

  if (is_not_null(weights)) {
    arg::arg_numeric(weights)
    arg::arg_gte(weights, 0)

    if (br) {
      weights <- .scale_br_weights(weights)
    }
  }

  #`survival::coxph.fit()` rejects weights of 0, but a unit with a weight of 0
  #contributes nothing to the weighted partial likelihood: neither its own term
  #nor any risk set denominator, which is weighted by the same weights. The fit
  #is therefore identical whether such units are retained or dropped, so they are
  #dropped here and scattered back afterward. The returned object must keep all
  #`n` rows because `.compute_vcov()` aligns it with the full-length weights of
  #the `weightit` object, and a unit with an outcome weight of 0 can still make a
  #nonzero contribution to the estimating equations for the weights.
  pos <- {
    if (is_null(weights)) rep.int(TRUE, n)
    else weights > 0
  }

  if (!any(pos)) {
    arg::err("all weights are 0; no units contribute to the model fit")
  }

  Xmeans <- colMeans(X[pos, , drop = FALSE])

  # Process offset
  offset <- as.vector(model.offset(mf))
  if (is_not_null(offset)) {
    arg::arg_numeric(offset)

    if (length(offset) != NROW(Y)) {
      arg::err("number of offsets is {length(offset)}; should equal {NROW(Y)} (number of observations)")
    }

    if (any(!is.finite(offset) | !is.finite(exp(offset)))) {
      arg::err("offsets must lead to a finite risk score")
    }

    meanoffset <- mean(offset)
    offset <- offset - meanoffset
  }
  else {
    meanoffset <- 0
    offset <- rep.int(0, nrow(mf))
  }

  temp <- c("(Intercept)", attr(Terms, "term.labels"))[attr(X, "assign") + 1L]
  assign <- split(seq(along.with = temp), factor(temp, levels = unique(temp)))

  # Fit Cox model, omitting any units with a weight of 0 (see above)
  fit <- {
    if (sum(Y[, ncol(Y)][pos]) == 0) {
      # No events, so nothing can be estimated; mirror `survival::coxph()` by
      # returning `NA` coefficients instead of failing. Shaped like the output
      # of `coxph.fit()` on the retained units so the processing below, including
      # the scatter back to the full sample, applies unchanged.
      ncoef <- ncol(X)

      list(coefficients = setNames(rep.int(NA_real_, ncoef), colnames(X)),
           var = sq_matrix(NA_real_, ncoef),
           loglik = c(0, 0), score = 0, iter = 0,
           linear.predictors = offset[pos],
           residuals = rep.int(0, sum(pos)),
           means = Xmeans,
           method = "breslow",
           class = "coxph")
    }
    else if (br && ncol(X) > 0L) {
      .coxph_firth.fit(x = X[pos, , drop = FALSE], y = Y[pos],
                       offset = offset[pos], weights = weights[pos],
                       control = control, rownames = row.names(mf)[pos])
    }
    else {
      survival::coxph.fit(x = X[pos, , drop = FALSE], y = Y[pos], strata = NULL,
                          offset = offset[pos], init = NULL,
                          control = control, weights = weights[pos],
                          method = "breslow", rownames = row.names(mf)[pos],
                          nocenter = c(-1, 0, 1))
    }
  }

  if (is.character(fit)) {
    fit <- list(fail = fit)
    class(fit) <- "coxph"
  }
  else {
    if (!all(pos)) {
      #Scatter the row-indexed components back to the full sample. The values for
      #units with a weight of 0 are arbitrary but must be finite; NAs there would
      #propagate through `0 * NA` into `psi` and the variance computation. The
      #linear predictor uses the same centering as `coxph.fit()` does.
      b <- fit$coefficients
      b[is.na(b)] <- 0

      fit$linear.predictors <- drop(X %*% b) + offset - sum(b * fit$means)

      res <- setNames(rep.int(0, n), row.names(mf))
      res[pos] <- fit$residuals
      fit$residuals <- res
    }

    fit$n <- data.n
    fit$nevent <- sum(Y[, ncol(Y)][pos])
    fit$terms <- Terms
    fit$assign <- assign

    class(fit) <- fit$class
    fit$class <- NULL

    fit$na.action <- attr(mf, "na.action")

    if (model) fit$model <- mf
    if (x) fit$x <- X
    if (y) fit$y <- Y

    fit$timefix <- control$timefix
  }

  if (is_not_null(weights) && !all_the_same(weights)) {
    fit$weights <- weights
  }

  names(fit$means) <- names(fit$coefficients)
  fit$formula <- formula(Terms)

  fit$xlevels <- .getXlevels(Terms, mf)
  fit$contrasts <- .attr(X, "contrasts")

  if (meanoffset != 0) {
    fit$linear.predictors <- fit$linear.predictors + meanoffset
  }

  if (x && !all(offset == 0)) {
    fit$offset <- offset
  }

  fit$call <- cal
  fit$br <- br

  fit
}

# Sums of the rows of `M` over each unit's risk set, i.e., over the units whose
# time is at least as late. `ranks` are the integer ranks of the times, with tied
# times sharing a rank.
.risk_set_sums <- function(M, ranks) {
  M <- as.matrix(M)

  agg <- rowsum(M, ranks, reorder = TRUE)
  k <- nrow(agg)

  cum <- apply(agg[k:1L, , drop = FALSE], 2L, cumsum)

  if (!is.matrix(cum)) {
    cum <- matrix(cum, nrow = k)
  }

  cum[k:1L, , drop = FALSE][ranks, , drop = FALSE]
}

# The Breslow log partial likelihood, its score, the information matrix, and
# Firth's (1993) adjustment to the score at `B`, as used by Heinze and Schemper
# (2001). The adjustment, half the trace of the inverse information times the
# derivative of the information with respect to each coefficient, is a sum over
# the events of the third central moments of the covariates in each risk set,
# contracted with the inverse information; `adj` splits it among the events, with
# a row of 0s for every other unit. `pen_loglik` is the log partial likelihood
# penalized by half the log determinant of the information, whose gradient is the
# adjusted score, and is `-Inf` where the information cannot be inverted.
.coxph_firth_parts <- function(B, X, status, weights, offset, ranks) {
  p <- ncol(X)

  #The largest linear predictor is subtracted to prevent overflow; it cancels out
  #of every quantity below
  eta <- drop(X %*% B) + offset
  eta <- eta - max(eta)
  r <- weights * exp(eta)

  XX <- X[, rep(seq_len(p), times = p), drop = FALSE] *
    X[, rep(seq_len(p), each = p), drop = FALSE]

  ev <- which(status > 0 & weights > 0)

  S <- .risk_set_sums(cbind(r, r * X, r * XX), ranks)[ev, , drop = FALSE]

  S0 <- S[, 1L]
  m <- S[, 1L + seq_len(p), drop = FALSE] / S0
  E2 <- S[, 1L + p + seq_len(p^2), drop = FALSE] / S0

  d <- weights[ev]

  loglik <- sum(d * (eta[ev] - log(S0)))
  score <- colSums(d * (X[ev, , drop = FALSE] - m))
  info <- matrix(colSums(d * E2), p, p) - crossprod(m, d * m)

  A <- try(solve(info), silent = TRUE)

  if (null_or_error(A)) {
    return(list(loglik = loglik, score = score, info = info,
                pen_loglik = -Inf))
  }

  #With q = x'Ax, the contracted third central moment for event i is
  #E[qx] - mE[q] - 2E[xx']Am + 2m(m'Am), where m is the mean of the covariates in
  #its risk set and the expectations are over that risk set
  q <- rowSums((X %*% A) * X)
  Q <- .risk_set_sums(cbind(r * q, r * q * X), ranks)[ev, , drop = FALSE] / S0

  Am <- m %*% A

  E2Am <- matrix(vapply(seq_len(p), function(j) {
    rowSums(E2[, (seq_len(p) - 1L) * p + j, drop = FALSE] * Am)
  }, numeric(length(ev))), nrow = length(ev))

  adj <- matrix(0, nrow = nrow(X), ncol = p)
  adj[ev, ] <- .5 * d * (Q[, -1L, drop = FALSE] - m * Q[, 1L] - 2 * E2Am +
                           2 * m * rowSums(Am * m))

  list(loglik = loglik,
       score = score,
       info = info,
       adj = adj,
       adj_score = score + colSums(adj),
       pen_loglik = loglik + .5 * as.numeric(determinant(info, logarithm = TRUE)$modulus))
}

# Fits a Cox model by maximizing the Breslow partial likelihood penalized by the
# Jeffreys invariant prior (Heinze & Schemper, 2001), returning the same
# components as `survival::coxph.fit()`. The penalized estimates are found by
# Newton-Raphson using the information matrix in place of the Hessian of the
# penalized log-likelihood, which, because the adjustment is evaluated afresh at
# each step, converges linearly rather than quadratically; the step is halved
# until the penalized log-likelihood does not decrease. Everything else is then
# computed by `coxph.fit()` itself, evaluated at those estimates without
# iterating, so that the residuals, linear predictors, and variance are exactly
# what an unpenalized fit with the same coefficients would have.
.coxph_firth.fit <- function(x, y, offset, weights, control, rownames) {
  if (is_null(weights)) {
    weights <- rep.int(1, nrow(x))
  }

  #Columns collinear with each other or with the baseline hazard are dropped, as
  #`coxph.fit()` drops them
  aliased <- colnames(x) %nin% colnames(make_full_rank(x, with.intercept = TRUE))

  if (all(aliased)) {
    return(survival::coxph.fit(x = x, y = y, strata = NULL, offset = offset,
                               init = NULL, control = control, weights = weights,
                               method = "breslow", rownames = rownames,
                               nocenter = c(-1, 0, 1)))
  }

  x_ <- x[, !aliased, drop = FALSE]

  status <- y[, "status"]

  ranks <- rank(y[, "time"]) |>
    factor() |>
    unclass()

  parts <- function(B) {
    .coxph_firth_parts(B, x_, status, weights, offset, ranks)
  }

  B <- rep.int(0, ncol(x_))
  pt <- parts(B)

  if (!is.finite(pt[["pen_loglik"]])) {
    .solve_info(pt[["info"]])
  }

  loglik0 <- pt[["loglik"]]

  converged <- FALSE
  iter <- 0L

  for (iter in seq_len(control$iter.max)) {
    step <- drop(solve(pt[["info"]], pt[["adj_score"]]))

    for (h in 0:30) {
      pt_new <- parts(B + step)

      if (pt_new[["pen_loglik"]] >= pt[["pen_loglik"]] - 1e-10 * (1 + abs(pt[["pen_loglik"]]))) {
        break
      }

      step <- step / 2
    }

    B <- B + step
    pt <- pt_new

    if (max(abs(step)) < control$eps) {
      converged <- TRUE
      break
    }
  }

  if (!converged) {
    arg::wrn("the penalized partial likelihood did not converge in {control$iter.max} iteration{?s}; estimates should not be trusted. Try increasing {.code iter.max} in {.arg control}")
  }

  control$iter.max <- 0L

  fit <- survival::coxph.fit(x = x_, y = y, strata = NULL, offset = offset,
                             init = B, control = control, weights = weights,
                             method = "breslow", rownames = rownames,
                             nocenter = c(-1, 0, 1))

  coefs <- setNames(rep.int(NA_real_, ncol(x)), colnames(x))
  coefs[!aliased] <- fit$coefficients

  V <- sq_matrix(0, n = ncol(x))
  V[!aliased, !aliased] <- fit$var

  means <- setNames(rep.int(0, ncol(x)), colnames(x))
  means[!aliased] <- fit$means

  fit$coefficients <- coefs
  fit$var <- V
  fit$means <- means
  fit$loglik[1L] <- loglik0
  fit$iter <- iter

  fit
}

.get_coxph_psi <- function(fit) {
  if (inherits(fit, "coxph.null")) {
    # A null model has no coefficients, so its score contributions form an
    # n x 0 matrix
    return(function(B, X, y, weights, offset = 0) {
      matrix(0, nrow = NROW(X), ncol = 0L)
    })
  }

  .y <- fit[["y"]] %or% model.response(model.frame(fit))

  ranks <- rank(.y[, "time"]) |>
    factor() |>
    unclass()

  #Units with a weight of 0 in the fit are excluded from the estimating equation
  #entirely: their own contribution is 0 and they are absent from every risk set.
  #The mask comes from the fit's weights rather than from the `weights` argument
  #so that `psi` describes the fitted estimating equation at whatever weights it
  #is evaluated, e.g., when the weights are perturbed to differentiate with
  #respect to the coefficients of the weighting model; a unit given a nonzero
  #weight there would otherwise return to the risk sets and perturb the psi of
  #the units that do contribute.
  .w <- fit[["weights"]]

  pos <- {
    if (is_null(.w)) TRUE
    else .w > 0
  }

  #For a bias-reduced fit, the weights are scaled as they were for the fit, and
  #each event's share of the adjustment to the score is added to its contribution
  br <- isTRUE(fit[["br"]])

  psi <- function(B, X, y, weights, offset = 0) {
    weights <- weights * pos

    if (br) {
      weights <- .scale_br_weights(weights)
    }

    time <- y[, "time"]
    status <- y[, "status"]

    p  <- exp(drop(offset + X %*% B))
    Wp <- weights * p

    # Map each unique Y_S value to an integer rank (ties get same rank)
    # This lets us treat the risk set condition as a comparison of integer ranks

    # Aggregate Wp and Wp*X by rank, then compute cumulative sums from the top
    # S0[i] = sum of Wp over all j where Y_S[j] >= Y_S[i]
    #       = sum of rank_Wp[r] for r >= ranks[i]  (cumsum from top)
    rank_Wp   <- rowsum(Wp,     ranks, reorder = TRUE)        # (max(ranks) x 1)
    rank_WpX  <- rowsum(Wp * X, ranks, reorder = TRUE)        # (max(ranks) x ncol(X))

    cum_Wp  <- rev(cumsum(rev(rank_Wp)))                      # cumsum from top rank down
    cum_WpX <- apply(rank_WpX, 2L, function(col) rev(cumsum(rev(col))))

    S0 <- cum_Wp[ranks] |> squish(lo = 1e-8, hi = Inf)        # (n x 1)
    S1 <- cum_WpX[ranks, , drop = FALSE]                      # (n x ncol(X))

    M <- status * (X - S1 / S0)

    WdS0     <- weights * status / S0
    WS1dS0sq <- weights * status * S1 / S0^2

    rank_WdS0     <- rowsum(WdS0,     ranks, reorder = TRUE)
    rank_WS1dS0sq <- rowsum(WS1dS0sq, ranks, reorder = TRUE)

    cum_WdS0     <- cumsum(rank_WdS0)                          # cumsum from bottom rank up
    cum_WS1dS0sq <- apply(rank_WS1dS0sq, 2L, cumsum)

    term1 <- cum_WdS0[ranks]
    term2 <- cum_WS1dS0sq[ranks, , drop = FALSE]

    M <- M - p * term1 * X + p * term2

    if (!br) {
      return(weights * M)
    }

    pt <- .coxph_firth_parts(B, X, status, weights, offset, ranks)

    if (is_null(pt[["adj"]])) {
      .solve_info(pt[["info"]])
    }

    weights * M + pt[["adj"]]
  }
}

.get_hess_coxph <- function(fit) {
  V <- fit[["naive.var"]] %or% fit[["var"]]

  if (is_null(V)) {
    return(NULL)
  }

  # `coxph.fit()` zeroes out the rows and columns of `var` belonging to aliased
  # (collinear) coefficients, which would make it singular
  aliased <- is.na(fit[["coefficients"]])

  if (any(aliased)) {
    V <- V[!aliased, !aliased, drop = FALSE]
  }

  # A degenerate fit can have a variance that isn't invertible; returning
  # `NULL` lets callers fall back to computing the Hessian numerically, which
  # routes the failure through `.solve_hessian()`'s more informative error
  tryCatch(solve(-V), error = function(e) NULL)
}

.removeDoubleColonSurv <- function(formula) {
  sname <- c("Surv", "strata", "cluster", "pspline", "tt",
             "frailty", "ridge", "frailty", "frailty.gaussian", "frailty.gamma",
             "frailty.t")

  cname <- paste0("survival::", sname[-1L])

  found1 <- found2 <- found3 <- NULL

  fix <- function(expr) {
    if (is.call(expr)) {
      if (!is.na(i <- match(deparse1(expr[[1L]]), sname)))
        found2 <<- c(found2, sname[i])
      else if (!is.na(i <- match(deparse1(expr[[1L]]), cname))) {
        found1 <<- c(found1, sname[i + 1])
        expr[[1]] <- str2lang(paste0(sname[i + 1L], "()"))[[1L]]
      }

      for (i in seq_along(expr)[-1L]) {
        if (is_not_null(expr[[i]])) {
          expr[[i]] <- fix(expr[[i]])
        }
      }
    }
    else if (is.name(expr) && !is.na(i <- match(as.character(expr), sname))) {
      found3 <<- c(found3, sname[i])
    }

    expr
  }

  newform <- fix(formula)
  found <- unique(c(found1, found2))

  if (is_not_null(found3)) {
    found <- found[!(found %in% found2)]
  }

  if (is_null(found)) {
    return(NULL)
  }

  list(formula = .addSurvFun(newform, found), newcall = FALSE)
}

.addSurvFun <- function(formula, found) {
  myenv <- new.env(parent = environment(formula))
  tt <- function(x) x
  for (i in found) {
    assign(i, get0(i, envir = rlang::current_env(), mode = "function", inherits = FALSE,
                   ifnotfound = get0(i, envir = asNamespace("survival"), mode = "function")),
           envir = myenv)
  }
  environment(formula) <- myenv
  formula
}

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
#' original `weightit()` call.
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
#' * Weights of 0 are allowed (`coxph()` rejects them). Such units are omitted from the model fit, which is exact rather than an approximation: a unit with a weight of 0 contributes neither its own term nor anything to any risk set denominator, both of which are weighted by the same weights. Missing values in the model variables are tolerated for these units, which makes it possible to fit an outcome model after censoring, where the event time is unascertained for censored units (see [`.cens()`][.cens]). Missing values in units with a nonzero weight produce an error.
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
#' * \pkgfun{survival}{coxph} for fitting Cox proportional hazards models without adjusting standard errors
#' for estimation of the weights.
#' * [glm_weightit()] for fitting generalized linear models that adjust for estimation of the weights.
#' * [ordinal_weightit()] and [multinom_weightit()] for fitting ordinal and multinomial regression models that adjust for estimation of the weights.
#'
#' @examples
#' # See `vignette("estimating-effects")` for an example

#' @export
coxph_weightit <- function(formula, data, weightit = NULL,
                           vcov = NULL, cluster, R = 500L,
                           control = list(...),
                           x = FALSE, y = TRUE,
                           fwb.args = list(), ...) {

  rlang::check_installed("survival")

  vcov <- .process_vcov(vcov, weightit, R, fwb.args,
                        m_est_supported = TRUE)

  if (missing(cluster)) {
    cluster <- NULL
  }

  model_call <- match.call()

  if (is_not_null(...get("weights"))) {
    arg::wrn("{.arg weights} is not an allowable argument to {.fun {rlang::call_name(model_call)}} and will be ignored. To fit a weighted model, supply a {.cls weightit} or {.cls weightitMSM} object to the {.arg weightit} argument")

    model_call[["weights"]] <- NULL
  }

  ##

  internal_model_call <- .build_internal_model_call(model = "coxph",
                                                    model_call = model_call,
                                                    weightit = weightit,
                                                    vcov = vcov)

  fit <- .eval_fit(internal_model_call,
                   errors = c("missing values in object" = "missing values are not allowed in the model variables"),
                   from = FALSE)

  #No `gradient` component is stored, so `.compute_vcov()` and `estfun()` evaluate
  #`psi` themselves. `residuals(fit, type = "score", weighted = TRUE)` would give
  #the same values (they agree to ~1e-14), but it returns one column per
  #model.matrix column, whereas `.compute_vcov()` reduces the model matrix to the
  #estimable columns; supplying it made rank-deficient fits fail with
  #"non-conformable arrays". It is also undefined for a fit whose row-indexed
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
                            x = FALSE, y = TRUE, contrasts = NULL, ...) {

  rlang::check_installed("survival")

  method <- "breslow"

  cal <- match.call()

  arg::arg_supplied(formula)
  arg::arg_formula(formula, one_sided = FALSE)
  arg::arg_flag(model)
  arg::arg_flag(x)
  arg::arg_flag(y)

  if (...length() > 0L) {
    controlargs <- names(formals(survival::coxph.control))
    indx <- pmatch(...names(), controlargs, nomatch = 0L)

    if (any(indx == 0L)) {
      bad_args <- ...names()[indx == 0L]
      arg::err("argument{?s} {.arg {bad_args}} not matched")
    }
  }

  arg::when_not_null(control, arg::arg_list)

  if (rlang::is_missing(control)) {
    control <- survival::coxph.control(...)
  }
  else {
    control <- do.call(survival::coxph.control, control)
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
    arg::err("No (non-missing) observations")
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
  #The mask must come from the fit's weights rather than from the `weights`
  #argument because `.compute_vcov()` evaluates `psi` at `weights = s.weights`
  #(i.e., with the estimated weights set to 1) to form the cross-derivative block,
  #which would otherwise return those units to the risk sets and perturb the psi
  #of the units that do contribute.
  .w <- fit[["weights"]]

  pos <- {
    if (is_null(.w)) TRUE
    else .w > 0
  }

  psi <- function(B, X, y, weights, offset = 0) {
    weights <- weights * pos

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

    weights * M
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
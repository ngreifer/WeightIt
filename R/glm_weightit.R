#' Fitting (Weighted) Generalized Linear Models
#'
#' @description
#' `glm_weightit()` is used to fit generalized linear models with a
#' covariance matrix that accounts for estimation of weights, if supplied.
#' `lm_weightit()` is a wrapper for `glm_weightit()` with the Gaussian family
#' and identity link (i.e., a linear model). By
#' default, these functions use M-estimation to construct a robust covariance
#' matrix using the estimating equations for the weighting model and the outcome
#' model when available.
#'
#' @param formula an object of class [`formula`] (or one that can be coerced to
#'   that class): a symbolic description of the model to be fitted.
#' @param data a data frame containing the variables in the model. If not found
#'   in data, the variables are taken from `environment(formula)`, typically the
#'   environment from which the function is called.
#' @param family a description of the error distribution and link function to be
#'   used in the model. This can be a character string naming a family function,
#'   a family function or the result of a call to a family function. See
#'   [family] for details of family functions.
#' @param weightit a `weightit` or `weightitMSM` object; the output of a call to
#'   [weightit()] or [weightitMSM()]. If not supplied, an unweighted model will
#'   be fit.
#' @param vcov string; the method used to compute the variance of the estimated
#'   parameters. Allowable options include `"asympt"`, which uses the
#'   asymptotically correct M-estimation-based method that accounts for
#'   estimation of the weights when available; `"const"`, which uses the usual
#'   maximum likelihood estimates (only available when `weightit` is not
#'   supplied); `"HC0"`, which computes the robust sandwich variance treating
#'   weights (if supplied) as fixed; `"BS"`, which uses the traditional
#'   bootstrap (including re-estimation of the weights, if supplied); `"FWB"`,
#'   which uses the fractional weighted bootstrap as implemented in
#'   \pkgfun{fwb}{fwb} (including re-estimation of the weights, if supplied);
#'   and `"none"` to omit calculation of a variance matrix. If `NULL` (the
#'   default), will use `"asympt"` if `weightit` is supplied and M-estimation is
#'   available and `"HC0"` otherwise. See the `vcov_type` component of the
#'   outcome object to see which was used.
#' @param cluster optional; for computing a cluster-robust variance matrix, a
#'   variable indicating the clustering of observations, a list (or data frame)
#'   thereof, or a one-sided formula specifying which variable(s) from the
#'   fitted model should be used. Note the cluster-robust variance matrix uses a
#'   correction for small samples, as is done in `sandwich::vcovCL()` by
#'   default. Cluster-robust variance calculations are available only when
#'   `vcov` is `"asympt"`, `"HC0"`, `"BS"`, or `"FWB"`.
#' @param R the number of bootstrap replications when `vcov` is `"BS"` or
#'   `"FWB"`. Default is 500. Ignored otherwise.
#' @param offset optional; a numeric vector containing the model offset. See
#'   [offset()]. An offset can also be present in the model formula.
#' @param start optional starting values for the coefficients.
#' @param etastart,mustart optional starting values for the linear predictor and
#'   vector of means. Passed to [glm()].
#' @param control a list of parameters for controlling the fitting process.
#' @param x,y logical values indicating whether the response vector and model
#'   matrix used in the fitting process should be returned as components of the
#'   returned value.
#' @param contrasts an optional list defining contrasts for factor variables.
#'   See [model.matrix()].
#' @param fwb.args an optional list of further arguments to supply to
#'   \pkgfun{fwb}{fwb} when `vcov = "FWB"`.
#' @param br `logical`; whether to use bias-reduced regression as implemented by
#'   \pkgfun{brglm2}{brglmFit} (including Firth logistic regression). If `TRUE`, arguments passed to `control` or
#'   \dots will be passed to \pkgfun{brglm2}{brglmControl}.
#' @param \dots arguments to be used to form the default control argument if it
#'   is not supplied directly.
#'
#' @returns
#' A `glm_weightit` object, which inherits from `glm`.
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
#' [glm_weightit()] is essentially a wrapper for [glm()] that
#' optionally computes a coefficient variance matrix that can be adjusted to
#' account for estimation of the weights if a `weightit` or `weightitMSM` object
#' is supplied to the `weightit` argument.
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
#' Functions in the \pkg{sandwich} package can be to compute standard errors
#' after fitting, regardless of how `vcov` was specified, though these will
#' ignore estimation of the weights, if any. When no adjustment is done for
#' estimation of the weights (i.e., because no `weightit` argument was supplied
#' or there was no `"Mparts"` component in the supplied object), the default
#' variance matrix produced by `glm_weightit()` should align with that from
#' `sandwich::vcovHC(. type = "HC0")` or `sandwich::vcovCL(., type = "HC0", cluster = cluster)`
#' when `cluster` is supplied. Not all types are available
#' for all models.
#'
#' @seealso
#' * [lm()] and [glm()] for fitting (generalized) linear models without
#' adjusting standard errors for estimation of the weights.
#' * [glm_weightit()] for fitting generalized linear models that adjust for estimation of the weights.
#' * [multinom_weightit()] for fitting multinomial regression models that adjust for estimation of the weights.
#' * [ordinal_weightit()] for fitting ordinal regression models that adjust for estimation of the weights.
#' * [coxph_weightit()] for fitting Cox proportional hazards models that adjust for estimation of the weights.
#'
#' @examples
#' data("lalonde", package = "cobalt")
#'
#' # Logistic regression ATT weights
#' w.out <- weightit(treat ~ age + educ + married + re74,
#'                   data = lalonde, method = "glm",
#'                   estimand = "ATT")
#'
#' # Linear regression outcome model that adjusts
#' # for estimation of weights
#' fit1 <- lm_weightit(re78 ~ treat, data = lalonde,
#'                     weightit = w.out)
#'
#' summary(fit1)
#'
#' # Linear regression outcome model that treats weights
#' # as fixed
#' fit2 <- lm_weightit(re78 ~ treat, data = lalonde,
#'                     weightit = w.out,
#'                     vcov = "HC0")
#'
#' summary(fit2)
#'
#' # Can also just call summary() with `vcov` option
#' summary(fit1, vcov = "HC0")
#' @examplesIf rlang::is_installed("fwb")
#' # Linear regression outcome model that bootstraps
#' # estimation of weights and outcome model fitting
#' # using fractional weighted bootstrap with "Mammen"
#' # weights
#' set.seed(123)
#' fit3 <- lm_weightit(re78 ~ treat, data = lalonde,
#'                     weightit = w.out,
#'                     vcov = "FWB",
#'                     R = 50, #should use way more
#'                     fwb.args = list(wtype = "mammen"))
#'
#' summary(fit3)

#' @export
#' @name glm_weightit
glm_weightit <- function(formula, data, family = gaussian, weightit = NULL,
                         vcov = NULL, cluster, R = 500L,
                         offset, start = NULL, etastart, mustart,
                         control = list(...),
                         x = FALSE, y = TRUE,
                         contrasts = NULL, fwb.args = list(),
                         br = FALSE, ...) {

  vcov <- .process_vcov(vcov, weightit, R, fwb.args)

  if (missing(cluster)) {
    cluster <- NULL
  }

  model_call <- match.call()

  if (identical(family, "multinomial")) {
    arg::wrn('using {.fun glm_weightit}  with {.code family = "multinomial"} is deprecated. Please use {.fun multinom_weightit} instead')
    model_call[[1L]] <- quote(WeightIt::multinom_weightit)
    model_call[["family"]] <- NULL

    return(eval.parent(model_call))
  }

  ###
  arg::arg_flag(br)

  internal_model_call <- .build_internal_model_call(model = "glm",
                                                    model_call = model_call,
                                                    weightit = weightit,
                                                    vcov = vcov,
                                                    br = br)

  fit <- .eval_fit(internal_model_call,
                   warnings = c("non-integer" = NA),
                   errors = c("missing values in object" = "missing values are not allowed in the model variables"))

  if (br && vcov %in% c("asympt", "HC0") && identical(fit$type, "correction")) {
    arg::err('{.code type = "correction"} cannot be used with the specified {.arg vcov}')
  }

  fit$psi <- .get_glm_psi(fit)
  fit$br <- br

  ###

  fit$vcov <- .compute_vcov(fit, weightit, vcov, cluster, model_call, internal_model_call)

  fit <- .process_fit(fit, weightit, vcov, model_call, x, y)

  class(fit) <- c("glm_weightit", class(fit))

  fit
}

#' @export
#' @rdname glm_weightit
lm_weightit <- function(formula, data, weightit = NULL,
                        vcov = NULL, cluster, R = 500L,
                        offset,
                        x = FALSE, y = TRUE,
                        contrasts = NULL, fwb.args = list(),
                        ...) {

  vcov <- .process_vcov(vcov, weightit, R, fwb.args)

  if (missing(cluster)) {
    cluster <- NULL
  }

  model_call <- match.call()

  ###
  if (is_not_null(...get("family"))) {
    arg::err(c("{.arg family} cannot be used with {.fun lm_weightit}.",
               "i" = "Did you mean to use {.fun glm_weightit} instead?"))
  }

  internal_model_call <- .build_internal_model_call(model = "lm",
                                                    model_call = model_call,
                                                    weightit = weightit,
                                                    vcov = vcov)

  fit <- .eval_fit(internal_model_call,
                   errors = c("missing values in object" = "missing values are not allowed in the model variables"))

  fit$psi <- .get_glm_psi(fit)

  fit$family <- gaussian()

  ###

  fit$vcov <- .compute_vcov(fit, weightit, vcov, cluster, model_call, internal_model_call)

  fit <- .process_fit(fit, weightit, vcov, model_call, x, y)

  class(fit) <- c("glm_weightit", class(fit))

  fit
}

.get_hess_glm <- function(fit) {
  X <- fit[["x"]] %or% model.matrix(fit)

  d1mus <- fit$family$mu.eta(fit$linear.predictors)
  varmus <- fit$family$variance(fit$fitted.values)

  crossprod(X, X * (-d1mus^2 * fit$prior.weights / varmus))
}
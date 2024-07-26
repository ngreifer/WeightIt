#' Generate Balancing Weights with Minimal Input Processing
#'
#' @description
#' `weightit.fit()` dispatches one of the weight estimation methods
#' determined by `method`. It is an internal function called by
#' [weightit()] and should probably not be used except in special cases. Unlike
#' `weightit()`, `weightit.fit()` does not accept a formula and data
#' frame interface and instead requires the covariates and treatment to be
#' supplied as a numeric matrix and atomic vector, respectively. In this way,
#' `weightit.fit()` is to `weightit()` what [lm.fit()] is to [lm()] -
#' a thinner, slightly faster interface that performs minimal argument
#' checking.
#'
#' @param covs a numeric matrix of covariates.
#' @param treat a vector of treatment statuses.
#' @param method a string of length 1 containing the name of the method that
#' will be used to estimate weights. See [weightit()] for allowable options.
#' The default is `"glm"` for propensity score weighting using a
#' generalized linear model to estimate the propensity score.
#' @param s.weights a numeric vector of sampling weights. See the individual
#' pages for each method for information on whether sampling weights can be
#' supplied.
#' @param by.factor a factor variable for which weighting is to be done within
#' levels. Corresponds to the `by` argument in [weightit()].
#' @param estimand the desired estimand. For binary and multi-category
#' treatments, can be "ATE", "ATT", "ATC", and, for some methods, "ATO", "ATM",
#' or "ATOS". The default for both is "ATE". This argument is ignored for
#' continuous treatments. See the individual pages for each method for more
#' information on which estimands are allowed with each method and what
#' literature to read to interpret these estimands.
#' @param stabilize `logical`; whether or not to stabilize the weights.
#' For the methods that involve estimating propensity scores, this involves
#' multiplying each unit's weight by the proportion of units in their treatment
#' group. Default is `FALSE`.
#' @param focal when multi-category treatments are used and ATT weights are
#' requested, which group to consider the "treated" or focal group. This group
#' will not be weighted, and the other groups will be weighted to be more like
#' the focal group. Must be non-`NULL` if `estimand = "ATT"` or
#' `"ATC"`.
#' @param ps a vector of propensity scores. If specified, `method` will be
#' ignored and set to `"glm"`.
#' @param moments,int,subclass arguments to customize the weight estimation.
#' See [weightit()] for details.
#' @param is.MSM.method see [weightitMSM()]. Typically can be ignored.
#' @param missing `character`; how missing data should be handled. The
#' options depend on the `method` used. If `NULL`, `covs` covs
#' will be checked for `NA` values, and if present, `missing` will be
#' set to `"ind"`. If `""`, `covs` covs will not be checked for
#' `NA` values; this can be faster when it is known there are none.
#' @param verbose whether to print additional information output by the fitting
#' function.
#' @param include.obj whether to include in the output any fit objects created
#' in the process of estimating the weights. For example, with `method = "glm"`, the `glm` objects containing the propensity score model will be
#' included. See the individual pages for each method for information on what
#' object will be included if `TRUE`.
#' @param ... other arguments for functions called by `weightit.fit()`
#' that control aspects of fitting that are not covered by the above arguments.
#'
#' @returns
#' A `weightit.fit` object with the following elements:
#' \item{weights}{The estimated weights, one for each unit.}
#' \item{treat}{The values of the treatment variable.}
#' \item{estimand}{The estimand requested. ATC is recoded as ATT.}
#' \item{method}{The weight estimation method specified.}
#' \item{ps}{The estimated or provided propensity scores. Estimated propensity scores are returned for binary treatments and only when `method` is `"glm"`, `"gbm"`, `"cbps"`, `"super"`, or `"bart"`.}
#' \item{s.weights}{The provided sampling weights.}
#' \item{focal}{The focal treatment level if the ATT or ATC was requested.}
#' \item{fit.obj}{When `include.obj = TRUE`, the fit object.}
#' \item{info}{Additional information about the fitting. See the individual methods pages for what is included.}
#'
#' The `weightit.fit` object does not have specialized `print()`,
#' `summary()`, or `plot()` methods. It is simply a list containing
#' the above components. Use [as.weightit()] to convert it to a `weightit` object, which does have these methods. See Examples.
#'
#' @details
#' `weightit.fit()` is called by [weightit()] after the arguments to
#' `weightit()` have been checked and processed. `weightit.fit()`
#' dispatches the function used to actually estimate the weights, passing on
#' the supplied arguments directly. `weightit.fit()` is not meant to be
#' used by anyone other than experienced users who have a specific use case in
#' mind. The returned object contains limited information about the supplied
#' arguments or details of the estimation method; all that is processed by
#' `weightit()`.
#'
#' Less argument checking or processing occurs in `weightit.fit()` than
#' does in `weightit()`, which means supplying incorrect arguments can
#' result in errors, crashes, and invalid weights, and error and warning
#' messages may not be helpful in diagnosing the problem. `weightit.fit()`
#' does check to make sure weights were actually estimated, though.
#'
#' `weightit.fit()` may be most useful in speeding up simulation
#' simulation studies that use `weightit()` because the covariates can be
#' supplied as a numeric matrix, which is often how they are generated in
#' simulations, without having to go through the potentially slow process of
#' extracting the covariates and treatment from a formula and data frame. If
#' the user is certain the arguments are valid (e.g., by ensuring the estimated
#' weights are consistent with those estimated from `weightit()` with the
#' same arguments), less time needs to be spent on processing the arguments.
#' Also, the returned object is much smaller than a `weightit` object
#' because the covariates are not returned alongside the weights.
#'
#' @seealso
#' [weightit()], which you should use for estimating weights unless
#' you know better.
#'
#' [as.weightit()] for converting a `weightit.fit` object to a `weightit` object.
#'
#' @examples
#'
#' library("cobalt")
#' data("lalonde", package = "cobalt")
#'
#' # Balancing covariates between treatment groups (binary)
#' covs <- lalonde[c("age", "educ", "race", "married",
#'                   "nodegree", "re74", "re75")]
#' ## Create covs matrix, splitting any factors using
#' ## cobalt::splitfactor()
#' covs_mat <- as.matrix(splitfactor(covs))
#'
#' WF1 <- weightit.fit(covs_mat, treat = lalonde$treat,
#'                     method = "glm", estimand = "ATT")
#' str(WF1)
#'
#' # Converting to a weightit object for use with
#' # summary() and bal.tab()
#' W1 <- as.weightit(WF1, covs = covs)
#' W1
#' summary(W1)
#' bal.tab(W1)
#'

#' @export
weightit.fit <- function(covs, treat, method = "glm", s.weights = NULL, by.factor = NULL,
                         estimand = "ATE", focal = NULL, stabilize = FALSE,
                         ps = NULL, moments = NULL, int = FALSE,
                         subclass = NULL, is.MSM.method = FALSE, missing = NULL,
                         verbose = FALSE, include.obj = FALSE, ...) {

  #Checks
  if (!check_if_call_from_fun(weightit) && !check_if_call_from_fun(weightitMSM)) {

    chk::chk_not_missing(covs, "`covs`")
    chk::chk_matrix(covs)
    chk::chk_numeric(covs)

    chk::chk_not_missing(treat, "`treat`")
    chk::chk_vector(treat)
    chk::chk_numeric(treat)
    chk::chk_not_any_na(treat)

    chk::chk_flag(stabilize)
    chk::chk_flag(verbose)
    chk::chk_flag(include.obj)

    if (length(treat) != nrow(covs)) {
      .err("`treat` and `'covs` must contain the same number of units")
    }

    if (!has_treat_type(treat)) treat <- assign_treat_type(treat)
    treat.type <- get_treat_type(treat)

    .check_acceptable_method(method, msm = FALSE, force = FALSE)

    if (is_not_null(ps)) {
      chk::chk_vector(ps)
      chk::chk_numeric(ps)

      if (length(ps) != length(treat)) {
        .err("`ps` and `treat` must be the same length")
      }

      method <- "glm"
    }

    if (is_null(s.weights)) {
      s.weights <- rep.int(1, length(treat))
    }
    else {
      chk::chk_vector(s.weights)
      chk::chk_numeric(s.weights)

      if (length(s.weights) != length(treat)) {
        .err("`s.weights` and `treat` must be the same length")
      }
    }

    if (is_null(by.factor)) {
      by.factor <- factor(rep.int(1, length(treat)), levels = 1)
    }
    else {
      chk::chk_factor(by.factor)

      if (length(by.factor) != length(treat)) {
        .err("`by.factor` and `treat` must be the same length")
      }
    }

    #Process estimand and focal
    estimand <- .process_estimand(estimand, method, treat.type)
    f.e.r <- .process_focal_and_estimand(focal, estimand, treat)
    focal <- f.e.r[["focal"]]
    estimand <- f.e.r[["estimand"]]

    .chk_null_or(missing, chk::chk_string)

    if (is_null(missing)) {
      if (anyNA(covs)) missing <- "ind"
      else missing <- ""
    }
    else if (missing != "" && anyNA(covs)) {
      missing <- .process_missing(missing, method, treat.type)
    }
    else missing <- ""

    #Check subclass
    if (is_not_null(subclass)) .check_subclass(method, treat.type)

    #Process moments and int
    moments.int <- .process_moments_int(moments, int, method)
    moments <- moments.int[["moments"]]
    int <- moments.int[["int"]]
  }
  else {
    if (!has_treat_type(treat)) treat <- assign_treat_type(treat)
    treat.type <- get_treat_type(treat)
  }

  out <- make_list(c("weights", "treat", "estimand", "method", "ps", "s.weights",
                     "focal", "fit.obj", "info"))
  out$weights <- out$ps <- rep.int(NA_real_, length(treat))

  if (include.obj) {
    fit.obj <- make_list(levels(by.factor))
  }
  info <- make_list(levels(by.factor))

  obj <- NULL

  if (!is.function(method)) {
    fun <- "weightit"

    if (is.MSM.method) fun <- paste0(fun, "MSM")

    if (is_null(ps))
      fun <- paste0(fun, "2", method)
    else
      fun <- paste0(fun, "2ps")

    if (!is.MSM.method) {
      fun <- switch(treat.type,
                    "multinomial" = paste.(fun, "multi"),
                    "continuous" = paste.(fun, "cont"),
                    fun)
    }

    ####
    # if (FORMULA.TEST) fun <- paste0(".", fun)
    ####
    # fun <- get1(fun, mode = "function")
  }

  for (i in levels(by.factor)) {
    #Run method
    if (!is.function(method)) {
      obj <- do.call(fun,
                     alist(covs = covs,
                           treat = treat,
                           s.weights = s.weights,
                           subset = by.factor == i,
                           estimand = estimand,
                           focal = focal,
                           stabilize = stabilize,
                           subclass = subclass,
                           ps = ps,
                           moments = moments,
                           int = int,
                           missing = missing,
                           verbose = verbose,
                           ...))
    }
    else if (is.MSM.method) {
      obj <- weightitMSM2user(Fun = method,
                              covs.list = covs,
                              treat.list = treat,
                              s.weights = s.weights,
                              subset = by.factor == i,
                              stabilize = stabilize,
                              missing = missing,
                              verbose = verbose,
                              ...)
    }
    else {
      obj <- weightit2user(Fun = method,
                           covs = covs,
                           treat = treat,
                           s.weights = s.weights,
                           subset = by.factor == i,
                           estimand = estimand,
                           focal = focal,
                           stabilize = stabilize,
                           subclass = subclass,
                           ps = ps,
                           missing = missing,
                           moments = moments,
                           int = int,
                           verbose = verbose,
                           ...)
    }

    #Extract weights
    if (is_null(obj)) {
      .err("no object was created. This is probably a bug,\n     and you should report it at https://github.com/ngreifer/WeightIt/issues")
    }
    if (is_null(obj$w) || all(is.na(obj$w))) {
      .wrn("no weights were estimated. This is probably a bug,\n     and you should report it at https://github.com/ngreifer/WeightIt/issues")
    }
    if (any(!is.finite(obj$w))) {
      .wrn("some weights were estimated as `NA`, which means a value was impossible to compute (e.g., Inf). Check for extreme values of the treatment or covariates and try removing them. Non-finite weights will be set to 0")
      obj$w[!is.finite(obj$w)] <- 0
    }
    # else if (any(!is.finite(obj$w))) probably.a.bug()

    out$weights[by.factor == i] <- obj$w
    if (is_not_null(obj$ps)) {
      out$ps[by.factor == i] <- obj$ps
    }

    if (include.obj) {
      fit.obj[[i]] <- obj$fit.obj
    }

    info[[i]] <- obj$info
  }

  if (include.obj) {
    if (nlevels(by.factor) == 1) {
      fit.obj <- fit.obj[[1]]
    }
    out$fit.obj <- fit.obj
  }

  if (nlevels(by.factor) == 1) {
    attr(out, "Mparts") <- obj$Mparts
  }

  if (is_not_null(info) && nlevels(by.factor) == 1) {
    info <- info[[1]]
  }
  out$info <- info

  if (all(is.na(out$ps))) {
    out$ps <- NULL
  }

  out$treat <- treat
  out$estimand <- estimand
  out$method <- method
  out$s.weights <- s.weights
  out$focal <- focal

  class(out) <- "weightit.fit"

  out
}

weightitMSM.fit <- function(covs.list, treat.list, method = "glm", s.weights = NULL, by.factor = NULL,
                            estimand = "ATE", focal = NULL, stabilize = FALSE,
                            moments = NULL, int = FALSE,
                            subclass = NULL, is.MSM.method = FALSE, missing = NULL,
                            verbose = FALSE, include.obj = FALSE, ...) {

  #Checks
  if (!check_if_call_from_fun(weightitMSM)) {

    chk::chk_not_missing(covs.list, "`covs.list`")
    chk::chk_list(covs.list)
    for (i in seq_along(covs.list)) {
      chk::chk_matrix(covs.list[[i]], sprintf("`covs.list[[%s]]", i))
      chk::chk_numeric(covs.list[[i]], sprintf("`covs.list[[%s]]", i))
    }

    chk::chk_not_missing(treat.list, "`treat.list`")
    chk::chk_list(treat.list)
    for (i in seq_along(treat.list)) {
      chk::chk_vector(treat.list[[i]], sprintf("`treat.list[[%s]]", i))
      chk::chk_numeric(treat.list[[i]], sprintf("`treat.list[[%s]]", i))
      chk::chk_not_any_na(treat.list[[i]], sprintf("`treat.list[[%s]]", i))
    }

    chk::chk_flag(stabilize)
    chk::chk_flag(verbose)
    chk::chk_flag(include.obj)

    for (i in seq_along(treat.list)) {
      n <- length(treat.list[[i]])

      if (nrow(covs.list[[i]]) != n) {
        .err("treatment and covariates must have the same number of units")
      }
      if (anyNA(treat.list[[i]])) {
        .err(sprintf("no missing values are allowed in the treatment variable. Missing values found in treat.list[[%s]]", i))
      }

      treat.list[[i]] <- assign_treat_type(treat.list[[i]])

      if (anyNA(covs.list[[i]])) {
        missing <- .process_missing(missing, method, get_treat_type(treat.list[[i]]))
      }
      else if (i == length(covs.list)) {
        missing <- ""
      }
    }

    .check_acceptable_method(method, msm = TRUE, force = FALSE)

    if (is_null(s.weights)) {
      s.weights <- rep.int(1, n)
    }
    else {
      chk::chk_vector(s.weights)
      chk::chk_numeric(s.weights)

      if (length(s.weights) != n) {
        .err("`s.weights` and `treat` must be the same length")
      }
    }

    if (is_null(by.factor)) {
      by.factor <- factor(rep.int(1, n), levels = 1)
    }
    else {
      chk::chk_factor(by.factor)

      if (length(by.factor) != n) {
        .err("`by.factor` and `treat` must be the same length")
      }
    }

    #Process moments and int
    moments.int <- .process_moments_int(moments, int, method)
    moments <- moments.int[["moments"]]
    int <- moments.int[["int"]]
  }
  else {
    for (i in seq_along(treat.list)) {
      if (!has_treat_type(treat.list[[i]])) treat.list[[i]] <- assign_treat_type(treat.list[[i]])
    }
  }

  out <- make_list(c("weights", "treat.list", "method", "s.weights",
                     "fit.obj", "info"))
  out$weights <- rep.int(NA_real_, length(treat.list[[1]]))

  if (include.obj) {
    fit.obj <- make_list(levels(by.factor))
  }
  info <- make_list(levels(by.factor))

  obj <- NULL

  if (!is.function(method)) {
    fun <- paste0("weightitMSM2", method)
  }

  for (i in levels(by.factor)) {
    #Run method
    if (!is.function(method)) {
      obj <- do.call(fun,
                     alist(covs.list = covs.list,
                           treat.list = treat.list,
                           s.weights = s.weights,
                           subset = by.factor == i,
                           estimand = estimand,
                           focal = focal,
                           stabilize = stabilize,
                           subclass = subclass,
                           moments = moments,
                           int = int,
                           missing = missing,
                           verbose = verbose,
                           ...))
    }
    else {
      obj <- weightitMSM2user(Fun = method,
                              covs.list = covs.list,
                              treat.list = treat.list,
                              s.weights = s.weights,
                              subset = by.factor == i,
                              stabilize = stabilize,
                              missing = missing,
                              verbose = verbose,
                              ...)
    }

    #Extract weights
    if (is_null(obj)) {
      .err("no object was created. This is probably a bug,\n     and you should report it at https://github.com/ngreifer/WeightIt/issues")
    }
    if (is_null(obj$w) || all(is.na(obj$w))) {
      .wrn("no weights were estimated. This is probably a bug,\n     and you should report it at https://github.com/ngreifer/WeightIt/issues")
    }
    if (any(!is.finite(obj$w))) {
      .wrn("some weights were estimated as `NA`, which means a value was impossible to compute (e.g., Inf). Check for extreme values of the treatment or covariates and try removing them. Non-finite weights will be set to 0")
      obj$w[!is.finite(obj$w)] <- 0
    }
    # else if (any(!is.finite(obj$w))) probably.a.bug()

    out$weights[by.factor == i] <- obj$w

    if (include.obj) {
      fit.obj[i] <- list(obj$fit.obj)
    }

    info[i] <- list(obj$info)
  }

  if (include.obj) {
    out$fit.obj <- {
      if (nlevels(by.factor) == 1) fit.obj[[1]]
      else fit.obj
    }
  }

  if (nlevels(by.factor) == 1) {
    attr(out, "Mparts") <- obj$Mparts
  }

  if (is_not_null(info) && nlevels(by.factor) == 1) {
    info <- info[[1]]
  }
  out$info <- info


  out$treat.list <- treat.list
  out$method <- method
  out$s.weights <- s.weights

  class(out) <- "weightitMSM.fit"

  out
}

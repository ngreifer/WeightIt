#' Generate Balancing Weights with Minimal Input Processing
#'
#' @description
#' `weightit.fit()` dispatches one of the weight estimation methods
#' determined by `method`. It is an internal function called by [weightit()] and
#' should probably not be used except in special cases. Unlike `weightit()`,
#' `weightit.fit()` does not accept a formula and data frame interface and
#' instead requires the covariates and treatment to be supplied as a numeric
#' matrix and atomic vector, respectively. In this way, `weightit.fit()` is to
#' `weightit()` what [lm.fit()] is to [lm()]: a thinner, slightly faster
#' interface that performs minimal argument checking.
#'
#' @inheritParams weightit
#' @param covs a numeric matrix of covariates.
#' @param treat a vector of treatment statuses.
#' @param method a string containing the name of the method that will be used to
#'   estimate weights. See [weightit()] for allowable options. The default is
#'   `"glm"` for propensity score weighting using a generalized linear model to
#'   estimate the propensity score.
#' @param s.weights a numeric vector of sampling weights. See the individual
#'   pages for each method for information on whether sampling weights can be
#'   supplied.
#' @param by.factor a factor variable for which weighting is to be done within
#'   levels. Corresponds to the `by` argument in [weightit()].
#' @param stabilize `logical`; whether or not to stabilize the weights. For the
#'   methods that involve estimating propensity scores, this involves
#'   multiplying each unit's weight by the proportion of units in their
#'   treatment group. Default is `FALSE`. Note this differs from its use with
#'   [weightit()].
#' @param focal when `estimand` is set to `"ATT"` or `"ATC"`, which group to
#'   consider the "treated" or "control" group. This group will not be weighted,
#'   and the other groups will be weighted to resemble the focal group. If
#'   specified, `estimand` will automatically be set to `"ATT"` (with a warning
#'   if `estimand` is not `"ATT"` or `"ATC"`). See section *`estimand` and `focal`* in Details at [weightit()].
#' @param ps a vector of propensity scores. If specified, `method` will be
#'   ignored and set to `"glm"`.
#' @param missing `character`; how missing data should be handled. The options
#'   depend on the `method` used. If `NULL`, `covs` will be checked for `NA`
#'   values, and if present, `missing` will be set to `"ind"`. If `""`, `covs`
#'   will not be checked for `NA` values; this can be faster when it is known
#'   there are none.
#' @param ... other arguments for functions called by `weightit.fit()` that
#'   control aspects of fitting that are not covered by the above arguments.
#'
#' @returns
#' A `weightit.fit` object with the following elements:
#'
#' \item{weights}{The estimated weights, one for each unit.}
#' \item{treat}{The values of the treatment variable.}
#' \item{estimand}{The estimand requested.}
#' \item{method}{The weight estimation method specified.}
#' \item{ps}{The estimated or provided propensity scores. Estimated propensity scores are returned for binary treatments and only when `method` is `"glm"`, `"gbm"`, `"cbps"`, `"ipt"`, `"super"`, or `"bart"`. The propensity score corresponds to the predicted probability of being treated; see section *`estimand` and `focal`* in Details at [weightit()] for how the treated group is determined.}
#' \item{s.weights}{The provided sampling weights.}
#' \item{focal}{The focal treatment level if the ATT or ATC was requested.}
#' \item{fit.obj}{When `include.obj = TRUE`, the fit object.}
#' \item{info}{Additional information about the fitting. See the individual methods pages for what is included.}
#'
#' The `weightit.fit` object does not have specialized `print()`, `summary()`,
#' or `plot()` methods. It is simply a list containing the above components. Use
#' [as.weightit()] to convert it to a `weightit` object, which does have these
#' methods. See Examples.
#'
#' @details
#' `weightit.fit()` is called by [weightit()] after the arguments to
#' `weightit()` have been checked and processed. `weightit.fit()` dispatches the
#' function used to actually estimate the weights, passing on the supplied
#' arguments directly. `weightit.fit()` is not meant to be used by anyone other
#' than experienced users who have a specific use case in mind. The returned
#' object contains limited information about the supplied arguments or details
#' of the estimation method; all that is processed by `weightit()`.
#'
#' Less argument checking or processing occurs in `weightit.fit()` than does in
#' `weightit()`, which means supplying incorrect arguments can result in errors,
#' crashes, and invalid weights, and error and warning messages may not be
#' helpful in diagnosing the problem. `weightit.fit()` does check to make sure
#' weights were actually estimated, though.
#'
#' `weightit.fit()` may be most useful in speeding up simulation simulation
#' studies that use `weightit()` because the covariates can be supplied as a
#' numeric matrix, which is often how they are generated in simulations, without
#' having to go through the potentially slow process of extracting the
#' covariates and treatment from a formula and data frame. If the user is
#' certain the arguments are valid (e.g., by ensuring the estimated weights are
#' consistent with those estimated from `weightit()` with the same arguments),
#' less time needs to be spent on processing the arguments. Also, the returned
#' object is much smaller than a `weightit` object because the covariates are
#' not returned alongside the weights.
#'
#' @seealso
#' [weightit()], which you should use for estimating weights unless you
#' know better.
#'
#' [as.weightit()] for converting a `weightit.fit` object to a `weightit`
#' object.
#'
#' @examples
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

#' @export
weightit.fit <- function(covs, treat, method = "glm", s.weights = NULL, by.factor = NULL,
                         estimand = "ATE", focal = NULL, stabilize = FALSE,
                         ps = NULL, missing = NULL, verbose = FALSE, include.obj = FALSE, ...) {

  A <- list(...)

  #Checks
  if (check_if_call_from_fun(weightit) || check_if_call_from_fun(weightitMSM)) {
    treat <- as.treat(treat)
    treat.type <- get_treat_type(treat)
  }
  else {
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
      .err("`treat` and `covs` must contain the same number of units")
    }

    treat <- as.treat(treat, process = TRUE)
    treat.type <- get_treat_type(treat)

    .check_acceptable_method(method, msm = FALSE, force = FALSE)

    .check_method_treat.type(method, treat.type)

    if (is_not_null(ps)) {
      chk::chk_vector(ps)
      chk::chk_numeric(ps)

      if (length(ps) != length(treat)) {
        .err("`ps` and `treat` must be the same length")
      }

      if (is_not_null(ps) && is_not_null(method) &&
          !is.function(method) &&
          !identical(method, "glm")) {
        .wrn("`ps` is supplied, so `method` will be ignored")
      }

      method <- NULL
    }

    if (is_null(s.weights)) {
      s.weights <- rep_with(1, treat)
    }
    else {
      .chk_basic_vector(s.weights)
      chk::chk_numeric(s.weights)

      .check_method_s.weights(method, s.weights)

      if (length(s.weights) != length(treat)) {
        .err("`s.weights` and `treat` must be the same length")
      }
    }

    if (is_null(by.factor)) {
      by.factor <- gl(1L, length(treat))
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

    if (treat.type == "binary") {
      attr(treat, "treated") <- f.e.r[["treated"]]
    }

    chk::chk_null_or(missing, vld = chk::vld_string)

    if (is_null(missing)) {
      if (is_not_null(method) && anyNA(covs)) missing <- "ind"
      else missing <- ""
    }
    else if (nzchar(missing) && anyNA(covs)) {
      missing <- .process_missing(missing, method)
    }
    else {
      missing <- ""
    }

    #Check subclass
    .check_subclass(...get("subclass"), method, treat.type)

    #Process moments and int
    m.i.q <- .process_moments_int_quantile(method = method, ...)

    A[c("moments", "int", "quantile")] <- m.i.q[c("moments", "int", "quantile")]
  }

  missing <- .process_missing2(missing, covs)

  out <- make_list(c("weights", "treat", "estimand", "method", "ps", "s.weights",
                     "focal", "missing", "fit.obj", "info"))
  out$weights <- out$ps <- rep_with(NA_real_, treat)

  if (include.obj) {
    fit.obj <- make_list(levels(by.factor))
  }

  info <- make_list(levels(by.factor))

  obj <- NULL

  .check_required_packages(method)

  if (is.function(method)) {
    fun <- "weightit2user"
  }
  else {
    fun <- {
      if (is_not_null(ps)) "weightit2ps"
      else if (is_null(method)) "weightit2null"
      else sprintf("weightit2%s", method)
    }

    fun <- switch(treat.type,
                  `multi-category` =,
                  multinomial = paste.(fun, "multi"),
                  continuous = paste.(fun, "cont"),
                  fun)
  }

  A["covs"] <- list(covs)
  A["treat"] <- list(treat)
  A["s.weights"] <- list(s.weights)
  A["estimand"] <- list(estimand)
  A["focal"] <- list(focal)
  A["stabilize"] <- list(stabilize)
  A["ps"] <- list(ps)
  A["missing"] <- list(missing)
  A["verbose"] <- list(verbose)

  for (i in levels(by.factor)) {
    A["subset"] <- list(by.factor == i)

    #Run method
    if (is.function(method)) {
      A["Fun"] <- list(method)
    }

    obj <- do.call(fun, A)

    #Extract weights
    if (is_null(obj)) {
      .err("no object was created. This is probably a bug, and you should report it at https://github.com/ngreifer/WeightIt/issues")
    }

    if (is_null(obj$w) || all(is.na(obj$w))) {
      .wrn("no weights were estimated. This is probably a bug, and you should report it at https://github.com/ngreifer/WeightIt/issues")
    }

    if (!all(is.finite(obj$w))) {
      .wrn("some weights were estimated as `NA`, which means a value was impossible to compute (e.g., Inf). Check for extreme values of the treatment or covariates and try removing them. Non-finite weights will be set to 0")
      obj$w[!is.finite(obj$w)] <- 0
    }
    # else if (any(!is.finite(obj$w))) probably_a_bug()

    out$weights[by.factor == i] <- obj$w
    if (is_not_null(obj$ps)) {
      out$ps[by.factor == i] <- obj$ps
    }

    if (include.obj && is_not_null(obj$fit.obj)) {
      fit.obj[[i]] <- obj$fit.obj
    }

    info[[i]] <- obj$info
  }

  if (include.obj) {
    if (any(lengths(fit.obj) > 0L)) {
      if (nlevels(by.factor) == 1L) {
        fit.obj <- fit.obj[[1L]]
      }
      out$fit.obj <- fit.obj
    }
    else {
      .wrn("`include.obj` was set to `TRUE`, but no fit object was returned by the fitting function")
    }
  }

  if (nlevels(by.factor) == 1L) {
    attr(out, "Mparts") <- obj$Mparts
  }

  if (is_not_null(info) && nlevels(by.factor) == 1L) {
    info <- info[[1L]]
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
  out$missing <- missing

  class(out) <- "weightit.fit"

  out
}

weightitMSM.fit <- function(covs.list, treat.list, method = "glm", s.weights = NULL, by.factor = NULL,
                            estimand = "ATE", focal = NULL, stabilize = FALSE,
                            is.MSM.method = FALSE, missing = NULL,
                            verbose = FALSE, include.obj = FALSE, ...) {

  A <- list(...)

  #Checks
  if (check_if_call_from_fun(weightitMSM)) {
    for (i in seq_along(treat.list)) {
      if (!has_treat_type(treat.list[[i]])) {
        treat.list[[i]] <- assign_treat_type(treat.list[[i]])
      }
    }
  }
  else {
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

      if (i > 1L && n != length(treat.list[[1L]])) {
        .err("the number of units must be the same for each time point")
      }

      if (nrow(covs.list[[i]]) != n) {
        .err("treatment and covariates must have the same number of units")
      }

      if (anyNA(treat.list[[i]])) {
        .err(sprintf("no missing values are allowed in the treatment variable. Missing values found in treat.list[[%s]]", i))
      }

      treat.list[[i]] <- assign_treat_type(treat.list[[i]])
    }

    .check_acceptable_method(method, msm = TRUE, force = FALSE)

    missing <- {
      if (is_null(method) || !any_apply(covs.list, anyNA)) ""
      else .process_missing(missing, method)
    }

    if (is_null(s.weights)) {
      s.weights <- rep.int(1, n)
    }
    else {
      chk::chk_vector(s.weights)
      chk::chk_numeric(s.weights)

      .check_method_s.weights(method, s.weights)

      if (length(s.weights) != n) {
        .err("`s.weights` must have length equal to the number of units")
      }
    }

    if (is_null(by.factor)) {
      by.factor <- gl(1L, n)
    }
    else {
      chk::chk_factor(by.factor)

      if (length(by.factor) != n) {
        .err("`by.factor` must have length equal to the number of units")
      }
    }

    #Process moments and int
    m.i.q <- .process_moments_int_quantile(method = method, ...)

    A[c("moments", "int", "quantile")] <- m.i.q[c("moments", "int", "quantile")]
  }

  out <- make_list(c("weights", "treat.list", "method", "s.weights", "missing",
                     "fit.obj", "info"))
  out$weights <- rep_with(NA_real_, treat.list[[1L]])

  if (include.obj) {
    fit.obj <- make_list(levels(by.factor))
  }
  info <- make_list(levels(by.factor))

  obj <- NULL

  .check_required_packages(method)

  fun <- {
    if (is_null(method)) "weightitMSM2null"
    else if (is.function(method)) "weightitMSM2user"
    else sprintf("weightitMSM2%s", method)
  }

  A["covs.list"] <- list(covs.list)
  A["treat.list"] <- list(treat.list)
  A["s.weights"] <- list(s.weights)
  A["estimand"] <- list(estimand)
  A["focal"] <- list(focal)
  A["stabilize"] <- list(stabilize)
  A["missing"] <- list(missing)
  A["verbose"] <- list(verbose)

  for (i in levels(by.factor)) {
    A["subset"] <- list(by.factor == i)

    #Run method
    if (is.function(method)) {
      A["Fun"] <- list(method)
    }

    obj <- do.call(fun, A)

    #Extract weights
    if (is_null(obj)) {
      .err("no object was created. This is probably a bug, and you should report it at https://github.com/ngreifer/WeightIt/issues")
    }

    if (is_null(obj$w) || all(is.na(obj$w))) {
      .wrn("no weights were estimated. This is probably a bug, and you should report it at https://github.com/ngreifer/WeightIt/issues")
    }

    if (!all(is.finite(obj$w))) {
      .wrn("some weights were estimated as `NA`, which means a value was impossible to compute (e.g., Inf). Check for extreme values of the treatment or covariates and try removing them. Non-finite weights will be set to 0")
      obj$w[!is.finite(obj$w)] <- 0
    }
    # else if (any(!is.finite(obj$w))) probably_a_bug()

    out$weights[by.factor == i] <- obj$w

    if (include.obj && is_not_null(obj$fit.obj)) {
      fit.obj[[i]] <- obj$fit.obj
    }

    info[i] <- list(obj$info)
  }

  if (include.obj) {
    if (any(lengths(fit.obj) > 0L)) {
      if (nlevels(by.factor) == 1L) {
        fit.obj <- fit.obj[[1L]]
      }
      out$fit.obj <- fit.obj
    }
    else {
      .wrn("`include.obj` was set to `TRUE`, but no fit object was returned by the fitting function")
    }
  }

  if (nlevels(by.factor) == 1L) {
    attr(out, "Mparts") <- obj$Mparts
  }

  if (is_not_null(info) && nlevels(by.factor) == 1L) {
    info <- info[[1L]]
  }
  out$info <- info

  out$treat.list <- treat.list
  out$method <- method
  out$s.weights <- s.weights
  out$missing <- missing

  class(out) <- "weightitMSM.fit"

  out
}

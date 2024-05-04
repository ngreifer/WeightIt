#' Estimate Balancing Weights
#'
#' @description
#' `weightit()` allows for the easy generation of balancing weights using
#' a variety of available methods for binary, continuous, and multi-category
#' treatments. Many of these methods exist in other packages, which
#' `weightit()` calls; these packages must be installed to use the desired
#' method.
#'
#' @param formula a formula with a treatment variable on the left hand side and
#' the covariates to be balanced on the right hand side. See [glm()] for more
#' details. Interactions and functions of covariates are allowed.
#' @param data an optional data set in the form of a data frame that contains
#' the variables in `formula`.
#' @param method a string of length 1 containing the name of the method that
#' will be used to estimate weights. See Details below for allowable options.
#' The default is `"glm"` for propensity score weighting using a
#' generalized linear model to estimate the propensity score.
#' @param estimand the desired estimand. For binary and multi-category
#' treatments, can be `"ATE"`, `"ATT"`, `"ATC"`, and, for some
#' methods, `"ATO"`, `"ATM"`, or `"ATOS"`. The default for both
#' is `"ATE"`. This argument is ignored for continuous treatments. See the
#' individual pages for each method for more information on which estimands are
#' allowed with each method and what literature to read to interpret these
#' estimands.
#' @param stabilize whether or not and how to stabilize the weights. If `TRUE`, each unit's weight will be multiplied by a standardization factor, which is the inverse of the unconditional probability (or density) of each unit's observed treatment value. If a formula, a generalized linear model will be fit with the included predictors, and the inverse of the corresponding weight will be used as the standardization factor. Can only be used with continuous treatments or when `estimand = "ATE"`. Default is `FALSE` for no standardization. See also the `num.formula` argument at [weightitMSM()]
#' @param focal when multi-category treatments are used and ATT weights are
#' requested, which group to consider the "treated" or focal group. This group
#' will not be weighted, and the other groups will be weighted to be more like
#' the focal group. If specified, `estimand` will automatically be set to
#' `"ATT"`.
#' @param by a string containing the name of the variable in `data` for
#' which weighting is to be done within categories or a one-sided formula with
#' the stratifying variable on the right-hand side. For example, if `by = "gender"` or `by = ~gender`, a separate propensity score model or optimization will occur within each level of the variable `"gender"`. Only one
#' `by` variable is allowed; to stratify by multiply variables
#' simultaneously, create a new variable that is a full cross of those
#' variables using [interaction()].
#' @param s.weights A vector of sampling weights or the name of a variable in
#' `data` that contains sampling weights. These can also be matching
#' weights if weighting is to be used on matched data. See the individual pages
#' for each method for information on whether sampling weights can be supplied.
#' @param ps A vector of propensity scores or the name of a variable in
#' `data` containing propensity scores. If not `NULL`, `method`
#' is ignored unless it is a user-supplied function, and the propensity scores will be used to create weights.
#' `formula` must include the treatment variable in `data`, but the
#' listed covariates will play no role in the weight estimation. Using
#' `ps` is similar to calling [get_w_from_ps()] directly, but produces a
#' full `weightit` object rather than just producing weights.
#' @param moments `numeric`; for some methods, the greatest power of each
#' covariate to be balanced. For example, if `moments = 3`, for each
#' non-categorical covariate, the covariate, its square, and its cube will be
#' balanced. This argument is ignored for other methods; to balance powers of
#' the covariates, appropriate functions must be entered in `formula`. See
#' the individual pages for each method for information on whether they accept
#' `moments`.
#' @param int `logical`; for some methods, whether first-order
#' interactions of the covariates are to be balanced. This argument is ignored
#' for other methods; to balance interactions between the variables,
#' appropriate functions must be entered in `formula`. See the individual
#' pages for each method for information on whether they accept `int`.
#' @param subclass `numeric`; the number of subclasses to use for
#' computing weights using marginal mean weighting with subclasses (MMWS). If
#' `NULL`, standard inverse probability weights (and their extensions)
#' will be computed; if a number greater than 1, subclasses will be formed and
#' weights will be computed based on subclass membership. Attempting to set a
#' non-`NULL` value for methods that don't compute a propensity score will
#' result in an error; see each method's help page for information on whether
#' MMWS weights are compatible with the method. See [get_w_from_ps()] for
#' details and references.
#' @param missing `character`; how missing data should be handled. The
#' options and defaults depend on the `method` used. Ignored if no missing
#' data is present. It should be noted that multiple imputation outperforms all
#' available missingness methods available in `weightit()` and should
#' probably be used instead. Consider the \CRANpkg{MatchThem}
#' package for the use of `weightit()` with multiply imputed data.
#' @param verbose `logical`; whether to print additional information
#' output by the fitting function.
#' @param include.obj `logical`; whether to include in the output any fit
#' objects created in the process of estimating the weights. For example, with
#' `method = "glm"`, the `glm` objects containing the propensity
#' score model will be included. See the individual pages for each method for
#' information on what object will be included if `TRUE`.
#' @param keep.mparts `logical`; whether to include in the output components necessary to estimate standard errors that account for estimation of the weights in [glm_weightit()]. Default is `TRUE` if such parts are present. See the individual pages for each method for whether these components are produced. Set to `FALSE` to keep the output object smaller, e.g., if standard errors will not be computed using `glm_weightit()`.
#' @param ... other arguments for functions called by `weightit()` that
#' control aspects of fitting that are not covered by the above arguments. See Details.
#'
#' @returns
#' A `weightit` object with the following elements:
#' \item{weights}{The estimated weights, one for each unit.}
#' \item{treat}{The values of the treatment variable.}
#' \item{covs}{The covariates used in the fitting. Only includes the raw covariates, which may have been altered in
#' the fitting process.}
#' \item{estimand}{The estimand requested.}
#' \item{method}{The weight estimation method specified.}
#' \item{ps}{The estimated or provided propensity scores. Estimated propensity scores are
#' returned for binary treatments and only when `method` is `"glm"`, `"gbm"`, `"cbps"`, `"ipt"`, `"super"`, or `"bart"`.}
#' \item{s.weights}{The provided sampling weights.}
#' \item{focal}{The focal treatment level if the ATT or ATC was requested.}
#' \item{by}{A data.frame containing the `by` variable when specified.}
#' \item{obj}{When `include.obj = TRUE`, the fit object.}
#' \item{info}{Additional information about the fitting. See the individual
#' methods pages for what is included.}
#'
#' When `keep.mparts` is `TRUE` (the default) and the chosen method is compatible with M-estimation, the components related to M-estimation for use in [glm_weightit()] are stored in the `"Mparts"` attribute. When `by` is specified, `keep.mparts` is set to `FALSE`.
#'
#' @details
#' The primary purpose of `weightit()` is as a dispatcher to functions
#' that perform the estimation of balancing weights using the requested
#' `method`. Below are the methods allowed and links to pages containing
#' more information about them, including additional arguments and outputs
#' (e.g., when `include.obj = TRUE`), how missing values are treated,
#' which estimands are allowed, and whether sampling weights are allowed.
#'
#' | [`"glm"`][method_glm] | Propensity score weighting using generalized linear models |
#' | :---- | :---- |
#' | [`"gbm"`][method_gbm] | Propensity score weighting using generalized boosted modeling |
#' | [`"cbps"`][method_cbps]| Covariate Balancing Propensity Score weighting |
#' | [`"npcbps"`][method_npcbps]| Non-parametric Covariate Balancing Propensity Score weighting |
#' | [`"ebal"`][method_ebal] | Entropy balancing |
#' | [`"ipt"`][method_ipt] | Inverse probability tilting |
#' | [`"optweight"`][method_optweight] | Optimization-based weighting |
#' | [`"super"`][method_super] | Propensity score weighting using SuperLearner |
#' | [`"bart"`][method_bart] | Propensity score weighting using Bayesian additive regression trees (BART) |
#' | [`"energy"`][method_energy] | Energy balancing |
#'
#' `method` can also be supplied as a user-defined function; see
#' [`method_user`] for instructions and examples.
#'
#' When using `weightit()`, please cite both the \pkg{WeightIt} package
#' (using `citation("WeightIt")`) and the paper(s) in the references
#' section of the method used.
#'
#' @seealso
#' [weightitMSM()] for estimating weights with sequential (i.e.,
#' longitudinal) treatments for use in estimating marginal structural models
#' (MSMs).
#'
#' [weightit.fit()], which is a lower-level dispatcher function that accepts a
#' matrix of covariates and a vector of treatment statuses rather than a
#' formula and data frame and performs minimal argument checking and
#' processing. It may be useful for speeding up simulation studies for which
#' the correct arguments are known. In general `weightit()` should be
#' used.
#'
#' [summary.weightit()] for summarizing the weights
#'
#' @examples
#'
#' library("cobalt")
#' data("lalonde", package = "cobalt")
#'
#' #Balancing covariates between treatment groups (binary)
#' (W1 <- weightit(treat ~ age + educ + married +
#'                   nodegree + re74, data = lalonde,
#'                 method = "glm", estimand = "ATT"))
#' summary(W1)
#' bal.tab(W1)
#'
#' #Balancing covariates with respect to race (multi-category)
#' (W2 <- weightit(race ~ age + educ + married +
#'                   nodegree + re74, data = lalonde,
#'                 method = "ebal", estimand = "ATE"))
#' summary(W2)
#' bal.tab(W2)
#'
#' #Balancing covariates with respect to re75 (continuous)
#' (W3 <- weightit(re75 ~ age + educ + married +
#'                   nodegree + re74, data = lalonde,
#'                 method = "cbps"))
#' summary(W3)
#' bal.tab(W3)

#' @export
weightit <- function(formula, data = NULL, method = "glm", estimand = "ATE", stabilize = FALSE, focal = NULL,
                     by = NULL, s.weights = NULL, ps = NULL, moments = NULL, int = FALSE, subclass = NULL,
                     missing = NULL, verbose = FALSE, include.obj = FALSE, keep.mparts = TRUE, ...) {

  ## Checks and processing ----

  A <- list(...)

  #Checks
  if (is_null(formula) || !rlang::is_formula(formula, lhs = TRUE)) {
    .err("`formula` must be a formula relating treatment to covariates")
  }

  #Process treat and covs from formula and data
  t.c <- get_covs_and_treat_from_formula(formula, data)
  reported.covs <- t.c[["reported.covs"]]
  covs <- t.c[["model.covs"]]
  treat <- t.c[["treat"]]
  # treat.name <- t.c[["treat.name"]]

  if (is_null(treat)) {
    .err("no treatment variable was specified")
  }
  if (length(treat) != nrow(covs)) {
    .err("the treatment and covariates must have the same number of units")
  }

  n <- length(treat)

  if (anyNA(treat)) {
    .err("no missing values are allowed in the treatment variable")
  }

  #Get treat type
  treat <- assign_treat_type(treat)
  treat.type <- get_treat_type(treat)

  chk::chk_flag(verbose)
  chk::chk_flag(include.obj)
  chk::chk_flag(keep.mparts)

  #Process ps
  ps <- .process_ps(ps, data, treat)
  if (is_not_null(ps) && !identical(method, "glm") && !is.function(method)) {
    .wrn("`ps` is supplied, so `method` will be ignored")
    method <- "glm"
  }

  ##Process method
  .check_acceptable_method(method, msm = FALSE, force = FALSE)

  if (is.character(method)) {
    method <- .method_to_proper_method(method)
    attr(method, "name") <- method
  }
  else if (is.function(method)) {
    method.name <- deparse1(substitute(method))
    .check_user_method(method)
    attr(method, "name") <- method.name
  }

  #Process estimand and focal
  estimand <- .process_estimand(estimand, method, treat.type)
  f.e.r <- .process_focal_and_estimand(focal, estimand, treat)
  focal <- f.e.r[["focal"]]
  # estimand <- f.e.r[["estimand"]]
  estimand <- f.e.r[["reported.estimand"]]
  reported.estimand <- f.e.r[["reported.estimand"]]

  #Process missing
  missing <- {
    if (!anyNA(reported.covs)) ""
    else .process_missing(missing, method, treat.type)
  }

  #Check subclass
  if (is_not_null(subclass)) .check_subclass(method, treat.type)

  #Process s.weights
  s.weights <- process.s.weights(s.weights, data)

  if (is_null(s.weights)) s.weights <- rep(1, n)

  #Process stabilize
  if (isFALSE(stabilize)) {
    stabilize <- NULL
  }
  else if (isTRUE(stabilize)) {
    stabilize <- as.formula("~ 1")
  }

  if (is_not_null(stabilize)) {
    if (treat.type != "continious" && estimand != "ATE") {
      .err("`stabilize` can only be supplied when `estimand = \"ATE\"`")
    }

    if (!rlang::is_formula(stabilize)) {
      .err("`stabilize` must be `TRUE`, `FALSE`, or a formula with the stabilization factors on the right hand side")
    }

    stabilize <- update(formula, update(stabilize, NULL ~ .))
  }

  ##Process by
  if (is_not_null(A[["exact"]])) {
    .msg("`by` has replaced `exact` in the `weightit()` syntax, but `exact` will always work")
    by <- A[["exact"]]
    by.arg <- "exact"
  }
  else by.arg <- "by"

  # processed.by <- .process_by(by.name, data = data, treat = treat)
  processed.by <- .process_by(by, data = data, treat = treat, by.arg = by.arg)

  #Process moments and int
  moments.int <- .process_moments_int(moments, int, method)

  call <- match.call()
  # args <- list(...)

  ## Running models ----

  #Returns weights (weights) and propensity score (ps)
  A[["treat"]] <- treat
  A[["covs"]] <- covs
  A[["s.weights"]] <- s.weights
  A[["by.factor"]] <- attr(processed.by, "by.factor")
  A[["estimand"]] <- estimand
  A[["focal"]] <- focal
  A[["stabilize"]] <- FALSE
  A[["method"]] <- method
  A[["moments"]] <- moments.int[["moments"]]
  A[["int"]] <- moments.int[["int"]]
  A[["subclass"]] <- subclass
  A[["ps"]] <- ps
  A[["missing"]] <- missing
  A[["verbose"]] <- verbose
  A[["include.obj"]] <- include.obj
  A[[".data"]] <- data
  A[[".formula"]] <- formula
  A[[".covs"]] <- reported.covs

  obj <- do.call("weightit.fit", A)

  if (is_not_null(stabilize)) {
    stab.t.c <- get_covs_and_treat_from_formula(stabilize, data)

    A[["treat"]] <- stab.t.c[["treat"]]
    A[["covs"]] <- stab.t.c[["model.covs"]]
    A[["method"]] <- "glm"
    A[["moments"]] <- NULL
    A[["int"]] <- NULL
    A[[".formula"]] <- stabilize
    A[[".covs"]] <- stab.t.c[["reported.covs"]]

    sw_obj <- do.call("weightit.fit", A)

    obj$weights <- obj$weights / sw_obj[["weights"]]
  }

  .check_estimated_weights(obj$weights, treat, treat.type, s.weights)

  ## Assemble output object----
  out <- list(weights = obj$weights,
              treat = treat,
              covs = reported.covs,
              estimand = if (treat.type == "continuous") NULL else reported.estimand,
              method = method,
              ps = if (is_null(obj$ps) || all(is.na(obj$ps))) NULL else obj$ps,
              s.weights = s.weights,
              #discarded = NULL,
              focal = if (reported.estimand %in% c("ATT", "ATC")) focal else NULL,
              by = processed.by,
              call = call,
              formula = formula,
              stabilize = stabilize,
              env = parent.frame(),
              info = obj$info,
              obj = obj$fit.obj)

  out <- clear_null(out)

  if (keep.mparts && is_not_null(attr(obj, "Mparts"))) {
    if (is_null(stabilize)) {
      attr(out, "Mparts") <- attr(obj, "Mparts")
    }
    else {
      attr(sw_obj, "Mparts")$wfun <- Invert(attr(sw_obj, "Mparts")$wfun)
      attr(out, "Mparts.list") <- list(attr(obj, "Mparts"), attr(sw_obj, "Mparts"))
    }
  }

  class(out) <- "weightit"

  out
}

#' @exportS3Method print weightit
print.weightit <- function(x, ...) {
  treat.type <- get_treat_type(x[["treat"]])
  trim <- attr(x[["weights"]], "trim")

  cat("A " %+% italic("weightit") %+% " object\n")

  if (is_not_null(x[["method"]])) {
    cat(sprintf(" - method: %s (%s)\n",
                add_quotes(attr(x[["method"]], "name")),
                .method_to_phrase(x[["method"]])))
  }

  cat(sprintf(" - number of obs.: %s\n",
              length(x[["weights"]])))

  cat(sprintf(" - sampling weights: %s\n",
              if (is_null(x[["s.weights"]]) || all_the_same(x[["s.weights"]])) "none" else "present"))

  cat(sprintf(" - treatment: %s\n",
              if (treat.type == "continuous") "continuous"
              else paste0(nunique(x[["treat"]]), "-category",
                          if (treat.type == "multinomial") paste0(" (", paste(levels(x[["treat"]]), collapse = ", "), ")")
                          else "")))

  if (is_not_null(x[["estimand"]])) {
    cat(paste0(" - estimand: ", x[["estimand"]],
               if (is_not_null(x[["focal"]])) sprintf(" (focal: %s)", x[["focal"]]) else "", "\n"))
  }

  if (is_not_null(x[["covs"]])) {
    cat(paste0(" - covariates: ", ifelse(length(names(x[["covs"]])) > 60, "too many to name", paste(names(x[["covs"]]), collapse = ", ")), "\n"))
  }

  if (is_not_null(x[["by"]])) {
    cat(sprintf(" - by: %s\n", paste(names(x[["by"]]), collapse = ", ")))
  }

  if (is_not_null(trim)) {
    if (trim < 1) {
      if (attr(x[["weights"]], "trim.lower")) trim <- c(1 - trim, trim)
      cat(paste(" - weights trimmed at", word_list(paste0(round(100*trim, 2), "%")), "\n"))
    }
    else {
      if (attr(x[["weights"]], "trim.lower")) t.b <- "top and bottom" else t.b <- "top"
      cat(paste(" - weights trimmed at the", t.b, trim, "\n"))
    }
  }
  invisible(x)
}

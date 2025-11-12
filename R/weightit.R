#' Estimate Balancing Weights
#'
#' @description
#' `weightit()` allows for the easy generation of balancing weights
#' using a variety of available methods for binary, continuous, and
#' multi-category treatments. Many of these methods exist in other packages,
#' which `weightit()` calls; these packages must be installed to use the desired
#' method.
#'
#' @param formula a formula with a treatment variable on the left hand side and
#'   the covariates to be balanced on the right hand side. See [glm()] for more
#'   details. Interactions and functions of covariates are allowed.
#' @param data an optional data set in the form of a data frame that contains
#'   the variables in `formula`.
#' @param method a string of length 1 containing the name of the method that
#'   will be used to estimate weights. See Details below for allowable options.
#'   The default is `"glm"` for propensity score weighting using a generalized
#'   linear model to estimate the propensity score.
#' @param estimand the desired estimand. For binary and multi-category
#'   treatments, can be `"ATE"`, `"ATT"`, `"ATC"`, and, for some methods,
#'   `"ATO"`, `"ATM"`, or `"ATOS"`. The default for both is `"ATE"`. This
#'   argument is ignored for continuous treatments. See the individual pages for
#'   each method for more information on which estimands are allowed with each
#'   method and what literature to read to interpret these estimands.
#' @param stabilize whether or not and how to stabilize the weights. If `TRUE`,
#'   each unit's weight will be multiplied by a standardization factor, which is
#'   the the unconditional probability (or density) of each unit's observed
#'   treatment value. If a formula, a generalized linear model will be fit with
#'   the included predictors, and the inverse of the corresponding weight will
#'   be used as the standardization factor. Can only be used with continuous
#'   treatments or when `estimand = "ATE"`. Default is `FALSE` for no
#'   standardization. See also the `num.formula` argument at [weightitMSM()].
#'   For continuous treatments, weights are already stabilized, so setting
#'   `stabilize = TRUE` will be ignored with a warning (supplying a formula
#'   still works).
#' @param focal when `estimand` is set to `"ATT"` or `"ATC"`, which group to
#'   consider the "treated" or "control" group. This group will not be weighted,
#'   and the other groups will be weighted to resemble the focal group. If
#'   specified, `estimand` will automatically be set to `"ATT"` (with a warning
#'   if `estimand` is not `"ATT"` or `"ATC"`). See section *`estimand` and `focal`* in Details below.
#' @param by a string containing the name of the variable in `data` for which
#'   weighting is to be done within categories or a one-sided formula with the
#'   stratifying variable on the right-hand side. For example, if `by = "gender"` or `by = ~gender`,
#'   a separate propensity score model or
#'   optimization will occur within each level of the variable `"gender"`. Only
#'   one `by` variable is allowed; to stratify by multiply variables
#'   simultaneously, create a new variable that is a full cross of those
#'   variables using [interaction()].
#' @param s.weights a vector of sampling weights or the name of a variable in
#'   `data` that contains sampling weights. These can also be matching weights
#'   if weighting is to be used on matched data. See the individual pages for
#'   each method for information on whether sampling weights can be supplied.
#' @param ps a vector of propensity scores or the name of a variable in `data`
#'   containing propensity scores. If not `NULL`, `method` is ignored unless it
#'   is a user-supplied function, and the propensity scores will be used to
#'   create weights. `formula` must include the treatment variable in `data`,
#'   but the listed covariates will play no role in the weight estimation. Using
#'   `ps` is similar to calling [get_w_from_ps()] directly, but produces a full
#'   `weightit` object rather than just producing weights.
#' @param missing `character`; how missing data should be handled. The options
#'   and defaults depend on the `method` used. Ignored if no missing data is
#'   present. It should be noted that multiple imputation outperforms all
#'   available missingness methods available in `weightit()` and should probably
#'   be used instead. Consider the \CRANpkg{MatchThem} package for the use of
#'   `weightit()` with multiply imputed data.
#' @param verbose `logical`; whether to print additional information output by
#'   the fitting function.
#' @param include.obj `logical`; whether to include in the output any fit
#'   objects created in the process of estimating the weights. For example, with
#'   `method = "glm"`, the `glm` objects containing the propensity score model
#'   will be included. See the individual pages for each method for information
#'   on what object will be included if `TRUE`.
#' @param keep.mparts `logical`; whether to include in the output components
#'   necessary to estimate standard errors that account for estimation of the
#'   weights in [glm_weightit()]. Default is `TRUE` if such parts are present.
#'   See the individual pages for each method for whether these components are
#'   produced. Set to `FALSE` to keep the output object smaller, e.g., if
#'   standard errors will not be computed using `glm_weightit()`.
#' @param ... other arguments for functions called by `weightit()` that control
#'   aspects of fitting that are not covered by the above arguments. See
#'   Details.
#'
#' @returns
#' A `weightit` object with the following elements:
#'
#' \item{weights}{The estimated weights, one for each unit.}
#' \item{treat}{The values of the treatment variable.}
#' \item{covs}{The covariates used in the fitting. Only includes the raw covariates, which may have been altered in the fitting process.}
#' \item{estimand}{The estimand requested.}
#' \item{method}{The weight estimation method specified.}
#' \item{ps}{The estimated or provided propensity scores. Estimated propensity scores are returned for binary treatments and only when `method` is `"glm"`, `"gbm"`, `"cbps"`, `"ipt"`, `"super"`, or `"bart"`. The propensity score corresponds to the predicted probability of being treated; see section *`estimand` and `focal`* in Details for how the treated group is determined.}
#' \item{s.weights}{The provided sampling weights.}
#' \item{focal}{The focal treatment level if the ATT or ATC was requested.}
#' \item{by}{A data.frame containing the `by` variable when specified.}
#' \item{obj}{When `include.obj = TRUE`, the fit object.}
#' \item{info}{Additional information about the fitting. See the individual methods pages for what is included.}
#'
#' When `keep.mparts` is `TRUE` (the default) and the chosen method is
#' compatible with M-estimation, the components related to M-estimation for use
#' in [glm_weightit()] are stored in the `"Mparts"` attribute. When `by` is
#' specified, `keep.mparts` is set to `FALSE`.
#'
#' @details
#' The primary purpose of `weightit()` is as a dispatcher to functions
#' that perform the estimation of balancing weights using the requested
#' `method`. Below are the methods allowed and links to pages containing more
#' information about them, including additional arguments and outputs (e.g.,
#' when `include.obj = TRUE`), how missing values are treated, which estimands
#' are allowed, and whether sampling weights are allowed.
#'
#' | [`"glm"`][method_glm] | Propensity score weighting using generalized linear models |
#' | :---- | :---- |
#' | [`"gbm"`][method_gbm] | Propensity score weighting using generalized boosted modeling |
#' | [`"cbps"`][method_cbps] | Covariate Balancing Propensity Score weighting |
#' | [`"npcbps"`][method_npcbps]| Non-parametric Covariate Balancing Propensity Score weighting |
#' | [`"ebal"`][method_ebal] | Entropy balancing |
#' | [`"ipt"`][method_ipt] | Inverse probability tilting |
#' | [`"optweight"`][method_optweight] | Stable balancing weights |
#' | [`"super"`][method_super] | Propensity score weighting using SuperLearner |
#' | [`"bart"`][method_bart] | Propensity score weighting using Bayesian additive regression trees (BART) |
#' | [`"energy"`][method_energy] | Energy balancing |
#'
#' `method` can also be supplied as a user-defined function; see [`method_user`]
#' for instructions and examples. Setting `method = NULL` computes unit weights.
#'
#' ## `estimand` and `focal`
#'
#' For binary and multi-category treatments, the
#' argument to `estimand` determines what distribution the weighted sample
#' should resemble. When set to `"ATE"`, this requests that each group resemble
#' the full sample. When set to `"ATO"`, `"ATM"`, or `"ATOS"` (for the methods
#' that allow them), this requests that each group resemble an "overlap" sample.
#' When set to `"ATT"` or `"ATC"`, this requests that each group resemble the
#' treated or control group, respectively (termed the "focal" group). Weights
#' are set to 1 for the focal group.
#'
#' How does `weightit()` decide which group is the treated and which group is
#' the control? For binary treatments, several heuristics are used. The first is
#' by checking whether a valid argument to `focal` was supplied containing the
#' name of the focal group, which is the treated group when `estimand = "ATT"`
#' and the control group when `estimand = "ATC"`. If `focal` is not supplied,
#' guesses are made using the following criteria, evaluated in order:
#'
#' * If the treatment variable is `logical`, `TRUE` is considered treated and `FALSE` control.
#' * If the treatment is numeric (or a string or factor with values that can be coerced to numeric values), if 0 is one of the values, it is considered the control, and otherwise, the lower value is considered the control (with the other considered treated).
#' * If exactly one of the treatment values is `"t"`, `"tr"`, `"treat"`, `"treated"`, or `"exposed"`, it is considered the treated (and the other control).
#' * If exactly one of the treatment values is `"c"`, `"co"`, `"ctrl"`, `"control"`, or `"unexposed"`, it is considered the control (and the other treated).
#' * If the treatment variable is a factor, the first level is considered control and the second treated.
#' * The lowest value after sorting with [sort()] is considered control and the other treated.
#'
#' To be safe, it is best to code your binary treatment variable as `0` for
#' control and `1` for treated. Otherwise, `focal` should be supplied when
#' requesting the ATT or ATC. For multi-category treatments, `focal` is required
#' when requesting the ATT or ATC; none of the heuristics above are used.
#'
#' ## Citing \pkg{WeightIt}
#'
#' When using `weightit()`, please cite both the
#' \pkg{WeightIt} package (using `citation("WeightIt")`) and the paper(s) in the
#' references section of the method used.
#'
#' @seealso
#' [weightitMSM()] for estimating weights with sequential (i.e.,
#' longitudinal) treatments for use in estimating marginal structural models
#' (MSMs).
#'
#' [weightit.fit()], which is a lower-level dispatcher function that accepts a
#' matrix of covariates and a vector of treatment statuses rather than a formula
#' and data frame and performs minimal argument checking and processing. It may
#' be useful for speeding up simulation studies for which the correct arguments
#' are known. In general, `weightit()` should be used.
#'
#' [summary.weightit()] for summarizing the weights
#'
#' @examples
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
                     by = NULL, s.weights = NULL, ps = NULL, missing = NULL, verbose = FALSE,
                     include.obj = FALSE, keep.mparts = TRUE, ...) {

  ## Checks and processing ----

  A <- list(...)

  #Checks
  if (is_null(formula) || !rlang::is_formula(formula, lhs = TRUE)) {
    .err("`formula` must be a formula relating treatment to covariates")
  }

  #Process treat and covs from formula and data
  t.c <- get_covs_and_treat_from_formula2(formula, data)
  reported.covs <- t.c[["reported.covs"]]
  simple.covs <- t.c[["simple.covs"]]

  covs <- t.c[["model.covs"]]
  treat <- t.c[["treat"]]

  if (is_null(treat)) {
    .err("no treatment variable was specified")
  }

  if (length(treat) != nrow(covs)) {
    .err("the treatment and covariates must have the same number of units")
  }

  n <- length(treat)

  if (anyNA(treat)) {
    .err("missing values are not allowed in the treatment variable")
  }

  #Get treat type
  treat <- as.treat(treat, process = TRUE)
  treat.type <- get_treat_type(treat)

  chk::chk_flag(verbose)
  chk::chk_flag(include.obj)
  chk::chk_flag(keep.mparts)

  #Process ps
  ps <- .process_ps(ps, data, treat)
  if (is_not_null(ps) && is_not_null(method) &&
      !is.function(method) &&
      !identical(method, "glm")) {
    .wrn("`ps` is supplied, so `method` will be ignored")
    method <- "glm"
  }

  ##Process method
  .check_acceptable_method(method, msm = FALSE, force = FALSE)

  if (is_null(method)) {
    method <- NULL
  }
  else if (is.character(method)) {
    method <- .method_to_proper_method(method)
    .check_method_treat.type(method, treat.type)
    attr(method, "name") <- method
  }
  else { #function
    method.name <- deparse1(substitute(method))
    .check_user_method(method)
    attr(method, "name") <- method.name
  }

  #Process estimand and focal
  estimand <- .process_estimand(estimand, method, treat.type)
  f.e.r <- .process_focal_and_estimand(focal, estimand, treat)
  focal <- f.e.r[["focal"]]
  estimand <- f.e.r[["reported.estimand"]]
  reported.estimand <- f.e.r[["reported.estimand"]]

  if (treat.type == "binary") {
    attr(treat, "treated") <- f.e.r[["treated"]]
  }

  #Process missing
  missing <- {
    if (is_null(method) || !anyNA(reported.covs)) ""
    else .process_missing(missing, method)
  }

  #Check subclass
  .check_subclass(...get("subclass"), method, treat.type)

  #Process s.weights
  s.weights <- .process.s.weights(s.weights, data)

  if (is_null(s.weights)) s.weights <- rep.int(1, n)
  else .check_method_s.weights(method, s.weights)

  #Process stabilize
  if (is_null(method) || isFALSE(stabilize)) {
    stabilize <- NULL
  }
  else if (isTRUE(stabilize)) {
    if (treat.type == "continuous") {
      .wrn('setting `stabilize = TRUE` does nothing for continuous treatments, so it will be set to `FALSE`. See `help("weightit")` for details')
      stabilize <- NULL
    }
    else {
      stabilize <- as.formula("~ 1")
    }
  }

  if (is_not_null(stabilize)) {
    if (is.character(method) && !.weightit_methods[[method]]$stabilize_ok) {
      .wrn(sprintf("`stabilize` cannot be used with %s and will be ignored",
                   .method_to_phrase(method)))
      stabilize <- NULL
    }
    else {
      if (treat.type != "continuous" && estimand != "ATE") {
        .err('`stabilize` can only be supplied when `estimand = "ATE"`')
      }

      if (!rlang::is_formula(stabilize)) {
        .err("`stabilize` must be `TRUE`, `FALSE`, or a formula with the stabilization factors on the right hand side")
      }

      stabilize <- update(formula, update(stabilize, NULL ~ .))
    }
  }

  ##Process by
  if (is_not_null(A[["exact"]])) {
    .wrn("`by` has replaced `exact` in the `weightit()` syntax, but `exact` will always work")
    by <- A[["exact"]]
    by.arg <- "exact"
  }
  else {
    by.arg <- "by"
  }

  processed.by <- .process_by(by, data = data, treat = treat, by.arg = by.arg)

  #Process moments and int
  m.i.q <- .process_moments_int_quantile(method = method, ...)

  call <- match.call()

  ## Running models ----

  A["treat"] <- list(treat)
  A["covs"] <- list(covs)
  A["s.weights"] <- list(s.weights)
  A["by.factor"] <- list(.attr(processed.by, "by.factor"))
  A["estimand"] <- list(estimand)
  A["focal"] <- list(focal)
  A["stabilize"] <- list(FALSE)
  A["method"] <- list(method)
  A[c("moments", "int", "quantile")] <- m.i.q[c("moments", "int", "quantile")]
  A["ps"] <- list(ps)
  A["missing"] <- list(missing)
  A["verbose"] <- list(verbose)
  A["include.obj"] <- list(include.obj)
  A[".data"] <- list(data)
  A[".formula"] <- list(formula)
  A[".covs"] <- list(reported.covs)

  #Returns weights (weights) and propensity score (ps)
  obj <- do.call("weightit.fit", A)

  if (is_not_null(stabilize)) {
    stab.t.c <- get_covs_and_treat_from_formula2(stabilize, data)

    A["covs"] <- list(stab.t.c[["model.covs"]])
    A["method"] <- list("glm")
    A["moments"] <- list(integer())
    A["int"] <- list(FALSE)
    A["quantile"] <- list(list())
    A[".formula"] <- list(stabilize)
    A[".covs"] <- list(stab.t.c[["reported.covs"]])

    sw_obj <- do.call("weightit.fit", A)

    obj$weights <- obj$weights / sw_obj[["weights"]]
  }

  if (is_not_null(method) && is_null(ps)) {
    .check_estimated_weights(obj$weights, treat, treat.type, s.weights)
  }

  ## Assemble output object----
  out <- list(weights = obj$weights,
              treat = treat,
              covs = simple.covs,
              estimand = if (treat.type == "continuous") NULL else reported.estimand,
              method = method,
              ps = if (is_null(obj$ps) || all(is.na(obj$ps))) NULL else obj$ps,
              s.weights = s.weights,
              focal = if (reported.estimand %in% c("ATT", "ATC")) focal else NULL,
              by = processed.by,
              call = call,
              formula = formula,
              stabilize = stabilize,
              missing = if (nzchar(missing)) missing else NULL,
              env = parent.frame(),
              info = obj$info,
              obj = obj$fit.obj) |>
    clear_null()

  if (keep.mparts && is_not_null(.attr(obj, "Mparts"))) {
    if (is_null(stabilize)) {
      attr(out, "Mparts") <- .attr(obj, "Mparts")
    }
    else {
      stab.Mparts <- .attr(sw_obj, "Mparts")
      #Invert wfun and compute derivative of inverted wfun
      .wfun <- stab.Mparts$wfun
      stab.Mparts$wfun <- Invert(.wfun)

      if (is_not_null(stab.Mparts$dw_dBtreat)) {
        .dw_dBtreat <- stab.Mparts$dw_dBtreat
        stab.Mparts$dw_dBtreat <- function(Btreat, Xtreat, A, SW) {
          -.dw_dBtreat(Btreat, Xtreat, A, SW) / .wfun(Btreat, Xtreat, A)^2
        }
      }

      attr(out, "Mparts.list") <- list(.attr(obj, "Mparts"), stab.Mparts)
    }
  }

  class(out) <- "weightit"

  out
}

#' @exportS3Method print weightit
print.weightit <- function(x, ...) {
  treat.type <- get_treat_type(x[["treat"]])

  cat(sprintf("A %s object\n", .it(class(x)[1L])))

  if (is_not_null(x[["method"]])) {
    method_name <- {
      if (is_not_null(.attr(x[["method"]], "name"))) add_quotes(.attr(x[["method"]], "name"))
      else if (is.character(x[["method"]])) add_quotes(x[["method"]])
      else "user-defined"
    }

    method_note <- {
      if (is_not_null(.attr(x[["method"]], "package")))
        sprintf(" (converted from %s)", .it(.attr(x[["method"]], "package")))
      else if (is_not_null(x[["method"]]))
        sprintf(" (%s)", .method_to_phrase(x[["method"]]))
      else
        ""
    }

    cat(sprintf(" - method: %s%s\n",
                method_name,
                method_note))
  }
  else if (all_the_same(x[["weights"]])) {
    cat(" - method: no weighting\n")
  }
  else if (is_not_null(x[["ps"]])) {
    cat(" - method: propensity score weighting\n")
  }

  cat(sprintf(" - number of obs.: %s\n",
              nobs(x)))

  cat(sprintf(" - sampling weights: %s\n",
              if (is_null(x[["s.weights"]]) || all_the_same(x[["s.weights"]])) "none" else "present"))

  cat(sprintf(" - treatment: %s\n",
              switch(treat.type,
                     continuous = "continuous",
                     `multi-category` =,
                     multinomial = sprintf("%s-category (%s)",
                                           nunique(x[["treat"]]),
                                           word_list(levels(x[["treat"]]), and.or = FALSE)),
                     binary = "2-category")))

  if (is_not_null(x[["estimand"]])) {
    cat(sprintf(" - estimand: %s\n",
                if (is_null(x[["focal"]])) x[["estimand"]]
                else sprintf("%s (focal: %s)", x[["estimand"]], x[["focal"]])))
  }

  if (is_not_null(x[["covs"]])) {
    cat(sprintf(" - covariates: %s\n",
                if (length(names(x[["covs"]])) > 60L) "too many to name"
                else word_list(names(x[["covs"]]), and.or = FALSE)))
  }

  if (is_not_null(x[["missing"]]) && !identical(x[["missing"]], "")) {
    cat(sprintf(" - missingness method: %s\n",
                .missing_to_phrase(x[["missing"]])))
  }

  if (is_not_null(x[["by"]])) {
    cat(sprintf(" - by: %s\n", word_list(names(x[["by"]]), and.or = FALSE)))
  }

  if (is_not_null(x[["moderator"]])) {
    nsubgroups <- nlevels(.attr(x[["moderator"]], "by.factor"))
    cat(sprintf(" - moderator: %s (%s subgroups)\n",
                word_list(names(x[["moderator"]]), and.or = FALSE),
                nsubgroups))
  }

  #trim
  if (is_not_null(.attr(x, "trim"))) {
    trim.at <- .attr(x, "trim")[["at"]]
    trim.lower <- .attr(x, "trim")[["lower"]]
    trim.drop <- .attr(x, "trim")[["drop"]]
  }
  else if (is_not_null(.attr(x[["weights"]], "trim"))) {
    trim.at <- .attr(x[["weights"]], "trim")
    trim.lower <- .attr(x[["weights"]], "trim.lower")
    trim.drop <- FALSE
  }
  else {
    trim.at <- NULL
  }

  if (is_not_null(trim.at) && chk::vld_number(trim.at)) {
    if (trim.at < 1) {
      if (trim.lower) {
        trim.at <- c(1 - trim.at, trim.at)
      }

      cat(sprintf(" - weights trimmed at %s%s\n",
                  word_list(paste0(round(100 * trim.at, 2L), "%")),
                  if (trim.drop) " and units dropped" else ""))
    }
    else {
      cat(sprintf(" - weights trimmed at the %s %s%s\n",
                  if (trim.lower) "top and bottom" else "top",
                  trim.at,
                  if (trim.drop) " and units dropped" else ""))
    }
  }

  invisible(x)
}

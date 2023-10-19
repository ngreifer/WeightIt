#' Covariate Balancing Propensity Score Weighting
#' @name method_cbps
#' @aliases method_cbps
#' @usage NULL
#'
#' @description
#'
#' This page explains the details of estimating weights from covariate balancing propensity scores by setting `method = "cbps"` in the call to [weightit()] or [weightitMSM()]. This method can be used with binary, multinomial, and continuous treatments.
#'
#' In general, this method relies on estimating propensity scores using generalized method of moments and then converting those propensity scores into weights using a formula that depends on the desired estimand. This method relies on \pkgfun{CBPS}{CBPS} from the \CRANpkg{CBPS} package.
#'
#' ## Binary Treatments
#'
#' For binary treatments, this method estimates the propensity scores and weights using \pkgfun{CBPS}{CBPS}. The following estimands are allowed: ATE, ATT, and ATC. The weights are taken from the output of the `CBPS` fit object. When the estimand is the ATE, the return propensity score is the probability of being in the "second" treatment group, i.e., `levels(factor(treat))[2]`; when the estimand is the ATC, the returned propensity score is the probability of being in the control (i.e., non-focal) group.
#'
#' ## Multinomial Treatments
#'
#' For multinomial treatments with three or four categories and when the estimand is the ATE, this method estimates the propensity scores and weights using one call to \pkgfun{CBPS}{CBPS}. For multinomial treatments with three or four categories or when the estimand is the ATT, this method estimates the propensity scores and weights using multiple calls to \pkgfun{CBPS}{CBPS}. The following estimands are allowed: ATE and ATT. The weights are taken from the output of the `CBPS` fit objects.
#'
#' ## Continuous Treatments
#'
#' For continuous treatments, the generalized propensity score and weights are estimated using \pkgfun{CBPS}{CBPS}.
#'
#' ## Longitudinal Treatments
#'
#' For longitudinal treatments, the weights are the product of the weights estimated at each time point. This is not how \pkgfun{CBPS}{CBMSM} in the \pkg{CBPS} package estimates weights for longitudinal treatments.
#'
#' ## Sampling Weights
#'
#' Sampling weights are supported through `s.weights` in all scenarios. See Note about sampling weights.
#'
#' ## Missing Data
#'
#' In the presence of missing data, the following value(s) for `missing` are allowed:
#'     \describe{
#'       \item{`"ind"` (default)}{
#'         First, for each variable with missingness, a new missingness indicator variable is created which takes the value 1 if the original covariate is `NA` and 0 otherwise. The missingness indicators are added to the model formula as main effects. The missing values in the covariates are then replaced with the covariate medians (this value is arbitrary and does not affect estimation). The weight estimation then proceeds with this new formula and set of covariates. The covariates output in the resulting `weightit` object will be the original covariates with the `NA`s.
#'       }
#'     }
#'
#' @section Additional Arguments:
#' All arguments to `CBPS()` can be passed through `weightit()` or `weightitMSM()`, with the following exceptions:
#'   * `method` in `CBPS()` is replaced with the argument `over` in `weightit()`. Setting `over = FALSE` in `weightit()` is the equivalent of setting `method = "exact"` in `CBPS()`.
#'   * `sample.weights` is ignored because sampling weights are passed using `s.weights`.
#'   * `standardize` is ignored.
#'
#' All arguments take on the defaults of those in `CBPS()`. It may be useful in many cases to set `over = FALSE`, especially with continuous treatments.
#'
#' @section Additional Outputs:
#' \describe{
#'   \item{`obj`}{When `include.obj = TRUE`, the CB(G)PS model fit. For binary treatments, multinomial treatments with `estimand = "ATE"` and four or fewer treatment levels, and continuous treatments, the output of the call to \pkgfun{CBPS}{CBPS}. For multinomial treatments with `estimand = "ATT"` or with more than four treatment levels, a list of `CBPS` fit objects.
#'   }
#' }
#'
#' @details
#' CBPS estimates the coefficients of a logistic regression model (for binary treatments), multinomial logistic regression model (form multinomial treatments), or linear regression model (for continuous treatments) that is used to compute (generalized) propensity scores, from which the weights are computed. It involves augmenting the standard regression score equations with the balance constraints in an over-identified generalized method of moments estimation. The idea is to nudge the estimation of the coefficients toward those that produce balance in the weighted sample. The just-identified version (with `exact = FALSE`) does away with the score equations for the coefficients so that only the balance constraints (and the score equation for the variance of the error with a continuous treatment) are used. The just-identified version will therefore produce superior balance on the means (i.e., corresponding to the balance constraints) for binary and multinomial treatments and linear terms for continuous treatments than will the over-identified version.
#'
#' Note that \pkg{WeightIt} provides less functionality than does the \pkg{CBPS} package in terms of the versions of CBPS available; for extensions to CBPS, the \pkg{CBPS} package may be preferred.
#'
#' @note
#' When sampling weights are used with `CBPS::CBPS()`, the estimated weights already incorporate the sampling weights. When `weightit()` is used with `method = "cbps"`, the estimated weights are separated from the sampling weights, as they are with all other methods.
#'
#' @seealso
#' [weightit()], [weightitMSM()]
#'
#' \pkgfun{CBPS}{CBPS} for the fitting function.
#'
#' @references
#' ## Binary treatments
#'
#' Imai, K., & Ratkovic, M. (2014). Covariate balancing propensity score. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 76(1), 243–263.
#'
#' ## Multinomial Treatments
#'
#' Imai, K., & Ratkovic, M. (2014). Covariate balancing propensity score. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 76(1), 243–263.
#'
#'
#' ## Continuous treatments
#'
#' Fong, C., Hazlett, C., & Imai, K. (2018). Covariate balancing propensity score for a continuous treatment: Application to the efficacy of political advertisements. The Annals of Applied Statistics, 12(1), 156–177. \doi{10.1214/17-AOAS1101}
#'
#' @examplesIf requireNamespace("CBPS", quietly = TRUE)
#' library("cobalt")
#' data("lalonde", package = "cobalt")
#'
#' #Balancing covariates between treatment groups (binary)
#' (W1 <- weightit(treat ~ age + educ + married +
#'                   nodegree + re74, data = lalonde,
#'                 method = "cbps", estimand = "ATT"))
#' summary(W1)
#' bal.tab(W1)
#'
#' \dontrun{
#'   #Balancing covariates with respect to race (multinomial)
#'   (W2 <- weightit(race ~ age + educ + married +
#'                     nodegree + re74, data = lalonde,
#'                   method = "cbps", estimand = "ATE"))
#'   summary(W2)
#'   bal.tab(W2)
#' }
#'
#' #Balancing covariates with respect to re75 (continuous)
#' (W3 <- weightit(re75 ~ age + educ + married +
#'                   nodegree + re74, data = lalonde,
#'                 method = "cbps", over = FALSE))
#' summary(W3)
#' bal.tab(W3)
NULL

weightit2cbps <- function(covs, treat, s.weights, estimand, focal, subset,
                          stabilize, subclass, missing, verbose, ...) {

  rlang::check_installed("CBPS")

  A <- list(...)

  covs <- covs[subset, , drop = FALSE]
  treat <- factor(treat[subset])
  s.weights <- s.weights[subset]

  if (missing == "ind") {
    covs <- add_missing_indicators(covs)
  }

  for (i in seq_col(covs)) covs[,i] <- make.closer.to.1(covs[,i])

  colinear.covs.to.remove <- colnames(covs)[colnames(covs) %nin% colnames(make_full_rank(covs))]
  covs <- covs[, colnames(covs) %nin% colinear.covs.to.remove, drop = FALSE]

  sw0 <- check_if_zero(s.weights)

  if (estimand == "ATT") {
    treat_ <- as.integer(treat == focal)
    new.data <- as.data.frame(cbind(treat_, covs))

    tryCatch({verbosely({
      fit <- CBPS::CBPS(formula(new.data),
                        data = new.data[!sw0,],
                        method = if (is_not_null(A$over) && A$over == FALSE) "exact" else "over",
                        standardize = FALSE,
                        sample.weights = s.weights[!sw0],
                        ATT = 1,
                        ...)
    }, verbose = verbose)},
    error = function(e) {
      e. <- conditionMessage(e)
      e. <- gsub("method = \"exact\"", "`over = FALSE`", e., fixed = TRUE)
      .err("(from `CBPS::CBPS()`) ", e., tidy = FALSE)
    })

    if (!any(sw0)) {
      ps <- fit[["fitted.values"]]
    }
    else {
      ps <- plogis(drop(cbind(1, covs) %*% fit[["coefficients"]]))
    }
  }
  else {
    new.data <- cbind(treat, as.data.frame(covs))

    tryCatch({verbosely({
      fit <- CBPS::CBPS(formula(new.data),
                             data = new.data[!sw0,],
                             method = if (isFALSE(A$over)) "exact" else "over",
                             standardize = FALSE,
                             sample.weights = s.weights[!sw0],
                             ATT = 0,
                             ...)
    }, verbose = verbose)},
    error = function(e) {
      e. <- conditionMessage(e)
      e. <- gsub("method = \"exact\"", "`over = FALSE`", e., fixed = TRUE)
      .err("(from `CBPS::CBPS()`) ", e., tidy = FALSE)
    })

    if (!any(sw0)) {
      ps <- fit[["fitted.values"]]
    }
    else {
      ps <- plogis(drop(cbind(1, covs) %*% fit[["coefficients"]]))
    }
  }

  w <- get_w_from_ps(ps, treat, estimand = estimand, subclass = subclass,
                     focal = focal, stabilize = stabilize)

  p.score <- {
    if (is_null(dim(ps)) || length(dim(ps)) != 2) ps
    else ps[[get_treated_level(treat)]]
  }

  list(w = w, ps = p.score, fit.obj = fit)
}

weightit2cbps.multi <- function(covs, treat, s.weights, estimand, focal, subset,
                                stabilize, subclass, missing, verbose, ...) {

  rlang::check_installed("CBPS")

  A <- list(...)

  covs <- covs[subset, , drop = FALSE]
  treat <- factor(treat[subset])
  s.weights <- s.weights[subset]

  if (missing == "ind") {
    covs <- add_missing_indicators(covs)
  }

  for (i in seq_col(covs)) covs[,i] <- make.closer.to.1(covs[,i])

  colinear.covs.to.remove <- colnames(covs)[colnames(covs) %nin% colnames(make_full_rank(covs))]
  covs <- covs[, colnames(covs) %nin% colinear.covs.to.remove, drop = FALSE]

  sw0 <- check_if_zero(s.weights)

  if (estimand == "ATT") {
    ps <- make_df(levels(treat), length(treat))

    control.levels <- levels(treat)[levels(treat) != focal]
    fit.list <- make_list(control.levels)

    for (i in control.levels) {
      treat.in.i.focal <- treat %in% c(focal, i)
      treat_ <- as.integer(treat[treat.in.i.focal] != i)
      covs_ <- covs[treat.in.i.focal, , drop = FALSE]
      new.data <- as.data.frame(cbind(treat_, covs_))

      tryCatch({verbosely({
        fit.list[[i]] <- CBPS::CBPS(formula(new.data),
                                    data = new.data[!sw0[treat.in.i.focal],],
                                    method = if (is_not_null(A$over) && A$over == FALSE) "exact" else "over",
                                    standardize = FALSE,
                                    sample.weights = s.weights[!sw0 & treat.in.i.focal],
                                    ATT = 1,
                                    ...)
      }, verbose = verbose)},
      error = function(e) {
        e. <- conditionMessage(e)
        e. <- gsub("method = \"exact\"", "`over = FALSE`", e., fixed = TRUE)
        .err("(from `CBPS::CBPS()`) ", e., tidy = FALSE)
      })

      if (!any(sw0[treat.in.i.focal])) {
        ps[[focal]][treat.in.i.focal] <- fit.list[[i]][["fitted.values"]]
        ps[[i]][treat.in.i.focal] <- 1 - ps[[focal]][treat.in.i.focal]
      }
      else {
        ps[[focal]][treat.in.i.focal] <- plogis(drop(cbind(1, covs_) %*% fit.list[[i]][["coefficients"]]))
        ps[[i]][treat.in.i.focal] <- 1 - ps[[focal]][treat.in.i.focal]
      }
    }
  }
  else {
    new.data <- cbind(treat, as.data.frame(covs))
    if (nlevels(treat) <= 4) {
      tryCatch({verbosely({
        fit.list <- CBPS::CBPS(formula(new.data),
                               data = new.data[!sw0,],
                               method = if (isFALSE(A$over)) "exact" else "over",
                               standardize = FALSE,
                               sample.weights = s.weights[!sw0],
                               ATT = 0,
                               ...)
      }, verbose = verbose)},
      error = function(e) {
        e. <- conditionMessage(e)
        e. <- gsub("method = \"exact\"", "`over = FALSE`", e., fixed = TRUE)
        .err("(from `CBPS::CBPS()`) ", e., tidy = FALSE)
      })

      if (!any(sw0)) {
        ps <- fit.list[["fitted.values"]]
      }
      else {
        ps <- make_df(levels(treat), length(treat))
        base.lvl <- setdiff(levels(treat), colnames(fit.list$coefficients))

        lin.pred <- cbind(1, covs) %*% fit.list[["coefficients"]]
        ps[, base.lvl] <- 1/rowSums(exp(lin.pred))
        ps[, colnames(fit.list$coefficients)] <- exp(lin.pred[, colnames(fit.list$coefficients)]) * ps[, base.lvl]
        ps <- ps/rowSums(ps)
      }
    }
    else {
      ps <- make_df(levels(treat), length(treat))

      fit.list <- make_list(levels(treat))

      for (i in levels(treat)) {
        new.data[[1]] <- as.integer(treat == i)

        tryCatch({verbosely({
          fit.list[[i]] <- CBPS::CBPS(formula(new.data), data = new.data[!sw0,],
                                      method = if (isFALSE(A$over)) "exact" else "over",
                                      standardize = FALSE,
                                      sample.weights = s.weights[!sw0],
                                      ATT = 0, ...)
        }, verbose = verbose)},
        error = function(e) {
          e. <- conditionMessage(e)
          e. <- gsub("method = \"exact\"", "`over = FALSE`", e., fixed = TRUE)
          .err("(from `CBPS::CBPS()`) ", e., tidy = FALSE)
        })

        if (!any(sw0)) {
          ps[, i] <- fit.list[[i]][["fitted.values"]]
        }
        else {
          ps[, i] <- plogis(drop(cbind(1, covs) %*% fit.list[[i]][["coefficients"]]))
        }
      }
    }
  }

  w <- get_w_from_ps(ps, treat, estimand = estimand, subclass = subclass,
                     focal = focal, stabilize = stabilize)

  p.score <- NULL

  list(w = w, ps = p.score, fit.obj = fit.list)
}

weightit2cbps.cont <- function(covs, treat, s.weights, subset, missing, verbose, ...) {
  rlang::check_installed("CBPS")

  A <- list(...)

  covs <- covs[subset, , drop = FALSE]
  treat <- treat[subset]
  s.weights <- s.weights[subset]

  if (missing == "ind") {
    covs <- add_missing_indicators(covs)
  }

  for (i in seq_col(covs)) covs[,i] <- make.closer.to.1(covs[,i])

  colinear.covs.to.remove <- colnames(covs)[colnames(covs) %nin% colnames(make_full_rank(covs))]
  covs <- covs[, colnames(covs) %nin% colinear.covs.to.remove, drop = FALSE]

  sw0 <- check_if_zero(s.weights)

  new.data <- data.frame(treat = treat, covs)

  w <- rep(0, length(treat))

  tryCatch({verbosely({
    fit <- CBPS::CBPS(formula(new.data),
                      data = new.data[!sw0,],
                      method = if (isFALSE(A$over)) "exact" else "over",
                      standardize = FALSE,
                      sample.weights = s.weights[!sw0],
                      ...)
  }, verbose = verbose)},
  error = function(e) {
    e. <- conditionMessage(e)
    e. <- gsub("method = \"exact\"", "`over = FALSE`", e., fixed = TRUE)
    .err("(from `CBPS::CBPS()`) ", e., tidy = FALSE)
  })

  w[!sw0] <- fit$weights / s.weights[!sw0]

  list(w = w, fit.obj = fit)
}

weightitMSM2cbps <- function(covs.list, treat.list, s.weights, subset, missing, verbose, ...) {
  .err("CBMSM doesn't work yet")
}
#' Inverse Probability Tilting
#' @name method_ipt
#' @aliases method_ipt
#' @usage NULL
#'
#' @description
#' This page explains the details of estimating weights using inverse probability tilting by setting `method = "ipt"` in the call to [weightit()] or [weightitMSM()]. This method can be used with binary, multi-category, and continuous treatments.
#'
#' In general, this method relies on estimating propensity scores using a modification of the usual regression score equations to enforce balance and the converting those propensity scores into weights using a formula that depends on the desired estimand. This method relies on code written for \pkg{WeightIt} using [optim()].
#'
#' ## Binary Treatments
#'
#' For binary treatments, this method estimates the weights using `optim()` using formulas described by Graham, Pinto, and Egel (2012). The following estimands are allowed: ATE, ATT, and ATC. When the ATE is requested, the optimization is run twice, once for each treatment group, and the return propensity score is the probability of being in the "second" treatment group, i.e., `levels(factor(treat))[2]`.
#'
#' ## Multi-Category Treatments
#'
#' For multi-category treatments, this method estimates the weights using `optim()`. The following estimands are allowed: ATE and ATT. When the ATE is requested, `optim()` is run once for each treatment group. When the ATT is requested, `optim()` is run once for each non-focal (i.e., control) group.
#'
#' ## Continuous Treatments
#'
#' For continuous treatments, this method estimates the weights using `optim()`.
#'
#' ## Longitudinal Treatments
#'
#' For longitudinal treatments, the weights are the product of the weights estimated at each time point. This method is not guaranteed to yield exact balance at each time point. NOTE: the use of inverse probability tilting with longitudinal treatments has not been validated!
#'
#' ## Sampling Weights
#'
#' Sampling weights are supported through `s.weights` in all scenarios.
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
#' `moments` and `int` are accepted. See [weightit()] for details.
#'
#' \describe{
#'   \item{`quantile`}{
#'     A named list of quantiles (values between 0 and 1) for each continuous covariate, which are used to create additional variables that when balanced ensure balance on the corresponding quantile of the variable. For example, setting `quantile = list(x1 = c(.25, .5. , .75))` ensures the 25th, 50th, and 75th percentiles of `x1` in each treatment group will be balanced in the weighted sample. Can also be a single number (e.g., `.5`) or an unnamed list of length 1 (e.g., `list(c(.25, .5, .75))`) to request the same quantile(s) for all continuous covariates, or a named vector (e.g., `c(x1 = .5, x2 = .75`) to request one quantile for each covariate. Only allowed with binary and multi-category treatments.
#'   }
#' }
#'
#' The arguments `maxit` and `reltol` can be supplied and are passed to the `control` argument of [optim()]. The `"BFGS"` method is used, so the defaults correspond to this.
#'
#' The `stabilize` argument is ignored.
#'
#' @section Additional Outputs:
#' \describe{
#'   \item{`obj`}{When `include.obj = TRUE`, the output of the call to [optim()], which contains the coefficient estimates and convergence information. For ATE fits or with multi-category treatments, a list of `optim()` outputs, one for each weighted group.
#'   }
#' }
#'
#' @details
#' Inverse probability tilting (IPT) involves specifying estimating equations that fit the parameters of one or more logistic regression models (for binary and multi-category treatments) or a linear regression model (for continuous treatments) with a modification that ensure exact balance on the covariates means. These estimating equations are solved, and the estimated parameters are used in the (generalized) propensity score, which is used to compute the weights. Conceptually and mathematically, IPT is very similar to entropy balancing and just-identified CBPS. For the ATT and ATC, entropy balancing, just-identified CBPS, and IPT will yield identical results. For the ATE, the three methods differ.
#'
#' Treatment effect estimates for binary treatments are consistent if the true propensity score is a logistic regression or the outcome model is linear in the covariates and their interaction with treatments. For entropy balancing, this is only true for the ATT, and fo just-identified CBPS, this is only true if there is no effect modification by covariates. In this way, IPT provides additional theoretical guarantees over the other two methods, though potentially with some cost in precision.
#'
#' IPT for continuous treatments has not been described in the literature; here we use the score equations for CBPS with continuous treatments, which correspond to the score equations for the residual variance of a linear regression for the generalized propensity score and the weighted covariance between the treatment and covariates. This method should be used with caution until it is better understood.
#'
#' @seealso
#' [weightit()], [weightitMSM()]
#'
#' [method_ebal] and [method_cbps] for entropy balancing and CBPS, which work similarly.
#'
#' @references
#' ## ATE
#'
#' Graham, B. S., De Xavier Pinto, C. C., & Egel, D. (2012). Inverse Probability Tilting for Moment Condition Models with Missing Data. *The Review of Economic Studies*, 79(3), 1053–1079. \doi{10.1093/restud/rdr047}
#'
#' ## ATT
#'
#' Sant’Anna, P. H. C., & Zhao, J. (2020). Doubly robust difference-in-differences estimators. *Journal of Econometrics*, 219(1), 101–122. \doi{10.1016/j.jeconom.2020.06.003}
#'
#' @examples
#' data("lalonde", package = "cobalt")
#'
#' #Balancing covariates between treatment groups (binary)
#' (W1 <- weightit(treat ~ age + educ + married +
#'                   nodegree + re74, data = lalonde,
#'                 method = "ipt", estimand = "ATT"))
#' summary(W1)
#' cobalt::bal.tab(W1)
#'
#' #Balancing covariates with respect to race (multi-category)
#' (W2 <- weightit(race ~ age + educ + married +
#'                   nodegree + re74, data = lalonde,
#'                 method = "ipt", estimand = "ATE"))
#' summary(W2)
#' cobalt::bal.tab(W2)
#'
#' #Balancing covariates and squares with respect to
#' #re75 (continuous)
#' (W3 <- weightit(re75 ~ age + educ + married +
#'                   nodegree + re74, data = lalonde,
#'                 method = "ipt"))
#' summary(W3)
#' cobalt::bal.tab(W3)
NULL

weightit2ipt <- function(covs, treat, s.weights, subset, estimand, focal,
                         stabilize, missing, moments, int, verbose, ...) {
  A <- list(...)

  covs <- covs[subset, , drop = FALSE]
  treat <- factor(treat[subset])
  s.weights <- s.weights[subset]

  if (missing == "ind") {
    covs <- add_missing_indicators(covs)
  }

  covs <- cbind(covs, int.poly.f(covs, poly = moments, int = int, center = TRUE))

  covs <- cbind(covs, quantile_f(covs, qu = A[["quantile"]], s.weights = s.weights,
                                 focal = focal, treat = treat))

  for (i in seq_col(covs)) covs[,i] <- make.closer.to.1(covs[,i])

  C <- cbind(1, covs)

  if (estimand == "ATT") {
    groups_to_weight <- levels(treat)[levels(treat) != focal]
    targets <- cobalt::col_w_mean(C, s.weights = s.weights, subset = treat == focal)
  }
  else if (estimand == "ATE") {
    groups_to_weight <- levels(treat)
    targets <- cobalt::col_w_mean(C, s.weights = s.weights)
  }

  reltol <- if_null_then(A[["reltol"]], sqrt(.Machine$double.eps))

  if (estimand == "ATT") {
    f <- function(theta, C, A, s) {
      p <- plogis(drop(C %*% theta))
      colMeans(C * s * (1 - A) * p / (1 - p)) - targets
    }

    start <- rep(0, ncol(C))

    # out <- rootSolve::multiroot(f, start = start,
    #                             A = as.numeric(treat == focal),
    #                             C = C,
    #                             s = s.weights)
    # par <- out$root
    verbosely({
      fit.list <- optim(par = start,
                        fn = function(b, ...) sum(f(b, ...)^2),
                        method = "BFGS",
                        control = list(trace = 1,
                                       reltol = reltol,
                                       maxit = if_null_then(A[["maxit"]], 200)),
                        A = as.numeric(treat == focal),
                        C = C,
                        s = s.weights)
    }, verbose = verbose)

    par <- fit.list$par

    ps <- plogis(drop(C %*% par))
  }
  else {
    ps <- make_df(levels(treat), length(treat))

    #Version that only balances groups together, not to target; equal to CBPS with over = F
    # f <- function(theta, C, A, s) {
    #   p <- plogis(drop(C %*% theta))
    #   colMeans(s * (A/p - (1 - A)/(1 - p)) * C)
    # }
    #
    # out <- rootSolve::multiroot(f, start = rep(0, ncol(C)),
    #                             A = as.numeric(treat == levels(treat)[2]),
    #                             C = C,
    #                             s = s.weights)
    #
    # ps <- plogis(drop(C %*% out$root))

    fit.list <- make_list(levels(treat))

    for (i in levels(treat)) {
      f <- function(theta, C, A, s) {
        p <- plogis(drop(C %*% theta))
        colMeans(C * s * (A/p)) - targets
      }

      start <- rep(0, ncol(C))

      # out <- rootSolve::multiroot(f, start = start,
      #                             A = as.numeric(treat == i),
      #                             C = C,
      #                             s = s.weights)
      # par <- out$root

      verbosely({
        fit.list[[i]] <- optim(par = start,
                               fn = function(b, ...) sum(f(b, ...)^2),
                               method = "BFGS",
                               control = list(trace = 1,
                                              reltol = reltol,
                                              maxit = if_null_then(A[["maxit"]], 200)),
                               A = as.numeric(treat == i),
                               C = C,
                               s = s.weights)
      }, verbose = verbose)

      par <- fit.list[[i]]$par

      ps[[i]] <- plogis(drop(C %*% par))
    }
  }

  w <- get_w_from_ps(ps, treat, estimand = estimand,
                     focal = focal)

  p.score <- {
    if (is_null(dim(ps)) || length(dim(ps)) != 2) ps
    else ps[[get_treated_level(treat)]]
  }

  list(w = w, ps = p.score, fit.obj = fit.list)
}

weightit2ipt.multi <- function(covs, treat, s.weights, subset, estimand, focal,
                               stabilize, missing, moments, int, verbose, ...) {
  A <- list(...)

  covs <- covs[subset, , drop = FALSE]
  treat <- factor(treat[subset])
  s.weights <- s.weights[subset]

  if (missing == "ind") {
    covs <- add_missing_indicators(covs)
  }

  covs <- cbind(covs, int.poly.f(covs, poly = moments, int = int, center = TRUE))

  covs <- cbind(covs, quantile_f(covs, qu = A[["quantile"]], s.weights = s.weights,
                                 focal = focal, treat = treat))

  for (i in seq_col(covs)) covs[,i] <- make.closer.to.1(covs[,i])

  C <- cbind(1, covs)

  if (estimand == "ATT") {
    groups_to_weight <- levels(treat)[levels(treat) != focal]
    targets <- cobalt::col_w_mean(C, s.weights = s.weights, subset = treat == focal)
  }
  else if (estimand == "ATE") {
    groups_to_weight <- levels(treat)
    targets <- cobalt::col_w_mean(C, s.weights = s.weights)
  }

  reltol <- if_null_then(A[["reltol"]], sqrt(.Machine$double.eps))

  w <- rep(1, length(treat))

  fit.list <- make_list(groups_to_weight)

  if (estimand == "ATT") {

    for (i in groups_to_weight) {
      f <- function(theta, C, A, s) {
        p <- plogis(drop(C %*% theta))
        colMeans(C * (s * (1 - A) * p / (1 - p))) - targets
      }

      start <- rep(0, ncol(C))

      # out <- rootSolve::multiroot(f, start = start,
      #                             A = as.numeric(treat[treat %in% c(i, focal)] == focal),
      #                             C = C[treat %in% c(i, focal),, drop = FALSE],
      #                             s = s.weights[treat %in% c(i, focal)])
      # par <- out$root
      verbosely({
        fit.list[[i]] <- optim(par = start,
                               fn = function(b, ...) sum(f(b, ...)^2),
                               A = as.numeric(treat[treat %in% c(i, focal)] == focal),
                               method = "BFGS",
                               control = list(trace = 1,
                                              reltol = reltol,
                                              maxit = if_null_then(A[["maxit"]], 200)),
                               C = C[treat %in% c(i, focal),, drop = FALSE],
                               s = s.weights[treat %in% c(i, focal)])
      }, verbose = verbose)

      par <- fit.list[[i]]$par

      ps_i <- plogis(drop(C[treat  == i,, drop = FALSE] %*% par))
      w[treat == i] <- ps_i / (1 - ps_i)
    }

  }
  else {
    #Version that only balances groups together, not to target; equal to CBPS with over = F
    # f <- function(theta, C, A, s) {
    #   p <- plogis(drop(C %*% theta))
    #   colMeans(s * (A/p - (1 - A)/(1 - p)) * C)
    # }
    #
    # out <- rootSolve::multiroot(f, start = rep(0, ncol(C)),
    #                             A = as.numeric(treat == levels(treat)[2]),
    #                             C = C,
    #                             s = s.weights)
    #
    # ps <- plogis(drop(C %*% out$root))

    for (i in groups_to_weight) {
      f <- function(theta, C, A, s) {
        p <- plogis(drop(C %*% theta))
        colMeans(C * (s * A / p)) - targets
      }

      start <- rep(0, ncol(C))

      # out <- rootSolve::multiroot(f, start = start,
      #                             A = as.numeric(treat == i),
      #                             C = C,
      #                             s = s.weights)
      # par <- out$root

      verbosely({
        fit.list[[i]] <- optim(par = start,
                               fn = function(b, ...) sum(f(b, ...)^2),
                               method = "BFGS",
                               control = list(trace = 1,
                                              reltol = reltol,
                                              maxit = if_null_then(A[["maxit"]], 200)),
                               A = as.numeric(treat == i),
                               C = C,
                               s = s.weights)
      }, verbose = verbose)

      par <- fit.list[[i]]$par

      ps_i <- plogis(drop(C[treat == i,, drop = FALSE] %*% par))

      w[treat == i] <- 1 / ps_i
    }
  }

  list(w = w, fit.obj = fit.list)
}

weightit2ipt.cont <- function(covs, treat, s.weights, subset, missing, moments, int, verbose, ...) {
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

  treat_c <- (treat - w.m(treat, s.weights)) / sqrt(col.w.v(treat, s.weights))
  covs_c <- cbind(1, scale(covs,
                           center = col.w.m(covs, s.weights),
                           scale = sqrt(col.w.v(covs, s.weights))))

  reltol <- if_null_then(A[["reltol"]], sqrt(.Machine$double.eps))

  #Stabilization - get dens.num
  dens.num <- dnorm(treat_c, log = TRUE)

  f <- function(theta, C, A, sw) {
    s2 <- exp(theta[1])
    b <- theta[1 + 1:ncol(C)]

    p <- drop(C %*% b)

    dens.denom <- dnorm(treat_c, p, sqrt(s2), log = TRUE)

    w <- exp(dens.num - dens.denom)

    c(mean(sw * (A - p)^2/s2 - 1),
      colMeans(w * sw * A * C))
  }

  start <- rep(0, 1 + ncol(covs_c))

  verbosely({
    fit.list <- optim(par = start,
                      fn = function(b, ...) sum(f(b, ...)^2),
                      method = "BFGS",
                      control = list(trace = 1,
                                     reltol = reltol,
                                     maxit = if_null_then(A[["maxit"]], 200)),
                      A = treat_c,
                      C = covs_c,
                      sw = s.weights)
  }, verbose = verbose)

  par <- fit.list$par

  # initfit <- lm.wfit(covs_c, treat_c, s.weights)
  # start <- c(log(var(initfit$residuals)), initfit$coefficients)
  #
  # out <- rootSolve::multiroot(f, start = start,
  #                             A = treat_c,
  #                             C = covs_c,
  #                             sw = s.weights)
  # par <- out$root

  s2 <- exp(par[1])
  b <- par[1 + 1:ncol(covs_c)]
  p <- drop(covs_c %*% b)
  dens.denom <- dnorm(treat_c, p, sqrt(s2), log = TRUE)

  w <- exp(dens.num - dens.denom)

  list(w = w/mean_fast(w), fit.obj = fit.list)
}

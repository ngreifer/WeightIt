#' Entropy Balancing
#' @name method_ebal
#' @aliases method_ebal
#' @usage NULL
#'
#' @description
#' This page explains the details of estimating weights using entropy balancing by setting `method = "ebal"` in the call to [weightit()] or [weightitMSM()]. This method can be used with binary, multi-category, and continuous treatments.
#'
#' In general, this method relies on estimating weights by minimizing the negative entropy of the weights subject to exact moment balancing constraints. This method relies on code written for \pkg{WeightIt} using [optim()].
#'
#' ## Binary Treatments
#'
#' For binary treatments, this method estimates the weights using `optim()` using formulas described by Hainmueller (2012). The following estimands are allowed: ATE, ATT, and ATC. When the ATE is requested, the optimization is run twice, once for each treatment group.
#'
#' ## Multi-Category Treatments
#'
#' For multi-category treatments, this method estimates the weights using `optim()`. The following estimands are allowed: ATE and ATT. When the ATE is requested, `optim()` is run once for each treatment group. When the ATT is requested, `optim()` is run once for each non-focal (i.e., control) group.
#'
#' ## Continuous Treatments
#'
#' For continuous treatments, this method estimates the weights using `optim()` using formulas described by Tübbicke (2022) and Vegetabile et al. (2021).
#'
#' ## Longitudinal Treatments
#'
#' For longitudinal treatments, the weights are the product of the weights estimated at each time point. This method is not guaranteed to yield exact balance at each time point. NOTE: the use of entropy balancing with longitudinal treatments has not been validated!
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
#' ## M-estimation
#'
#' M-estimation is supported for all scenarios. See [glm_weightit()] and `vignette("estimating-effects")` for details.
#'
#' @section Additional Arguments:
#' `moments` and `int` are accepted. See [weightit()] for details.
#'
#' \describe{
#'   \item{`base.weights`}{
#'     A vector of base weights, one for each unit. These correspond to the base weights $q$ in Hainmueller (2012). The estimated weights minimize the Kullback entropy divergence from the base weights, defined as \eqn{\sum w \log(w/q)}, subject to exact balance constraints. These can be used to supply previously estimated weights so that the newly estimated weights retain the some of the properties of the original weights while ensuring the balance constraints are met. Sampling weights should not be passed to `base.weights` but can be included in a `weightit()` call that includes `s.weights`.
#'   }
#'   \item{`quantile`}{
#'     A named list of quantiles (values between 0 and 1) for each continuous covariate, which are used to create additional variables that when balanced ensure balance on the corresponding quantile of the variable. For example, setting `quantile = list(x1 = c(.25, .5. , .75))` ensures the 25th, 50th, and 75th percentiles of `x1` in each treatment group will be balanced in the weighted sample. Can also be a single number (e.g., `.5`) or an unnamed list of length 1 (e.g., `list(c(.25, .5, .75))`) to request the same quantile(s) for all continuous covariates, or a named vector (e.g., `c(x1 = .5, x2 = .75`) to request one quantile for each covariate. Only allowed with binary and multi-category treatments.
#'   }
#'   \item{`d.moments`}{
#'     With continuous treatments, the number of moments of the treatment and covariate distributions that are constrained to be the same in the weighted sample as in the original sample. For example, setting `d.moments = 3` ensures that the mean, variance, and skew of the treatment and covariates are the same in the weighted sample as in the unweighted sample. `d.moments` should be greater than or equal to `moments` and will be automatically set accordingly if not (or if not specified). Vegetabile et al. (2021) recommend setting `d.moments = 3`, even if `moments` is less than 3. This argument corresponds to the tuning parameters $r$ and $s$ in Vegetabile et al. (2021) (which here are set to be equal). Ignored for binary and multi-category treatments.
#'   }
#' }
#'
#' The arguments `maxit` and `reltol` can be supplied and are passed to the `control` argument of [optim()]. The `"BFGS"` method is used, so the defaults correspond to this.
#'
#' The `stabilize` argument is ignored; in the past it would reduce the variability of the weights through an iterative process. If you want to minimize the variance of the weights subject to balance constraints, use `method = "optweight"`.
#'
#' @section Additional Outputs:
#' \describe{
#'   \item{`obj`}{When `include.obj = TRUE`, the output of the call to [optim()], which contains the dual variables and convergence information. For ATE fits or with multi-category treatments, a list of `optim()` outputs, one for each weighted group.
#'   }
#' }
#'
#' @details
#' Entropy balancing involves the specification of an optimization problem, the solution to which is then used to compute the weights. The constraints of the primal optimization problem correspond to covariate balance on the means (for binary and multi-category treatments) or treatment-covariate covariances (for continuous treatments), positivity of the weights, and that the weights sum to a certain value. It turns out that the dual optimization problem is much easier to solve because it is over only as many variables as there are balance constraints rather than over the weights for each unit and it is unconstrained. Zhao and Percival (2017) found that entropy balancing for the ATT of a binary treatment actually involves the estimation of the coefficients of a logistic regression propensity score model but using a specialized loss function different from that optimized with maximum likelihood. Entropy balancing is doubly robust (for the ATT) in the sense that it is consistent either when the true propensity score model is a logistic regression of the treatment on the covariates or when the true outcome model for the control units is a linear regression of the outcome on the covariates, and it attains a semi-parametric efficiency bound when both are true. Entropy balancing will always yield exact mean balance on the included terms.
#'
#' @seealso
#' [weightit()], [weightitMSM()]
#'
#' [method_ipt] and [method_cbps] for inverse probability tilting and CBPS, which work similarly.
#'
#' @references
#' ## Binary Treatments
#'
#' ### `estimand = "ATT"`
#' Hainmueller, J. (2012). Entropy Balancing for Causal Effects: A Multivariate Reweighting Method to Produce Balanced Samples in Observational Studies. *Political Analysis*, 20(1), 25–46. \doi{10.1093/pan/mpr025}
#'
#' Zhao, Q., & Percival, D. (2017). Entropy balancing is doubly robust. Journal of Causal Inference, 5(1). \doi{10.1515/jci-2016-0010}
#'
#' ### `estimand = "ATE"`
#'
#' Källberg, D., & Waernbaum, I. (2023). Large Sample Properties of Entropy Balancing Estimators of Average Causal Effects. *Econometrics and Statistics*. \doi{10.1016/j.ecosta.2023.11.004}
#'
#' ## Continuous Treatments
#'
#' Tübbicke, S. (2022). Entropy Balancing for Continuous Treatments. *Journal of Econometric Methods*, 11(1), 71–89. \doi{10.1515/jem-2021-0002}
#'
#' Vegetabile, B. G., Griffin, B. A., Coffman, D. L., Cefalu, M., Robbins, M. W., & McCaffrey, D. F. (2021). Nonparametric estimation of population average dose-response curves using entropy balancing weights for continuous exposures. *Health Services and Outcomes Research Methodology*, 21(1), 69–110. \doi{10.1007/s10742-020-00236-2}
#'
#' @examples
#' data("lalonde", package = "cobalt")
#'
#' #Balancing covariates between treatment groups (binary)
#' (W1 <- weightit(treat ~ age + educ + married +
#'                   nodegree + re74, data = lalonde,
#'                 method = "ebal", estimand = "ATT"))
#' summary(W1)
#' cobalt::bal.tab(W1)
#'
#' #Balancing covariates with respect to race (multi-category)
#' (W2 <- weightit(race ~ age + educ + married +
#'                   nodegree + re74, data = lalonde,
#'                 method = "ebal", estimand = "ATE"))
#' summary(W2)
#' cobalt::bal.tab(W2)
#'
#' #Balancing covariates and squares with respect to
#' #re75 (continuous), maintaining 3 moments of the
#' #covariate and treatment distributions
#' (W3 <- weightit(re75 ~ age + educ + married +
#'                   nodegree + re74, data = lalonde,
#'                 method = "ebal", moments = 2,
#'                 d.moments = 3))
#' summary(W3)
#' cobalt::bal.tab(W3)
NULL

weightit2ebal <- function(covs, treat, s.weights, subset, estimand, focal,
                          stabilize, missing, moments, int, verbose, ...) {
  A <- list(...)

  covs <- covs[subset, , drop = FALSE]
  treat <- factor(treat[subset])
  s.weights <- s.weights[subset]

  if (missing == "ind") {
    covs <- add_missing_indicators(covs)
  }

  covs <- cbind(covs, .int_poly_f(covs, poly = moments, int = int, center = TRUE))

  covs <- cbind(covs, .quantile_f(covs, qu = A[["quantile"]], s.weights = s.weights,
                                 focal = focal, treat = treat))

  for (i in seq_col(covs)) covs[,i] <- .make_closer_to_1(covs[,i])

  colinear.covs.to.remove <- setdiff(colnames(covs), colnames(make_full_rank(covs)))
  covs <- covs[, colnames(covs) %nin% colinear.covs.to.remove, drop = FALSE]

  if (is_not_null(A[["base.weights"]])) A[["base.weight"]] <- A[["base.weights"]]
  if (is_null(A[["base.weight"]])) {
    bw <- rep(1, length(treat))
  }
  else {
    if (!is.numeric(A[["base.weight"]]) || length(A[["base.weight"]]) != length(treat)) {
      .err("the argument to `base.weight` must be a numeric vector with length equal to the number of units")
    }

    bw <- A[["base.weight"]]
  }

  reltol <- if (is_null(A$reltol)) (.Machine$double.eps) else A$reltol
  chk::chk_number(reltol)

  maxit <- if (is_null(A$maxit)) 1e4 else A$maxit
  chk::chk_count(maxit)

  eb <- function(C, s.weights_t, Q) {
    n <- nrow(C)

    W <- function(Z, Q, C) {
      drop(Q * exp(-C %*% Z))
    }

    objective.EB <- function(Z, Q, C) {
      log(sum(W(Z, Q, C)))
    }

    gradient.EB <- function(Z, Q, C) {
      w <- W(Z, Q, C)
      drop(-(w %*% C)/sum(w))
    }

    opt.out <- optim(par = rep(0, ncol(C)),
                     fn = objective.EB,
                     gr = gradient.EB,
                     method = "BFGS",
                     control = list(trace = 1,
                                    reltol = reltol,
                                    maxit = maxit),
                     C = C, Q = Q)

    w <- W(opt.out$par, Q, C)
    opt.out$gradient <- gradient.EB(opt.out$par, Q, C)

    if (opt.out$convergence == 1) {
      .wrn("the optimization failed to converge in the alotted number of iterations. Try increasing `maxit`")
    }
    else if (any(abs(opt.out$gradient) > 1e-3)) {
      .wrn("the estimated weights do not balance the covariates, indicating the optimization arrived at a degenerate solution. Try decreasing the number of variables supplied to the optimization")
    }

    if (sum(w) > n*.Machine$double.eps) w <- w*n/sum(w)

    list(Z = setNames(opt.out$par, colnames(C)),
         w = w/s.weights_t,
         opt.out = opt.out)
  }

  w <- rep(1, length(treat))
  sw0 <- check_if_zero(s.weights)

  if (estimand == "ATE") {
    groups_to_weight <- levels(treat)
    targets <- cobalt::col_w_mean(covs, s.weights = s.weights)
  }
  else {
    groups_to_weight <- setdiff(levels(treat), focal)
    targets <- cobalt::col_w_mean(covs, s.weights = s.weights, subset = treat == focal)
  }

  covs <- sweep(covs, 2, targets, check.margin = FALSE)

  fit.list <- make_list(groups_to_weight)
  for (i in groups_to_weight) {
    verbosely({
      fit.list[[i]] <- eb(C = covs[treat == i & !sw0,,drop = FALSE],
                          s.weights_t = s.weights[treat == i & !sw0],
                          Q = bw[treat == i & !sw0])
    }, verbose = verbose)

    w[treat == i & !sw0] <- fit.list[[i]]$w
  }

  Mparts <- list(
    psi_treat = function(Btreat, A, Xtreat, SW) {
      coef_ind <- setNames(lapply(seq_along(groups_to_weight), function(i) {
        (i - 1) * ncol(Xtreat) + seq_len(ncol(Xtreat))
      }), groups_to_weight)

      sw0 <- check_if_zero(SW)

      do.call("cbind", lapply(groups_to_weight, function(i) {
        m <- matrix(0, nrow = length(A), ncol = length(Btreat[coef_ind[[i]]]))

        C <- Xtreat[A == i & !sw0,,drop = FALSE]
        w <- drop(bw[A == i & !sw0] * exp(-C %*% Btreat[coef_ind[[i]]]))

        m[A == i,] <- -(w * C)/sum(w)
        m
      }))
    },
    wfun = function(Btreat, Xtreat, A) {
      coef_ind <- setNames(lapply(seq_along(groups_to_weight), function(i) {
        (i - 1) * ncol(Xtreat) + seq_len(ncol(Xtreat))
      }), groups_to_weight)

      sw0 <- check_if_zero(s.weights)
      w <- rep(1, length(A))

      for (i in groups_to_weight) {
        C <- Xtreat[A == i & !sw0,,drop = FALSE]
        n <- nrow(C)
        w[A == i & !sw0] <- drop(bw[A == i & !sw0] * exp(-C %*% Btreat[coef_ind[[i]]]))
        if (sum(w[A == i & !sw0]) > n * .Machine$double.eps) {
          w[A == i & !sw0] <- w[A == i & !sw0] * n / sum(w[A == i & !sw0])
        }
      }

      w[!sw0] <- w[!sw0] / s.weights[!sw0]

      w
    },
    Xtreat = covs,
    A = treat,
    btreat = unlist(grab(fit.list, "Z"))
  )

  list(w = w, fit.obj = grab(fit.list, "opt.out"),
       Mparts = Mparts)
}

weightit2ebal.multi <- weightit2ebal

weightit2ebal.cont <- function(covs, treat, s.weights, subset, missing, moments, int, verbose, ...) {
  A <- list(...)

  covs <- covs[subset, , drop = FALSE]
  treat <- treat[subset]
  s.weights <- s.weights[subset]

  if (missing == "ind") {
    covs <- add_missing_indicators(covs)
  }

  if (is_not_null(A[["base.weights"]])) A[["base.weight"]] <- A[["base.weights"]]
  if (is_null(A[["base.weight"]])) {
    bw <- rep(1, length(treat))
  }
  else {
    if (!is.numeric(A[["base.weight"]]) || length(A[["base.weight"]]) != length(treat)) {
      .err("the argument to `base.weight` must be a numeric vector with length equal to the number of units")
    }

    bw <- A[["base.weight"]]
  }

  bw <- bw[subset]

  s.weights <- s.weights / mean_fast(s.weights)

  reltol <- if (is_null(A$reltol)) (.Machine$double.eps) else A$reltol
  chk::chk_number(reltol)

  maxit <- if (is_null(A$maxit)) 1e4 else A$maxit
  chk::chk_count(maxit)

  d.moments <- max(if_null_then(A[["d.moments"]], 1), moments)
  chk::chk_count(d.moments)

  k <- ncol(covs)

  poly.covs <- .int_poly_f(covs, poly = moments)
  int.covs <- .int_poly_f(covs, int = int)

  treat <- .make_closer_to_1(treat)
  for (i in seq_col(poly.covs)) poly.covs[,i] <- .make_closer_to_1(poly.covs[,i])
  for (i in seq_col(int.covs)) int.covs[,i] <- .make_closer_to_1(int.covs[,i])
  if (d.moments == moments) {
    d.poly.covs <- poly.covs
  }
  else {
    d.poly.covs <- .int_poly_f(covs, poly = d.moments)
    for (i in seq_col(d.poly.covs)) d.poly.covs[,i] <- .make_closer_to_1(d.poly.covs[,i])
  }
  for (i in seq_col(covs)) covs[,i] <- .make_closer_to_1(covs[,i])

  covs <- cbind(covs, poly.covs, int.covs, d.poly.covs)
  # colinear.covs.to.remove <- colnames(covs)[colnames(covs) %nin% colnames(make_full_rank(covs))]
  # covs <- covs[, colnames(covs) %nin% colinear.covs.to.remove, drop = FALSE]

  t.mat <- matrix(treat, ncol = 1, dimnames = list(NULL, "treat"))
  if (d.moments > 1) t.mat <- cbind(t.mat, .int_poly_f(t.mat, poly = d.moments))

  treat_c <- sweep(t.mat, 2, cobalt::col_w_mean(t.mat, s.weights))
  covs_c <- sweep(covs, 2, cobalt::col_w_mean(covs, s.weights))

  covs.ind <- seq_len(k)
  poly.covs.ind <- k + seq_col(poly.covs)
  int.covs.ind <- k + ncol(poly.covs) + seq_col(int.covs)
  d.poly.covs.ind <- k + ncol(poly.covs) + ncol(int.covs) + seq_col(d.poly.covs)

  C <- cbind(treat_c, covs_c[, c(covs.ind, int.covs.ind, d.poly.covs.ind)],
             treat_c[,1] * covs_c[, c(covs.ind, int.covs.ind, poly.covs.ind)])

  colnames(C) <- c(paste(colnames(treat_c), "(mean)"),
                   paste(colnames(covs_c)[c(covs.ind, int.covs.ind, d.poly.covs.ind)], "(mean)"),
                   colnames(covs_c)[c(covs.ind, int.covs.ind, poly.covs.ind)])

  colinear.covs.to.remove <- setdiff(colnames(C), colnames(make_full_rank(C)))
  C <- C[, colnames(C) %nin% colinear.covs.to.remove, drop = FALSE]

  eb <- function(C, s.weights, Q) {
    n <- nrow(C)

    W <- function(Z, Q, C) {
      drop(Q * exp(-C %*% Z))
    }

    objective.EB <- function(Z, Q, C) {
      log(sum(W(Z, Q, C)))
    }

    gradient.EB <- function(Z, Q, C) {
      w <- W(Z, Q, C)
      drop(-(w %*% C)/sum(w))
    }

    opt.out <- optim(par = rep(0, ncol(C)),
                     fn = objective.EB,
                     gr = gradient.EB,
                     method = "BFGS",
                     control = list(trace = 0,
                                    reltol = reltol,
                                    maxit = maxit),
                     Q = Q, C = C)

    w <- W(opt.out$par, Q, C)
    opt.out$gradient <- gradient.EB(opt.out$par, Q, C)

    if (opt.out$convergence == 1) {
      .wrn("the optimization failed to converge in the alotted number of iterations. Try increasing `maxit`")
    }
    else if (any(abs(opt.out$gradient) > 1e-3)) {
      .wrn("the estimated weights do not balance the covariates, indicating the optimization arrived at a degenerate solution. Try decreasing the number of variables supplied to the optimization")
    }

    if (sum(w) > n*.Machine$double.eps) w <- w*n/sum(w)

    list(Z = setNames(opt.out$par, colnames(C)),
         w = w/s.weights,
         opt.out = opt.out)
  }

  w <- rep(1, length(treat))
  sw0 <- check_if_zero(s.weights)

  verbosely({
    fit <- eb(C[!sw0,, drop = FALSE], s.weights[!sw0], bw[!sw0])
  }, verbose = verbose)

  w[!sw0] <- fit$w

  Mparts <- list(
    psi_treat = function(Btreat, A, Xtreat, SW) {
      sw0 <- check_if_zero(SW)

      C <- Xtreat[!sw0,,drop = FALSE]
      w <- drop(bw[!sw0] * exp(-C %*% Btreat))

      -(w * C)/sum(w)

    },
    wfun = function(Btreat, Xtreat, A) {
      sw0 <- check_if_zero(s.weights)
      w <- rep(1, length(A))

      C <- Xtreat[!sw0,,drop = FALSE]
      n <- nrow(C)
      w[!sw0] <- drop(bw[!sw0] * exp(-C %*% Btreat))

      if (sum(w[!sw0]) > n * .Machine$double.eps) {
        w[!sw0] <- w[!sw0] * n / sum(w[!sw0])
      }

      w[!sw0] <- w[!sw0] / s.weights[!sw0]

      w
    },
    Xtreat = C,
    A = treat,
    btreat = fit$Z
  )

  list(w = w, fit.obj = fit$opt.out,
       Mparts = Mparts)
}

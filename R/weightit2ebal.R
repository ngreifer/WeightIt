#' Entropy Balancing
#' @name method_ebal
#' @aliases method_entropy
#' @usage NULL
#'
#' @description
#' This page explains the details of estimating weights using
#' entropy balancing by setting `method = "ebal"` in the call to [weightit()] or
#' [weightitMSM()]. This method can be used with binary, multi-category, and
#' continuous treatments.
#'
#' In general, this method relies on estimating weights by minimizing the
#' negative entropy of the weights subject to exact moment balancing
#' constraints. This method relies on code written for \pkg{WeightIt} using
#' [optim()].
#'
#' ## Binary Treatments
#'
#' For binary treatments, this method estimates the weights using `optim()`
#' using formulas described by Hainmueller (2012). The following estimands are
#' allowed: ATE, ATT, and ATC. When the ATE is requested, the optimization is
#' run twice, once for each treatment group.
#'
#' ## Multi-Category Treatments
#'
#' For multi-category treatments, this method estimates the weights using
#' `optim()`. The following estimands are allowed: ATE and ATT. When the ATE is
#' requested, the optimization is run once for each treatment group. When the ATT is
#' requested, the optimization is run once for each non-focal (i.e., control) group.
#'
#' ## Continuous Treatments
#'
#' For continuous treatments, this method estimates the weights using `optim()`
#' using formulas described by Tübbicke (2022) and Vegetabile et al. (2021).
#'
#' ## Longitudinal Treatments
#'
#' For longitudinal treatments, the weights are the product of the weights
#' estimated at each time point. This method is not guaranteed to yield exact
#' balance at each time point. **NOTE: the use of entropy balancing with
#' longitudinal treatments has not been validated and should not be done!**
#'
#' ## Sampling Weights
#'
#' Sampling weights are supported through `s.weights` in all scenarios.
#'
#' ## Missing Data
#'
#' In the presence of missing data, the following value(s) for `missing` are
#' allowed:
#'     \describe{
#'       \item{`"ind"` (default)}{
#'         First, for each variable with missingness, a new missingness indicator variable is created which takes the value 1 if the original covariate is `NA` and 0 otherwise. The missingness indicators are added to the model formula as main effects. The missing values in the covariates are then replaced with the covariate medians (this value is arbitrary and does not affect estimation). The weight estimation then proceeds with this new formula and set of covariates. The covariates output in the resulting `weightit` object will be the original covariates with the `NA`s.
#'       }
#'     }
#'
#' ## M-estimation
#'
#' M-estimation is supported for all treatment types and estimands except when `tols` is greater than 0. See [glm_weightit()] and
#' `vignette("estimating-effects")` for details.
#'
#' @section Additional Arguments:
#'
#' \describe{
#'   \item{`base.weights`}{a vector of base weights, one for each unit. These correspond to the base weights $q$ in Hainmueller (2012). The estimated weights minimize the Kullback entropy divergence from the base weights, defined as \eqn{\sum w \log(w/q)}, subject to exact balance constraints. These can be used to supply previously estimated weights so that the newly estimated weights retain the some of the properties of the original weights while ensuring the balance constraints are met. Sampling weights should not be passed to `base.weights` but can be included in a `weightit()` call that includes `s.weights`.}
#'   \item{`reltol`}{the relative tolerance for convergence of the optimization. Passed to the `control` argument of `optim()`. Default is `1e-10`.}
#'   \item{`maxit`}{the maximum number of iterations for convergence of the optimization. Passed to the `control` argument of `optim()`. Default is 1000 for binary and multi-category treatments and 10000 for continuous and longitudinal treatments.}
#'   \item{`solver`}{the solver to use to estimate the parameters. Allowable options include `"multiroot"` to use \pkgfun{rootSolve}{multiroot} and `"optim"` to use [stats::optim()]. `"multiroot"` is the default when \pkg{rootSolve} is installed, as it tends to be much faster and more accurate; otherwise, `"optim"` is the default and requires no dependencies. Regardless of `solver`, the output of `optim()` is returned when `include.obj = TRUE` (see below).}
#'   \item{`moments`}{`integer`; the highest power of each covariate to be balanced. For example, if `moments = 3`, each covariate, its square, and its cube will be balanced. Can also be a named vector with a value for each covariate (e.g., `moments = c(x1 = 2, x2 = 4)`). Values greater than 1 for categorical covariates are ignored. Default is 1 to balance covariate means.
#'     }
#'     \item{`int`}{`logical`; whether first-order interactions of the covariates are to be balanced. Default is `FALSE`.
#'     }
#'     \item{`quantile`}{a named list of quantiles (values between 0 and 1) for each continuous covariate, which are used to create additional variables that when balanced ensure balance on the corresponding quantile of the variable. For example, setting `quantile = list(x1 = c(.25, .5. , .75))` ensures the 25th, 50th, and 75th percentiles of `x1` in each treatment group will be balanced in the weighted sample. Can also be a single number (e.g., `.5`) or a vector (e.g., `c(.25, .5, .75)`) to request the same quantile(s) for all continuous covariates. Only allowed with binary and multi-category treatments.
#'     }
#'     \item{`tols`}{a number corresponding to the maximum allowed standardized mean difference (for binary and multi-category treatments) or treatment-covariate correlation (for continuous treatments) allowed. Default is 0 for exact balance.}
#'   \item{`d.moments`}{`integer`; with continuous treatments, the number of moments of the treatment and covariate distributions that are constrained to be the same in the weighted sample as in the original sample. For example, setting `d.moments = 3` ensures that the mean, variance, and skew of the treatment and covariates are the same in the weighted sample as in the unweighted sample. `d.moments` should be greater than or equal to `moments` and will be automatically set accordingly if not (or if not specified). Vegetabile et al. (2021) recommend setting `d.moments = 3`, even if `moments` is less than 3. This argument corresponds to the tuning parameters \eqn{r} and \eqn{s} in Vegetabile et al. (2021) (which here are set to be equal). Ignored for binary and multi-category treatments.}
#' }
#'
#' The `stabilize` argument is ignored; in the past it would reduce the
#' variability of the weights through an iterative process. If you want to
#' minimize the variance of the weights subject to balance constraints, use
#' `method = "optweight"`.
#'
#' @section Additional Outputs:
#' \describe{
#'   \item{`obj`}{When `include.obj = TRUE`, the output of the call to [optim()], which contains the dual variables and convergence information. For ATE fits or with multi-category treatments, a list of `optim()` outputs, one for each weighted group.}
#' }
#'
#' @details
#' Entropy balancing involves the specification of an optimization
#' problem, the solution to which is then used to compute the weights. The
#' constraints of the primal optimization problem correspond to covariate
#' balance on the means (for binary and multi-category treatments) or
#' treatment-covariate covariances (for continuous treatments), positivity of
#' the weights, and that the weights sum to a certain value. It turns out that
#' the dual optimization problem is much easier to solve because it is over only
#' as many variables as there are balance constraints rather than over the
#' weights for each unit, and it is unconstrained.
#'
#' Zhao and Percival (2017) found
#' that entropy balancing for the ATT of a binary treatment actually involves
#' the estimation of the coefficients of a logistic regression propensity score
#' model but using a specialized loss function different from that optimized
#' with maximum likelihood. Entropy balancing is doubly robust (for the ATT) in
#' the sense that it is consistent either when the true propensity score model
#' is a logistic regression of the treatment on the covariates or when the true
#' outcome model for the control units is a linear regression of the outcome on
#' the covariates, and it attains a semi-parametric efficiency bound when both
#' are true. Entropy balancing will always yield exact mean balance on the
#' included terms unless an imbalance tolerance is requested with `tols`.
#'
#' When `tols` is greater than 0, inexact balance is allowed on the covariates. This can improve precision while allowing a small amount of bias in. The optimization problem is an L1 regularization problem and is solved using the Fast Iterative Shrinkage-Thresholding Algorithm (FISTA).
#'
#' @seealso [weightit()], [weightitMSM()]
#'
#' [method_ipt] and [method_cbps] for inverse probability tilting and CBPS,
#' which work similarly. [method_optweight] for another implementation of entropy balancing (by setting `"norm = entropy"`).
#'
#' @references
#'
#' ## Binary Treatments
#'
#' ### `estimand = "ATT"`
#'
#' Hainmueller, J. (2012). Entropy Balancing for Causal Effects: A Multivariate Reweighting Method to Produce Balanced Samples in Observational Studies. *Political Analysis*, 20(1), 25–46. \doi{10.1093/pan/mpr025}
#'
#' Zhao, Q., & Percival, D. (2017). Entropy balancing is doubly robust. *Journal of Causal Inference*, 5(1). \doi{10.1515/jci-2016-0010}
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
#'
#' summary(W1)
#'
#' cobalt::bal.tab(W1)
#'
#' #Balancing covariates with respect to race (multi-category)
#' (W2 <- weightit(race ~ age + educ + married +
#'                   nodegree + re74, data = lalonde,
#'                 method = "ebal", estimand = "ATE"))
#'
#' summary(W2)
#'
#' cobalt::bal.tab(W2)
#'
#' #Balancing covariates and squares with respect to
#' #re75 (continuous), maintaining 3 moments of the
#' #covariate and treatment distributions
#' (W3 <- weightit(re75 ~ age + educ + married +
#'                   nodegree + re74, data = lalonde,
#'                 method = "ebal", moments = 2,
#'                 d.moments = 3))
#'
#' summary(W3)
#'
#' cobalt::bal.tab(W3, poly = 2,
#'                 stats = c("c", "m"))
#'
#' #Balancing covariates between treatment groups (binary),
#' #allowing for inexact balance
#' (W1b <- weightit(treat ~ age + educ + married +
#'                   nodegree + re74, data = lalonde,
#'                 method = "ebal", estimand = "ATT",
#'                 tols = .02))
#'
#' summary(W1b)
#'
#' cobalt::bal.tab(W1, weights = list(inexact = W1b))

NULL

weightit2ebal <- function(covs, treat, s.weights, subset, estimand, focal,
                          stabilize, missing, verbose, ...) {

  bw <- .process_b.weights(..., treat = treat)

  treat <- factor(treat[subset])
  covs <- covs[subset, , drop = FALSE]
  s.weights <- s.weights[subset]
  bw <- bw[subset]

  missing <- .process_missing2(missing, covs)

  if (missing == "ind") {
    covs <- add_missing_indicators(covs)
  }

  covs <- covs |>
    .apply_moments_int_quantile(moments = ...get("moments"),
                                int = ...get("int"),
                                quantile = ...get("quantile"),
                                s.weights = s.weights, focal = focal,
                                treat = treat) |>
    .make_covs_closer_to_1() |>
    .make_covs_full_rank()

  tols <- ...get("tols", 0)
  arg::arg_numeric(tols)
  arg::arg_length(tols, c(1L, ncol(covs)))
  tols <- abs(tols)

  reltol <- ...get("reltol", 1e-10)
  arg::arg_number(reltol)

  maxit <- ...get("maxit", 1e4L)
  arg::arg_count(maxit)

  solver <- ...get("solver", NULL)
  if (is_null(solver)) {
    solver <- {
      if (any(tols > 0)) "FISTA"
      else if (rlang::is_installed("rootSolve")) "multiroot"
      else "optim"
    }
  }
  else if (any(tols > 0)) {
    solver <- arg::match_arg(solver, c("optim", "FISTA"))
  }
  else {
    solver <- arg::match_arg(solver, c("optim", "multiroot"))

    if (solver == "multiroot") {
      rlang::check_installed("rootSolve")
    }
  }

  if (tols > 0) {
    sds <- rep.int(1, ncol(covs))

    bin.vars <- is_binary_col(covs)

    if (!all(bin.vars)) {
      if (estimand == "ATE") {
        sds[!bin.vars] <- sqrt(colMeans(do.call("rbind",
                                                lapply(levels(treat), function(t) {
                                                  col.w.v(covs[treat == t, !bin.vars, drop = FALSE], w = s.weights[treat == t])
                                                }))))
      }
      else {
        sds[!bin.vars] <- cobalt::col_w_sd(covs[, !bin.vars, drop = FALSE],
                                           s.weights = s.weights, subset = treat == focal)
      }
    }

    tols <- c(0, tols * sds)

    if (estimand == "ATE") {
      tols <- tols / 2
    }
  }

  sw0 <- check_if_zero(s.weights)

  covs <- cbind(`(Intercept)` = 1, covs)

  if (estimand == "ATE") {
    groups_to_weight <- levels(treat)
    N <- sum(!sw0)
    s.weights <- s.weights / mean(s.weights)
    targets <- cobalt::col_w_mean(covs, s.weights = s.weights)
  }
  else {
    groups_to_weight <- setdiff(levels(treat), focal)
    N <- sum(!sw0[treat == focal])
    s.weights <-  s.weights / mean(s.weights[treat == focal])
    targets <- cobalt::col_w_mean(covs, s.weights = s.weights, subset = treat == focal)
  }

  eb <- function(C, s.weights_t, Q, M) {
    W <- function(Z, Q, C) {
      Q * exp(-drop(C %*% Z))
    }

    objective.EB <- function(Z, S, Q, C) {
      sum(S * W(Z, Q, C)) / N + sum(Z * M) + sum(tols * abs(Z))
    }

    gradient.EB <- function(Z, S, Q, C) {
      w <- S * W(Z, Q, C)
      -drop(w %*% C) / N + M + tols * sign1(Z)
    }

    coef_start <- setNames(rep.int(0, ncol(C)),
                           colnames(C))

    if (solver == "FISTA") {
      opt.out <- .FISTA_ebal(coef_start, S = s.weights_t,
                            C = C, Q = Q, N = N,
                            M = M, tols = tols,
                            reltol = reltol, maxit = maxit)
    }
    else if (solver == "multiroot") {
      out <- suppressWarnings({
        try(verbosely({
          rootSolve::multiroot(f = gradient.EB,
                               start = coef_start,
                               S = s.weights_t, C = C, Q = Q,
                               rtol = reltol,
                               atol = reltol,
                               ctol = reltol)
        }, verbose = verbose), silent = TRUE)
      })

      if (!null_or_error(out) && utils::hasName(out, "root") &&
          utils::hasName(out, "estim.precis") &&
          is_number(out[["estim.precis"]]) &&
          out[["estim.precis"]] < 1e-5) {
        coef_start <- out$root
      }

      opt.out <- optim(par = coef_start,
                       fn = objective.EB,
                       method = "BFGS",
                       control = list(trace = 1,
                                      reltol = reltol,
                                      maxit = if (all(coef_start == 0)) maxit else 0),
                       S = s.weights_t, C = C, Q = Q)

      opt.out$gradient <- gradient.EB(opt.out$par, s.weights_t, Q, C)
    }
    else {
      opt.out <- optim(par = coef_start,
                       fn = objective.EB,
                       gr = gradient.EB,
                       method = "BFGS",
                       control = list(trace = 1,
                                      reltol = reltol,
                                      maxit = maxit),
                       S = s.weights_t, C = C, Q = Q)

      opt.out$gradient <- gradient.EB(opt.out$par, s.weights_t, Q, C)
    }

    w <- W(opt.out$par, Q, C)

    if (opt.out$convergence != 0) {
      arg::wrn("the optimization failed to converge in the alotted number of iterations. Try increasing {.arg maxit}")
    }
    else if (any(abs(opt.out$gradient) > tols + 1e-4)) {
      arg::wrn("the estimated weights do not balance the covariates, indicating the optimization arrived at a degenerate solution. Try decreasing the number of variables supplied to the optimization")
    }

    list(Z = setNames(opt.out$par, colnames(C)),
         w = w,
         opt.out = opt.out)
  }

  w <- rep_with(1, treat)
  fit.list <- make_list(groups_to_weight)

  for (i in groups_to_weight) {
    in_i <- which(treat == i & !sw0)

    verbosely({
      fit.list[[i]] <- eb(C = covs[in_i, , drop = FALSE],
                          s.weights_t = s.weights[in_i],
                          Q = bw[in_i], M = targets)
    }, verbose = verbose)

    w[in_i] <- fit.list[[i]]$w
  }

  beta <- unlist(grab(fit.list, "Z"))

  Mparts <- if (all(tols == 0)) {list(
    psi_treat = function(Btreat, Xtreat, A, SW) {
      coef_ind <- lapply(seq_along(groups_to_weight), function(i) {
        (i - 1L) * ncol(Xtreat) + seq_col(Xtreat)
      }) |>
        setNames(groups_to_weight)

      Atarget  <- {
        if (is_null(focal)) rep_with(1.0, A)
        else 1.0 * (A == focal)
      }

      w <- rep_with(1.0, A)

      for (i in groups_to_weight) {
        in_i <- which(A == i)
        w[in_i] <- bw[in_i] * exp(-drop(Xtreat[in_i, , drop = FALSE] %*% Btreat[coef_ind[[i]]]))
      }

      M <- do.call("cbind", lapply(groups_to_weight, function(i) {
        (SW * Xtreat * (Atarget - (A == i) * w))
      }))

      if (any(tols > 0)) {
        M <- M + tcrossprod(Atarget, tols * sign1(Btreat))
      }

      M
    },
    wfun = function(Btreat, Xtreat, A) {
      coef_ind <- lapply(seq_along(groups_to_weight), function(i) {
        (i - 1L) * ncol(Xtreat) + seq_col(Xtreat)
      }) |>
        setNames(groups_to_weight)

      w <- rep_with(1.0, A)

      for (i in groups_to_weight) {
        in_i <- which(A == i)
        w[in_i] <- bw[in_i] * exp(-drop(Xtreat[in_i, , drop = FALSE] %*% Btreat[coef_ind[[i]]]))
      }

      w
    },
    Xtreat = covs,
    A = treat,
    btreat = beta
  )}

  list(w = w,
       fit.obj = fit.list,
       Mparts = Mparts)
}

weightit2ebal.multi <- weightit2ebal

weightit2ebal.cont <- function(covs, treat, s.weights, subset, missing, verbose, ...) {

  bw <- .process_b.weights(..., treat = treat)

  treat <- treat[subset]
  covs <- covs[subset, , drop = FALSE]
  s.weights <- s.weights[subset]
  bw <- bw[subset]

  missing <- .process_missing2(missing, covs)

  if (missing == "ind") {
    covs <- add_missing_indicators(covs)
  }

  tols <- ...get("tols", 0)
  arg::arg_numeric(tols)
  arg::arg_length(tols, c(1L, ncol(covs)))
  tols <- abs(tols)

  reltol <- ...get("reltol", 1e-10)
  arg::arg_number(reltol)

  maxit <- ...get("maxit", 1e4)
  arg::arg_count(maxit)

  solver <- ...get("solver", NULL)
  if (is_null(solver)) {
    solver <- {
      if (any(tols > 0)) "FISTA"
      else if (rlang::is_installed("rootSolve")) "multiroot"
      else "optim"
    }
  }
  else if (any(tols > 0)) {
    solver <- arg::match_arg(solver, c("optim", "FISTA"))
  }
  else {
    solver <- arg::match_arg(solver, c("optim", "multiroot"))

    if (solver == "multiroot") {
      rlang::check_installed("rootSolve")
    }
  }

  moments <- ...get("moments", 1L)

  d.moments <- ...get("d.moments", 1L)
  arg::arg_count(d.moments)

  treat <- .make_closer_to_1(treat)

  t.mat <- matrix(treat, ncol = 1L, dimnames = list(NULL, "treat")) |>
    .apply_moments_int_quantile(moments = d.moments) |>
    scale_w(s.weights)

  bal.covs <- covs |>
    .apply_moments_int_quantile(moments = moments,
                                int = ...get("int")) |>
    scale_w(s.weights) |>
    .make_covs_full_rank()

  if (all(d.moments <= moments)) {
    C <- cbind(t.mat,
               bal.covs,
               t.mat[, 1L] * bal.covs)

    colnames(C) <- c(paste(colnames(t.mat), "(mean)"),
                     paste(colnames(bal.covs), "(mean)"),
                     colnames(bal.covs))
  }
  else {
    d.covs <- covs |>
      .apply_moments_int_quantile(moments = pmax(d.moments, moments),
                                  int = ...get("int")) |>
      scale_w(s.weights) |>
      .make_covs_full_rank()

    C <- cbind(t.mat,
               d.covs,
               t.mat[, 1L] * bal.covs)

    colnames(C) <- c(paste(colnames(t.mat), "(mean)"),
                     paste(colnames(d.covs), "(mean)"),
                     colnames(bal.covs))
  }

  dist_vars <- colnames(C)[seq_len(ncol(C) - ncol(bal.covs))]

  sw0 <- check_if_zero(s.weights)

  N <- sum(!sw0)

  s.weights <- s.weights / mean(s.weights)

  C <- cbind(`(Intercept)` = 1, C)
  dist_vars <- c(colnames(C)[1L], dist_vars)

  colinear.covs.to.remove <- setdiff(colnames(C),
                                     colnames(make_full_rank(C, with.intercept = FALSE)))

  C <- C[, colnames(C) %nin% colinear.covs.to.remove, drop = FALSE]
  dist_vars <- setdiff(dist_vars, colinear.covs.to.remove)
  dist_ind <- match(dist_vars, colnames(C))

  targets <- rep.int(0, ncol(C))
  targets[dist_ind] <- cobalt::col_w_mean(C[, dist_ind, drop = FALSE], s.weights = s.weights)

  tols_ <- rep.int(0, ncol(C))
  tols_[-dist_ind] <- tols * (1 - 1 / N)
  tols <- tols_

  eb <- function(C, s.weights_t, Q, M) {
    W <- function(Z, Q, C) {
      Q * exp(-drop(C %*% Z))
    }

    objective.EB <- function(Z, S, Q, C) {
      sum(S * W(Z, Q, C)) / N + sum(Z * M) + sum(tols * abs(Z))
    }

    gradient.EB <- function(Z, S, Q, C) {
      w <- S * W(Z, Q, C)
      -drop(w %*% C) / N + M + tols * sign1(Z)
    }

    coef_start <- setNames(rep.int(0, ncol(C)),
                           colnames(C))

    if (solver == "FISTA") {
      opt.out <- .FISTA_ebal(coef_start, S = s.weights_t,
                            C = C, Q = Q, N = N,
                            M = M, tols = tols,
                            reltol = reltol, maxit = maxit)
    }
    else if (solver == "multiroot") {
      out <- suppressWarnings({
        try(verbosely({
          rootSolve::multiroot(f = gradient.EB,
                               start = coef_start,
                               S = s.weights_t, C = C, Q = Q,
                               rtol = reltol,
                               atol = reltol,
                               ctol = reltol)
        }, verbose = verbose), silent = TRUE)
      })

      if (!null_or_error(out) && utils::hasName(out, "root") &&
          utils::hasName(out, "estim.precis") &&
          is_number(out[["estim.precis"]]) &&
          out[["estim.precis"]] < 1e-5) {
        coef_start <- out$root
      }

      opt.out <- optim(par = coef_start,
                       fn = objective.EB,
                       method = "BFGS",
                       control = list(trace = 1,
                                      reltol = reltol,
                                      maxit = if (all(coef_start == 0)) maxit else 0),
                       S = s.weights_t, C = C, Q = Q)

      opt.out$gradient <- gradient.EB(opt.out$par, s.weights_t, Q, C)
    }
    else {
      opt.out <- optim(par = coef_start,
                       fn = objective.EB,
                       gr = gradient.EB,
                       method = "BFGS",
                       control = list(trace = 1,
                                      reltol = reltol,
                                      maxit = maxit),
                       S = s.weights_t, C = C, Q = Q)

      opt.out$gradient <- gradient.EB(opt.out$par, s.weights_t, Q, C)
    }

    w <- W(opt.out$par, Q, C)

    if (opt.out$convergence != 0) {
      arg::wrn("the optimization failed to converge in the alotted number of iterations. Try increasing {.arg maxit}")
    }
    else if (any(abs(opt.out$gradient) > tols + 1e-4)) {
      arg::wrn("the estimated weights do not balance the covariates, indicating the optimization arrived at a degenerate solution. Try decreasing the number of variables supplied to the optimization")
    }

    list(Z = setNames(opt.out$par, colnames(C)),
         w = w,
         opt.out = opt.out)
  }

  w <- rep_with(1, treat)

  verbosely({
    fit <- eb(C[!sw0, , drop = FALSE], s.weights[!sw0], bw[!sw0],
              M = targets)
  }, verbose = verbose)

  w[!sw0] <- fit$w

  Mparts <- if (all(tols == 0)) {list(
    psi_treat = function(Btreat, Xtreat, A, SW) {
      w <- bw * exp(-drop(Xtreat %*% Btreat))

      cbind(
        SW * (w - 1) * Xtreat[, dist_ind, drop = FALSE], #weighted mean
        SW * w * Xtreat[, -dist_ind, drop = FALSE] #weighted covariance
      )
    },
    wfun = function(Btreat, Xtreat, A) {
      bw * exp(-drop(Xtreat %*% Btreat))
    },
    dw_dBtreat = function(Btreat, Xtreat, A, SW) {
      w <- bw * exp(-drop(Xtreat %*% Btreat))

      -w * Xtreat
    },
    Xtreat = C,
    A = treat,
    btreat = fit$Z
  )}

  list(w = w, fit.obj = fit$opt.out,
       Mparts = Mparts)
}

.FISTA_ebal <- function(start, S, Q, C, N, M, tols, reltol = 1e-10, maxit = 10000) {

  W <- function(Z) {
    S * Q * exp(-drop(C %*% Z))
  }

  loss_EB_smooth <- function(Z, w = NULL) {
    if (is_null(w)) w <- W(Z)
    sum(w) / N + sum(Z * M)
  }

  grad_EB_smooth <- function(Z, w = NULL) {
    if (is_null(w)) w <- W(Z)
    -drop(w %*% C) / N + M
  }

  # Gradient of smooth part only (no L1 term)

  # Soft threshold operator
  soft_thresh <- function(Z, lambda) {
    sign(Z) * squish(abs(Z) - lambda, lo = 0, hi = Inf)
  }

  # Estimate Lipschitz constant for step size if not provided
  # L <= S_max * ||C||^2 / N  (crude but safe bound)
  L <- 1
  Z <- start
  Y <- Z          # FISTA momentum term
  t_old <- 1

  w <- W(Y)
  lossY <- loss_EB_smooth(Y, w)
  gY    <- grad_EB_smooth(Y, w)

  loss_new <- lossY + sum(tols * abs(Z))

  for (i in seq_len(maxit)) {
    Z_old <- Z
    loss_old <- loss_new

    # Backtracking line search
    for (j in seq_len(10L)) {
      Z_new <- soft_thresh(Y - gY / L, tols / L)
      diff  <- Z_new - Y
      w <- W(Z_new)

      if (loss_EB_smooth(Z_new, w) <= lossY + sum(gY * diff) + (L / 2) * sum(diff^2)) {
        break
      }

      L <- L * 2
    }

    Z <- Z_new

    # FISTA momentum update
    t_new <- (1 + sqrt(1 + 4 * t_old^2)) / 2
    Y     <- Z + ((t_old - 1) / t_new) * (Z - Z_old)
    t_old <- t_new

    # Recompute at new Y for next iteration
    w <- W(Y)
    lossY <- loss_EB_smooth(Y, w)

    # Termination check on full loss (including L1)
    loss_new <- lossY + sum(tols * abs(Z))

    if (i > 1 && abs(loss_new - loss_old) < reltol * (abs(loss_new) + reltol)) {
      break
    }

    gY <- grad_EB_smooth(Y, w)

    loss_old <- loss_new
  }

  w <- W(Z)

  list(
    par = Z,
    value = loss_EB_smooth(Z, w) + sum(abs(Z) * tols),
    counts = i,
    convergence = as.integer(i == maxit),
    gradient = grad_EB_smooth(Z, w) + sign1(Z) * tols
  )
}
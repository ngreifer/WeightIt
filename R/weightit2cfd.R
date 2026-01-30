#' Characteristic Function Distance Balancing
#' @name method_cfd
#' @usage NULL
#'
#' @description
#' This page explains the details of estimating weights using
#' characteristic function distance (CFD) balancing by setting `method = "cfd"` in the call to [weightit()]
#' or [weightitMSM()]. This method can be used with binary, multi-category, and
#' continuous treatments.
#'
#' In general, this method relies on estimating weights by minimizing a scalar measure of covariate balance, the CFD. The CFD is related to the maximum mean discrepancy and captures covariate balance for the joint covariate distribution as determined by a specific choice of kernel. This
#' method relies on code written for \pkg{WeightIt} using \pkgfun{osqp}{osqp}
#' from the \CRANpkg{osqp} package to perform the optimization. This method may
#' be slow or memory-intensive for large datasets.
#'
#' ## Binary Treatments
#'
#' For binary treatments, this method estimates the weights using `osqp()` using
#' formulas described by Santra, Chen, and Park (2026). The following estimands are
#' allowed: ATE, ATT, and ATC.
#'
#' ## Multi-Category Treatments
#'
#' For multi-category treatments, this method estimates the weights using
#' `osqp()` using formulas described by Santra, Chen, and Park (2026). The following
#' estimands are allowed: ATE and ATT.
#'
#' ## Continuous Treatments
#'
#' CFD balancing is not compatible with continuous treatments.
#'
#' ## Longitudinal Treatments
#'
#' For longitudinal treatments, the weights are the product of the weights
#' estimated at each time point. This method is not guaranteed to yield optimal
#' balance at each time point. **NOTE: the use of CFD balancing with longitudinal treatments has not been validated!**
#'
#' ## Sampling Weights
#'
#' Sampling weights are supported through `s.weights` in all scenarios. In some
#' cases, sampling weights will cause the optimization to fail due to lack of
#' convexity or infeasible constraints.
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
#' M-estimation is not supported.
#'
#' @section Additional Arguments:
#'
#' The following following additional arguments can be specified:
#'
#' \describe{
#'   \item{`kernel`}{the name of the kernel used to characterize the CFD. Allowable optiosn include `"gaussian"` for the multivariate Gaussian kernel (the default), `"matern"` for the multivariate Matern kernel, `"energy"` for the energy distance kernel, `"laplace"` for the univariate Laplacian kernel, and `"t"` for the univariate t-dsitribution kernel.}
#'   \item{`nu`}{for `kernel = "matern"`, the \eqn{\nu} parameter used to control smoothness. The default value is 3/2. For any values other than 1/2, 3/2, and 5/2, the \pkg{GPBayes} package is required to compute the Matern kernel. For `kernel = "t"`, the degrees of freedom for the univariate t-distributions used in the kernel. The default value is 5. Ignored for other kernels.}
#'   \item{`nsim`}{for `kernel = "t"`, the number of simulations to use to compute the t-distribution kernel. Default is 5000. Greater is better but takes longer and uses more memory.}
#'   \item{`lambda`}{a positive numeric scalar used to penalize the square of the weights. This value divided by the square of the total sample size is added to the diagonal of the quadratic part of the loss function. Higher values favor weights with less variability. Default is .0001, which is essentially 0.}
#'   \item{`moments`}{`integer`; the highest power of each covariate to be balanced. For example, if `moments = 3`, each covariate, its square, and its cube will be balanced. Can also be a named vector with a value for each covariate (e.g., `moments = c(x1 = 2, x2 = 4)`). Values greater than 1 for categorical covariates are ignored. Default is 0 to impose no constraint on balance.}
#'   \item{`int`}{`logical`; whether first-order interactions of the covariates are to be balanced. Default is `FALSE`.}
#'   \item{`tols`}{when `moments` is positive, a number corresponding to the maximum allowed standardized mean difference (for binary and multi-category treatments) or treatment-covariate correlation (for continuous treatments) allowed. Default is 0. Ignored when `moments = 0`.}
#'   \item{`min.w`}{the minimum allowable weight. Negative values (including `-Inf`) are allowed. Default is `1e-8`.}
#' }
#'
#' For binary and multi-category treatments, the following additional arguments can be specified:
#'   \describe{
#'     \item{`improved`}{`logical`; whether to include an additional term in the CFD objective function to minimize the distance between pairs of groups when `estimand = "ATE"`. Default is `TRUE`.}
#'     \item{`quantile`}{a named list of quantiles (values between 0 and 1) for each continuous covariate, which are used to create additional variables that when balanced ensure balance on the corresponding quantile of the variable. For example, setting `quantile = list(x1 = c(.25, .5. , .75))` ensures the 25th, 50th, and 75th percentiles of `x1` in each treatment group will be balanced in the weighted sample. Can also be a single number (e.g., `.5`) or a vector (e.g., `c(.25, .5, .75)`) to request the same quantile(s) for all continuous covariates.
#'     }
#'   }
#'
#' The `moments` argument functions differently for `method = "cfd"` from
#' how it does with some other methods. When unspecified or set to zero, CFD
#' balancing weights are estimated as described by Santra et al. (2026) for
#' binary and multi-category treatments. When `moments` is set to an integer larger than 0,
#' additional balance constraints on the requested moments of the covariates
#' are also included, guaranteeing exact moment balance on these covariates
#' while minimizing the CFD of the weighted sample. This involves exact balance on the means of the
#' entered covariates. The constraint on exact balance can be relaxed using the `tols` argument.
#'
#' Any other arguments will be passed to \pkgfun{osqp}{osqpSettings}. Some defaults differ from those in `osqpSettings()`; see *Reproducibility* section.
#'
#' @section Additional Outputs:
#'
#' \describe{
#'   \item{`obj`}{When `include.obj = TRUE`, the output of the call to \pkgfun{osqp}{solve_osqp}, which contains the dual variables and convergence information.
#'   }
#' }
#'
#' @details
#' CFD balancing is a method of estimating weights using
#' optimization without a propensity score. The weights are the solution to a
#' constrained quadratic optimization problem where the objective function
#' concerns covariate balance as measured by the CFD between groups.
#'
#' CFD balancing for binary and multi-category treatments involves minimizing
#' the CFD between the treatment groups and between each treatment
#' group and a target group (e.g., the full sample for the ATE). The CFD is a scalar measure of the difference between two multivariate
#' distributions. The performance of CFD balance depends on the choice of kernel, controlled by the `kernel` argument. Each kernel corresponds to different assumptions about the form of the true outcome model. See Santra et al. (2026) for a comparison of these different kernels. Setting `kernel = "energy"` is equivalent to entropy balancing, which can also be requested by using [`method = "energy"`][method_energy].
#'
#' The primary benefit of CFD balancing is that all features of the covariate
#' distribution are balanced, not just means, as with other optimization-based
#' methods like entropy balancing. Still, it is possible to add additional
#' balance constraints to require balance on individual terms using the
#' `moments` argument, just like with entropy balancing. CFD balancing can
#' sometimes yield weights with high variability; the `lambda` argument can be
#' supplied to penalize highly variable weights to increase the effective sample
#' size at the expense of balance.
#'
#' ## Reproducibility
#'
#' Although there are no stochastic components to the optimization, a feature
#' turned off by default is to update the optimization based on how long the
#' optimization has been running, which will vary across runs even when a seed
#' is set and no parameters have been changed. See the discussion
#' [here](https://github.com/osqp/osqp-r/issues/19) for more details. To ensure
#' reproducibility by default, `adaptive_rho_interval` is set to 10. See
#' \pkgfun{osqp}{osqpSettings} for details.
#'
#' @note
#' Sometimes the optimization can fail to converge because the problem is
#' not convex. A warning will be displayed if so. In these cases, try simply
#' re-fitting the weights without changing anything (but see the
#' *Reproducibility* section above). If the method repeatedly fails, you should
#' try another method or change the supplied parameters (though this is
#' uncommon). Increasing `max_iter` or changing `adaptive_rho_interval` might
#' help.
#'
#' If it seems like the weights are balancing the covariates but you still get a
#' failure to converge, this usually indicates that more iterations are needs to
#' find the optimal solutions. This can occur when `moments` or `int` are
#' specified. `max_iter` should be increased, and setting `verbose = TRUE`
#' allows you to monitor the process and examine if the optimization is
#' approaching convergence.
#'
#' If `min.w` is positive and you still get a warning about the presence of negative weights, try setting `eps` to a smaller number (e.g., to `1e-12`).
#'
#' As of version 1.5.0, `polish` is now set to `TRUE` by default. This should yield slightly improved solutions but may be a little slower.
#'
#' @seealso [weightit()], [weightitMSM()]
#'
#' @references
#' Santra, D., Chen, G., & Park, C. (2026). Distributional Balancing for Causal Inference: A Unified Framework via Characteristic Function Distance (arXiv:2601.15449). arXiv. \doi{10.48550/arXiv.2601.15449}
#'
#' @examplesIf rlang::is_installed("osqp")
#' data("lalonde", package = "cobalt")
#'
#' #Balancing covariates between treatment groups (binary)
#' (W1 <- weightit(treat ~ age + educ + married +
#'                   nodegree + re74, data = lalonde,
#'                 method = "cfd", estimand = "ATE"))
#'
#' summary(W1)
#'
#' cobalt::bal.tab(W1)
#'
#' #Using a different kernel:
#' (W1b <- weightit(treat ~ age + educ + married +
#'                   nodegree + re74, data = lalonde,
#'                 method = "cfd", estimand = "ATE",
#'                 kernel = "matern", nu = 5/2))
#'
#' summary(W1b)
#'
#' cobalt::bal.tab(W1b)
#'
#' #Balancing covariates with respect to race (multi-category)
#' (W2 <- weightit(race ~ age + educ + married +
#'                   nodegree + re74, data = lalonde,
#'                 method = "cfd", estimand = "ATT",
#'                 focal = "black"))
#'
#' summary(W2)
#'
#' cobalt::bal.tab(W2)
NULL

weightit2cfd <- function(covs, treat, s.weights, subset, estimand, focal,
                         missing, verbose, ...) {

  missing <- .process_missing2(missing, covs)

  if (missing == "ind") {
    covs <- add_missing_indicators(covs)
  }

  d <- .process_kernel(covs, s.weights, subset,
                       kernel = ...get("kernel", "gaussian"),
                       bw_scale = ...get("bw_scale", 1),
                       nu = ...get("nu"),
                       nsim = ...get("nsim", 5000L))

  d <- unname(d[subset, subset, drop = FALSE])

  covs <- covs[subset, , drop = FALSE]
  treat <- treat[subset]
  s.weights <- s.weights[subset]

  t.lev <- get_treated_level(treat, estimand, focal)
  treat <- binarize(treat, one = t.lev)

  n <- length(treat)

  sw0 <- check_if_zero(s.weights)

  diagn <- diag(n)

  min.w <- ...get("min.w", 1e-8)
  chk::chk_number(min.w)

  lambda <- ...get("lambda", 1e-4)
  chk::chk_number(lambda)

  moments <- ...get("moments", 0)
  int <- isTRUE(...get("int", FALSE))
  quantile <- ...get("quantile")
  add_constraints <- any(moments > 0) || int || is_not_null(quantile)

  if (add_constraints) {
    covs <- scale(covs)

    tols <- ...get("tols", 0)
    chk::chk_number(tols)
    tols <- abs(tols)
  }

  A <- ...mget(names(formals(osqp::osqpSettings)))

  eps <- ...get("eps", squish(min.w, lo = 1e-12, hi = 1e-8))
  if (is_not_null(eps)) {
    chk::chk_number(eps)
    if (is_null(A[["eps_abs"]])) A[["eps_abs"]] <- eps
    if (is_null(A[["eps_rel"]])) A[["eps_rel"]] <- eps
  }

  if (is_null(A[["max_iter"]])) A[["max_iter"]] <- 5e4L
  chk::chk_count(A[["max_iter"]], "`max_iter`")
  chk::chk_lt(A[["max_iter"]], Inf, "`max_iter`")
  if (is_null(A[["eps_abs"]])) A[["eps_abs"]] <- 1e-8
  chk::chk_number(A[["eps_abs"]], "`eps_abs`")
  if (is_null(A[["eps_rel"]])) A[["eps_rel"]] <- 1e-6
  chk::chk_number(A[["eps_rel"]], "`eps_rel`")
  if (is_null(A[["time_limit"]])) A[["time_limit"]] <- 0
  chk::chk_number(A[["time_limit"]], "`time_limit`")
  if (is_null(A[["adaptive_rho_interval"]])) A[["adaptive_rho_interval"]] <- 10L
  chk::chk_count(A[["adaptive_rho_interval"]], "`adaptive_rho_interval`")
  if (is_null(A[["polish"]])) A[["polish"]] <- TRUE
  chk::chk_flag(A[["polish"]], "`polish`")
  if (is_null(A[["polish_refine_iter"]])) A[["polish_refine_iter"]] <- 20L
  chk::chk_count(A[["polish_refine_iter"]], "`polish_refine_iter`")
  A[["verbose"]] <- verbose

  options.list <- do.call(osqp::osqpSettings, A)

  t0 <- which(treat == 0)
  t1 <- which(treat == 1)

  s.weights[t0] <- s.weights[t0] / mean_fast(s.weights[t0])
  s.weights[t1] <- s.weights[t1] / mean_fast(s.weights[t1])

  n0 <- length(t0)
  n1 <- length(t1)

  s.weights_n_0 <- s.weights_n_1 <- rep.int(0, n)
  s.weights_n_0[t0] <- s.weights[t0] / n0
  s.weights_n_1[t1] <- s.weights[t1] / n1

  if (estimand == "ATE") {
    improved <- ...get("improved", TRUE)
    chk::chk_flag(improved)

    nn <- tcrossprod(cbind(s.weights_n_0, s.weights_n_1))

    if (improved) {
      nn <- nn + tcrossprod(s.weights_n_0 - s.weights_n_1)
    }

    P <- d * nn

    q <- -((s.weights * 2 / n) %*% d) * (s.weights_n_0 + s.weights_n_1)

    #Constraints for positivity and sum of weights
    Amat <- cbind(diagn, s.weights_n_0, s.weights_n_1)
    lvec <- c(ifelse(sw0, 1, min.w), 1, 1)
    uvec <- c(ifelse(sw0, 1, Inf), 1, 1)

    unbounded <- lvec == -Inf & uvec == Inf

    if (any(unbounded)) {
      Amat <- Amat[, !unbounded, drop = FALSE]
      lvec <- lvec[!unbounded]
      uvec <- uvec[!unbounded]
    }

    if (add_constraints) {
      #Exactly balance moments, interactions, and/or quantiles
      covs <- .apply_moments_int_quantile(covs,
                                          moments = moments,
                                          int = int,
                                          quantile = quantile,
                                          s.weights = s.weights)

      targets <- col.w.m(covs, s.weights)

      sds <- rep.int(1, ncol(covs))

      if (tols > 0) {
        bin.vars <- is_binary_col(covs)

        if (!all(bin.vars)) {
          sds[!bin.vars] <- sqrt(colMeans(rbind(col.w.v(covs[t1, !bin.vars, drop = FALSE], w = s.weights[t1]),
                                                col.w.v(covs[t0, !bin.vars, drop = FALSE], w = s.weights[t0]))))
        }
      }

      Amat <- cbind(Amat, covs * s.weights_n_0, covs * s.weights_n_1)
      lvec <- c(lvec, targets - sds * tols / 2, targets - sds * tols / 2)
      uvec <- c(uvec, targets + sds * tols / 2, targets + sds * tols / 2)
    }
  }
  else if (estimand == "ATT") {
    nn <- tcrossprod(s.weights_n_0[t0])

    P <- d[t0, t0, drop = FALSE] * nn

    q <- -2 * (s.weights_n_1[t1] %*% d[t1, t0, drop = FALSE]) * s.weights_n_0[t0]

    Amat <- cbind(diag(n0), s.weights_n_0[t0])
    lvec <- c(ifelse(sw0[t0], 1, min.w), 1)
    uvec <- c(ifelse(sw0[t0], 1, Inf), 1)

    unbounded <- lvec == -Inf & uvec == Inf

    if (any(unbounded)) {
      Amat <- Amat[, !unbounded, drop = FALSE]
      lvec <- lvec[!unbounded]
      uvec <- uvec[!unbounded]
    }

    if (add_constraints) {
      #Exactly balance moments, interactions, and/or quantiles
      covs <- .apply_moments_int_quantile(covs,
                                          moments = moments,
                                          int = int,
                                          quantile = quantile,
                                          s.weights = s.weights, focal = 1,
                                          treat = treat)

      targets <- col.w.m(covs[t1, , drop = FALSE], s.weights[t1])

      sds <- rep.int(1, ncol(covs))

      if (tols > 0) {
        bin.vars <- is_binary_col(covs)

        if (!all(bin.vars)) {
          sds[!bin.vars] <- sqrt(col.w.v(covs[t1, !bin.vars, drop = FALSE], w = s.weights[t1]))
        }
      }

      Amat <- cbind(Amat, covs[t0, , drop = FALSE] * s.weights_n_0[t0])
      lvec <- c(lvec, targets - sds * tols)
      uvec <- c(uvec, targets + sds * tols)
    }
  }
  else if (estimand == "ATC") {
    nn <- tcrossprod(s.weights_n_1[t1])

    P <- d[t1, t1, drop = FALSE] * nn

    q <- -2 * (s.weights_n_0[t0] %*% d[t0, t1, drop = FALSE]) * s.weights_n_1[t1]

    Amat <- cbind(diag(n1), s.weights_n_1[t1])
    lvec <- c(ifelse(sw0[t1], 1, min.w), 1)
    uvec <- c(ifelse(sw0[t1], 1, Inf), 1)

    unbounded <- lvec == -Inf & uvec == Inf

    if (any(unbounded)) {
      Amat <- Amat[, !unbounded, drop = FALSE]
      lvec <- lvec[!unbounded]
      uvec <- uvec[!unbounded]
    }

    if (add_constraints) {
      #Exactly balance moments, interactions, and/or quantiles
      covs <- .apply_moments_int_quantile(covs,
                                          moments = moments,
                                          int = int,
                                          quantile = quantile,
                                          s.weights = s.weights, focal = 0,
                                          treat = treat)

      targets <- col.w.m(covs[t0, , drop = FALSE], s.weights[t0])

      if (tols > 0) {
        bin.vars <- is_binary_col(covs)

        if (!all(bin.vars)) {
          sds[!bin.vars] <- sqrt(col.w.v(covs[t0, !bin.vars, drop = FALSE], w = s.weights[t0]))
        }
      }

      Amat <- cbind(Amat, covs[t1, , drop = FALSE] * s.weights_n_1[t1])
      lvec <- c(lvec, targets - sds * tols)
      uvec <- c(uvec, targets + sds * tols)
    }
  }

  #Add weight penalty
  if (lambda < 0) {
    # es <- eigen(P, symmetric = TRUE)
    # esv <- es$values
    #
    # tol <- nrow(P) * max(abs(esv)) * .Machine$double.eps
    #
    # if (!all(esv > tol)) {
    #   tau <- pmax(0, 2 * tol - esv)
    #   P <- P + tcrossprod(es$vectors %*% diag(sqrt(tau), nrow(P)))
    # }

    #Find lambda to make P PSD
    e <- eigen(P, symmetric = TRUE, only.values = TRUE)
    e.min <- min(e$values)

    if (e.min < 0) {
      diag(P) <- diag(P) - e.min + .Machine$double.eps
    }
  }
  else if (lambda != 0) {
    diag(P) <- switch(estimand,
                      ATE = diag(P) + lambda * (s.weights_n_0 + s.weights_n_1)^2 / 2,
                      ATT = diag(P) + lambda * s.weights_n_0[t0]^2 / 2,
                      ATC = diag(P) + lambda * s.weights_n_1[t1]^2 / 2)
  }


  verbosely({
    opt.out <- osqp::solve_osqp(P = 2 * P, q = q,
                                A = t(Amat), l = lvec, u = uvec,
                                pars = options.list)
  }, verbose = verbose)

  if (identical(opt.out$info$status, "maximum iterations reached")) {
    .wrn(sprintf("the optimization failed to converge. Try increasing `max_iter` (current value: %s)",
                 A[["max_iter"]]))
  }
  else if (identical(opt.out$info$status, "run time limit reached")) {
    .wrn(sprintf("the optimization failed to converge. Try increasing `time_limit` (current value: %s)",
                 A[["time_limit"]]))
  }
  else if (!startsWith(opt.out$info$status, "solved")) {
    .wrn("no feasible solution could be found that satisfies all constraints. Relax any constraints supplied")
  }

  if (estimand == "ATT") {
    w <- rep.int(1, n)
    w[t0] <- opt.out$x
  }
  else if (estimand == "ATC") {
    w <- rep.int(1, n)
    w[t1] <- opt.out$x
  }
  else {
    w <- opt.out$x
  }

  # Shrink tiny weights to 0
  if (abs(min.w) <= 1e-10) {
    w[abs(w) <= 1e-10] <- 0
  }

  opt.out$lambda <- lambda

  list(w = w, fit.obj = opt.out)
}

weightit2cfd.multi <- function(covs, treat, s.weights, subset, estimand, focal,
                               missing, verbose, ...) {

  missing <- .process_missing2(missing, covs)

  if (missing == "ind") {
    covs <- add_missing_indicators(covs)
  }

  d <- .process_kernel(covs, s.weights, subset,
                       kernel = ...get("kernel", "gaussian"),
                       bw_scale = ...get("bw_scale", 1),
                       nu = ...get("nu"),
                       nsim = ...get("nsim", 5000L))

  d <- unname(d[subset, subset, drop = FALSE])

  covs <- covs[subset, , drop = FALSE]
  treat <- factor(treat[subset])
  s.weights <- s.weights[subset]

  n <- length(treat)

  sw0 <- check_if_zero(s.weights)

  diagn <- diag(n)

  levels_treat <- levels(treat)

  min.w <- ...get("min.w", 1e-8)
  chk::chk_number(min.w)

  lambda <- ...get("lambda", 1e-4)
  chk::chk_number(lambda)

  moments <- ...get("moments", 0)
  int <- isTRUE(...get("int", FALSE))
  quantile <- ...get("quantile")
  add_constraints <- any(moments > 0) || int || is_not_null(quantile)

  if (add_constraints) {
    covs <- scale(covs)

    tols <- ...get("tols", 0)
    chk::chk_number(tols)
    tols <- abs(tols)
  }

  A <- ...mget(names(formals(osqp::osqpSettings)))

  eps <- ...get("eps", squish(min.w, lo = 1e-12, hi = 1e-8))
  if (is_not_null(eps)) {
    chk::chk_number(eps)
    if (is_null(A[["eps_abs"]])) A[["eps_abs"]] <- eps
    if (is_null(A[["eps_rel"]])) A[["eps_rel"]] <- eps
  }

  if (is_null(A[["max_iter"]])) A[["max_iter"]] <- 5e4L
  chk::chk_count(A[["max_iter"]], "`max_iter`")
  chk::chk_lt(A[["max_iter"]], Inf, "`max_iter`")
  if (is_null(A[["eps_abs"]])) A[["eps_abs"]] <- 1e-8
  chk::chk_number(A[["eps_abs"]], "`eps_abs`")
  if (is_null(A[["eps_rel"]])) A[["eps_rel"]] <- 1e-6
  chk::chk_number(A[["eps_rel"]], "`eps_rel`")
  if (is_null(A[["time_limit"]])) A[["time_limit"]] <- 0
  chk::chk_number(A[["time_limit"]], "`time_limit`")
  if (is_null(A[["adaptive_rho_interval"]])) A[["adaptive_rho_interval"]] <- 10L
  chk::chk_count(A[["adaptive_rho_interval"]], "`adaptive_rho_interval`")
  if (is_null(A[["polish"]])) A[["polish"]] <- TRUE
  chk::chk_flag(A[["polish"]], "`polish`")
  if (is_null(A[["polish_refine_iter"]])) A[["polish_refine_iter"]] <- 20L
  chk::chk_count(A[["polish_refine_iter"]], "`polish_refine_iter`")
  A[["verbose"]] <- verbose

  options.list <- do.call(osqp::osqpSettings, A)

  treat_t <- matrix(0, nrow = n, ncol = length(levels_treat),
                    dimnames = list(NULL, levels_treat))

  for (t in levels_treat) {
    in_t <- which(treat == t)
    s.weights[in_t] <- s.weights[in_t] / mean_fast(s.weights[in_t])
    treat_t[in_t, t] <- 1
  }

  n_t <- colSums(treat_t)

  s.weights_n_t <- vapply(levels_treat, function(t) treat_t[, t] * s.weights / n_t[t],
                          numeric(n))

  if (estimand == "ATE") {
    improved <- ...get("improved", TRUE)
    chk::chk_flag(improved)

    nn <- tcrossprod(s.weights_n_t)

    if (improved) {
      .col_diff <- function(x) x[, 1L] - x[, 2L]
      all_pairs <- utils::combn(levels_treat, 2L, simplify = FALSE)
      nn <- nn + tcrossprod(vapply(all_pairs, function(p) .col_diff(s.weights_n_t[, p, drop = FALSE]),
                                   numeric(n)))
    }

    P <- d * nn

    q <- -((s.weights * 2 / n) %*% d) * rowSums(s.weights_n_t)

    #Constraints for positivity and sum of weights
    Amat <- cbind(diagn, s.weights_n_t)
    lvec <- c(ifelse(sw0, 1, min.w), rep.int(1, length(levels_treat)))
    uvec <- c(ifelse(sw0, 1, Inf), rep.int(1, length(levels_treat)))

    unbounded <- lvec == -Inf & uvec == Inf

    if (any(unbounded)) {
      Amat <- Amat[, !unbounded, drop = FALSE]
      lvec <- lvec[!unbounded]
      uvec <- uvec[!unbounded]
    }

    if (add_constraints) {
      #Exactly balance moments, interactions, and/or quantiles
      covs <- .apply_moments_int_quantile(covs,
                                          moments = moments,
                                          int = int,
                                          quantile = quantile,
                                          s.weights = s.weights)

      targets <- col.w.m(covs, s.weights)

      sds <- rep.int(1, ncol(covs))

      if (tols > 0) {
        bin.vars <- is_binary_col(covs)

        if (!all(bin.vars)) {
          sds[!bin.vars] <- sqrt(colMeans(do.call("rbind",
                                                  lapply(levels_treat, function(t) {
                                                    col.w.v(covs[treat == t, !bin.vars, drop = FALSE], w = s.weights[treat == t])
                                                  }))))
        }
      }

      Amat <- cbind(Amat, do.call("cbind", apply(s.weights_n_t, 2L, function(x) covs * x, simplify = FALSE)))
      lvec <- c(lvec, rep.int(targets - sds * tols / 2, length(levels_treat)))
      uvec <- c(uvec, rep.int(targets + sds * tols / 2, length(levels_treat)))
    }
  }
  else {
    non_focal <- setdiff(levels_treat, focal)
    in_focal <- treat == focal

    nn <- tcrossprod(s.weights_n_t[!in_focal, non_focal, drop = FALSE])

    P <- d[!in_focal, !in_focal, drop = FALSE] * nn

    q <- -2 * (s.weights_n_t[in_focal, focal] %*% d[in_focal, !in_focal, drop = FALSE]) *
      rowSums(s.weights_n_t[!in_focal, non_focal, drop = FALSE])

    Amat <- cbind(diag(sum(!in_focal)), s.weights_n_t[!in_focal, non_focal])
    lvec <- c(ifelse(sw0[!in_focal], 1, min.w), rep.int(1, length(non_focal)))
    uvec <- c(ifelse(sw0[!in_focal], 1, Inf), rep.int(1, length(non_focal)))

    unbounded <- lvec == -Inf & uvec == Inf

    if (any(unbounded)) {
      Amat <- Amat[, !unbounded, drop = FALSE]
      lvec <- lvec[!unbounded]
      uvec <- uvec[!unbounded]
    }

    if (add_constraints) {
      #Exactly balance moments, interactions, and/or quantiles
      covs <- .apply_moments_int_quantile(covs,
                                          moments = moments,
                                          int = int,
                                          quantile = quantile,
                                          s.weights = s.weights, focal = focal,
                                          treat = treat)

      targets <- col.w.m(covs[in_focal, , drop = FALSE], s.weights[in_focal])

      sds <- rep.int(1, ncol(covs))

      if (tols > 0) {
        bin.vars <- is_binary_col(covs)

        if (!all(bin.vars)) {
          sds[!bin.vars] <- sqrt(col.w.v(covs[in_focal, !bin.vars, drop = FALSE], w = s.weights[in_focal]))
        }
      }

      Amat <- cbind(Amat, do.call("cbind", apply(s.weights_n_t[!in_focal, non_focal, drop = FALSE], 2L,
                                                 function(x) covs[!in_focal, , drop = FALSE] * x,
                                                 simplify = FALSE)))
      lvec <- c(lvec, rep.int(targets - sds * tols, length(non_focal)))
      uvec <- c(uvec, rep.int(targets + sds * tols, length(non_focal)))
    }
  }

  #Add weight penalty
  if (lambda < 0) {
    #Find lambda to make P PSD
    e <- eigen(P, symmetric = TRUE, only.values = TRUE)
    e.min <- min(e$values)
    if (e.min < 0) {
      lambda <- -e.min * n^2
    }
  }

  diag(P) <- diag(P) + lambda / n^2

  verbosely({
    opt.out <- osqp::solve_osqp(P = 2 * P, q = q, A = t(Amat), l = lvec, u = uvec,
                                pars = options.list)
  }, verbose = verbose)

  if (identical(opt.out$info$status, "maximum iterations reached")) {
    .wrn(sprintf("the optimization failed to converge. Try increasing `max_iter` (current value: %s)",
                 A[["max_iter"]]))
  }
  else if (identical(opt.out$info$status, "run time limit reached")) {
    .wrn(sprintf("the optimization failed to converge. Try increasing `time_limit` (current value: %s)",
                 A[["time_limit"]]))
  }
  else if (!startsWith(opt.out$info$status, "solved")) {
    .wrn("no feasible solution could be found that satisfies all constraints. Relax any constraints supplied")
  }

  if (estimand == "ATE") {
    w <- opt.out$x
  }
  else {
    w <- rep.int(1, n)
    w[treat != focal] <- opt.out$x
  }

  # Shrink tiny weights to 0
  if (abs(min.w) <= 1e-10) {
    w[abs(w) <= 1e-10] <- 0
  }

  opt.out$lambda <- lambda

  list(w = w, fit.obj = opt.out)
}

.process_kernel <- function(X, s.weights, subset, kernel = "gaussian",
                            bw_scale = 1, nu = NULL, nsim = 5000) {
  chk::chk_string(kernel)

  kernel <- match_arg(kernel, c("gaussian", "matern", "energy", "laplace", "t"))

  transform_method <- {
    if (kernel %in% c("gaussian", "matern")) "mahalanobis"
    else "scaled_euclidean"
  }

  X <- transform_covariates(data = X,
                            method = transform_method,
                            s.weights = s.weights,
                            discarded = !subset)

  D <- {
    if (kernel %in% c("laplace", "t")) as.matrix(dist(X, "manhattan"))
    else eucdist_internal(X)
  }

  if (kernel == "energy") {
    return(-D)
  }

  chk::chk_number(bw_scale)
  chk::chk_gt(bw_scale, 0)

  bw <- median(D[lower.tri(D)]) * bw_scale

  if (kernel == "gaussian") {
    return(exp(-D^2 / bw^2))
  }

  if (kernel == "matern") {
    if (is_null(nu)) {
      nu <- 3/2
    }

    chk::chk_number(nu)
    chk::chk_gt(nu, 0)
    chk::chk_lte(nu, 21/2)

    if (nu == 1/2) {
      return(exp(-D / bw))
    }

    if (nu == 3/2) {
      return((1 + D / bw) * exp(-D / bw))
    }

    if (nu == 5/2) {
      return((1 + D / bw + D^2 / (3 * bw^2)) * exp(-D / bw))
    }

    rlang::check_installed("GPBayes")
    return(GPBayes::matern(D, bw, nu))
  }

  if (kernel == "laplace") {
    return(exp(-D / bw))
  }

  # t-distribution
  # Monte carlo simulation for separable t-distribution
  if (is_null(nu)) {
    nu <- 5
  }

  chk::chk_number(nu)
  chk::chk_gt(nu, 0)

  chk::chk_count(nsim)
  chk::chk_gte(nsim, 10)

  V.random <- matrix(rt(ncol(X) * nsim, df = nu) / bw,
                     nrow = nsim)

  PX <- t(tcrossprod(X, V.random))

  D <- do.call("rbind", apply(PX, 2L, function(Pxii) {
    cos(PX - Pxii) |> colMeans()
  }, simplify = FALSE))

  D
}
#' Energy Balancing
#' @name method_energy
#' @aliases method_energy
#' @usage NULL
#'
#' @description This page explains the details of estimating weights using
#' energy balancing by setting `method = "energy"` in the call to [weightit()]
#' or [weightitMSM()]. This method can be used with binary, multi-category, and
#' continuous treatments.
#'
#' In general, this method relies on estimating weights by minimizing an energy
#' statistic related to covariate balance. For binary and multi-category
#' treatments, this is the energy distance, which is a multivariate distance
#' between distributions, between treatment groups. For continuous treatments,
#' this is the sum of the distance covariance between the treatment variable and
#' the covariates and the energy distances between the treatment and covariates
#' in the weighted sample and their distributions in the original sample. This
#' method relies on code written for \pkg{WeightIt} using \pkgfun{osqp}{osqp}
#' from the \CRANpkg{osqp} package to perform the optimization. This method may
#' be slow or memory-intensive for large datasets.
#'
#' ## Binary Treatments
#'
#' For binary treatments, this method estimates the weights using `osqp()` using
#' formulas described by Huling and Mak (2024). The following estimands are
#' allowed: ATE, ATT, and ATC.
#'
#' ## Multi-Category Treatments
#'
#' For multi-category treatments, this method estimates the weights using
#' `osqp()` using formulas described by Huling and Mak (2024). The following
#' estimands are allowed: ATE and ATT.
#'
#' ## Continuous Treatments
#'
#' For continuous treatments, this method estimates the weights using `osqp()`
#' using formulas described by Huling, Greifer, and Chen (2023).
#'
#' ## Longitudinal Treatments
#'
#' For longitudinal treatments, the weights are the product of the weights
#' estimated at each time point. This method is not guaranteed to yield optimal
#' balance at each time point. NOTE: the use of energy balancing with
#' longitudinal treatments has not been validated!
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
#' @section Additional Arguments: The following following additional arguments
#'   can be specified:
#' \describe{
#'   \item{`dist.mat`}{the name of the method used to compute the distance matrix of the covariates or the numeric distance matrix itself. Allowable options include `"scaled_euclidean"` for the Euclidean (L2) distance on the scaled covariates (the default), `"mahalanobis"` for the Mahalanobis distance, and `"euclidean"` for the raw Euclidean distance. Abbreviations allowed. Note that some user-supplied distance matrices can cause the R session to abort due to a bug within \pkg{osqp}, so this argument should be used with caution. A distance matrix must be a square, symmetric, numeric matrix with zeros along the diagonal and a row and column for each unit. Can also be supplied as the output of a call to [dist()].
#'   }
#'   \item{`lambda`}{a positive numeric scalar used to penalize the square of the weights. This value divided by the square of the total sample size is added to the diagonal of the quadratic part of the loss function. Higher values favor weights with less variability. Note this is distinct from the lambda value described in Huling and Mak (2024), which penalizes the complexity of individual treatment rules rather than the weights, but does correspond to lambda from Huling et al. (2023). Default is .0001, which is essentially 0.
#'   }
#' }
#'   For binary and multi-category treatments, the following additional
#'   arguments can be specified:
#'   \describe{
#'     \item{`improved`}{`logical`; whether to use the improved energy balancing weights as described by Huling and Mak (2024) when `estimand = "ATE"`. This involves optimizing balance not only between each treatment group and the overall sample, but also between each pair of treatment groups. Huling and Mak (2024) found that the improved energy balancing weights generally outperformed standard energy balancing. Default is `TRUE`; set to `FALSE` to use the standard energy balancing weights instead (not recommended).
#'     }
#'   \item{`quantile`}{
#'     A named list of quantiles (values between 0 and 1) for each continuous covariate, which are used to create additional variables that when balanced ensure balance on the corresponding quantile of the variable. For example, setting `quantile = list(x1 = c(.25, .5. , .75))` ensures the 25th, 50th, and 75th percentiles of `x1` in each treatment group will be balanced in the weighted sample. Can also be a single number (e.g., `.5`) or an unnamed list of length 1 (e.g., `list(c(.25, .5, .75))`) to request the same quantile(s) for all continuous covariates, or a named vector (e.g., `c(x1 = .5, x2 = .75)` to request one quantile for each covariate.
#'   }
#'   }
#'
#'   For continuous treatments, the following additional arguments can be
#'   specified:
#'   \describe{
#'     \item{`d.moments`}{
#'       The number of moments of the treatment and covariate distributions that are constrained to be the same in the weighted sample as in the original sample. For example, setting `d.moments = 3` ensures that the mean, variance, and skew of the treatment and covariates are the same in the weighted sample as in the unweighted sample. `d.moments` should be greater than or equal to `moments` and will be automatically set accordingly if not (or if not specified).
#'     }
#'     \item{`dimension.adj`}{
#'       `logical`; whether to include the dimensionality adjustment described by Huling et al. (2023). If `TRUE`, the default, the energy distance for the covariates is weighted \eqn{\sqrt{p}} times as much as the energy distance for the treatment, where \eqn{p} is the number of covariates. If `FALSE`, the two energy distances are given equal weights. Default is `TRUE`.
#'     }
#'   }
#'
#'   The `moments` argument functions differently for `method = "energy"` from
#'   how it does with other methods. When unspecified or set to zero, energy
#'   balancing weights are estimated as described by Huling and Mak (2024) for
#'   binary and multi-category treatments or by Huling et al. (2023) for
#'   continuous treatments. When `moments` is set to an integer larger than 0,
#'   additional balance constraints on the requested moments of the covariates
#'   are also included, guaranteeing exact moment balance on these covariates
#'   while minimizing the energy distance of the weighted sample. For binary and
#'   multi-category treatments, this involves exact balance on the means of the
#'   entered covariates; for continuous treatments, this involves exact balance
#'   on the treatment-covariate correlations of the entered covariates.
#'
#'   Any other arguments will be passed to \pkgfun{osqp}{osqpSettings}. Some
#'   defaults differ from those in `osqpSettings()`; see *Reproducibility*
#'   below.
#'
#' @section Additional Outputs:
#' \describe{
#'   \item{`obj`}{When `include.obj = TRUE`, the output of the call to \pkgfun{osqp}{solve_osqp}, which contains the dual variables and convergence information.
#'   }
#' }
#'
#' @details Energy balancing is a method of estimating weights using
#' optimization without a propensity score. The weights are the solution to a
#' constrain quadratic optimization problem where the objective function
#' concerns covariate balance as measured by the energy distance and (for
#' continuous treatments) the distance covariance.
#'
#' Energy balancing for binary and multi-category treatments involves minimizing
#' the energy distance between the treatment groups and between each treatment
#' group and a target group (e.g., the full sample for the ATE). The energy
#' distance is a scalar measure of the difference between two multivariate
#' distributions and is equal to 0 when the two distributions are identical.
#'
#' Energy balancing for continuous treatments involves minimizing the distance
#' covariance between the treatment and the covariates; the distance covariance
#' is a scalar measure of the association between two (possibly multivariate)
#' distributions that is equal to 0 when the two distributions are independent.
#' In addition, the energy distances between the treatment and covariate
#' distributions in the weighted sample and the treatment and covariate
#' distributions in the original sample are minimized.
#'
#' The primary benefit of energy balancing is that all features of the covariate
#' distribution are balanced, not just means, as with other optimization-based
#' methods like entropy balancing. Still, it is possible to add additional
#' balance constraints to require balance on individual terms using the
#' `moments` argument, just like with entropy balancing. Energy balancing can
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
#' @note Sometimes the optimization can fail to converge because the problem is
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
#' @author Noah Greifer, using code from Jared Huling's
#' \CRANpkg{independenceWeights} package for continuous treatments.
#'
#' @seealso [weightit()], [weightitMSM()]
#'
#' @references ## Binary and multi-category treatments
#'
#' Huling, J. D., & Mak, S. (2024). Energy balancing of covariate distributions.
#' *Journal of Causal Inference*, 12(1). \doi{10.1515/jci-2022-0029}
#'
#' ## Continuous treatments
#'
#' Huling, J. D., Greifer, N., & Chen, G. (2023). Independence weights for
#' causal inference with continuous treatments. *Journal of the American
#' Statistical Association*, 0(ja), 1–25. \doi{10.1080/01621459.2023.2213485}
#'
#' @examplesIf requireNamespace("osqp", quietly = TRUE)
#' library("cobalt")
#' data("lalonde", package = "cobalt")
#'
#' #Balancing covariates between treatment groups (binary)
#' (W1 <- weightit(treat ~ age + educ + married +
#'                   nodegree + re74, data = lalonde,
#'                 method = "energy", estimand = "ATE"))
#' summary(W1)
#' bal.tab(W1)
#'
#' #Balancing covariates with respect to race (multi-category)
#' (W2 <- weightit(race ~ age + educ + married +
#'                   nodegree + re74, data = lalonde,
#'                 method = "energy", estimand = "ATT",
#'                 focal = "black"))
#' summary(W2)
#' bal.tab(W2)
#' \donttest{
#'   #Balancing covariates with respect to re75 (continuous)
#'   (W3 <- weightit(re75 ~ age + educ + married +
#'                     nodegree + re74, data = lalonde,
#'                   method = "energy", moments = 1))
#'   summary(W3)
#'   bal.tab(W3, poly = 2)
#' }
NULL

weightit2energy <- function(covs, treat, s.weights, subset, estimand, focal,
                            missing, moments, int, verbose, ...) {

  missing <- .process_missing2(missing, covs)

  if (missing == "ind") {
    covs <- add_missing_indicators(covs)
  }

  d <- ...get("dist.mat", "scaled_euclidean")

  if (chk::vld_string(d)) {
    dist.covs <- transform_covariates(data = covs, method = d,
                                      s.weights = s.weights, discarded = !subset)
    d <- unname(eucdist_internal(dist.covs))
  }
  else {
    if (inherits(d, "dist")) d <- as.matrix(d)

    if (!is.matrix(d) || !all(dim(d) == length(treat)) ||
        !all(check_if_zero(diag(d))) ||
        any(d < 0) ||
        !isSymmetric(unname(d))) {
      .err(sprintf("`dist.mat` must be one of %s or a square, symmetric distance matrix with a value for all pairs of units",
                   word_list(weightit_distances(), "or", quotes = TRUE)))
    }
  }

  d <- unname(d[subset, subset])

  covs <- covs[subset, , drop = FALSE]
  treat <- treat[subset]
  s.weights <- s.weights[subset]

  t.lev <- get_treated_level(treat, estimand, focal)
  treat <- binarize(treat, one = t.lev)

  n <- length(treat)
  diagn <- diag(n)

  covs <- scale(covs)

  min.w <- ...get("min.w", 1e-8)
  chk::chk_number(min.w)

  lambda <- ...get("lambda", 1e-4)
  chk::chk_number(lambda)

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

    P <- -d * nn

    q <- ((s.weights * 2 / n) %*% d) * (s.weights_n_0 + s.weights_n_1)

    #Constraints for positivity and sum of weights
    Amat <- cbind(diagn, s.weights_n_0, s.weights_n_1)
    lvec <- c(rep.int(min.w, n), 1, 1)
    uvec <- c(ifelse(check_if_zero(s.weights), min.w, Inf), 1, 1)

    if (moments > 0 || int || is_not_null(...get("quantile"))) {
      #Exactly balance moments, interactions, and/or quantiles
      covs <- cbind(.int_poly_f(covs, poly = moments, int = int, center = TRUE),
                    .quantile_f(covs, qu = ...get("quantile"), s.weights = s.weights))

      targets <- col.w.m(covs, s.weights)

      Amat <- cbind(Amat, covs * s.weights_n_0, covs * s.weights_n_1)
      lvec <- c(lvec, targets, targets)
      uvec <- c(uvec, targets, targets)
    }
  }
  else if (estimand == "ATT") {
    nn <- tcrossprod(s.weights_n_0[t0])

    P <- -d[t0, t0] * nn

    q <- 2 * (s.weights_n_1[t1] %*% d[t1, t0]) * s.weights_n_0[t0]

    Amat <- cbind(diag(n0), s.weights_n_0[t0])
    lvec <- c(rep.int(min.w, n0), 1)
    uvec <- c(ifelse(check_if_zero(s.weights[t0]), min.w, Inf), 1)

    if (moments > 0 || int || is_not_null(...get("quantile"))) {
      #Exactly balance moments, interactions, and/or quantiles
      covs <- cbind(.int_poly_f(covs, poly = moments, int = int, center = TRUE),
                    .quantile_f(covs, qu = ...get("quantile"), s.weights = s.weights,
                                focal = 1, treat = treat))

      targets <- col.w.m(covs[t1, , drop = FALSE], s.weights[t1])

      Amat <- cbind(Amat, covs[t0, , drop = FALSE] * s.weights_n_0[t0])

      lvec <- c(lvec, targets)
      uvec <- c(uvec, targets)
    }
  }
  else if (estimand == "ATC") {
    nn <- tcrossprod(s.weights_n_1[t1])

    P <- -d[t1, t1] * nn

    q <- 2 * (s.weights_n_0[t0] %*% d[t0, t1]) * s.weights_n_1[t1]

    Amat <- cbind(diag(n1), s.weights_n_1[t1])
    lvec <- c(rep.int(min.w, n1), 1)
    uvec <- c(ifelse(check_if_zero(s.weights[t1]), min.w, Inf), 1)

    if (moments > 0 || int || is_not_null(...get("quantile"))) {
      #Exactly balance moments, interactions, and/or quantiles
      covs <- cbind(.int_poly_f(covs, poly = moments, int = int, center = TRUE),
                    .quantile_f(covs, qu = ...get("quantile"), s.weights = s.weights,
                                focal = 0, treat = treat))

      targets <- col.w.m(covs[t0, , drop = FALSE], s.weights[t0])

      Amat <- cbind(Amat, covs[t1, , drop = FALSE] * s.weights_n_1[t1])

      lvec <- c(lvec, targets)
      uvec <- c(uvec, targets)
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

  A <- ...mget(names(formals(osqp::osqpSettings)))

  if (is_not_null(...get("eps"))) {
    chk::chk_number(...get("eps"), "`eps`")
    if (is_null(A[["eps_abs"]])) A[["eps_abs"]] <- ...get("eps")
    if (is_null(A[["eps_rel"]])) A[["eps_rel"]] <- ...get("eps")
  }

  if (is_null(A[["max_iter"]])) A[["max_iter"]] <- 4e3L
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
  A[["verbose"]] <- TRUE

  options.list <- do.call(osqp::osqpSettings, A)

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

  w[w <= min.w] <- min.w

  opt.out$lambda <- lambda

  list(w = w, fit.obj = opt.out)
}

weightit2energy.multi <- function(covs, treat, s.weights, subset, estimand, focal,
                                  missing, moments, int, verbose, ...) {

  missing <- .process_missing2(missing, covs)

  if (missing == "ind") {
    covs <- add_missing_indicators(covs)
  }

  d <- ...get("dist.mat", "scaled_euclidean")

  if (chk::vld_string(d)) {
    dist.covs <- transform_covariates(data = covs, method = d,
                                      s.weights = s.weights, discarded = !subset)
    d <- unname(eucdist_internal(dist.covs))
  }
  else {
    if (inherits(d, "dist")) d <- as.matrix(d)

    if (!is.matrix(d) || !all(dim(d) == length(treat)) ||
        !all(check_if_zero(diag(d))) || any(d < 0) ||
        !isSymmetric(unname(d))) {
      .err(sprintf("`dist.mat` must be one of %s or a square, symmetric distance matrix with a value for all pairs of units",
                   word_list(weightit_distances(), "or", quotes = TRUE)))
    }

  }

  d <- unname(d[subset, subset])

  covs <- covs[subset, , drop = FALSE]
  treat <- factor(treat[subset])
  s.weights <- s.weights[subset]

  n <- length(treat)
  levels_treat <- levels(treat)
  diagn <- diag(n)

  covs <- scale(covs)

  min.w <- ...get("min.w", 1e-8)
  chk::chk_number(min.w)

  lambda <- ...get("lambda", 1e-4)
  chk::chk_number(lambda)

  for (t in levels_treat) {
    s.weights[treat == t] <- s.weights[treat == t] / mean_fast(s.weights[treat == t])
  }

  treat_t <- vapply(levels_treat, function(t) treat == t, logical(n))
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

    P <- -d * nn

    q <- ((s.weights * 2 / n) %*% d) * rowSums(s.weights_n_t)

    #Constraints for positivity and sum of weights
    Amat <- cbind(diagn, s.weights_n_t)
    lvec <- c(rep.int(min.w, n), rep.int(1, length(levels_treat)))
    uvec <- c(ifelse(check_if_zero(s.weights), min.w, Inf), rep.int(1, length(levels_treat)))

    if (moments > 0 || int || is_not_null(...get("quantile"))) {
      #Exactly balance moments, interactions, and/or quantiles
      covs <- cbind(.int_poly_f(covs, poly = moments, int = int, center = TRUE),
                    .quantile_f(covs, qu = ...get("quantile"), s.weights = s.weights))

      targets <- col.w.m(covs, s.weights)

      Amat <- cbind(Amat, do.call("cbind", apply(s.weights_n_t, 2L, function(x) covs * x, simplify = FALSE)))
      lvec <- c(lvec, rep.int(targets, length(levels_treat)))
      uvec <- c(uvec, rep.int(targets, length(levels_treat)))
    }
  }
  else {
    non_focal <- setdiff(levels_treat, focal)
    in_focal <- treat == focal

    nn <- tcrossprod(s.weights_n_t[!in_focal, non_focal, drop = FALSE])

    P <- -d[!in_focal, !in_focal] * nn

    q <- 2 * (s.weights_n_t[in_focal, focal] %*% d[in_focal, !in_focal]) *
      rowSums(s.weights_n_t[!in_focal, non_focal, drop = FALSE])

    Amat <- cbind(diag(sum(!in_focal)), s.weights_n_t[!in_focal, non_focal])
    lvec <- c(rep.int(min.w, sum(!in_focal)), rep.int(1, length(non_focal)))
    uvec <- c(ifelse(check_if_zero(s.weights[!in_focal]), min.w, Inf), rep.int(1, length(non_focal)))

    if (moments > 0 || int || is_not_null(...get("quantile"))) {
      #Exactly balance moments, interactions, and/or quantiles
      covs <- cbind(.int_poly_f(covs, poly = moments, int = int, center = TRUE),
                    .quantile_f(covs, qu = ...get("quantile"), s.weights = s.weights,
                                focal = focal, treat = treat))

      targets <- col.w.m(covs[in_focal, , drop = FALSE], s.weights[in_focal])

      Amat <- cbind(Amat, do.call("cbind", apply(s.weights_n_t[!in_focal, non_focal, drop = FALSE], 2L,
                                                 function(x) covs[!in_focal, , drop = FALSE] * x,
                                                 simplify = FALSE)))
      lvec <- c(lvec, rep.int(targets, length(non_focal)))
      uvec <- c(uvec, rep.int(targets, length(non_focal)))
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

  A <- ...mget(names(formals(osqp::osqpSettings)))

  if (is_not_null(...get("eps"))) {
    chk::chk_number(...get("eps"), "`eps`")
    if (is_null(A[["eps_abs"]])) A[["eps_abs"]] <- ...get("eps")
    if (is_null(A[["eps_rel"]])) A[["eps_rel"]] <- ...get("eps")
  }

  if (is_null(A[["max_iter"]])) A[["max_iter"]] <- 4e3L
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
  A[["verbose"]] <- TRUE

  options.list <- do.call(osqp::osqpSettings, A)

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

  w[w <= min.w] <- min.w

  opt.out$lambda <- lambda

  list(w = w, fit.obj = opt.out)
}

weightit2energy.cont <- function(covs, treat, s.weights, subset, missing, moments, int, verbose, ...) {

  missing <- .process_missing2(missing, covs)

  if (missing == "ind") {
    covs <- add_missing_indicators(covs)
  }

  Xdist <- ...get("dist.mat", "scaled_euclidean")

  if (chk::vld_string(Xdist)) {
    dist.covs <- transform_covariates(data = covs, method = Xdist,
                                      s.weights = s.weights, discarded = !subset)
    Xdist <- unname(eucdist_internal(X = dist.covs))
  }
  else {
    if (inherits(Xdist, "dist")) Xdist <- as.matrix(Xdist)

    if (!is.matrix(Xdist) || !all(dim(Xdist) == length(treat)) ||
        !all(check_if_zero(diag(Xdist))) || any(Xdist < 0) ||
        !isSymmetric(unname(Xdist))) {
      .err(sprintf("`dist.mat` must be one of %s or a square, symmetric distance matrix with a value for all pairs of units",
                   word_list(weightit_distances(), "or", quotes = TRUE)))
    }
  }

  Xdist <- unname(Xdist[subset, subset])

  if (is_null(...get("treat.dist.mat"))) {
    Adist <- eucdist_internal(X = treat / sqrt(col.w.v(treat, s.weights)))
  }
  else {
    Adist <- ...get("treat.dist.mat")

    if (inherits(Adist, "dist")) Adist <- as.matrix(Adist)

    if (!is.matrix(Adist) || !all(dim(Adist) == length(treat)) ||
        !all(check_if_zero(diag(Adist))) || any(Adist < 0) ||
        !isSymmetric(unname(Adist))) {
      .err("`treat.dist.mat` must be a square, symmetric distance matrix with a value for all pairs of units")
    }
  }

  Adist <- unname(Adist[subset, subset])

  covs <- covs[subset, , drop = FALSE]
  treat <- treat[subset]
  s.weights <- s.weights[subset]

  n <- length(treat)

  sw0 <- check_if_zero(s.weights)

  s.weights <- n * s.weights / sum(s.weights)

  min.w <- ...get("min.w", 1e-8)
  chk::chk_number(min.w)

  lambda <- ...get("lambda", 1e-4)
  chk::chk_number(lambda)

  d.moments <- max(...get("d.moments", 0), moments)
  chk::chk_count(d.moments)

  dimension.adj <- ...get("dimension.adj", TRUE)
  chk::chk_flag(dimension.adj)

  Xmeans <- colMeans(Xdist)
  Xgrand_mean <- mean(Xmeans)
  XA <- Xdist + Xgrand_mean - outer(Xmeans, Xmeans, "+")

  Ameans <- colMeans(Adist)
  Agrand_mean <- mean(Ameans)
  AA <- Adist + Agrand_mean - outer(Ameans, Ameans, "+")

  Pdcow <- XA * AA / n^2
  PebA <- -Adist / n^2
  PebX <- -Xdist / n^2

  qebA <- drop(s.weights %*% Adist) * 2 / n^2
  qebX <- drop(s.weights %*% Xdist) * 2 / n^2

  Q_energy_A_adj <- {
    if (dimension.adj) 1 / (1 + sqrt(ncol(covs)))
    else .5
  }

  Q_energy_X_adj <- 1 - Q_energy_A_adj

  PebA <- PebA * Q_energy_A_adj
  PebX <- PebX * Q_energy_X_adj

  qebA <- qebA * Q_energy_A_adj
  qebX <- qebX * Q_energy_X_adj

  P <- Pdcow + PebA + PebX
  q <- qebA + qebX

  P <- P * tcrossprod(s.weights)
  q <- q * s.weights

  Amat <- cbind(diag(n), s.weights)
  lvec <- c(rep.int(min.w, n), n)
  uvec <- c(ifelse(sw0, min.w, Inf), n)

  if (d.moments > 0) {
    d.covs <- .int_poly_f(covs, poly = d.moments)
    d.treat <- cbind(poly(treat, degree = d.moments))

    d.covs <- center(d.covs, col.w.m(d.covs, s.weights))
    d.treat <- center(d.treat, col.w.m(d.treat, s.weights))

    Amat <- cbind(Amat, d.covs * s.weights, d.treat * s.weights)
    lvec <- c(lvec, rep.int(0, ncol(d.covs)), rep.int(0, ncol(d.treat)))
    uvec <- c(uvec, rep.int(0, ncol(d.covs)), rep.int(0, ncol(d.treat)))
  }

  if (moments > 0 || int) {
    covs <- .int_poly_f(covs, poly = moments, int = int)

    X.means <- col.w.m(covs, s.weights)
    A.mean <- w.m(treat, s.weights)

    covs <- center(covs, X.means)
    treat <- treat - A.mean

    Amat <- cbind(Amat, covs * treat * s.weights)

    lvec <- c(lvec, rep.int(0, ncol(covs)))
    uvec <- c(uvec, rep.int(0, ncol(covs)))
  }

  #Add weight penalty
  if (lambda < 0) {
    #Find lambda to make P PSD
    e <- eigen(P, symmetric = TRUE, only.values = TRUE)
    e.min <- min(e$values)

    lambda <- -e.min * n^2
  }

  diag(P) <- diag(P) + lambda / n^2

  A <- ...mget(names(formals(osqp::osqpSettings)))

  if (is_not_null(...get("eps"))) {
    chk::chk_number(...get("eps"), "`eps`")
    if (is_null(A[["eps_abs"]])) A[["eps_abs"]] <- ...get("eps")
    if (is_null(A[["eps_rel"]])) A[["eps_rel"]] <- ...get("eps")
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
  A[["verbose"]] <- TRUE

  options.list <- do.call(osqp::osqpSettings, A)

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

  w <- opt.out$x
  w[w <= min.w] <- min.w

  opt.out$lambda <- lambda

  list(w = w, fit.obj = opt.out)
}

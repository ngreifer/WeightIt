#Kernel balancing using Gaussian kernel; currently only works for ATT, but works
data("lalonde")

d <- MatchIt::scaled_euclidean_dist(~ age + educ+ race+married +nodegree +re74 + re75, data = lalonde)

sigma <- median(d[lower.tri(d)])

d <- exp(-d/sigma)
# diag(d) <- 0

We <- weightit(treat~ age + educ+ race+married +nodegree +re74 + re75, data = lalonde,
              dist.mat = "scaled_euclidean", method = "energy", estimand = "ATT")

weightit2kb <- function(covs, treat, s.weights, subset, estimand, focal,
                        missing, moments, int, verbose, ...) {

  rlang::check_installed("osqp")

  A <- list(...)

  if (missing == "ind") {
    covs <- add_missing_indicators(covs)
  }

  d <- if_null_then(A[["dist.mat"]], "scaled_euclidean")
  A[["dist.mat"]] <- NULL

  if (is.character(d) && length(d) == 1L) {
    dist.covs <- transform_covariates(data = covs, method = d,
                                      s.weights = s.weights, discarded = !subset)
    d <- unname(eucdist_internal(dist.covs))
  }
  else {
    if (inherits(d, "dist")) d <- as.matrix(d)

    if (!is.matrix(d) || !all(dim(d) == length(treat)) ||
        # !all(check_if_zero(diag(d))) ||
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

  t.lev <- get_treated_level(treat)
  treat <- binarize(treat, one = t.lev)

  n <- length(treat)
  diagn <- diag(n)

  covs <- scale(covs)

  min.w <- if_null_then(A[["min.w"]], 1e-8)
  chk::chk_number(min.w)

  lambda <- if_null_then(A[["lambda"]], 1e-4)
  chk::chk_number(lambda)

  t0 <- which(treat == 0)
  t1 <- which(treat == 1)


  n0 <- length(t0)
  n1 <- length(t1)

  n <- n0 + n1

  sigma2 <- median(d[lower.tri(d)]^2)^2

  K <- exp(-d/(sigma2))

  K[t0, t1] <- -K[t0, t1]
  K[t1, t0] <- -K[t1, t0]

  #Constraints for positivity and sum of weights
  if (estimand == "ATE") {
    # Need to fix q to make ATE valid; currently balances toward overlap
    # q <-
    #   #t1 vs all
    #   diag(treat) %*% K %*% rep(1, n) -
    #   #t0 vs all
    #   diag(1 - treat) %*% K %*% rep(1, n)
    q <- NULL

    P <- K

    Amat <- cbind(diagn, treat, 1 - treat)
    lvec <- c(rep(min.w, n), 1, 1)
    uvec <- c(rep(Inf, n), 1, 1)
  }
  else if (estimand == "ATT") {
    P <- K
    q <- NULL

    Amat <- cbind(diagn, treat, 1 - treat)
    lvec <- c(ifelse(treat == 1, 1 / n1, min.w), 1, 1)
    uvec <- c(ifelse(treat == 1, 1 / n1, Inf), 1, 1)
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

  if (is_not_null(A[["eps"]])) {
    chk::chk_number(A[["eps"]], "`eps`")
    if (is_null(A[["eps_abs"]])) A[["eps_abs"]] <- A[["eps"]]
    if (is_null(A[["eps_rel"]])) A[["eps_rel"]] <- A[["eps"]]
  }
  A[names(A) %nin% names(formals(osqp::osqpSettings))] <- NULL
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


  w <- opt.out$x

  for (i in 0:1) {
    w[treat == i] <- w[treat == i] / mean(w[treat == i])
  }

  w[w <= min.w] <- min.w

  opt.out$lambda <- lambda

  list(w = w, fit.obj = opt.out)
}

W <- weightit(treat~ age + educ+ race+married +nodegree +re74 + re75, data = lalonde,
              method = weightit2kb, estimand = "ATE")

bal.tab(W, weights = We)

#Optimization-based weights using CVXR for arbitrary objective functions, with
#support for targets, tolerance, and base weights. It works! (Able to reproduce
#results from entropy balancing and SBW)
#Try to figure out balancing subject to contraints on ESS?

ow.cvxr <- function(formula, data, estimand = "ATE", targets = NULL, focal = NULL,
                    tols = 0, obj = "l2",
                    s.weights = NULL, b.weights = NULL, min.w = 0,
                    verbose = FALSE) {

  #Process treat and covs from formula and data
  t.c <- get_covs_and_treat_from_formula(formula, data)
  reported.covs <- t.c[["reported.covs"]]
  covs <- t.c[["model.covs"]]
  treat <- t.c[["treat"]]

  treat <- assign_treat_type(treat)
  treat.type <- get_treat_type(treat)

  n <- length(treat)

  #Process s.weights
  s.weights <- .process.s.weights(s.weights, data)

  if (is_null(s.weights)) s.weights <- rep.int(1, n)

  #Process estimand and focal and targets
  if (is.null(targets)) {
    estimand <- .process_estimand(estimand, "optweight", treat.type)
    f.e.r <- .process_focal_and_estimand(focal, estimand, treat)
    focal <- f.e.r[["focal"]]
    # estimand <- f.e.r[["estimand"]]
    estimand <- f.e.r[["reported.estimand"]]
    reported.estimand <- f.e.r[["reported.estimand"]]

    targets <- {
      if (is_null(focal)) col.w.m(covs, s.weights)
      else col.w.m(covs[treat == focal,,drop = FALSE], s.weights[treat == focal])
    }
  }

  #Process tols
  if (length(tols) == 1) {
    tols <- rep(tols, ncol(covs))
  }

  #Standardize
  if (is_null(focal)) {
    sds <- sqrt(colMeans(do.call("rbind", lapply(unique(treat), function(i) {
      col.w.v(covs[treat == i,,drop = FALSE], s.weights[treat == i])
    }))))

    tols <- tols / 2
  }
  else {
    sds <- sqrt(col.w.v(covs[treat == focal,,drop = FALSE],
                        s.weights[treat == focal]))
  }

  covs <- sweep(covs, 2, sds, "/")
  targets <- targets / sds

  out <- ow.cvxr.fit(treat, covs, tols, targets, obj, focal, s.weights, b.weights,
                     min.w, verbose)

  out
}

ow.cvxr.fit <- function(A, X, tols = 0, targets, obj = "l2", focal = NULL,
                        s.weights = NULL, b.weights = NULL, min.w = 0,
                        verbose = FALSE, ...) {

  s.weights <- rep(1, length(A))
  for (t in unique(A)) {
    s.weights[A == t] <- s.weights[A == t]/mean(s.weights[A == t])
  }

  # Groups to weight
  g2w <- if (is_null(focal)) unique(A) else setdiff(unique(A), focal)

  A_ <- A[A %in% g2w]

  N <- length(A_)

  if (is_null(b.weights)) {
    b.weights <- rep(1, N)
    for (i in g2w) {
      b.weights[A == i] <- 1/sum(A == i)
    }
  }

  B <- b.weights[A %in% g2w]
  s.weights <- s.weights[A %in% g2w]

  W <- CVXR::Variable(length(A_))

  #Sum constraints
  const_sum <- lapply(g2w, function(i) {
    sum(s.weights[A_ == i] * W[A_ == i]) == 1
  })

  #Positivity Constraints
  const_pos <- list(W >= min.w)

  #Target constraints
  X <- sweep(X, 2, targets)
  const_target <- lapply(g2w, function(i) {
    abs(t(X[A == i,]) %*% (s.weights[A_ == i] * W[A_ == i])) <= tols
  })

  #Objective
  objective <- switch(obj,
                      "l2" = CVXR::Minimize(CVXR::p_norm(W - B, 2)),
                      "l1" = CVXR::Minimize(CVXR::p_norm(W - B, 1)),
                      "linf" = CVXR::Minimize(CVXR::p_norm(W - B, Inf)),
                      "entropy" = CVXR::Minimize(sum(CVXR::kl_div(W, B))),
                      "log" = CVXR::Minimize(-sum(log(W))), #npCBPS
                      "logit" = CVXR::Minimize(sum(((B - exp(-B) - W)))),
                      "DS1" = CVXR::Minimize(sum((W - B)^2/B)/2),
                      "DS2" = CVXR::Minimize(sum(-CVXR::entr(W) - W * log(B) - W + B)),
                      "DS3" = CVXR::Minimize(sum(W - 2 * sqrt(W * B) + B)),
                      "DS4" = CVXR::Minimize(sum(-B * log(W) -CVXR::entr(B) + W - B)),
                      "DS5" = CVXR::Minimize(sum(W - 2*B + (B^2) * CVXR::inv_pos(W))))

  prob <- CVXR::Problem(objective, c(const_sum,
                                     const_pos,
                                     const_target))

  res <- CVXR::psolve(prob, ignore_dcp = TRUE)

  w_ <- drop(res$getValue(W))

  duals <- lapply(const_target, res$getDualValue)

  if (is_null(focal)) {
    w <- w_
  }
  else {
    w <- rep(1, length(A))
    w[A %in% g2w] <- w_
  }

  for (i in g2w) {
    w[A == i] <- w[A == i] / mean(w[A == i])
  }

  opt_out <- list(w = w,
                  duals = duals,
                  out = res)

  opt_out
}

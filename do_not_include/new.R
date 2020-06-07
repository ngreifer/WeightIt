#New methods and functions

#------Preliminary template----
weightit2XXX <- function(covs, treat...) {
  stop("method = \"XXX\" isn't ready to use yet.", call. = FALSE)
}
#------Template----
weightit2XXX <- function(covs, treat, s.weights, subset, estimand, focal, missing, moments, int, ...) {
  A <- list(...)

  covs <- covs[subset, , drop = FALSE]
  treat <- factor(treat[subset])

  if (missing == "ind") {
    missing.ind <- apply(covs[, apply(covs, 2, anyNA), drop = FALSE], 2, function(x) as.numeric(is.na(x)))
    covs[is.na(covs)] <- 0
    covs <- cbind(covs, missing.ind)
  }

  covs <- cbind(covs, int.poly.f(covs, poly = moments, int = int, center = TRUE))
  covs <- apply(covs, 2, make.closer.to.1)

  new.data <- data.frame(treat, covs)
  new.formula <- formula(new.data)

  for (f in names(formals(PACKAGE::FUNCTION))) {
    if (is_null(A[[f]])) A[[f]] <- formals(PACKAGE::FUNCTION)[[f]]
  }
  A[names(A) %in% names(formals(weightit2XXX))] <- NULL

  A[["formula"]] <- new.formula
  A[["data"]] <- new.data
  A[["estimand"]] <- estimand
  A[["s.weights"]] <- s.weights[subset]
  A[["focal"]] <- focal
  A[["verbose"]] <- TRUE

  if (check.package("optweight")) {
    out <- do.call(PACKAGE::FUNCTION, A, quote = TRUE)
    obj <- list(w = out[["weights"]], fit.obj = out)
    return(obj)
  }
}

#------Under construction----


#------Ready for use, but not ready for CRAN----
#KBAL
weightit2kbal <- function(covs, treat, s.weights, subset, estimand, focal, ...) {
  A <- list(...)

  covs <- covs[subset, , drop = FALSE]
  treat <- factor(treat)[subset]

  covs <- apply(covs, 2, make.closer.to.1)

  if (any(vars.w.missing <- apply(covs, 2, function(x) anyNA(x)))) {
    missing.ind <- apply(covs[, vars.w.missing, drop = FALSE], 2, function(x) as.numeric(is.na(x)))
    covs[is.na(covs)] <- 0
    covs <- cbind(covs, missing.ind)
  }

  if ("kbal.method" %in% names(A)) {
    names(A)[names(A) == "kbal.method"] <- "method"
  }
  for (f in names(formals(KBAL::kbal))) {
    if (is_null(A[[f]])) A[[f]] <- formals(KBAL::kbal)[[f]]
  }
  A[names(A) %nin% setdiff(names(formals(KBAL::kbal)), c("X", "D"))] <- NULL

  if (check.package("KBAL")) {
    if (hasName(A, "method")) {
      if (A[["method"]] == "el") check.package(c("glmc", "emplik"))
    }

    if (estimand == "ATT") {
      w <- rep(1, length(treat))
      control.levels <- levels(treat)[levels(treat) != focal]
      fit.list <- setNames(vector("list", length(control.levels)), control.levels)

      covs[treat == focal,] <- covs[treat == focal, , drop = FALSE] * s.weights[subset][treat == focal] * sum(treat == focal)/sum(s.weights[subset][treat == focal])

      for (i in control.levels) {
        treat.in.i.focal <- treat %in% c(focal, i)
        treat_ <- ifelse(treat[treat.in.i.focal] == i, 0L, 1L)
        covs_ <- covs[treat.in.i.focal, , drop = FALSE]

        colinear.covs.to.remove <- colnames(covs_)[colnames(covs_) %nin% colnames(make_full_rank(covs_[treat_ == 0, , drop = FALSE]))]

        covs_ <- covs_[, colnames(covs_) %nin% colinear.covs.to.remove, drop = FALSE]

        kbal.out <- do.call(KBAL::kbal, c(list(X = covs_, D = treat_), args))

        w[treat == i] <- (kbal.out$w / s.weights[subset])[treat_ == 0L]
        fit.list[[i]] <- kbal.out
      }
    }
    else if (estimand == "ATE") {
      w <- rep(1, length(treat))
      fit.list <- setNames(vector("list", nlevels(treat)), levels(treat))

      for (i in levels(treat)) {
        covs_i <- rbind(covs, covs[treat==i, , drop = FALSE])
        treat_i <- c(rep(1, nrow(covs)), rep(0, sum(treat==i)))

        colinear.covs.to.remove <- colnames(covs_i)[colnames(covs_i) %nin% colnames(make_full_rank(covs_i[treat_i == 0, , drop = FALSE]))]

        covs_i <- covs_i[, colnames(covs_i) %nin% colinear.covs.to.remove, drop = FALSE]

        covs_i[treat_i == 1,] <- covs_i[treat_i == 1,] * s.weights[subset] * sum(treat_i == 1) / sum(s.weights[subset])

        kbal.out_i <- do.call(KBAL::kbal, c(list(X = covs_i, D = treat_i), args))

        w[treat == i] <- kbal.out_i$w[treat_i == 0] / s.weights[subset][treat == i]
        fit.list[[i]] <- kbal.out_i
      }
    }
  }

  obj <- list(w = w)
  return(obj)

}

#Energy balancing
weightit2energy.cont <- function(covs, treat, s.weights, subset, missing, moments, int, ...) {
  A <- list(...)

  covs <- covs[subset, , drop = FALSE]
  treat <- treat[subset]
  sw <- s.weights[subset]

  sw <- sw/mean(sw)

  if (missing == "ind") {
    missing.ind <- apply(covs[, apply(covs, 2, anyNA), drop = FALSE], 2, function(x) as.numeric(is.na(x)))
    covs[is.na(covs)] <- 0
    covs <- cbind(covs, missing.ind)
  }

  treat <- (treat - w.m(treat, sw)) / sqrt(col.w.v(treat, sw))
  covs <- mat_div(center(covs, at = col.w.m(covs, sw)),
                  sqrt(col.w.v(covs, sw)))
  p <- ncol(covs)

  if (check.package("osqp")) {

    covs_dist <- as.matrix(dist(covs))
    covs_means <- colMeans(covs_dist)
    covs_grand_mean <- mean(covs_means)
    covs_A <- covs_dist + covs_grand_mean - outer(covs_means, covs_means, "+")

    treat_dist <- as.matrix(dist(treat))
    treat_means <- colMeans(treat_dist)
    treat_grand_mean <- mean(treat_means)
    treat_A <- treat_dist + treat_grand_mean - outer(treat_means, treat_means, "+")

    n <- length(treat)

    min.w <- if_null_then(A[["min.w"]], 1e-8)
    if (!is.numeric(min.w) || length(min.w) != 1 || min.w < 0) {
      warning("'min.w' must be a nonnegative number. Setting min.w = 1e-8.", call. = FALSE)
      min.w <- 1e-8
    }

    Pmat <- (covs_A * treat_A) %*% diag((sw/n)^2)

    Amat <- rbind(diag(n), sw, t(covs * sw)/n, treat * sw/n)
    lvec <- c(rep(min.w, n), n, rep(0, p), 0)
    uvec <- c(ifelse(check_if_zero(sw), min.w, Inf), n, rep(0, p), 0)

    if ((is_not_null(moments) && moments != 0) || int) {
      #Exactly balance correlations of moments and/or interactions
      covs <- cbind(covs, int.poly.f(covs, poly = moments, int = int))
      covs <- center(covs, at = col.w.m(covs, sw))

      Amat <- do.call("rbind", list(Amat,
                                    t(covs * treat * sw / n),
                                    if (moments > 1) t(covs[,-seq_len(p), drop = FALSE] * sw)/n))
      lvec <- do.call("c", list(lvec,
                                rep(0, ncol(covs)),
                                if (moments > 1) rep(0, ncol(covs)-p)))
      uvec <- do.call("c", list(uvec,
                                rep(0, ncol(covs)),
                                if (moments > 1) rep(0, ncol(covs)-p)))
    }

    if (is_not_null(A[["eps"]])) {
      if (is_null(A[["eps_abs"]])) A[["eps_abs"]] <- A[["eps"]]
      if (is_null(A[["eps_rel"]])) A[["eps_rel"]] <- A[["eps"]]
    }
    A[names(A) %nin% names(formals(osqp::osqpSettings))] <- NULL
    if (is_null(A[["max_iter"]])) A[["max_iter"]] <- 2E5L
    if (is_null(A[["eps_abs"]])) A[["eps_abs"]] <- 1E-8
    if (is_null(A[["eps_rel"]])) A[["eps_rel"]] <- 1E-8
    A[["verbose"]] <- TRUE

    options.list <- do.call(osqp::osqpSettings, A)

    opt.out <- do.call(osqp::solve_osqp, list(2*Pmat, A = Amat, l = lvec, u = uvec,
                                              pars = options.list),
                       quote = TRUE)

    w <- opt.out$x

    w[w <= min.w] <- min.w

    obj <- list(w = w, fit.obj = opt.out)
    return(obj)
  }
}

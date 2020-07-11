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

weightit2enet <- function(covs, treat, s.weights, subset, estimand, focal, stabilize, subclass, missing, moments, int, ...) {
  A <- list(...)

  covs <- covs[subset, , drop = FALSE]
  treat <- factor(treat[subset])
  s.weights <- s.weights[subset]

  if (anyNA(covs) && missing == "ind") {
    missing.ind <- apply(covs[, apply(covs, 2, anyNA), drop = FALSE], 2, function(x) as.numeric(is.na(x)))
    covs[is.na(covs)] <- 0
    covs <- cbind(covs, missing.ind)
  }

  covs <- apply(covs, 2, make.closer.to.1)
  model.covs <- cbind(covs, int.poly.f(covs, int = int, poly = moments))
  model.covs <- apply(model.covs, 2, make.closer.to.1)

  if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
  treat.type <- get.treat.type(treat)

  if (is_null(A[["stop.method"]])) {
    warning("No stop.method was provided. Using \"es.mean\".",
            call. = FALSE, immediate. = TRUE)
    A[["stop.method"]] <- "es.mean"
  }
  else if (length(A[["stop.method"]]) > 1) {
    warning("Only one stop.method is allowed at a time. Using just the first stop.method.",
            call. = FALSE, immediate. = TRUE)
    A[["stop.method"]] <- A[["stop.method"]][1]
  }

  cv <- 0
  available.stop.methods <- bal_criterion(treat.type, list = TRUE)
  s.m.matches <- charmatch(A[["stop.method"]], available.stop.methods)
  if (is.na(s.m.matches) || s.m.matches == 0L) {
    if (startsWith(A[["stop.method"]], "cv") && can_str2num(numcv <- substr(A[["stop.method"]], 3, nchar(A[["stop.method"]])))) {
      cv <- round(str2num(numcv))
      if (cv < 3) stop("At least 3 CV-folds must be specified in stop.method.", call. = FALSE)
    }
    else stop(paste0("'stop.method' must be one of ", word_list(c(available.stop.methods, "cv{#}"), "or", quotes = TRUE), "."), call. = FALSE)
  }
  else stop.method <- available.stop.methods[s.m.matches]

  tunable <- c("alpha", "relax", "type.multinomial", "reg.covs")

  trim.at <- if_null_then(A[["trim.at"]], 0)
  if (is_null(A[["alpha"]])) A[["alpha"]] <- 1 - .0001
  if (is_null(A[["thresh"]])) A[["thresh"]] <- 1e-7
  if (is_null(A[["maxit"]])) A[["maxit"]] <- 10^5
  if (is_null(A[["relax"]])) A[["relax"]] <- FALSE
  if (is_null(A[["reg.covs"]])) A[["reg.covs"]] <- TRUE
  gamma <- if (isTRUE(is_null(A[["relax"]]))) 0 else 1
  nlambda <- if_null_then(A[["nlambda"]], 5000)

  if (moments == 1 && !int && any(!A[["reg.covs"]])) {
    stop("If moments = 1 and int = FALSE (the default), 'reg.covs' cannot be FALSE.", call. = FALSE)
  }

  if (is_null(A[["lambda"]])) {
    lambda <- c(exp(seq(log(1/ncol(model.covs)), log(1/ncol(model.covs)/nlambda), length.out = nlambda - 1)), 0)
  }
  else {
    if (is.numeric(A[["lambda"]])) {
      lambda <- sort.int(A[["lambda"]], decreasing = TRUE)
      nlambda <- length(lambda)
    }
    else {
      stop("'lambda' must be a numeric vector.")
    }
  }

  if (treat.type == "binary")  {
    family <- "binomial"
    treat <- binarize(treat, one = focal)
    if (is_not_null(focal)) focal <- "1"
    A[["type.multinomial"]] <- NULL
  }
  else {
    family <- "multinomial"
    A[["type.multinomial"]] <- if_null_then(A[["type.multinomial"]], "ungrouped")
  }

  tune <- do.call("expand.grid", c(A[names(A) %in% tunable],
                                   list(stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)))
  if (cv == 0) {
    start.lambda <- if_null_then(A[["start.lambda"]], 1)
    if (is_null(A[["n.grid"]])) {
      n.grid <- round(1 + sqrt(2*(nlambda-start.lambda+1)))
    }
    else if (!is_(A[["n.grid"]], "numeric") || length(A[["n.grid"]]) > 1 ||
             !between(A[["n.grid"]], c(2, nlambda))) {
      stop(paste0("'n.grid' must be a numeric value between 2 and ", nlambda, "."), call. = FALSE)
    }
    else n.grid <- round(A[["n.grid"]])

    if (n.grid >= nlambda/3) n.grid <- nlambda
  }
  else {
    foldid <- sample(rep(seq_len(cv), length = length(treat)))
    type.measure <- if_null_then(A[["type.measure"]], "default")
    A[["type.measure"]] <-  match_arg(type.measure, formals(glmnet::cv.glmnet)[["type.measure"]])
  }

  current.best.loss <- Inf

  for (i in seq_len(nrow(tune))) {

    A[["penalty.factor"]] <- rep(1, ncol(model.covs))
    if (!tune[["reg.covs"]][i]) A[["penalty.factor"]][seq_len(ncol(covs))] <- 0

    if (cv == 0) {
      fit <- do.call(glmnet::glmnet, list(model.covs, treat, family = family, standardize = FALSE,
                                          lambda = lambda, alpha = tune[["alpha"]][i], thresh = A[["thresh"]],
                                          maxit = A[["maxit"]], relax = tune[["relax"]][i], weights = s.weights,
                                          penalty.factor = A[["penalty.factor"]],
                                          type.multinomial = tune[["type.multinomial"]][i]))

      if (treat.type == "binary") {
        treat <- binarize(treat, one = focal)
        if (is_not_null(focal)) focal <- "1"
      }

      crit <- bal_criterion(treat.type, stop.method)
      init <- crit$init(covs, treat, estimand = estimand, s.weights = s.weights, focal = focal)

      iters <- seq_along(fit$lambda)
      iters.grid <- round(seq(1, length(fit$lambda), length.out = n.grid))

      ps <- predict(fit, newx = model.covs, type = "response", s = fit$lambda[iters.grid], gamma = gamma)
      w <- get.w.from.ps(ps, treat = treat, estimand = estimand, focal = focal, stabilize = stabilize, subclass = subclass)
      if (trim.at != 0) w <- suppressMessages(apply(w, 2, trim, at = trim.at, treat = treat))

      iter.grid.balance <- apply(w, 2, function(w_) {
        crit$fun(init = init, weights = w_)
      })

      if (n.grid == nlambda) {
        best.lambda.index <- which.min(iter.grid.balance)
        best.loss <- iter.grid.balance[best.lambda.index]
        best.lambda <- lambda[best.lambda.index]
        lambda.val <- setNames(data.frame(fit$lambda,
                                          iter.grid.balance),
                               c("lambda", stop.method))
      }
      else {
        it <- which.min(iter.grid.balance) + c(-1, 1)
        it[1] <- iters.grid[max(1, it[1])]
        it[2] <- iters.grid[min(length(iters.grid), it[2])]
        iters.to.check <- iters[between(iters, iters[it])]

        if (is_null(iters.to.check) || anyNA(iters.to.check) || any(iters.to.check > nlambda)) stop("A problem has occurred")

        ps <- predict(fit, newx = model.covs, type = "response", s = fit$lambda[iters.to.check], gamma = gamma)
        w <- get.w.from.ps(ps, treat = treat, estimand = estimand, focal = focal, stabilize = stabilize, subclass = subclass)
        if (trim.at != 0) w <- suppressMessages(apply(w, 2, trim, at = trim.at, treat = treat))

        iter.grid.balance.fine <- apply(w, 2, function(w_) {
          crit$fun(init = init, weights = w_)
        })

        best.lambda.index <- which.min(iter.grid.balance.fine)
        best.loss <- iter.grid.balance.fine[best.lambda.index]
        best.lambda <- lambda[iters.to.check[best.lambda.index]]
        lambda.val <- setNames(data.frame(c(fit$lambda[iters.grid], fit$lambda[iters.to.check]),
                                          c(iter.grid.balance, iter.grid.balance.fine)),
                               c("lambda", stop.method))
      }

      lambda.val <- unique(lambda.val[order(lambda.val$lambda),])
      w <- w[,best.lambda.index]
      ps <- if (treat.type == "binary") ps[,best.lambda.index] else NULL

      tune[[paste.("best", stop.method)]][i] <- best.loss
      tune[["best.lambda"]][i] <- best.lambda

      if (best.loss < current.best.loss) {
        best.fit <- fit
        best.w <- w
        best.ps <- ps
        current.best.loss <- best.loss
        best.tune.index <- i

        info <- list(best.lambda = best.lambda,
                     lambda.val = lambda.val,
                     coef = predict(fit, type = "coef", s = best.lambda))
      }

    }
    else {
      fit <- do.call(glmnet::cv.glmnet, list(model.covs, treat, family = family, standardize = FALSE,
                                             lambda = lambda, alpha = tune[["alpha"]][i], thresh = A[["thresh"]],
                                             maxit = A[["maxit"]], relax = tune[["relax"]][i], weights = s.weights,
                                             penalty.factor = A[["penalty.factor"]],
                                             type.multinomial = tune[["type.multinomial"]][i],
                                             foldid = foldid))

      best.lambda.index <- which.min(fit$cvm)
      best.lambda <- fit$lambda[best.lambda.index]
      best.loss <- fit$cvm[best.lambda.index]

      tune[[paste.("best", names(fit$name))]][i] <- best.loss
      tune[["best.lambda"]][i] <- best.lambda

      if (best.loss < current.best.loss) {
        best.fit <- fit
        best.ps <- predict(fit, newx = model.covs, type = "response", s = best.lambda, gamma = gamma)
        best.w <- drop(get.w.from.ps(best.ps, treat = treat, estimand = estimand, focal = focal, stabilize = stabilize, subclass = subclass))
        current.best.loss <- best.loss
        best.tune.index <- i

        lambda.val <- setNames(data.frame(fit$lambda,
                                          fit$cvm),
                               c("lambda", names(fit$name)))

        info <- list(best.lambda = best.lambda,
                     lambda.val = lambda.val,
                     coef = predict(fit, type = "coef", s = best.lambda))

        if (treat.type == "multinomial") best.ps <- NULL
      }
    }
  }

  tune[tunable[vapply(tunable, function(x) length(A[[x]]) == 1, logical(1L))]] <- NULL

  if (ncol(tune) > 1) {
    info[["tune"]] <- tune
    info[["best.tune"]] <- tune[best.tune.index,]
  }

  obj <- list(w = best.w, ps = best.ps, info = info, fit.obj = best.fit)
  return(obj)

}

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

#New methods and functions

#------Preliminary template----
weightit2XXX <- function(covs, treat...) {
  stop("method = \"XXX\" isn't ready to use yet.", call. = FALSE)
}
#------Template----
weightit2XXX <- function(covs, treat, s.weights, subset, estimand, focal, moments, int, ...) {
  A <- list(...)

  covs <- covs[subset, , drop = FALSE]
  treat <- factor(treat[subset])

  covs <- cbind(covs, int.poly.f(covs, poly = moments, int = int))
  covs <- apply(covs, 2, make.closer.to.1)

  if (any(vars.w.missing <- apply(covs, 2, function(x) anyNA(x)))) {
    missing.ind <- apply(covs[, vars.w.missing, drop = FALSE], 2, function(x) as.numeric(is.na(x)))
    covs[is.na(covs)] <- 0
    covs <- cbind(covs, missing.ind)
  }

  new.data <- data.frame(treat, covs)
  new.formula <- formula(new.data)

  for (f in names(formals(PACKAGE::FUNCTION))) {
    if (is_null(A[[f]])) A[[f]] <- formals(PACKAGE::FUNCTION)[[f]]
  }
  A[names(A) %in% names(formals(weightit2optweight))] <- NULL

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

#Subgroup Balancing PS
weightit2sbps <- function(covs, treat, s.weights, subset, estimand, focal, stabilize, ...) {
  A <- list(...)

  fit.obj <- NULL

  covs <- covs[subset, , drop = FALSE]
  t <- factor(treat[subset])

  if (!is_binary(t)) stop("Subgroup balancing propensity score weighting is not yet compatible with non-binary treatments.", call. = FALSE)

  if (any(vars.w.missing <- apply(covs, 2, function(x) anyNA(x)))) {
    missing.ind <- apply(covs[, vars.w.missing, drop = FALSE], 2, function(x) as.numeric(is.na(x)))
    covs[is.na(covs)] <- 0
    covs <- cbind(covs, missing.ind)
  }
  covs <- apply(covs, 2, make.closer.to.1)

  smd <- function(x, t, w, estimand, std = TRUE) {
    m <- vapply(levels(t), function(t.lev) w.m(x[t==t.lev], w = w[t==t.lev]), numeric(1L))
    mdiff <- abs(diff(m))

    if (check_if_zero(mdiff)) return(0)
    else {
      if (!std) sd <- 1
      else if (estimand == "ATT") sd <- sd(x[t==1])
      else if (estimand == "ATC") sd <- sd(x[t==0])
      else sd <- sqrt(.5 * (var(x[t==1]) + var(x[t==0])))
      return(mdiff/sd)
    }
  }

  loss <- A[["loss"]]
  loss <- match_arg(loss, c("weighting", "matching"))

  if (loss == "matching") {
    F_ <- function(covs, sub, t, w) {
      #Overall Balance of covs
      Mk <- apply(covs, 2, function(x) smd(x, t, w, estimand))
      #Subgroup Balance
      Mkr <- unlist(lapply(levels(sub), function(s) {apply(covs, 2,
                                                           function(x) smd(x[sub==s], t[sub==s], w[sub==s], estimand))}))
      return(sum(c(Mk, Mkr) ^ 2))
    }
  }
  else if (loss == "weighting") {
    F_ <- function(covs, sub, t, w) {
      #Overall Balance of covs
      Mk <- apply(covs, 2, function(x) smd(x, t, w, estimand))
      #Overall balance of subgroups
      Mr <- vapply(levels(sub), function(s) {smd(as.numeric(sub == s), t, w, std = FALSE)}, numeric(1L))
      #Subgroup Balance
      Mkr <- unlist(lapply(levels(sub), function(s) {apply(covs, 2,
                                                           function(x) smd(x[sub==s], t[sub==s], w[sub==s], estimand))}))
      return(sum(c(Mk, Mr, Mkr) ^ 2))
    }
  }
  else stop()

  #Process subgroup
  subgroup <- process.by(by = A[["subgroup"]], data = covs, treat = t, by.arg = "subgroup")$by.factor

  overall.weights <- subgroup.weights <- NULL
  if (is_not_null(A[["overall.ps"]])) {
    if ((is.matrix(A[["overall.ps"]]) || is.data.frame(A[["overall.ps"]])) &&
        ncol(A[["overall.ps"]]) == nlevels(t) && all(colnames(A[["overall.ps"]] %in% levels(t)))) {
      ps.mat <- A[["overall.ps"]]
    }
    else if (is.numeric(A[["overall.ps"]])) {
      ps.mat <- matrix(NA_real_, nrow = length(t), ncol = nlevels(t), dimnames = list(NULL, levels(t)))
      ps.mat[, 2] <- A[["overall.ps"]]
      ps.mat[, 1] <- 1 - A[["overall.ps"]]
    }
    else {
      stop()
    }
    overall.weights <- get_w_from_ps(ps.mat, t, estimand, focal)
  }
  if (is_not_null(A[["subgroup.ps"]])) {
    if ((is.matrix(A[["subgroup.ps"]]) || is.data.frame(A[["subgroup.ps"]])) &&
        ncol(A[["subgroup.ps"]]) == nlevels(t) && all(colnames(A[["subgroup.ps"]] %in% levels(t)))) {
      ps.mat <- A[["subgroup.ps"]]
    }
    else if (is.numeric(A[["subgroup.ps"]])) {
      ps.mat <- matrix(NA_real_, nrow = length(t), ncol = nlevels(t), dimnames = list(NULL, levels(t)))
      ps.mat[, 2] <- A[["subgroup.ps"]]
      ps.mat[, 1] <- 1 - A[["subgroup.ps"]]
    }
    else {
      stop()
    }
    subgroup.weights <- get_w_from_ps(ps.mat, t, estimand, focal)
  }

  if (is_not_null(A[["overall.weights"]])) {
    if (!is.numeric(A[["overall.weights"]])) {
      stop()
    }
    overall.weights <- A[["overall.weights"]]
  }
  if (is_not_null(A[["subgroup.weights"]])) {
    if (!is.numeric(A[["subgroup.weights"]])) {
      stop()
    }
    subgroup.weights <- A[["subgroup.weights"]]
  }

  if (is_null(overall.weights) || is_null(subgroup.weights)) {
    #Process w.method
    w.method <- A[["w.method"]]
    check.acceptable.method(w.method, msm = FALSE, force = FALSE)

    if (is.character(w.method)) {
      w.method <- method.to.proper.method(w.method)
      attr(w.method, "name") <- w.method
    }
    else if (is.function(w.method)) {
      w.method.name <- paste(deparse(substitute(w.method)))
      check.user.method(w.method)
      attr(w.method, "name") <- w.method.name
    }

    if (loss == "matching") {
      t.bin <- binarize(t)
      overall.fit <- weightit.fit(covs = covs, treat = t, method = "ps",
                                  treat.type = "binary", s.weights = s.weights,
                                  by.factor = factor(rep(1, length(t))), estimand = estimand,
                                  focal = focal, stabilize = stabilize,
                                  ps = NULL, moments = 1, int = FALSE)
      overall.ps <- overall.fit$ps
      overall.match <- Matching::Match(Tr = t.bin, X = matrix(c(overall.ps, as.numeric(subgroup)), ncol = 2),
                                       estimand = estimand, caliper = .25,
                                       M = 1, replace = FALSE, exact = c(FALSE, TRUE), ties = TRUE)
      overall.weights <- cobalt::get.w(overall.match)

      subgroup.fit <- weightit.fit(covs = covs, treat = t, method = "ps",
                                   treat.type = "binary", s.weights = s.weights,
                                   by.factor = subgroup, estimand = estimand,
                                   focal = focal, stabilize = stabilize,
                                   ps = NULL, moments = 1, int = FALSE)
      subgroup.ps <- subgroup.fit$ps
      subgroup.match <- Matching::Match(Tr = t.bin, X = matrix(c(subgroup.ps, as.numeric(subgroup)), ncol = 2),
                                        estimand = estimand, caliper = .25,
                                        M = 1, replace = FALSE, exact = c(FALSE, TRUE), ties = TRUE)
      subgroup.weights <- cobalt::get.w(subgroup.match)
    }
    if (loss == "weighting") {
      #Estimate overall weights
      overall.fit <- weightit.fit(covs = covs, treat = t, method = w.method,
                                  treat.type = "binary", s.weights = s.weights,
                                  by.factor = factor(rep(1, length(t))), estimand = estimand,
                                  focal = focal, stabilize = stabilize,
                                  ps = NULL, moments = 1, int = FALSE)
      overall.weights <- overall.fit$w
      #Estimate subgroup weights
      subgroup.fit <- weightit.fit(covs = covs, treat = t, method = w.method,
                                   treat.type = "binary", s.weights = s.weights,
                                   by.factor = subgroup, estimand = estimand,
                                   focal = focal, stabilize = stabilize,
                                   ps = NULL, moments = 1, int = FALSE)
      subgroup.weights <- subgroup.fit$w
    }


  }

  #Find combinations that minimize loss
  n.subgroups <- nunique(subgroup)
  if (n.subgroups > 8) {
    #Stochastic search
    L1 <- 10
    L2 <- 5
    S_ <- setNames(rep("overall", nlevels(subgroup)),
                   levels(subgroup))
    rep <- 0
    no.change.streak <- 0
    current.loss <- Inf
    while (rep <= L1 && no.change.streak <= L2) {
      rep <- rep + 1
      if (is_null(get0("last.loss"))) last.loss <- Inf
      else last.loss <- current.loss

      rand.subs <- sample(levels(subgroup))
      S__ <- setNames(sample(c("overall", "subgroup"), length(S_), replace = TRUE), rand.subs)

      for (i in 1:length(S__)) {
        S__[i] <- "overall"
        to.overall <- subgroup %in% rand.subs[S__[rand.subs] == "overall"]
        w_ <- subgroup.weights
        w_[to.overall] <- overall.weights[to.overall]
        loss.o <- F_(covs, subgroup, t, w_)

        S__[i] <- "subgroup"
        to.overall <- subgroup %in% rand.subs[S__[rand.subs] == "overall"]
        w_ <- subgroup.weights
        w_[to.overall] <- overall.weights[to.overall]
        loss.s <- F_(covs, subgroup, t, w_)

        if (loss.o < loss.s) {
          S__[i] <- "overall"
          if (loss.o < current.loss) {
            current.loss <- loss.o
            attr(current.loss, "S") <- S__
          }
        }
        else {
          S__[i] <- "subgroup"
          if (loss.s < current.loss) {
            current.loss <- loss.s
            attr(current.loss, "S") <- S__
          }
        }



      }

      to.overall <- subgroup %in% rand.subs[S__[rand.subs] == "overall"]
      w_ <- subgroup.weights
      w_[to.overall] <- overall.weights[to.overall]
      current.loss <- F_(covs, subgroup, t, w_)
      if (check_if_zero(current.loss - last.loss)) no.change.streak <- no.change.streak + 1
      print(current.loss)
      print(S__)
    }

    best.S <- attr(current.loss, "S")
    to.overall <- subgroup %in% rand.subs[best.S[rand.subs] == "overall"]
    w <- subgroup.weights
    w[to.overall] <- overall.weights[to.overall]

  }
  else {
    S <- setNames(do.call("expand.grid", lapply(integer(n.subgroups), function(x) (c("overall", "subgroup")))),
                  levels(subgroup))
    print(S)
    w.list <<- lapply(seq_len(nrow(S)), function(i) {
      to.overall <- subgroup %in% levels(subgroup)[S[i, levels(subgroup)] == "overall"]
      w_ <- subgroup.weights
      w_[to.overall] <- overall.weights[to.overall]
      return(w_)
    })

    loss.val <- vapply(w.list, function(w_) F_(covs, subgroup, t, w_), numeric(1L))
    best.loss <- which.min(loss.val)
    w <- w.list[[best.loss]]

    if (is_not_null(overall.fit$ps)) {
      to.overall <- subgroup %in% levels(subgroup)[S[best.loss, levels(subgroup)] == "overall"]
      p.score <- subgroup.fit$ps
      p.score[to.overall] <- overall.fit$ps[to.overall]
    }
    else p.score <- NULL

  }

  obj <- list(w = w
              , ps = p.score
              #, fit.obj = fit.obj
  )
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

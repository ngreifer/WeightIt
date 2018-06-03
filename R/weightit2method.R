#Propensity score estimation with regression
weightit2ps <- function(covs, treat, s.weights, subset, estimand, focal, stabilize, ps, ...) {
  A <- list(...)

  if (is_null(ps)) {

    covs <- covs[subset, , drop = FALSE]
    t <- factor(treat)[subset]

    if (any(vars.w.missing <- apply(covs, 2, function(x) any(is.na(x))))) {
      missing.ind <- apply(covs[, vars.w.missing, drop = FALSE], 2, function(x) as.numeric(is.na(x)))
      covs[is.na(covs)] <- 0
      covs <- cbind(covs, missing.ind)
    }
    covs <- apply(covs, 2, make.closer.to.1)
    colinear.covs.to.remove <- colnames(covs)[colnames(covs) %nin% colnames(remove.collinearity(covs))]

    covs <- covs[, colnames(covs) %nin% colinear.covs.to.remove, drop = FALSE]

    if (is_null(A$link)) A$link <- "logit"
    else {
      if (nunique.gt(t, 2)) {
        acceptable.links <- c("logit", "probit", "bayes.probit")
        which.link <- acceptable.links[pmatch(A$link, acceptable.links, nomatch = 0)][1]
        if (is.na(which.link)) {
          A$link <- "logit"
          warning("Only \"logit\",\"probit\" and \"bayes.probit\" are allowed as links for multinomial treatments. Using link = \"logit\".",
                  call. = FALSE)
        }
        else A$link <- which.link
      }
      else {
        acceptable.links <- c("logit", "probit", "cloglog")
        which.link <- acceptable.links[pmatch(A$link, acceptable.links, nomatch = 0)][1]
        if (is.na(which.link)) {
          A$link <- "logit"
          warning("Only \"logit\",\"probit\" and \"cloglog\" are allowed as links for binary treatments. Using link = \"logit\".",
                  call. = FALSE)
        }
        else A$link <- which.link
      }

    }

    if (!nunique.gt(t, 2)) {
      data <- data.frame(t, covs)
      formula <- formula(data)

      ps <- setNames(as.data.frame(matrix(NA_real_, ncol = nunique(t), nrow = length(t))),
                     levels(t))

      fit <- glm(formula, data = data,
                 weights = s.weights[subset],
                 family = binomial(link = A$link),
                 control = list(),
                 ...)
      ps[[2]] <- p.score <- fit$fitted.values
      ps[[1]] <- 1 - ps[[2]]

    }
    else {
      if (A$link %in% c("logit", "probit")) {
        if (check.package("mlogit", alternative = TRUE) && (is_null(A$use.mlogit) || A$use.mlogit == TRUE)) {
          message(paste0("Using multinomial ", A$link, " regression."))
          data <- data.frame(t = t , s.weights = s.weights[subset], covs)
          covnames <- names(data)[-c(1,2)]
          mult <- mlogit::mlogit.data(data, varying = NULL, shape = "wide", sep = "", choice = "t")
          tryCatch({fit <- mlogit::mlogit(as.formula(paste0("t ~ 1 | ", paste(covnames, collapse = " + "),
                                                            " | 1")), data = mult, estimate = TRUE,
                                          probit = ifelse(A$link[1] == "probit", TRUE, FALSE),
                                          weights = s.weights, ...)},
                   error = function(e) {stop(paste0("There was a problem fitting the multinomial ", A$link, " regressions with mlogit().\n       Try again with use.mlogit = FALSE."), call. = FALSE)}
          )
          ps <- fitted(fit, outcome = FALSE)
        }
        else {
          message(paste0("Using a series of ", nunique(t), " binomial ", A$link, " regressions."))
          ps <- setNames(as.data.frame(matrix(NA_real_, ncol = nunique(t), nrow = length(t))),
                         levels(t))

          for (i in levels(t)) {
            t_i <- rep(0, length(t)); t_i[t == i] <- 1
            data_i <- data.frame(t_i, covs)
            fit_i <- glm(formula(data_i), data = data_i,
                         family = binomial(link = A$link),
                         weights = s.weights[subset])
            ps[[i]] <- fit_i$fitted.values
          }
        }
      }
      else if (A$link == "bayes.probit") {
        check.package("MNP")
        data <- data.frame(t, covs)
        formula <- formula(data)
        tryCatch({fit <- MNP::mnp(formula, data, verbose = TRUE)},
                 error = function(e) stop("There was a problem with the Bayes probit regression. Try a different link.", call. = FALSE))
        ps <- MNP::predict.mnp(fit, type = "prob")$p
      }
      else {
        stop('link must be "logit", "probit", or "bayes.probit".', call. = FALSE)
      }
      p.score <- NULL
    }
  }
  else {
    if (ps) {
      ps <- setNames(as.data.frame(matrix(c(1-ps, ps), ncol = 2)),
                     levels(t))
      p.score <- ps
    }
  }

  #ps should be matrix of probs for each treat
  #Computing weights
  w <- rep(0, nrow(ps))
  for (i in seq_len(nunique(t))) {
    w[t == levels(t)[i]] <- 1/ps[t == levels(t)[i], i]
  }

  if (toupper(estimand) == "ATE") {
    w <- w
  }
  else if (toupper(estimand) == "ATT") {
    w <- w*ps[, levels(t) == focal]
  }
  else if (toupper(estimand) == "ATO") {
    w <- w*apply(ps, 1, prod)
  }
  else if (toupper(estimand) == "ATM") {
    w <- w*apply(ps, 1, min)
  }
  else w <- NULL

  if (stabilize) {
    w <- w * sapply(t, function(x) sum(t==x) / sum(1*(t==x)*w))
  }

  obj <- list(w = w, ps = p.score)
  return(obj)
}
weightit2ps.cont <- function(covs, treat, s.weights, subset, stabilize, ps, ...) {
  A <- list(...)
  if (is_null(A$link)) A$link <- "identity"
  if (is_null(A$family)) A$family <- gaussian(link = A$link)

  covs <- covs[subset, , drop = FALSE]
  t <- treat[subset]

  if (any(vars.w.missing <- apply(covs, 2, function(x) any(is.na(x))))) {
    missing.ind <- apply(covs[, vars.w.missing, drop = FALSE], 2, function(x) as.numeric(is.na(x)))
    covs[is.na(covs)] <- 0
    covs <- cbind(covs, missing.ind)
  }
  covs <- apply(covs, 2, make.closer.to.1)
  data <- data.frame(t, covs)
  formula <- formula(data)

  stabilize <- TRUE

  if (is_null(ps)) {
    fit <- glm(formula, data = data,
               weights = s.weights[subset],
               family = A$family,
               control = list(),
               ...)
    p.denom <- fit$fitted.values

    if (isTRUE(A[["use.kernel"]])) {
      if (is_null(A[["bw"]])) A[["bw"]] <- "nrd0"
      if (is_null(A[["adjust"]])) A[["adjust"]] <- 1
      if (is_null(A[["kernel"]])) A[["kernel"]] <- "gaussian"
      if (is_null(A[["n"]])) A[["n"]] <- 10*length(t)

      d.d <- density(t - p.denom, n = A[["n"]],
                     weights = s.weights[subset]/sum(s.weights[subset]), give.Rkern = FALSE,
                     bw = A[["bw"]], adjust = A[["adjust"]], kernel = A[["kernel"]])
      if (isTRUE(A[["plot"]])) plot(d.d, main = "Denominator density")
      dens.denom <- with(d.d, approxfun(x = x, y = y))(t)
    }
    else {
      dens.denom <- dnorm(t, mean = p.denom, sd = sqrt(summary(fit)$dispersion))
    }

    if (stabilize) {
      num.fit <- glm(t ~ 1,
                     data = data.frame(t = t),
                     weights = s.weights[subset],
                     family = A$family,
                     control = list(), ...)
      p.num <- num.fit$fitted.values

      if (isTRUE(A[["use.kernel"]])) {
        d.n <- density(t - p.num, n = A[["n"]],
                       weights = s.weights[subset]/sum(s.weights[subset]), give.Rkern = FALSE,
                       bw = A[["bw"]], adjust = A[["adjust"]], kernel = A[["kernel"]])
        if (isTRUE(A[["plot"]])) plot(d.n, main = "Numerator density")
        dens.num <- with(d.n, approxfun(x = x, y = y))(t)
      }
      else {
        dens.num <- dnorm(t, p.num, sqrt(summary(num.fit)$dispersion))
      }
      w <- dens.num/dens.denom
    }
    else {
      w <- 1/dens.denom
    }
  }

  obj <- list(ps = p.denom,
              w = w)
  return(obj)
}

#Generalized boosted modeling with twang
weightit2gbm <- function(covs, treat, s.weights, estimand, focal, subset, stabilize, ...) {
  A <- list(...)
  if (is_null(A$stop.method)) {
    warning("No stop.method was provided. Using \"es.mean\".",
            call. = FALSE, immediate. = TRUE)
    A$stop.method <- "es.mean"
  }
  else if (length(A$stop.method) > 1) {
    warning("Only one stop.method is allowed at a time. Using just the first stop.method.",
            call. = FALSE, immediate. = TRUE)
    A$stop.method <- A$stop.method[1]
  }

  available.stop.methods <- c("ks.mean", "es.mean", "ks.max", "es.max")
  s.m.matches <- charmatch(A[["stop.method"]], available.stop.methods)
  if (is.na(s.m.matches) || s.m.matches == 0L) {stop(paste0("stop.method must be one of ", word.list(available.stop.methods, "or", quotes = TRUE), "."), call. = FALSE)}
  else A[["stop.method"]] <- available.stop.methods[s.m.matches]

  for (f in names(formals(twang::ps))) {
    if (is_null(A[[f]])) A[[f]] <- formals(twang::ps)[[f]]
  }

  covs <- covs[subset, , drop = FALSE]
  treat <- factor(treat[subset])

  covs <- apply(covs, 2, make.closer.to.1)

  if (check.package("twang")) {
    if (estimand == "ATT") {
      w <- rep(1, length(treat))
      control.levels <- levels(treat)[levels(treat) != focal]

      for (i in control.levels) {
        treat.in.i.focal <- treat %in% c(focal, i)
        treat_ <- ifelse(treat[treat.in.i.focal] == i, 0, 1)
        covs_ <- covs[treat.in.i.focal, , drop = FALSE]
        new.data <- data.frame(treat_, covs_)

        fit <- twang::ps(formula(new.data),
                         data = new.data,
                         estimand = "ATT",
                         stop.method = A[["stop.method"]],
                         sampw = s.weights[subset][treat.in.i.focal],
                         verbose = TRUE,
                         print.level = 2,
                         n.trees = A[["n.trees"]],
                         interaction.depth = A[["interaction.depth"]],
                         shrinkage = A[["shrinkage"]],
                         bag.fraction = A[["bag.fraction"]],
                         perm.test.iters = A[["perm.test.iters"]],
                         iterlim = A[["iterlim"]],
                         multinom = FALSE)

        s <- fit$stopMethods[1]

        w[treat == i] <- cobalt::get.w(fit, stop.method = s)[treat_ == 0]

      }
    }
    else if (estimand == "ATE") {
      w <- rep(1, length(treat))

      for (i in levels(treat)) {
        #Mimicking twang
        #Seeks balance between weighted treat group and all others combined
        #Note: Gives different answer than twang for some reason; has better balance though.
        treat_i <- ifelse(treat == i, 0, 1)
        new.data <- data.frame(treat_i, covs)

        fit <- twang::ps(formula(new.data),
                         data = new.data,
                         estimand = "ATE",
                         stop.method = A[["stop.method"]],
                         sampw = s.weights[subset],
                         verbose = TRUE,
                         print.level = 2,
                         n.trees = A[["n.trees"]],
                         interaction.depth = A[["interaction.depth"]],
                         shrinkage = A[["shrinkage"]],
                         bag.fraction = A[["bag.fraction"]],
                         perm.test.iters = A[["perm.test.iters"]],
                         iterlim = A[["iterlim"]],
                         multinom = FALSE)

        s <- fit$stopMethods[1]

        if (nlevels(treat) == 2) {
          w <- cobalt::get.w(fit, stop.method = s)
          break
        }
        else {
          w[treat == i] <- cobalt::get.w(fit, stop.method = s)[treat == i]
        }
      }
    }

    if (stabilize) {
      w <- w * sapply(t, function(x) sum(t==x) / sum(1*(t==x)*w))
    }
  }
  out <- list(w = w)
}
weightit2gbm.cont <- function(covs, treat, s.weights, subset, stabilize, ...) {
  A <- list(...)
  if (is_null(A$stop.method)) {
    warning("No stop.method was provided. Using \"s.mean.z\".",
            call. = FALSE, immediate. = TRUE)
    A$stop.method <- "s.mean.z"
  }
  else if (length(A$stop.method) > 1) {
    warning("Only one stop.method is allowed at a time. Using just the first stop.method.",
            call. = FALSE, immediate. = TRUE)
    A$stop.method <- A$stop.method[1]
  }

  for (f in names(formals(ps.cont))) {
    if (is_null(A[[f]])) A[[f]] <- formals(ps.cont)[[f]]
  }

  covs <- covs[subset, , drop = FALSE]
  treat <- treat[subset]

  covs <- apply(covs, 2, make.closer.to.1)

  new.data <- data.frame(treat, covs)

  if (check.package("wCorr") && check.package("gbm")) {
    fit <- ps.cont(formula(new.data), data = new.data,
                   n.trees = A[["n.trees"]],
                   interaction.depth = A[["interaction.depth"]],
                   shrinkage = A[["shrinkage"]],
                   bag.fraction = A[["bag.fraction"]],
                   stop.method = A[["stop.method"]],
                   use.optimize = A[["use.optimize"]],
                   sampw = s.weights[subset],
                   verbose = TRUE)
    w <- cobalt::get.w(fit, stop.method = A[["stop.method"]])
  }

  #ps <- fit[["ps"]][[A[["stop.method"]]]]

  out <- list(w = w)
}

#CBPS
weightit2cbps <- function(covs, treat, s.weights, estimand, focal, subset, stabilize, ...) {
  A <- list(...)

  covs <- covs[subset, , drop = FALSE]
  treat <- factor(treat)[subset]

  if (any(vars.w.missing <- apply(covs, 2, function(x) any(is.na(x))))) {
    missing.ind <- apply(covs[, vars.w.missing, drop = FALSE], 2, function(x) as.numeric(is.na(x)))
    covs[is.na(covs)] <- 0
    covs <- cbind(covs, missing.ind)
  }
  covs <- apply(covs, 2, make.closer.to.1)

  if (check.package("CBPS")) {
    if (estimand == "ATT") {
      w <- rep(1, length(treat))
      control.levels <- levels(treat)[levels(treat) != focal]

      for (i in control.levels) {
        treat.in.i.focal <- treat %in% c(focal, i)
        treat_ <- ifelse(treat[treat.in.i.focal] == i, 0, 1)
        covs_ <- covs[treat.in.i.focal, , drop = FALSE]
        new.data <- data.frame(treat_, covs_)

        tryCatch({fit <- CBPS::CBPS(formula(new.data),
                                    data = new.data,
                                    method = if (is_null(A$over) || A$over == TRUE) "over" else "exact",
                                    standardize = FALSE,
                                    sample.weights = s.weights[subset][treat.in.i.focal],
                                    ATT = 1,
                                    ...)},
                 error = function(e) {
                   e. <- conditionMessage(e)
                   e. <- gsub("method = \"exact\"", "over = FALSE", e., fixed = TRUE)
                   stop(e., call. = FALSE)
                 }
        )

        w[treat == i] <- cobalt::get.w(fit, estimand = "ATT")[treat_ == 0] / s.weights[subset][treat.in.i.focal][treat_ == 0]

      }
    }
    else {
      new.data <- data.frame(treat, covs)
      if (nunique(treat) <= 4) {
        tryCatch({fit <- CBPS::CBPS(formula(new.data),
                                    data = new.data,
                                    method = if (is_null(A$over) || A$over == TRUE) "over" else "exact",
                                    standardize = FALSE,
                                    sample.weights = s.weights[subset],
                                    ATT = 0,
                                    ...)},
                 error = function(e) {
                   e. <- conditionMessage(e)
                   e. <- gsub("method = \"exact\"", "over = FALSE", e., fixed = TRUE)
                   stop(e., call. = FALSE)
                 }
        )

        w <- cobalt::get.w(fit, estimand = "ATE") / s.weights[subset]
      }
      else {
        w <- rep(1, length(treat))
        for (i in levels(treat)) {
          new.data[[1]] <- ifelse(treat == i, 1, 0)
          fit <- CBPS::CBPS(formula(new.data), data = new.data,
                            method = if (is_null(A$over) || A$over == TRUE) "over" else "exact",
                            standardize = FALSE,
                            sample.weights = s.weights[subset],
                            ATT = 0, ...)

          w[treat==i] <- cobalt::get.w(fit, estimand = "ATE")[treat==i] / s.weights[subset][treat==i]
        }
      }
    }

  }
  if (stabilize) {
    w <- w * sapply(t, function(x) sum(t==x) / sum(1*(t==x)*w))
  }

  obj <- list(w = w)
}
weightit2cbps.cont <- function(covs, treat, s.weights, subset, ...) {
  A <- list(...)

  covs <- covs[subset, , drop = FALSE]
  treat <- treat[subset]

  if (any(vars.w.missing <- apply(covs, 2, function(x) any(is.na(x))))) {
    missing.ind <- apply(covs[, vars.w.missing, drop = FALSE], 2, function(x) as.numeric(is.na(x)))
    covs[is.na(covs)] <- 0
    covs <- cbind(covs, missing.ind)
  }
  covs <- apply(covs, 2, make.closer.to.1)

  new.data <- data.frame(treat = treat, covs)

  if (check.package("CBPS")) {
    tryCatch({fit <- CBPS::CBPS(formula(new.data),
                                data = new.data,
                                method = ifelse(is_null(A$over) || A$over == TRUE, "over", "exact"),
                                standardize = FALSE,
                                sample.weights = s.weights[subset],
                                ...)},
             error = function(e) {
               e. <- conditionMessage(e)
               e. <- gsub("method = \"exact\"", "over = FALSE", e., fixed = TRUE)
               stop(e., call. = FALSE)
             }
    )
  }
  w <- cobalt::get.w(fit) #/ s.weights[subset]

  obj <- list(w = w)
}
weightit2npcbps <- function(covs, treat, s.weights, subset, ...) {
  A <- list(...)
  if (!all_the_same(s.weights)) stop(paste0("Sampling weights cannot be used with method = \"npcbps\"."),
                                     call. = FALSE)

  covs <- covs[subset, , drop = FALSE]
  treat <- factor(treat)[subset]

  if (any(vars.w.missing <- apply(covs, 2, function(x) any(is.na(x))))) {
    missing.ind <- apply(covs[, vars.w.missing, drop = FALSE], 2, function(x) as.numeric(is.na(x)))
    covs[is.na(covs)] <- 0
    covs <- cbind(covs, missing.ind)
  }
  covs <- apply(covs, 2, make.closer.to.1)
  colinear.covs.to.remove <- colnames(covs)[colnames(covs) %nin% colnames(remove.collinearity(covs))]
  covs <- covs[, colnames(covs) %nin% colinear.covs.to.remove, drop = FALSE]

  new.data <- data.frame(treat = treat, covs)

  if (check.package("CBPS")) {
    fit <- CBPS::npCBPS(formula(new.data), data = new.data, print.level = 1, ...)
  }
  w <- cobalt::get.w(fit)

  obj <- list(w = w)

  return(obj)
}
weightit2npcbps.cont <- function(covs, treat, s.weights, subset, estimand, ...) {
  A <- list(...)

  if (!all_the_same(s.weights)) stop(paste0("Sampling weights cannot be used with method = \"npcbps\"."),
                                     call. = FALSE)

  covs <- covs[subset, , drop = FALSE]
  treat <- treat[subset]

  if (any(vars.w.missing <- apply(covs, 2, function(x) any(is.na(x))))) {
    missing.ind <- apply(covs[, vars.w.missing, drop = FALSE], 2, function(x) as.numeric(is.na(x)))
    covs[is.na(covs)] <- 0
    covs <- cbind(covs, missing.ind)
  }
  covs <- apply(covs, 2, make.closer.to.1)
  colinear.covs.to.remove <- colnames(covs)[colnames(covs) %nin% colnames(remove.collinearity(covs))]
  covs <- covs[, colnames(covs) %nin% colinear.covs.to.remove, drop = FALSE]

  new.data <- data.frame(treat = treat, covs)

  if (check.package("CBPS")) {
    fit <- CBPS::npCBPS(formula(new.data), data = new.data, print.level = 1, ...)
  }
  w <- cobalt::get.w(fit)

  obj <- list(w = w)

  return(obj)
}

#Entropy balancing with ebal
weightit2ebal <- function(covs, treat, s.weights, subset, estimand, focal, stabilize, moments, int, ...) {
  A <- list(...)

  covs <- covs[subset, , drop = FALSE]
  treat <- factor(treat)[subset]

  covs <- cbind(covs, int.poly.f(covs, poly = moments, int = int))
  covs <- apply(covs, 2, make.closer.to.1)

  if (any(vars.w.missing <- apply(covs, 2, function(x) any(is.na(x))))) {
    missing.ind <- apply(covs[, vars.w.missing, drop = FALSE], 2, function(x) as.numeric(is.na(x)))
    covs[is.na(covs)] <- 0
    covs <- cbind(covs, missing.ind)
  }

  for (f in names(formals(ebal::ebalance))) {
    if (is_null(A[[f]])) A[[f]] <- formals(ebal::ebalance)[[f]]
  }
  if (stabilize) {
    for (f in names(formals(ebal::ebalance.trim))) {
      if (is_null(A[[f]])) A[[f]] <- formals(ebal::ebalance.trim)[[f]]
    }
  }

  if (check.package("ebal")) {
    if (estimand == "ATT") {
      w <- rep(1, length(treat))
      control.levels <- levels(treat)[levels(treat) != focal]

      covs[treat == focal,] <- covs[treat == focal, , drop = FALSE] * s.weights[subset][treat == focal] * sum(treat == focal)/sum(s.weights[subset][treat == focal])

      for (i in control.levels) {
        treat.in.i.focal <- treat %in% c(focal, i)
        treat_ <- ifelse(treat[treat.in.i.focal] == i, 0, 1)
        covs_ <- covs[treat.in.i.focal, , drop = FALSE]

        colinear.covs.to.remove <- colnames(covs_)[colnames(covs_) %nin% colnames(remove.collinearity(covs_[treat_ == 0, , drop = FALSE]))]

        covs_ <- covs_[, colnames(covs_) %nin% colinear.covs.to.remove, drop = FALSE]

        ebal.out <- ebal::ebalance(Treatment = treat_, X = covs_,
                                   base.weight = A[["base.weight"]],
                                   norm.constant = A[["norm.constant"]],
                                   coefs = A[["coefs"]],
                                   max.iterations = A[["max.iterations"]],
                                   constraint.tolerance = A[["constraint.tolerance"]],
                                   print.level = 3)
        if (stabilize) ebal.out <- ebal::ebalance.trim(ebalanceobj = ebal.out,
                                                       max.weight = A[["max.weight"]],
                                                       min.weight = A[["min.weight"]],
                                                       max.trim.iterations = A[["max.trim.iterations"]],
                                                       max.weight.increment = A[["max.weight.increment"]],
                                                       min.weight.increment = A[["min.weight.increment"]],
                                                       print.level = 3)
        w[treat == i] <- ebal.out$w / s.weights[subset][treat.in.i.focal][treat_ == 0]

      }
    }
    else if (estimand == "ATE") {
      w <- rep(1, length(treat))

      for (i in levels(treat)) {
        covs_i <- rbind(covs, covs[treat==i, , drop = FALSE])
        treat_i <- c(rep(1, nrow(covs)), rep(0, sum(treat==i)))

        colinear.covs.to.remove <- colnames(covs_i)[colnames(covs_i) %nin% colnames(remove.collinearity(covs_i[treat_i == 0, , drop = FALSE]))]

        covs_i <- covs_i[, colnames(covs_i) %nin% colinear.covs.to.remove, drop = FALSE]

        covs_i[treat_i == 1,] <- covs_i[treat_i == 1,] * s.weights[subset] * sum(treat_i == 1) / sum(s.weights[subset])

        ebal.out_i <- ebal::ebalance(Treatment = treat_i, X = covs_i,
                                     base.weight = A[["base.weight"]],
                                     norm.constant = A[["norm.constant"]],
                                     coefs = A[["coefs"]],
                                     max.iterations = A[["max.iterations"]],
                                     constraint.tolerance = A[["constraint.tolerance"]],
                                     print.level = 3)

        if (stabilize) ebal.out_i <- ebal::ebalance.trim(ebalanceobj = ebal.out_i,
                                                         max.weight = A[["max.weight"]],
                                                         min.weight = A[["min.weight"]],
                                                         max.trim.iterations = A[["max.trim.iterations"]],
                                                         max.weight.increment = A[["max.weight.increment"]],
                                                         min.weight.increment = A[["min.weight.increment"]],
                                                         print.level = 3)

        w[treat == i] <- ebal.out_i$w / s.weights[subset][treat == i]
      }
    }
  }

  obj <- list(w = w)
  return(obj)
}

#Empirical Balancing Calibration weights with ATE
weightit2ebcw <- function(covs, treat, s.weights, subset, estimand, focal, moments, int, ...) {
  A <- list(...)

  covs <- covs[subset, , drop = FALSE]
  treat <- factor(treat)[subset]

  covs <- cbind(covs, int.poly.f(covs, poly = moments, int = int))
  covs <- apply(covs, 2, make.closer.to.1)

  if (any(vars.w.missing <- apply(covs, 2, function(x) any(is.na(x))))) {
    missing.ind <- apply(covs[, vars.w.missing, drop = FALSE], 2, function(x) as.numeric(is.na(x)))
    covs[is.na(covs)] <- 0
    covs <- cbind(covs, missing.ind)
  }

  for (f in names(formals(ATE::ATE))) {
    if (is_null(A[[f]])) A[[f]] <- formals(ATE::ATE)[[f]]
  }

  if (check.package("ATE")) {
    if (estimand == "ATT") {
      w <- rep(1, length(treat))
      control.levels <- levels(treat)[levels(treat) != focal]

      for (i in control.levels) {
        treat.in.i.focal <- treat %in% c(focal, i)
        treat_ <- ifelse(treat[treat.in.i.focal] == i, 0, 1)
        covs_ <- covs[treat.in.i.focal, , drop = FALSE]

        colinear.covs.to.remove <- colnames(covs_)[colnames(covs_) %nin% colnames(remove.collinearity(covs_[treat_ == 0, , drop = FALSE]))]

        covs_ <- covs_[, colnames(covs_) %nin% colinear.covs.to.remove, drop = FALSE]

        covs_[treat_ == 1,] <- covs_[treat_ == 1,] * s.weights[subset][treat == focal] * sum(treat == focal)/ sum(s.weights[subset][treat == focal])

        Y <- rep(0, length(treat_))

        ate.out <- ATE::ATE(Y = Y, Ti = treat_, X = covs_,
                            ATT = TRUE,
                            theta = A[["theta"]],
                            verbose = TRUE,
                            max.iter = A[["max.iter"]],
                            tol = A[["tol"]],
                            initial.values = A[["initial.values"]],
                            backtrack = A[["backtrack"]],
                            backtrack.alpha = A[["backtrack.alpha"]],
                            backtrack.beta = A[["backtrack.beta"]])
        w[treat == i] <- ate.out$weights.q[treat_ == 0] / s.weights[subset][treat == i]

      }
    }
    else if (estimand == "ATE") {
      w <- rep(1, length(treat))

      for (i in levels(treat)) {
        covs_i <- rbind(covs, covs[treat==i, , drop = FALSE])
        treat_i <- c(rep(1, nrow(covs)), rep(0, sum(treat==i)))

        colinear.covs.to.remove <- colnames(covs_i)[colnames(covs_i) %nin% colnames(remove.collinearity(covs_i[treat_i == 0, , drop = FALSE]))]

        covs_i <- covs_i[, colnames(covs_i) %nin% colinear.covs.to.remove, drop = FALSE]

        covs_i[treat_i == 1,] <- covs_i[treat_i == 1,] * s.weights[subset] * sum(treat_i == 1) / sum(s.weights[subset])

        Y <- rep(0, length(treat_i))

        ate.out <- ATE::ATE(Y = Y, Ti = treat_i, X = covs_i,
                            ATT = TRUE,
                            theta = A[["theta"]],
                            verbose = TRUE,
                            max.iter = A[["max.iter"]],
                            tol = A[["tol"]],
                            initial.values = A[["initial.values"]],
                            backtrack = A[["backtrack"]],
                            backtrack.alpha = A[["backtrack.alpha"]],
                            backtrack.beta = A[["backtrack.beta"]])
        w[treat == i] <- ate.out$weights.q[treat_i == 0] / s.weights[subset][treat == i]
      }
    }
    obj <- list(w = w)
    return(obj)
  }
}

#Stable balancing weights with sbw
# weightit2sbw <- function(...) {
#   stop("Stable balancing weights are not currently supported. Please choose another method.\n        The github version of WeightIt may allow stable balancing weights.\n        Install it with devtools::install_github(\"ngreifer/WeightIt\").", call. = FALSE)
# }

#SBW--------
weightit2sbw <- function(covs, treat, s.weights, subset, estimand, focal, moments, int, ...) {
  A <- list(...)

  if (check.package("sbw")) {
    check.package("slam")
    if (!"package:slam" %in% search()) {
      need.to.detach.slam <- TRUE
      attachNamespace("slam")
    } else {
      need.to.detach.slam <- FALSE
    }

    if (is_null(A$l_norm)) A$l_norm <- "l_2"
    if (is_null(A$solver)) A$solver <- "quadprog"
    if (is_null(A$max_iter)) A$max_iter <- 100000
    if (is_null(A$rel_tol)) A$rel_tol <- 1e-4
    if (is_null(A$abs_tol)) A$abs_tol <- 1e-4
    if (is_null(A$gap_stop)) A$gap_stop <- TRUE
    if (is_null(A$adaptive_rho)) A$adaptive_rho <- TRUE

    check.package(A$solver)
    if (!paste0("package:", A$solver) %in% search()) {
      need.to.detach.solver <- TRUE
      attachNamespace(A$solver)
    } else {
      need.to.detach.solver <- FALSE
    }

    covs <- covs[subset, , drop = FALSE]
    treat <- factor(treat)[subset]

    if (any(is.na(covs))) {
      stop("Stable balancing weights are not compatible with missing values.", call. = FALSE)
    }
    covs <- cbind(covs, int.poly.f(covs, poly = moments, int = int))
    covs <- apply(covs, 2, make.closer.to.1)

    #new.data <- setNames(data.frame(treat, covs), as.character(seq_len(1+ncol(covs))))

    binary.vars <- apply(covs, 2, function(x) !nunique.gt(x, 2))

    if (is_null(A$bal_tols)) bal_tols <- .0001
    else {
      bal_tols <- A$bal_tols
      if (length(bal_tols) != 1 && length(bal_tols) != ncol(covs)) {
        stop(paste0("bal_tols needs to be of length 1 or equal to the number of covariates (", ncol(covs),
                    ").\nThe covariates (in order) are:\n   ", paste0(colnames(covs), collapse = " ")), call.= FALSE)
      }
    }
    if (is_null(A$bal_tols_sd)) bal_tols_sd <- TRUE
    else bal_tols_sd <- A$bal_tols_sd

    if (is_null(A$bal_tols) && is_null(A$bal_tols_sd)) {
      message("Using bal_tols = 0.0001 and bal_tols_sd = TRUE.")
    }
    else if (is_null(A$bal_tols)) {
      message("Using bal_tols = 0.0001.")
    }
    else if (is_null(A$bal_tols_sd)) {
      message("Using bal_tols_sd = TRUE.")
    }

    if (estimand == "ATT") {
      control.levels <- levels(treat)[levels(treat) != focal]

      w <- rep(1, length(treat))

      if (bal_tols_sd) {
        bal_tols <- bal_tols * sapply(1:ncol(covs), function(x) {if (binary.vars[x]) 1 else sqrt(cov.wt(covs[treat == focal, x, drop = FALSE], s.weights[subset][treat == focal])$cov[1,1])})
      }

      covs[treat == focal,] <- covs[treat == focal,] * s.weights[subset][treat == focal] * sum(treat == focal)/sum(s.weights[subset][treat == focal])

      for (i in control.levels) {

        treat.in.i.focal <- treat %in% c(focal, i)
        treat_ <- ifelse(treat[treat.in.i.focal] == i, 0, 1)
        covs_ <- covs[treat.in.i.focal, , drop = FALSE]

        new.data_ <- data.frame(treat_, covs_)
        t_ind <- names(new.data_)[1]
        bal_covs = names(new.data_)[-1]

        sbw.fit <- sbw::sbw(new.data_,
                            t_ind = t_ind,
                            bal_covs = bal_covs,
                            bal_tols = bal_tols,
                            bal_tols_sd = FALSE,
                            target = "treated",
                            l_norm = A[["l_norm"]],
                            w_min = 0,
                            normalize = TRUE,
                            solver = A[["solver"]],
                            display = 1,
                            max_iter = A[["max_iter"]],
                            rel_tol = A[["rel_tol"]],
                            abs_tol = A[["abs_tol"]],
                            gap_stop = A[["gap_stop"]],
                            adaptive_rho = A[["adaptive_rho"]])

        w[treat==i] <- sbw.fit$data_frame_weights$weights[sbw.fit$data_frame_weights[[t_ind]] == 0]*sum(treat == i) / s.weights[subset][treat == i]
      }
      w[w < 0] <- 0
    }
    else if (estimand == "ATE") {
      if (bal_tols_sd) {
        bal_tols <- bal_tols * sapply(1:ncol(covs), function(x) {if (binary.vars[x]) 1 else sqrt(mean(sapply(unique(treat), function(t) cov.wt(covs[treat == t, x, drop = FALSE], s.weights[subset][treat == t])$cov[1,1])))})
      }

      bal_tols <- bal_tols/nunique(treat)

      w <- rep(1, length(treat))

      for (i in levels(treat)) {
        covs_i <- rbind(covs, covs[treat==i, , drop = FALSE])
        treat_i <- c(rep(1, nrow(covs)), rep(0, sum(treat==i)))

        covs_i[treat_i == 1,] <- covs_i[treat_i == 1, , drop = FALSE] * s.weights[subset] * sum(treat_i == 1) / sum(s.weights[subset])

        new.data_i <- data.frame(treat_i, covs_i)
        t_ind <- names(new.data_i)[1]
        bal_covs = names(new.data_i)[-1]

        sbw.fit_i <- sbw::sbw(new.data_i, t_ind = t_ind,
                              bal_covs = bal_covs,
                              bal_tols = bal_tols,
                              bal_tols_sd = FALSE,
                              target = "treated",
                              l_norm = A$l_norm,
                              w_min = 0,
                              normalize = TRUE,
                              solver = A[["solver"]],
                              display = 1,
                              max_iter = A[["max_iter"]],
                              rel_tol = A[["rel_tol"]],
                              abs_tol = A[["abs_tol"]],
                              gap_stop = A[["gap_stop"]],
                              adaptive_rho = A[["adaptive_rho"]])

        w[treat==i] <- sbw.fit_i$data_frame_weights$weights[sbw.fit_i$data_frame_weights[[t_ind]] == 0]*sum(treat == i) / s.weights[subset][treat == i]

      }
      w[w < 0] <- 0
    }
  }

  if (need.to.detach.slam) detach("package:slam", character.only = TRUE)
  if (need.to.detach.solver) detach(paste0("package:", A$solver), character.only = TRUE)

  obj <- list(w = w)
  return(obj)

}

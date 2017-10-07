#Propensity score estimation with regression
weightit2ps <- function(formula, data, s.weights, estimand, subset, stabilize, ps, ...) {
  #exact.factor should be a factor vector of values where each uniue value gets its own weights
  A <- list(...)
  if (length(A$link) == 0) A$link <- "logit"
  if (length(A$family) == 0) A$family <- binomial(link = A$link)

  if (length(ps) > 0) {
    mf <- model.frame(formula, data[subset,])
    t <- model.response(mf)
  }
  else {
    fit <- glm(formula, data = data[subset,],
               #weights = s.weights,
               family = A$family,
               ...)
    ps <- fit$fitted.values
    t <- fit$y
  }

  #Computing weights
  if (toupper(estimand) == "ATE") {
    w <- t/ps + (1-t)/(1-ps)
  }
  else if (toupper(estimand) == "ATT") {
    w <- t + (1-t)*ps/(1-ps)
  }
  else if (toupper(estimand) == "ATC") {
    w <- (1-t) + t*(1-ps)/ps
  }
  else if (toupper(estimand) == "ATO") {
    w <- t*(1-ps) + (1-t)*ps
  }
  if (stabilize) w <- w*(t*mean(t) + (1-t)*(1-mean(t)))

  obj <- list(ps = ps,
              t = t,
              w = w)
  return(obj)

}
weightit2ps.multi <- function(formula, data, s.weights, subset, estimand, focal, stabilize, ps, ...) {

}
weightit2ps.cont <- function(formula, data, s.weights, subset, stabilize, ps, ...) {

}

#Generalized boosted modeling with twang
weightit2gbm <- function(formula, data, s.weights, estimand, subset, stabilize, verbose, ...) {
  A <- list(...)

  tt <- terms(formula)
  mf <- model.frame(tt, data[subset,], drop.unused.levels = TRUE)
  treat <- model.response(mf)
  covs <- model.matrix(tt, data=mf)[,-1]

  if (estimand == "ATC") {
    treat <- 1 - treat
    estimand <- "ATT"
  }

  fit <- twang::ps(cobalt::f.build("treat", covs), data = data.frame(treat = treat, covs),
                   estimand = estimand, sampw = s.weights,
                   verbose = verbose, print.level = 2*verbose, ...)

  s <- names(fit$ps)[1]

  w <- cobalt::get.w(fit, stop.method = s)

  obj <- list(w = w,
              ps = fit$ps[,s])

  return(obj)
}
weightit2gbm.multi <- function(formula, data, s.weights, estimand, focal, subset, stabilize, verbose, ...) {
  A <- list(...)

  tt <- terms(formula)
  mf <- model.frame(tt, data[subset,], drop.unused.levels = TRUE)
  treat <- model.response(mf)
  covs <- model.matrix(tt, data=mf)[,-1]

  fit <- twang::mnps(data = data.frame(treat = treat, covs),
                     estimand = "ATE", sampw = s.weights,
                     verbose = verbose, print.level = 2*verbose,
                     treatATT = focal, ...)

  s <- names(fit$ps)[1]

  w <- cobalt::get.w(fit, stop.method = s)

  out <- list(w = w)
}

#CBPS
weightit2cbps <- function(formula, data, subset, estimand, verbose, ...) {
  A <- list(...)

  if (verbose) {
    fit <- CBPS::CBPS(formula, data = data[subset, ], ATT = switch(estimand, ATT = 1, ATC = 2, ATE = 0),
                                                method = ifelse(length(A$over) == 0 || isTRUE(A$over), "over", "exact"), ...)
  }
  else {
    capture.output(fit <- CBPS::CBPS(formula, data = data[subset, ], ATT = switch(estimand, ATT = 1, ATC = 2, ATE = 0),
                                                method = ifelse(length(A$over) == 0 || isTRUE(A$over), "over", "exact"), ...))
  }

  w <- cobalt::get.w(fit, estimand = estimand)

  obj <- list(w = w,
              ps = fit$fitted.values)

  return(obj)
}
weightit2cbps.multi <- function(formula, data, subset, verbose, ...) {
  A <- list(...)
  if (verbose) {
    fit <- CBPS::CBPS(formula, data = data[subset, ],
                      method = ifelse(length(A$over) == 0 || isTRUE(A$over), "over", "exact"), ...)
  }
  else {
    capture.output(fit <- CBPS::CBPS(formula, data = data[subset, ],
                                     method = ifelse(length(A$over) == 0 || isTRUE(A$over), "over", "exact"), ...))
  }

  w <- cobalt::get.w(fit)

  obj <- list(w = w)
}
weightit2cbps.cont <- weightit2cbps.multi
weightit2nbcbps <- function(formula, data, subset, verbose, ...) {
  A <- list(...)

  tt <- terms(formula)
  mf <- model.frame(tt, data[subset,], drop.unused.levels = TRUE)
  treat <- model.response(mf)
  covs <- model.matrix(tt, data=mf)[,-1]

  treat <- factor(treat)
  fit <- CBPS::nbCBPS(treat = treat, X = covs, print.level = verbose, ...)

  w <- cobalt::get.w(fit)

  obj <- list(w = w)

  return(obj)
}
weightit2nbcbps.multi <- function(formula, data, subset, estimand, verbose, ...) {
  A <- list(...)

  tt <- terms(formula)
  mf <- model.frame(tt, data[subset,], drop.unused.levels = TRUE)
  treat <- model.response(mf)
  covs <- model.matrix(tt, data=mf)[,-1]

  treat <- factor(treat)
  fit <- CBPS::nbCBPS(treat = treat, X = covs, print.level = verbose, ...)

  w <- cobalt::get.w(fit)

  obj <- list(w = w)

  return(obj)
}
weightit2nbcbps.cont <- weightit2nbcbps.multi

#Entropy balancing with ebal
weightit2ebal <- function(formula, data, s.weights, subset, estimand, stabilize, ...) {
  A <- list(...)

  tt <- terms(formula)
  mf <- model.frame(tt, data[subset,], drop.unused.levels = TRUE)
  treat <- model.response(mf)
  covs <- model.matrix(tt, data=mf)[,-1]

  covs <- covs * replicate(ncol(covs), s.weights)

  if (estimand %in% c("ATT", "ATC")) {
    if (estimand == "ATC" ) treat_ <- 1 - treat
    else treat_ <- treat

    ebal.out <- ebalance(Treatment = treat_, X = covs, ...)
    if (stabilize) ebal.out <- ebalance.trim(ebal.out, ...)

    w <- cobalt::get.w(ebal.out, treat = treat_)
  }
  else if (estimand == "ATE") {
    w <- rep(1, length(treat))

    #Reweight controls to be like total (need treated to look like total)
    covs1 <- rbind(covs, covs[treat==0,])
    treat1 <- c(rep(1, nrow(covs)), treat[treat==0])

    ebal.out1 <- ebalance(Treatment = treat1, X = covs1, ...)
    if (stabilize) ebal.out1 <- ebalance.trim(ebal.out1, ...)

    w[treat == 0] <- ebal.out1$w

    #Reweight treated to be like total

    covs0 <- rbind(covs, covs[treat==1,])
    treat0 <- c(rep(1, nrow(covs)), 1 - treat[treat==1])

    ebal.out0 <- ebalance(Treatment = treat0, X = covs0, ...)
    if (stabilize) ebal.out0 <- ebalance.trim(ebal.out0, ...)

    w[treat == 1] <- ebal.out0$w
  }

  obj <- list(w = w)
  return(obj)

}
weightit2ebal.multi <- function(formula, data, s.weights, subset, estimand, stabilize, ...) {
  A <- list(...)

  tt <- terms(formula)
  mf <- model.frame(tt, data[subset,], drop.unused.levels = TRUE)
  treat <- model.response(mf)
  covs <- model.matrix(tt, data=mf)[,-1]

  covs <- covs * replicate(ncol(covs), s.weights)

  if (estimand %in% c("ATT")) {
    w <- rep(1, length(treat))
    control.levels <- levels(treat)[levels(treat) != focal]

    for (i in control.levels) {
      treat_ <- ifelse(treat[treat %in% c(focal, i)] == i, 0, 1)
      covs_ <- covs[treat %in% c(focal, i),]
      ebal.out <- ebalance(Treatment = treat_, X = covs_, ...)
      if (stabilize) ebal.out <- ebalance.trim(ebal.out, ...)
      w[treat == i] <- e$w
    }
  }
  else if (estimand == "ATE") {
    w <- rep(1, length(treat))

    for (i in levels(treat)) {
      covs_i <- rbind(covs, covs[treat==i,])
      treat_i <- c(rep(1, nrow(covs)), rep(0, sum(treat==i)))

      ebal.out_i <- ebalance(Treatment = treat_i, X = covs_i, ...)
      if (stabilize) ebal.out_i <- ebalance.trim(ebal.out_i, ...)

      w[treat == i] <- ebal.out_i$w
    }
  }

  obj <- list(w = w)
  return(obj)
}

#Stable balancing weights with sbw
weightit2sbw <- function() {

}

#Empirical Balancing Calibration weights with ATE
weightit2ebcw <- function(formula, data, subset, estimand, ...) {
  A <- list(...)

  tt <- terms(formula)
  #attr(tt, "intercept") <- 0
  mf <- model.frame(tt, data[subset,], drop.unused.levels = TRUE)
  treat <- model.response(mf)
  covs <- model.matrix(tt, data=mf)[,-1]
  Y <- rep(0, length(treat))

  if (estimand %in% c("ATT", "ATC")) {
    if (estimand == "ATC" ) treat_ <- 1 - treat
    else treat_ <- treat

    ate.out <- ATE(Y = Y, Ti = treat_, X = covs,
                   ATT = TRUE, ...)
    w <- ate.out$weights.q
    w[treat_ == 1] <- 1

  }
  else if (estimand == "ATE") {
    ate.out <- ATE(Y = Y, Ti = treat, X = covs,
                   ATT = FALSE, ...)
    w <- ate.out$weights.q + ate.out$weights.p
  }

  obj <- list(w = w)
  return(obj)
}
weightit2ebcw.multi <- function(formula, data, subset, estimand, ...) {
  A <- list(...)

  tt <- terms(formula)
  mf <- model.frame(tt, data[subset,], drop.unused.levels = TRUE)
  treat <- as.numeric(model.response(mf)) - 1
  covs <- model.matrix(tt, data=mf)[,-1]
  Y <- rep(0, length(treat))

  ate.out <- ATE(Y = Y, Ti = treat, X = covs,
                 ATT = FALSE, ...)

  w <- apply(ate.out$weights.mat, 2, sum)

  obj <- list(w = w)
  return(obj)
}

#Propensity score estimation with ipw (probably won't use)
weightit2ipw <- function(formula, data, truncate.q, ...) {
  A <- list(...)

  tt <- terms(formula)
  #attr(tt, "intercept") <- 0
  mf <- model.frame(tt, data[subset,], drop.unused.levels = TRUE)
  treat <- model.response(mf)
  covs <- model.matrix(tt, data=mf)[,-1]

  out <- ipwpoint(exposure = treat, family = A$family, link = A$link, numerator = A$numerator,
                  denominator = f.build("", covs), data = data,
                  trunc = truncate.q, ...)
  if (length(out$weights.trunc) == 0) w <- obj$ipw.weights
  else w <- obj$weights.trunc
  obj <- list(w = w)
}
weightit2ipw.multi <- function(formula, data, truncate.q, ...) {

}
weightit2ipw.cont <- function(formula, data, truncate.q, ...) {

}

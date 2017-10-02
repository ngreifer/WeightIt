weightit2ps <- function(formula, data, s.weights, subset, ...) {
  #exact.factor should be a factor vector of values where each uniue value gets its own weights
  A <- list(...)
  if (length(A$family) == 0) A$family <- binomial(link = "logit")

  fit <- glm(formula, data = data[subset,],
             #weights = s.weights,
             family = A$family,
             ...)
  ps <- fit$fitted.values

  obj <- list(ps = ps,
              t = data[subset, all.vars(fit$terms)[1]])
  return(obj)

}
weightit2gbm <- function() {

}
weightit2cbps <- function(formula, data, subset, estimand, ...) {
print(estimand)
  A <- list(...)
  fit <- CBPS(formula, data = data[subset, ], ATT = switch(estimand, ATT = 1, ATC = 2, ATE = 0),
              method = ifelse(length(A$over) == 0 || isTRUE(A$over), "over", "exact"), ...)

  obj <- list(w = fit$weights,
              ps = fit$fitted.values)

  return(obj)
}
weightit2ebal <- function(formula, data, s.weights, subset, estimand, stabilize, ...) {
  A <- list(...)

  tt <- terms(formula)
  #attr(tt, "intercept") <- 0
  mf <- model.frame(tt, data[subset,], drop.unused.levels = TRUE)
  treat <- model.response(mf)
  covs <- model.matrix(tt, data=mf)[,-1]

  if (estimand %in% c("ATT", "ATC")) {
    if (estimand == "ATC" ) treat_ <- 1 - treat
    else treat_ <- treat

    ebal.out <- ebalance(Treatment = treat_, X = covs, ...)
    if (stabilize) ebal.out <- ebalance.trim(ebal.out, ...)

    w <- rep(1, length(treat_))
    if (length(ebal.out$w) != sum(treat_ == 0)) {
      stop("There are more control units in treat than weights in the ebalance object.", call. = FALSE)
    }
    w[treat_ == 0] <- ebal.out$w
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
    treat0 <- c(rep(1, nrow(covs)), 1- treat[treat==1])

    ebal.out0 <- ebalance(Treatment = treat0, X = covs0, ...)
    if (stabilize) ebal.out0 <- ebalance.trim(ebal.out0, ...)

    w[treat == 1] <- ebal.out0$w
  }

  obj <- list(w = w)
  return(obj)

}
weightit2sbw <- function() {

}
weightit2ecbw <- function(formula, data, subset, estimand, ...) {
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

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
    fit <- glm(formula, data = data.frame(data, .s.weights = s.weights)[subset,],
               weights = .s.weights,
               family = A$family,
               control = list(),
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
  else w <- NULL

  if (stabilize) {
    w <- w * sapply(t, function(x) sum(t==x) / sum(1*(t==x)*w))
  }

  obj <- list(ps = ps,
              w = w)
  return(obj)

}
weightit2ps.multi <- function(formula, data, s.weights, subset, estimand, focal, stabilize, ps, ...) {
  A <- list(...)
  if (length(A$link) == 0) A$link <- "logit"
  else {
    acceptable.links <- c("logit", "probit", "bayes.probit")
    which.link <- acceptable.links[pmatch(A$link, acceptable.links, nomatch = 0)][1]
    if (length(which.link) == 0) {
      A$link <- "logit"
      warning("Only  \"logit\",\"probit\" and \"bayes.probit\" are allowed as links for multinomial treatments. Using link = \"logit\".",
              call. = FALSE)
    }
    else A$link <- which.link
  }

  tt <- terms(formula)
  mf <- model.frame(tt, data[subset,], drop.unused.levels = TRUE)
  t <- model.response(mf)

  if (length(ps) > 0) {
  }
  else {
    if (A$link %in% c("logit", "probit")) {
      if (check.package("mlogit", alternative = TRUE) && (length(A$use.mlogit) == 0 ||
                                                          (length(A$use.mlogit) > 0 && A$use.mlogit != FALSE))) {
        message(paste0("Using multinomial ", A$link, " regression."))
        covs <- model.matrix(tt, data=mf)[,-1, drop = FALSE]
        mult <- mlogit::mlogit.data(data.frame(.t = t, covs, .s.weights = s.weights[subset]), varying = NULL, shape = "wide", sep = "", choice = ".t")
        tryCatch({fit <- mlogit::mlogit(as.formula(paste0(".t ~ 1 | ", paste(colnames(covs), collapse = " + "),
                                                " | 1")), data = mult, estimate = TRUE,
                              probit = ifelse(A$link[1] == "probit", TRUE, FALSE),
                              weights = .s.weights, ...)},
                 error = function(e) {stop(paste0("There was a problem fitting the multinomial ", A$link, " regressions with mlogit().\n       Try again with use.mlogit = FALSE."), call. = FALSE)}
        )
        ps <- fitted(fit, outcome = FALSE)
      }
      else {
        message(paste0("Using a series of ", nunique(t), " binomial ", A$link, " regressions."))
        covs <- apply(model.matrix(tt, data=mf)[,-1, drop = FALSE], 2, make.closer.to.1)
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
      fit <- MNP::mnp(formula, data[subset,], verbose = TRUE)
      ps <- MNP::predict.mnp(fit, type = "prob")$p
    }
    else {
      stop('link must be "logit", "probit", or "bayes.probit".', call. = FALSE)
    }
  }
  #ps should be matrix of probs for each treat
  #Computing weights
  if (toupper(estimand) == "ATE") {
    w <- rep(0, nrow(ps))
    for (i in seq_len(nunique(t))) {
      w[t == levels(t)[i]] <- 1/ps[t == levels(t)[i], i]
    }
  }
  else if (toupper(estimand) == "ATT") {
    w <- rep(0, nrow(ps))
    for (i in seq_len(nunique(t))) {
      w[t == levels(t)[i]] <- 1/ps[t == levels(t)[i], i]
    }
    w <- w*ps[, which(levels(t) == focal)]
  }
  else if (toupper(estimand) == "ATO") {
    w <- rep(0, nrow(ps))
    for (i in seq_len(nunique(t))) {
      w[t == levels(t)[i]] <- 1/ps[t == levels(t)[i], i]
    }
    w <- w*apply(ps, 1, prod)
  }
  else w <- NULL

  #ps <- sapply(seq_along(t), function(i) ps[i, which(levels(t) == t[i])])

  if (stabilize) {
    w <- w * sapply(t, function(x) sum(t==x) / sum(1*(t==x)*w))
  }

  obj <- list(w = w)
  return(obj)
}
weightit2ps.cont <- function(formula, data, s.weights, subset, stabilize, ps, ...) {
  A <- list(...)
  if (length(A$link) == 0) A$link <- "identity"
  if (length(A$family) == 0) A$family <- gaussian(link = A$link)

  mf <- model.frame(formula, data[subset,])
  t <- model.response(mf)

  stabilize <- TRUE

  if (length(ps) == 0) {
    fit <- glm(formula, data = data.frame(data, .s.weights = s.weights)[subset,],
               weights = .s.weights,
               family = A$family,
               control = list(),
               ...)
    p.denom <- fit$fitted.values
    den.denom <- dnorm(t, p.denom, sqrt(summary(fit)$dispersion))

    if (stabilize) {
      if (length(A$num.formula) == 0) A$num.formula <- ~ 1
      num.fit <- glm(update.formula(A$num.formula, .t ~ .),
                     data = data.frame(.t = t, data, .s.weights = s.weights)[subset,],
                     weights = .s.weights,
                     family = A$family,
                     control = list(), ...)
      p.num <- num.fit$fitted.values
      den.num <- dnorm(t, p.num, sqrt(summary(num.fit)$dispersion))
      w <- den.num/den.denom
    }
    else {
      w <- 1/den.denom
    }
  }

  obj <- list(ps = p.denom,
              w = w)
  return(obj)
}

#Generalized boosted modeling with twang
weightit2gbm <- function(formula, data, s.weights, estimand, subset, stabilize, ...) {
  A <- list(...)
  if (length(A$stop.method) == 0) {
    warning("No stop.method was provided. Using \"es.mean\".",
            call. = FALSE, immediate. = TRUE)
    A$stop.method <- "es.mean"
  }
  else if (length(A$stop.method) > 1) {
    warning("Only one stop.method is allowed at a time. Using just the first stop.method (\"", A$stop.method[1],"\").",
            call. = FALSE, immediate. = TRUE)
    A$stop.method <- A$stop.method[1]
  }

  tt <- terms(formula)
  mf <- model.frame(tt, data[subset,], drop.unused.levels = TRUE)
  treat <- model.response(mf)
  covs <- apply(model.matrix(tt, data=mf)[,-1, drop = FALSE], 2, make.closer.to.1)

  if (estimand == "ATC") {
    treat <- 1 - treat
    estimand <- "ATT"
  }

  new.data <- data.frame(treat = treat, covs)

  check.package("twang")
  fit <- do.call(twang::ps, c(list(formula = formula(new.data),
                                   data = new.data,
                                   estimand = estimand, sampw = s.weights[subset],
                                   verbose = TRUE, print.level = 2), A))

  s <- names(fit$ps)[1]
  w <- cobalt::get.w(fit, stop.method = s)

  if (stabilize) {
    w <- w * sapply(t, function(x) sum(t==x) / sum(1*(t==x)*w))
  }

  obj <- list(w = w,
              ps = fit$ps[[s]])

  return(obj)
}
weightit2gbm.multi <- function(formula, data, s.weights, estimand, focal, subset, stabilize, ...) {
  A <- list(...)
  if (length(A$stop.method) == 0) {
    warning("No stop.method was provided. Using \"es.mean\".",
            call. = FALSE, immediate. = TRUE)
    A$stop.method <- "es.mean"
  }
  else if (length(A$stop.method) > 1) {
    warning("Only one stop.method is allowed at a time. Using just the first stop.method.",
            call. = FALSE, immediate. = TRUE)
    A$stop.method <- A$stop.method[1]
  }

  tt <- terms(formula)
  mf <- model.frame(tt, data[subset,], drop.unused.levels = TRUE)
  treat <- model.response(mf)
  covs <- apply(model.matrix(tt, data=mf)[,-1, drop = FALSE], 2, make.closer.to.1)
  new.data <- data.frame(treat = treat, covs)

  check.package("twang")
  fit <- do.call(twang::mnps, c(list(formula = formula(new.data),
                                     data = new.data,
                                     estimand = estimand, sampw = s.weights[subset],
                                     verbose = TRUE, print.level = 2,
                                     treatATT = focal), A))

  s <- fit$stopMethods[1]

  w <- cobalt::get.w(fit, stop.method = s)[[1]]

  if (stabilize) {
    w <- w * sapply(t, function(x) sum(t==x) / sum(1*(t==x)*w))
  }

  out <- list(w = w)
}

#CBPS
weightit2cbps <- function(formula, data, subset, estimand, s.weights, stabilize, ...) {
  A <- list(...)

  tt <- terms(formula)
  mf <- model.frame(tt, data[subset,], drop.unused.levels = TRUE)
  treat <- model.response(mf)
  covs <- apply(model.matrix(tt, data=mf)[,-1, drop = FALSE], 2, make.closer.to.1)
  new.data <- data.frame(treat = treat, covs)

  check.package("CBPS")
  fit <- CBPS::CBPS(formula(new.data), data = new.data, ATT = switch(estimand, ATT = 1, ATC = 2, ATE = 0),
                    method = ifelse(length(A$over) == 0 || isTRUE(A$over), "over", "exact"),
                    standardize = FALSE,
                    sample.weights = s.weights[subset],
                    ...)

  w <- cobalt::get.w(fit, estimand = switch(estimand, ATE = "ate", "att"))

  if (stabilize) {
    w <- w * sapply(t, function(x) sum(t==x) / sum(1*(t==x)*w))
  }

  obj <- list(w = w,
              ps = fit$fitted.values)

  return(obj)
}
weightit2cbps.multi <- function(formula, data, subset, s.weights, stabilize, ...) {
  A <- list(...)

  tt <- terms(formula)
  mf <- model.frame(tt, data[subset,], drop.unused.levels = TRUE)
  treat <- model.response(mf)
  covs <- apply(model.matrix(tt, data=mf)[,-1, drop = FALSE], 2, make.closer.to.1)
  new.data <- data.frame(treat = treat, covs)

  check.package("CBPS")
  if (nunique(treat) <= 4) {
    fit <- CBPS::CBPS(formula(new.data),
                      data = new.data,
                      method = ifelse(length(A$over) == 0 || isTRUE(A$over), "over", "exact"),
                      standardize = FALSE,
                      sample.weights = s.weights[subset],
                      ...)


    w <- cobalt::get.w(fit)
  }
  else {
    w <- rep(1, length(treat))
    for (i in levels(treat)) {
      new.data[[1]] <- ifelse(treat == i, 1, 0)
      fit <- CBPS::CBPS(formula(new.data), data = new.data,
                        method = ifelse(length(A$over) == 0 || isTRUE(A$over), "over", "exact"),
                        standardize = FALSE,
                        sample.weights = s.weights[subset],
                        ATT = 0, ...)

      w[treat==i] <- cobalt::get.w(fit)[treat==i]
    }
  }

  if (stabilize) {
    w <- w * sapply(t, function(x) sum(t==x) / sum(1*(t==x)*w))
  }

  obj <- list(w = w)
}
weightit2cbps.cont <- function(formula, data, subset, s.weights, ...) {
  A <- list(...)

  tt <- terms(formula)
  mf <- model.frame(tt, data[subset,], drop.unused.levels = TRUE)
  treat <- model.response(mf)
  covs <- apply(model.matrix(tt, data=mf)[,-1, drop = FALSE], 2, make.closer.to.1)
  new.data <- data.frame(treat = treat, covs)

  check.package("CBPS")
  fit <- CBPS::CBPS(formula(new.data),
                    data = new.data,
                    method = ifelse(length(A$over) == 0 || isTRUE(A$over), "over", "exact"),
                    standardize = FALSE,
                    sample.weights = s.weights[subset],
                    ...)

  w <- cobalt::get.w(fit)

  obj <- list(w = w)
}
weightit2npcbps <- function(formula, data, subset, ...) {
  A <- list(...)

  tt <- terms(formula)
  mf <- model.frame(tt, data[subset,], drop.unused.levels = TRUE)
  treat <- factor(model.response(mf))
  covs <- apply(model.matrix(tt, data=mf)[,-1, drop = FALSE], 2, make.closer.to.1)
  new.data <- data.frame(treat = treat, covs)

  check.package("CBPS")
  fit <- CBPS::npCBPS(formula(new.data), data = new.data, print.level = 1, ...)

  w <- cobalt::get.w(fit)

  obj <- list(w = w)

  return(obj)
}
weightit2npcbps.multi <- function(formula, data, subset, estimand, ...) {
  A <- list(...)

  tt <- terms(formula)
  mf <- model.frame(tt, data[subset,], drop.unused.levels = TRUE)
  treat <- model.response(mf)
  covs <- apply(model.matrix(tt, data=mf)[,-1, drop = FALSE], 2, make.closer.to.1)
  new.data <- data.frame(.t = factor(treat), covs)

  check.package("CBPS")
  fit <- CBPS::npCBPS(formula(new.data), data = new.data, print.level = 1, ...)

  w <- cobalt::get.w(fit)

  obj <- list(w = w)

  return(obj)
}
weightit2npcbps.cont <- function(formula, data, subset, estimand, ...) {
  A <- list(...)

  tt <- terms(formula)
  mf <- model.frame(tt, data[subset,], drop.unused.levels = TRUE)
  treat <- model.response(mf)
  covs <- apply(model.matrix(tt, data=mf)[,-1, drop = FALSE], 2, make.closer.to.1)
  new.data <- data.frame(treat = treat, covs)

  check.package("CBPS")
  fit <- CBPS::npCBPS(formula(new.data), data = new.data, print.level = 1, ...)

  w <- cobalt::get.w(fit)

  obj <- list(w = w)

  return(obj)
}

#Entropy balancing with ebal
weightit2ebal <- function(formula, data, s.weights, subset, estimand, stabilize, ...) {
  A <- list(...)

  tt <- terms(formula)
  mf <- model.frame(tt, data[subset, , drop = FALSE], drop.unused.levels = TRUE)
  treat <- model.response(mf)
  covs <- apply(model.matrix(tt, data=mf)[, -1, drop = FALSE], 2, make.closer.to.1)
  #covs <- model.matrix(tt, data=mf)[, -1, drop = FALSE]

  #covs <- covs * replicate(ncol(covs), s.weights)

  check.package("ebal")
  if (estimand %in% c("ATT", "ATC")) {
    if (estimand == "ATC" ) treat_ <- 1 - treat
    else treat_ <- treat

    covs <- covs[, Reduce("intersect", lapply(unique(treat_, nmax = 2), function(j) colnames(remove.collinearity(covs[treat_ == j, , drop = FALSE])))), drop = FALSE]

    ebal.out <- ebal::ebalance(Treatment = treat_, X = covs,
                               print.level = 3,
                               base.weight = s.weights[subset][treat_ == 0], ...)
    if (stabilize) ebal.out <- ebal::ebalance.trim(ebal.out,
                                                   print.level = 3,
                                                   ...)

    w <- cobalt::get.w(ebal.out, treat = treat_)
  }
  else if (estimand == "ATE") {
    w <- rep(1, length(treat))

    for (i in unique(treat)) {
      #Reweight controls to be like total (need treated to look like total)

      covs_i <- rbind(covs, covs[treat==i, , drop = FALSE])
      treat_i <- c(rep(1, nrow(covs)), rep(0, sum(treat==i)))

      covs_i <- covs_i[, Reduce("intersect", lapply(unique(treat_i, nmax = 2), function(j) colnames(remove.collinearity(covs_i[treat_i == j, , drop = FALSE]))))]

      ebal.out_i <- ebal::ebalance(Treatment = treat_i, X = covs_i,
                                   base.weight = s.weights[subset][treat==i],
                                   print.level = 3, ...)
      if (stabilize) ebal.out_i <- ebal::ebalance.trim(ebal.out_i, ...)

      w[treat == i] <- ebal.out_i$w
    }
  }

  #w <- w/s.weights
  obj <- list(w = w)
  return(obj)

}
weightit2ebal.multi <- function(formula, data, s.weights, subset, estimand, focal, stabilize, ...) {
  A <- list(...)

  tt <- terms(formula)
  mf <- model.frame(tt, data[subset,], drop.unused.levels = TRUE)
  treat <- model.response(mf)
  covs <- apply(model.matrix(tt, data=mf)[,-1, drop = FALSE], 2, make.closer.to.1)

  check.package("ebal")
  if (estimand %in% c("ATT")) {
    w <- rep(1, length(treat))
    control.levels <- levels(treat)[levels(treat) != focal]

    for (i in control.levels) {
      treat_ <- ifelse(treat[treat %in% c(focal, i)] == i, 0, 1)
      covs_ <- covs[treat %in% c(focal, i), , drop = FALSE]

      covs_ <- covs_[, Reduce("intersect", lapply(unique(treat_, nmax = 2), function(j) colnames(remove.collinearity(covs_[treat_ == j, , drop = FALSE]))))]

      ebal.out <- ebal::ebalance(Treatment = treat_, X = covs_,
                                 base.weight = s.weights[subset][treat == i],
                                 print.level = 3, ...)
      if (stabilize) ebal.out <- ebal::ebalance.trim(ebal.out,
                                                     print.level = 3, ...)
      w[treat == i] <- ebal.out$w
    }
  }
  else if (estimand == "ATE") {
    w <- rep(1, length(treat))

    for (i in levels(treat)) {
      covs_i <- rbind(covs, covs[treat==i, , drop = FALSE])
      treat_i <- c(rep(1, nrow(covs)), rep(0, sum(treat==i)))

      covs_i <- covs_i[, Reduce("intersect", lapply(unique(treat_i, nmax = 2), function(j) colnames(remove.collinearity(covs_i[treat_i == j, , drop = FALSE]))))]

      ebal.out_i <- ebal::ebalance(Treatment = treat_i, X = covs_i,
                                   base.weight = s.weights[subset][treat == i],
                                   print.level = 3, ...)
      if (stabilize) ebal.out_i <- ebal::ebalance.trim(ebal.out_i,
                                                       print.level = 3,
                                                       ...)

      w[treat == i] <- ebal.out_i$w
    }
  }

  obj <- list(w = w)
  return(obj)
}

#Stable balancing weights with sbw
weightit2sbw <- function(formula, data, s.weights, subset, estimand, ...) {
  A <- list(...)

  check.package("sbw")
  check.package("slam"); requireNamespace("slam")

  if (length(A$l_norm) == 0) A$l_norm <- "l_2"
  if (length(A$solver) == 0) A$solver <- "quadprog"
  if (length(A$max_iter) == 0) A$max_iter <- 100000
  if (length(A$rel_tol) == 0) A$rel_tol <- 1e-4
  if (length(A$abs_tol) == 0) A$abs_tol <- 1e-4
  if (length(A$gap_stop) == 0) A$gap_stop <- TRUE
  if (length(A$adaptive_rho) == 0) A$adaptive_rho <- TRUE

  check.package(A$solver); require(A$solver, character.only = TRUE)

  tt <- terms(formula)
  mf <- model.frame(tt, data[subset,], drop.unused.levels = TRUE)
  treat <- model.response(mf)
  covs0 <- model.matrix(tt, data=mf,
                        contrasts.arg = lapply(mf[sapply(mf, is.factor)],
                                               contrasts, contrasts=FALSE))[,-1, drop = FALSE]
  covs <- apply(covs0, 2, make.closer.to.1)
  new.data <- setNames(data.frame(treat, covs), as.character(seq_len(1+ncol(covs))))

  if (length(A$bal_tols) == 0) bal_tols <- rep(.0001, ncol(covs))
  else {
    bal_tols <- A$bal_tols
    if (length(bal_tols) != 1 && length(bal_tols) != ncol(covs)) {
      stop(paste0("bal_tols needs to be of length 1 or equal to the number of covariates (", ncol(covs),
                  ").\nThe covariates (in order) are:\n   ", paste0(colnames(covs), collapse = " ")), call.= FALSE)
    }
  }
  if (length(A$bal_tols_sd) == 0) bal_tols_sd <- TRUE
  else bal_tols_sd <- A$bal_tols_sd

  if (length(A$bal_tols) == 0 && length(A$bal_tols_sd) == 0) {
    message("Using bal_tols = 0.0001 and bal_tols_sd = TRUE.")
  }
  else if (length(A$bal_tols) == 0) {
    message("Using bal_tols = 0.0001.")
  }
  else if (length(A$bal_tols_sd) == 0) {
    message("Using bal_tols_sd = TRUE.")
  }

  if (estimand %in% c("ATT", "ATC")) {
    new.data_ <- new.data
    if (estimand == "ATC") {
      new.data_[["1"]] <- 1 - new.data_[["1"]]
    }
    if (bal_tols_sd) {
      binary.vars <- apply(covs, 2, function(x) nunique(x) <= 2)
      bal_tols <- bal_tols * sapply(1:ncol(covs), function(x) {if (binary.vars[x]) 1 else sqrt(mean((covs[new.data_[["1"]] == 1,x] - mean(covs[new.data_[["1"]] == 1,x]))^2))})
    }
    sbw.fit <- sbw::sbw(new.data_, "1", names(new.data_)[-1],
                        bal_tols = bal_tols,
                        bal_tols_sd = FALSE,
                        target = "treated",
                        l_norm = A$l_norm,
                        w_min = 0,
                        normalize = TRUE,
                        solver = A$solver,
                        display = 1,
                        max_iter = A$max_iter,
                        rel_tol = A$rel_tol,
                        abs_tol = A$abs_tol,
                        gap_stop = A$gap_stop,
                        adaptive_rho = A$adaptive_rho)
    w <- rep(0, length(treat))
    w[new.data_[["1"]]==0] <- sbw.fit$data_frame_weights$weights[sbw.fit$data_frame_weights[["1"]] == 0]*sum(new.data_[["1"]] == 0)
    w[new.data_[["1"]]==1] <- 1
    w[w < 0] <- 0
  }
  else {
    if (bal_tols_sd) {
      binary.vars <- apply(covs, 2, function(x) nunique(x) <= 2)
      bal_tols <- bal_tols * sapply(1:ncol(covs), function(x) {if (binary.vars[x]) 1 else sqrt(mean((covs[,x] - mean(covs[,x]))^2))})
    }
    sbw.fit0 <- sbw::sbw(new.data, "1", names(new.data)[-1],
                         bal_tols = bal_tols,
                         bal_tols_sd = FALSE,
                         target = "all",
                         l_norm = A$l_norm,
                         w_min = 0,
                         normalize = TRUE,
                         solver = A$solver, display = 1,
                         max_iter = A$max_iter,
                         rel_tol = A$rel_tol,
                         abs_tol = A$abs_tol,
                         gap_stop = A$gap_stop,
                         adaptive_rho = A$adaptive_rho)
    new.data_ <- new.data; new.data_[["1"]] <- 1 - new.data_[["1"]]
    sbw.fit1 <- sbw::sbw(new.data_, "1", names(new.data)[-1],
                         bal_tols = bal_tols,
                         bal_tols_sd = FALSE,
                         target = "all",
                         l_norm = A$l_norm,
                         w_min = 0,
                         normalize = TRUE,
                         solver = A$solver, display = 1,
                         max_iter = A$max_iter,
                         rel_tol = A$rel_tol,
                         abs_tol = A$abs_tol,
                         gap_stop = A$gap_stop,
                         adaptive_rho = A$adaptive_rho)
    w <- rep(0, length(treat))
    w[treat==0] <- sbw.fit0$data_frame_weights$weights[sbw.fit0$data_frame_weights[["1"]]==0]
    w[treat==1] <- sbw.fit1$data_frame_weights$weights[sbw.fit1$data_frame_weights[["1"]]==0]
    w[w < 0] <- 0
  }

  obj <- list(w = w)
  return(obj)

}
weightit2sbw.multi <- function(formula, data, s.weights, subset, estimand, focal, ...) {
  A <- list(...)

  check.package("sbw")
  check.package("slam"); requireNamespace("slam")

  if (length(A$l_norm) == 0) A$l_norm <- "l_2"
  if (length(A$solver) == 0) A$solver <- "quadprog"
  if (length(A$max_iter) == 0) A$max_iter <- 100000
  if (length(A$rel_tol) == 0) A$rel_tol <- 1e-4
  if (length(A$abs_tol) == 0) A$abs_tol <- 1e-4
  if (length(A$gap_stop) == 0) A$gap_stop <- TRUE
  if (length(A$adaptive_rho) == 0) A$adaptive_rho <- TRUE

  tt <- terms(formula)
  mf <- model.frame(tt, data[subset,], drop.unused.levels = TRUE)
  treat <- model.response(mf)
  covs0 <- model.matrix(tt, data=mf,
                        contrasts.arg = lapply(mf[sapply(mf, is.factor)],
                                               contrasts, contrasts=FALSE))[,-1, drop = FALSE]
  covs <- apply(covs0, 2, make.closer.to.1)
  new.data <- setNames(data.frame(treat, covs), as.character(seq_len(1+ncol(covs))))

  if (length(A$bal_tols) == 0) bal_tols <- .0001
  else {
    bal_tols <- A$bal_tols
    if (length(bal_tols) != 1 && length(bal_tols) != ncol(covs)) {
      stop(paste0("bal_tols needs to be of length 1 or equal to the number of covariates (", ncol(covs),
                  ").\nThe covariates (in order) are:\n   ", paste0(colnames(covs), collapse = " ")), call.= FALSE)
    }
  }
  if (length(A$bal_tols_sd) == 0) bal_tols_sd <- TRUE
  else bal_tols_sd <- A$bal_tols_sd

  if (length(A$bal_tols) == 0 && length(A$bal_tols_sd) == 0) {
    message("Using bal_tols = 0.0001 and bal_tols_sd = TRUE.")
  }
  else if (length(A$bal_tols) == 0) {
    message("Using bal_tols = 0.0001.")
  }
  else if (length(A$bal_tols_sd) == 0) {
    message("Using bal_tols_sd = TRUE.")
  }

  if (estimand == "ATT") {
    control.levels <- levels(treat)[levels(treat) != focal]
    w <- rep(1, length(treat))
    for (i in control.levels) {
      if (bal_tols_sd) {
        binary.vars <- apply(covs, 2, function(x) nunique(x) <= 2)
        bal_tols <- bal_tols * sapply(1:ncol(covs), function(x) {if (binary.vars[x]) 1 else sqrt(mean((covs[treat == i, x] - mean(covs[treat == i, x]))^2))})
      }
      treat_ <- ifelse(treat[treat %in% c(focal, i)] == i, 0, 1)
      covs_ <- covs[treat %in% c(focal, i),]
      new.data_ <- setNames(data.frame(treat_, covs_), as.character(seq_len(1+ncol(covs))))

      sbw.fit <- sbw::sbw(new.data_, "1", names(new.data_)[-1],
                          bal_tols = bal_tols,
                          bal_tols_sd = bal_tols_sd,
                          target = "treated",
                          l_norm = A$l_norm,
                          w_min = 0,
                          normalize = TRUE,
                          solver = A$solver, display = 1,
                          max_iter = A$max_iter,
                          rel_tol = A$rel_tol,
                          abs_tol = A$abs_tol,
                          gap_stop = A$gap_stop,
                          adaptive_rho = A$adaptive_rho)

      w[treat==i] <- sbw.fit$data_frame_weights$weights[sbw.fit$data_frame_weights[["1"]] == 0]*sum(treat == i)
    }
    w[w < 0] <- 0
  }
  else {
    w <- rep(1, length(treat))
    for (i in levels(treat)) {
      new.data[["1"]] <- ifelse(treat == i, 0, 1)

      sbw.fit_i <- sbw::sbw(new.data, "1", names(new.data)[-1],
                            bal_tols = bal_tols,
                            bal_tols_sd = bal_tols_sd,
                            target = "all",
                            l_norm = A$l_norm,
                            w_min = 0,
                            normalize = TRUE,
                            solver = A$solver,
                            display = 1,
                            max_iter = A$max_iter,
                            rel_tol = A$rel_tol,
                            abs_tol = A$abs_tol,
                            gap_stop = A$gap_stop,
                            adaptive_rho = A$adaptive_rho)
      w[treat==i] <- sbw.fit_i$data_frame_weights$weights[sbw.fit_i$data_frame_weights[["1"]]==0]*sum(treat == i)
    }
    w[w < 0] <- 0
  }

  obj <- list(w = w)
  return(obj)

}

#Empirical Balancing Calibration weights with ATE
weightit2ebcw <- function(formula, data, subset, estimand, ...) {
  A <- list(...)

  tt <- terms(formula)
  mf <- model.frame(tt, data[subset,], drop.unused.levels = TRUE)
  treat <- model.response(mf)
  covs <- apply(model.matrix(tt, data=mf)[,-1, drop = FALSE], 2, make.closer.to.1)
  Y <- rep(0, length(treat))

  check.package("ATE")
  if (estimand %in% c("ATT", "ATC")) {
    if (estimand == "ATC" ) treat_ <- 1 - treat
    else treat_ <- treat

    ate.out <- ATE::ATE(Y = Y, Ti = treat_, X = covs,
                        ATT = TRUE, ...)
    w <- ate.out$weights.q
    w[treat_ == 1] <- 1

  }
  else if (estimand == "ATE") {
    ate.out <- ATE::ATE(Y = Y, Ti = treat, X = covs,
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
  covs <- apply(model.matrix(tt, data=mf)[,-1, drop = FALSE], 2, make.closer.to.1)
  Y <- rep(0, length(treat))

  check.package("ATE")
  ate.out <- ATE::ATE(Y = Y, Ti = treat, X = covs,
                      ATT = FALSE, ...)

  w <- apply(ate.out$weights.mat, 2, sum)

  obj <- list(w = w)
  return(obj)
}

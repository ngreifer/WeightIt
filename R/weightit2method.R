#User-defined weighting function
weightit2user <- function(Fun, covs, treat, s.weights, subset, estimand, focal, stabilize, ps, moments, int, ...) {
  A <- list(...)
  if (is_not_null(covs)) {
    covs <- covs[subset, , drop = FALSE]
  }
  if (is_not_null(treat)) {
    treat <- treat[subset]
  }
  if (is_not_null(s.weights)) {
    s.weights <- s.weights[subset]
  }
  if (is_not_null(ps)) {
    ps <- ps[subset]
  }

  #Get a list of function args for the user-defined function Fun
  Fun_formal <- as.list(formals(Fun))
  if (has_dots <- any(names(Fun_formal) == "...")) {
    Fun_formal[["..."]] <- NULL
  }

  fun_args <- Fun_formal
  for (i in names(fun_args)) {
    if (is_not_null(get0(i))) fun_args[[i]] <- get0(i)
    else if (is_not_null(A[[i]])) {
      fun_args[[i]] <- A[[i]]
      A[[i]] <- NULL
    }
    #else just use Fun default
  }

  if (has_dots) fun_args <- c(fun_args, A)
  else {
    if (is_not_null(A)) {
      Anames <- names(A)
      unnamedAnames <- Anames[Anames == ""]
      namedAnames <- Anames[Anames != ""]
      if (length(unnamedAnames) == 1) Anames <- c(namedAnames, "an unnamed argument")
      else if (length(unnamedAnames) > 1) Anames <- c(namedAnames, paste(length(unnamedAnames), "unnamed arguments"))

      if (length(Anames) > 1) warning(paste0("The following arguments were specified but are not suitable arguments to the provided function:\n\t", word.list(Anames)), call. = FALSE)
      else warning(paste0("The following argument was specified but is not a suitable argument to the provided function:\n\t", Anames), call. = FALSE)
    }
  }

  obj <- do.call(Fun, fun_args)

  if (is.numeric(obj)) {
    obj <- list(w = obj)
  }
  else if (!is.list(obj) || !any(c("w", "weights") %nin% names(obj))) {
    stop("The output of the user-provided function must be a list with an entry named \"w\" or \"weights\" containing the estimated weights.", call. = FALSE)
  }
  else {
    names(obj)[names(obj) == "weights"] <- "w"
  }
  if (is_null(obj[["w"]])) stop("No weights were estimated.", call. = FALSE)
  if (!is.vector(obj[["w"]], mode = "numeric")) stop("The \"w\" or \"weights\" entry of the output of the user-provided function must be a numeric vector of weights.", call. = FALSE)
  if (all(is.na(obj[["w"]]))) stop("All weights were generated as NA.", call = FALSE)
  if (length(obj[["w"]]) != length(treat)) {
    stop(paste(length(obj[["w"]]), "weights were estimated, but there are", length(treat), "units."), call. = FALSE)
  }

  return(obj)
}
weightitMSM2user <- function(Fun, covs.list, treat.list, s.weights, subset, stabilize, moments, int, ...) {
  A <- list(...)
  if (is_not_null(covs.list)) {
    covs.list <- covs <- lapply(covs.list, function(c) c[subset, , drop = FALSE])
  }
  if (is_not_null(treat.list)) {
    treat.list <- treat <- lapply(treat.list, function(t) t[subset])
  }
  if (is_not_null(s.weights)) {
    s.weights <- s.weights[subset]
  }

  #Get a list of function args for the user-defined function Fun
  Fun_formal <- as.list(formals(Fun))
  if (has_dots <- any(names(Fun_formal) == "...")) {
    Fun_formal[["..."]] <- NULL
  }

  fun_args <- Fun_formal
  for (i in names(fun_args)) {
    if (is_not_null(get0(i))) fun_args[[i]] <- get0(i)
    else if (is_not_null(A[[i]])) {
      fun_args[[i]] <- A[[i]]
      A[[i]] <- NULL
    }
    #else just use Fun default
  }

  if (has_dots) fun_args <- c(fun_args, A)
  else {
    if (is_not_null(A)) {
      Anames <- names(A)
      unnamedAnames <- Anames[Anames == ""]
      namedAnames <- Anames[Anames != ""]
      if (length(unnamedAnames) == 1) Anames <- c(namedAnames, "an unnamed argument")
      else if (length(unnamedAnames) > 1) Anames <- c(namedAnames, paste(length(unnamedAnames), "unnamed arguments"))

      if (length(Anames) > 1) warning(paste0("The following arguments were specified but are not suitable arguments to the provided function:\n\t", word.list(Anames)), call. = FALSE)
      else warning(paste0("The following argument was specified but is not a suitable argument to the provided function:\n\t", Anames), call. = FALSE)
    }
  }

  obj <- do.call(Fun, fun_args)

  if (is.numeric(obj)) {
    obj <- list(w = obj)
  }
  else if (!is.list(obj) || !any(c("w", "weights") %nin% names(obj))) {
    stop("The output of the user-provided function must be a list with an entry named \"w\" or \"weights\" containing the estimated weights.", call. = FALSE)
  }
  else {
    names(obj)[names(obj) == "weights"] <- "w"
  }
  if (is_null(obj[["w"]])) stop("No weights were estimated.", call. = FALSE)
  if (!is.vector(obj[["w"]], mode = "numeric")) stop("The \"w\" or \"weights\" entry of the output of the user-provided function must be a numeric vector of weights.", call. = FALSE)
  if (all(is.na(obj[["w"]]))) stop("All weights were generated as NA.", call = FALSE)
  if (length(obj[["w"]]) != length(treat.list[[1]])) {
    stop(paste(length(obj[["w"]]), "weights were estimated, but there are", length(treat.list[[1]]), "units."), call. = FALSE)
  }

  return(obj)
}

#Propensity score estimation with regression
weightit2ps <- function(covs, treat, s.weights, subset, estimand, focal, stabilize, ps, ...) {
  A <- list(...)

  fit.obj <- NULL

  if (is_null(ps)) {

    covs <- covs[subset, , drop = FALSE]
    t <- factor(treat[subset])

    if (any(vars.w.missing <- apply(covs, 2, function(x) any(is.na(x))))) {
      missing.ind <- apply(covs[, vars.w.missing, drop = FALSE], 2, function(x) as.numeric(is.na(x)))
      covs[is.na(covs)] <- 0
      covs <- cbind(covs, missing.ind)
    }
    covs <- apply(covs, 2, make.closer.to.1)

    if (ncol(covs) > 1) {
      colinear.covs.to.remove <- colnames(covs)[colnames(covs) %nin% colnames(make.full.rank(covs))]
      covs <- covs[, colnames(covs) %nin% colinear.covs.to.remove, drop = FALSE]
    }

    if (is_null(A$link)) A$link <- "logit"
    else {
      if (bin.treat <- is_binary(t)) acceptable.links <- c("logit", "probit", "cloglog", "identity", "log", "cauchit")
      else acceptable.links <- c("logit", "probit", "bayes.probit")
      which.link <- acceptable.links[pmatch(A$link, acceptable.links, nomatch = 0)][1]
      if (is.na(which.link)) {
        A$link <- acceptable.links[1]
        warning(paste0("Only ", word.list(acceptable.links, quotes = TRUE), " are allowed as links for ", if (bin.treat) "binary" else "multinomial", " treatments. Using link = ", word.list(acceptable.links[1], quotes = TRUE), "."),
                call. = FALSE)
      }
      else A$link <- which.link
    }

    if (is_binary(t)) {
      data <- data.frame(t, covs)
      formula <- formula(data)

      ps <- setNames(as.data.frame(matrix(NA_real_, ncol = nunique(t), nrow = length(t))),
                     levels(t))

      fit <- do.call("glm", list(formula, data = data,
                                 weights = s.weights[subset],
                                 family = binomial(link = A$link),
                                 control = as.list(A$control)), quote = TRUE)
      ps[[2]] <- p.score <- fit$fitted.values
      ps[[1]] <- 1 - ps[[2]]

      fit.obj <- fit
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
          fit.obj <- fit
        }
        else {
          message(paste0("Using a series of ", nunique(t), " binomial ", A$link, " regressions."))
          ps <- setNames(as.data.frame(matrix(NA_real_, ncol = nunique(t), nrow = length(t))),
                         levels(t))

          fit.list <- setNames(vector("list", nlevels(t)), levels(t))
          for (i in levels(t)) {
            t_i <- rep(0, length(t)); t_i[t == i] <- 1
            data_i <- data.frame(t_i, covs)
            fit.list[[i]] <- glm(formula(data_i), data = data_i,
                                 family = binomial(link = A$link),
                                 weights = s.weights[subset])
            ps[[i]] <- fit.list[[i]]$fitted.values
          }
          fit.obj <- fit.list
        }
      }
      else if (A$link == "bayes.probit") {
        check.package("MNP")
        data <- data.frame(t, covs)
        formula <- formula(data)
        tryCatch({fit <- MNP::mnp(formula, data, verbose = TRUE)},
                 error = function(e) stop("There was a problem with the Bayes probit regression. Try a different link.", call. = FALSE))
        ps <- MNP::predict.mnp(fit, type = "prob")$p
        fit.obj <- fit
      }
      else {
        stop('link must be "logit", "probit", or "bayes.probit".', call. = FALSE)
      }
      p.score <- NULL
    }
  }
  else {
      p.score <- ps
      ps <- setNames(as.data.frame(matrix(c(1-ps, ps), ncol = 2)),
                     levels(t))

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
    w <- w/apply(ps, 1, function(x) sum(1/x)) #Li & Li (2018)
  }
  else if (toupper(estimand) == "ATM") {
    w <- w*apply(ps, 1, min)
  }
  else w <- NULL

  if (stabilize) {
    tab <- vapply(levels(t), function(x) mean(t == x), numeric(1L))
    w <- w * tab[t]
  }

  obj <- list(w = w, ps = p.score, fit.obj = fit.obj)
  return(obj)
}
weightit2ps.cont <- function(covs, treat, s.weights, subset, stabilize, ps, ...) {
  A <- list(...)

  fit.obj <- NULL

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
    if (is_null(A$link)) A$link <- "identity"
    fit <- do.call("glm", c(list(formula, data = data,
                                 weights = s.weights[subset],
                                 family = gaussian(link = A$link),
                                 control = as.list(A$control))),
                   quote = TRUE)
    p.denom <- t - fit$fitted.values

    fit.obj <- fit

    if (isTRUE(A[["use.kernel"]])) {
      if (is_null(A[["bw"]])) A[["bw"]] <- "nrd0"
      if (is_null(A[["adjust"]])) A[["adjust"]] <- 1
      if (is_null(A[["kernel"]])) A[["kernel"]] <- "gaussian"
      if (is_null(A[["n"]])) A[["n"]] <- 10*length(t)

      d.d <- density(p.denom, n = A[["n"]],
                     weights = s.weights[subset]/sum(s.weights[subset]), give.Rkern = FALSE,
                     bw = A[["bw"]], adjust = A[["adjust"]], kernel = A[["kernel"]])
      if (isTRUE(A[["plot"]])) {
        par(mfrow=c(2,1))
        plot(d.d, main = "Denominator density")
      }
      dens.denom <- with(d.d, approxfun(x = x, y = y))(p.denom)
    }
    else {
      dens.denom <- dnorm(p.denom, 0, sd = sqrt(summary(fit)$dispersion))
    }

    if (stabilize) {

      num.fit <- do.call("glm", c(list(t ~ 1, data = data.frame(t = t),
                                       weights = s.weights[subset],
                                       family = gaussian(link = A$link),
                                       control = as.list(A$control))),
                         quote = TRUE)

      p.num <- t - num.fit$fitted.values

      if (isTRUE(A[["use.kernel"]])) {
        d.n <- density(p.num, n = A[["n"]],
                       weights = s.weights[subset]/sum(s.weights[subset]), give.Rkern = FALSE,
                       bw = A[["bw"]], adjust = A[["adjust"]], kernel = A[["kernel"]])
        if (isTRUE(A[["plot"]])) plot(d.n, main = "Numerator density")
        dens.num <- with(d.n, approxfun(x = x, y = y))(p.num)
      }
      else {
        dens.num <- dnorm(p.num, 0, sqrt(summary(num.fit)$dispersion))
      }
      w <- dens.num/dens.denom
    }
    else {
      w <- 1/dens.denom
    }
  }

  obj <- list(w = w, fit.obj = fit.obj)
  return(obj)
}

#MABW with optweight
weightit2optweight <- function(covs, treat, s.weights, subset, estimand, focal, moments, int, ...) {
  A <- list(...)

  covs <- covs[subset, , drop = FALSE]
  treat <- factor(treat[subset])

  covs <- cbind(covs, int.poly.f(covs, poly = moments, int = int))
  covs <- apply(covs, 2, make.closer.to.1)

  if (any(vars.w.missing <- apply(covs, 2, function(x) any(is.na(x))))) {
    missing.ind <- apply(covs[, vars.w.missing, drop = FALSE], 2, function(x) as.numeric(is.na(x)))
    covs[is.na(covs)] <- 0
    covs <- cbind(covs, missing.ind)
  }

  new.data <- data.frame(treat, covs)
  new.formula <- formula(new.data)

  for (f in names(formals(optweight::optweight))) {
    if (is_null(A[[f]])) A[[f]] <- formals(optweight::optweight)[[f]]
  }
  A[names(A) %in% names(formals(weightit2optweight))] <- NULL

  if ("tols" %in% names(A)) A[["tols"]] <- optweight::check.tols(new.formula, new.data, A[["tols"]], stop = TRUE)
  if ("targets" %in% names(A)) {
    warning("targets cannot be used through WeightIt and will be ignored.", call. = FALSE)
    A[["targets"]] <- NULL
  }

  A[["formula"]] <- new.formula
  A[["data"]] <- new.data
  A[["estimand"]] <- estimand
  A[["s.weights"]] <- s.weights[subset]
  A[["focal"]] <- focal
  A[["verbose"]] <- TRUE

  if (check.package("optweight")) {
    out <- do.call(optweight::optweight, A, quote = TRUE)
    obj <- list(w = out[["weights"]], fit.obj = out)
    return(obj)
  }
}
weightit2optweight.cont <- function(covs, treat, s.weights, subset, moments, int, ...) {
  A <- list(...)
  check.package("optweight")

  covs <- covs[subset, , drop = FALSE]
  treat <- treat[subset]

  covs <- cbind(covs, int.poly.f(covs, poly = moments, int = int))
  covs <- apply(covs, 2, make.closer.to.1)

  if (any(vars.w.missing <- apply(covs, 2, function(x) any(is.na(x))))) {
    missing.ind <- apply(covs[, vars.w.missing, drop = FALSE], 2, function(x) as.numeric(is.na(x)))
    covs[is.na(covs)] <- 0
    covs <- cbind(covs, missing.ind)
  }

  new.data <- data.frame(treat, covs)
  new.formula <- formula(new.data)

  for (f in names(formals(optweight::optweight))) {
    if (is_null(A[[f]])) A[[f]] <- formals(optweight::optweight)[[f]]
  }
  A[names(A) %in% names(formals(weightit2optweight.cont))] <- NULL

  if ("tols" %in% names(A)) A[["tols"]] <- optweight::check.tols(new.formula, new.data, A[["tols"]], stop = TRUE)
  if ("targets" %in% names(A)) {
    warning("targets cannot be used through WeightIt and will be ignored.", call. = FALSE)
    A[["targets"]] <- NULL
  }

  A[["formula"]] <- new.formula
  A[["data"]] <- new.data
  A[["s.weights"]] <- s.weights[subset]
  A[["verbose"]] <- TRUE


  out <- do.call(optweight::optweight, A, quote = TRUE)
  obj <- list(w = out[["weights"]], fit.obj = out)
  return(obj)

}
weightit2optweight.msm <- function(covs.list, treat.list, s.weights, subset, moments, int, ...) {
  A <- list(...)
  check.package("optweight")
  if (is_not_null(covs.list)) {
    covs.list <- lapply(covs.list, function(c) {
      covs <- c[subset, , drop = FALSE]
      covs <- cbind(covs, int.poly.f(covs, poly = moments, int = int))
      covs <- apply(covs, 2, make.closer.to.1)

      if (any(vars.w.missing <- apply(covs, 2, function(x) any(is.na(x))))) {
        missing.ind <- apply(covs[, vars.w.missing, drop = FALSE], 2, function(x) as.numeric(is.na(x)))
        covs[is.na(covs)] <- 0
        covs <- cbind(covs, missing.ind)
      }
      return(covs)
    })
  }
  if (is_not_null(treat.list)) {
    treat.list <- lapply(treat.list, function(t) {
      treat <- t[subset]
      if (attr(t, "treat.type") != "continuous") treat <- factor(treat)
      return(treat)
    })
  }
  if (is_not_null(s.weights)) {
    s.weights <- s.weights[subset]
  }

  baseline.data <- cbind(treat.list[[1]], covs.list[[1]])
  baseline.formula <- formula(baseline.data)
  if ("tols" %in% names(A)) A[["tols"]] <- optweight::check.tols(baseline.formula, baseline.data, A[["tols"]], stop = TRUE)
  if ("targets" %in% names(A)) {
    warning("targets cannot be used through WeightIt and will be ignored.", call. = FALSE)
    A[["targets"]] <- NULL
  }

  out <- do.call(optweight::optweight.fit, c(list(treat = treat.list,
                                                  covs = covs.list,
                                                  s.weights = s.weights,
                                                  verbose = TRUE),
                                             A), quote = TRUE)
  obj <- list(w = out$w, fit.obj = out)
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

  w <- rep(1, length(treat))
  ps <- NULL

  if (check.package("twang")) {
    if (estimand == "ATT") {

      control.levels <- levels(treat)[levels(treat) != focal]
      fit.list <- setNames(vector("list", length(control.levels)), control.levels)

      for (i in control.levels) {
        treat.in.i.focal <- treat %in% c(focal, i)
        treat_ <- ifelse(treat[treat.in.i.focal] == i, 0, 1)
        covs_ <- covs[treat.in.i.focal, , drop = FALSE]
        new.data <- data.frame(treat_, covs_)

        fit.list[[i]] <- twang::ps(formula(new.data),
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

        s <- fit.list[[i]]$stopMethods[1]

        w[treat == i] <- cobalt::get.w(fit.list[[i]], stop.method = s)[treat_ == 0]

        if (nunique(treat) == 2) {
          ps <- fit.list[[i]][["ps"]][[1]]
          fit.list <- fit.list[[i]]
        }

      }
    }
    else if (estimand == "ATE") {
      fit.list <- setNames(vector("list", nlevels(treat)), levels(treat))

      for (i in levels(treat)) {
        #Mimicking twang
        #Seeks balance between weighted treat group and all others combined
        #Note: Gives different answer than twang for some reason; has better balance though.
        treat_i <- ifelse(treat == i, 0, 1)
        new.data <- data.frame(treat_i, covs)

        fit.list[[i]] <- twang::ps(formula(new.data),
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

        s <- fit.list[[i]]$stopMethods[1]

        if (nunique(treat) == 2) {
          w <- cobalt::get.w(fit.list[[i]], stop.method = s)
          ps <- fit.list[[i]][["ps"]][[1]]
          fit.list <- fit.list[[i]]
          break
        }
        else {
          w[treat == i] <- cobalt::get.w(fit.list[[i]], stop.method = s)[treat == i]
        }
      }
    }

    if (stabilize) {
      tab <- vapply(levels(t), function(x) mean(t == x), numeric(1L))
      w <- w * tab[t]
    }
  }
  obj <- list(w = w, ps = ps, fit.obj = fit.list)
  return(obj)
}
weightit2gbm.cont <- function(covs, treat, s.weights, subset, stabilize, ...) {
  A <- list(...)
  A[c("formula", "data", "sampw", "verbose")] <- NULL
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

  covs <- covs[subset, , drop = FALSE]
  treat <- treat[subset]

  covs <- apply(covs, 2, make.closer.to.1)

  new.data <- data.frame(treat, covs)

  if (check.package("wCorr") && check.package("gbm")) {
    fit <- do.call("ps.cont", c(list(formula(new.data),
                                     data = new.data,
                                     sampw = s.weights[subset],
                                     verbose = TRUE), A))
    w <- cobalt::get.w(fit, stop.method = A[["stop.method"]])
  }

  #ps <- fit[["ps"]][[A[["stop.method"]]]]

  obj <- list(w = w, fit.obj = fit)
  return(obj)
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
  ps <- NULL

  if (check.package("CBPS")) {
    if (estimand == "ATT") {
      w <- rep(1, length(treat))
      control.levels <- levels(treat)[levels(treat) != focal]
      fit.list <- setNames(vector("list", length(control.levels)), control.levels)

      for (i in control.levels) {
        treat.in.i.focal <- treat %in% c(focal, i)
        treat_ <- ifelse(treat[treat.in.i.focal] == i, 0, 1)
        covs_ <- covs[treat.in.i.focal, , drop = FALSE]
        new.data <- data.frame(treat_, covs_)

        tryCatch({fit.list[[i]] <- CBPS::CBPS(formula(new.data),
                                              data = new.data,
                                              method = if (is_not_null(A$over) && A$over == FALSE) "exact" else "over",
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

        w[treat == i] <- cobalt::get.w(fit.list[[i]], estimand = "ATT")[treat_ == 0] / s.weights[subset][treat.in.i.focal][treat_ == 0]

        if (nlevels(treat) == 2) {
          ps <- fit.list[[i]][["fitted.values"]]
          fit.list <- fit.list[[1]]
        }
      }
    }
    else {
      new.data <- data.frame(treat, covs)
      if (!nunique.gt(treat, 4)) {
        tryCatch({fit.list <- CBPS::CBPS(formula(new.data),
                                         data = new.data,
                                         method = if (is_not_null(A$over) && A$over == FALSE) "exact" else "over",
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

        w <- cobalt::get.w(fit.list, estimand = "ATE") / s.weights[subset]
        if (nunique(treat) == 2) {
          ps <- fit.list[["fitted.values"]]
          fit.list <- fit.list[[1]]
        }
      }
      else {
        w <- rep(1, length(treat))
        fit.list <- setNames(vector("list", nlevels(treat)), levels(treat))

        for (i in levels(treat)) {
          new.data[[1]] <- ifelse(treat == i, 1, 0)
          fit.list[[i]] <- CBPS::CBPS(formula(new.data), data = new.data,
                                      method = if (is_null(A$over) || A$over == TRUE) "over" else "exact",
                                      standardize = FALSE,
                                      sample.weights = s.weights[subset],
                                      ATT = 0, ...)

          w[treat==i] <- cobalt::get.w(fit.list[[i]], estimand = "ATE")[treat==i] / s.weights[subset][treat==i]
        }
      }
    }

  }
  if (stabilize) {
    tab <- vapply(levels(t), function(x) mean(t == x), numeric(1L))
    w <- w * tab[t]
  }

  obj <- list(w = w, ps = ps, fit.obj = fit.list)
  return(obj)
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

  obj <- list(w = w, fit.obj = fit)
  return(obj)
}
weightit2cbps.msm <- function(covs.list, treat.list, s.weights, subset, ...) {
  stop("CBMSM doesn't work yet.")
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
  colinear.covs.to.remove <- colnames(covs)[colnames(covs) %nin% colnames(make.full.rank(covs))]
  covs <- covs[, colnames(covs) %nin% colinear.covs.to.remove, drop = FALSE]

  new.data <- data.frame(treat = treat, covs)

  if (check.package("CBPS")) {
    fit <- do.call(CBPS::npCBPS, c(list(formula(new.data), data = new.data, print.level = 1), A),
                   quote = TRUE)
  }
  w <- cobalt::get.w(fit)

  obj <- list(w = w, fit.obj = fit)

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
  colinear.covs.to.remove <- colnames(covs)[colnames(covs) %nin% colnames(make.full.rank(covs))]
  covs <- covs[, colnames(covs) %nin% colinear.covs.to.remove, drop = FALSE]

  new.data <- data.frame(treat = treat, covs)

  if (check.package("CBPS")) {
    fit <- do.call(CBPS::npCBPS, c(list(formula(new.data), data = new.data, print.level = 1), A),
                   quote = TRUE)
  }
  w <- cobalt::get.w(fit)

  obj <- list(w = w, fit.obj = fit)

  return(obj)
}

#Entropy balancing with ebal
weightit2ebal <- function(covs, treat, s.weights, subset, estimand, focal, stabilize, moments, int, ...) {
  A <- list(...)

  covs <- covs[subset, , drop = FALSE]
  treat <- factor(treat)[subset]

  covs <- cbind(covs, int.poly.f(covs, poly = moments, int = int, center = TRUE))
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
      fit.list <- setNames(vector("list", length(control.levels)), control.levels)

      covs[treat == focal,] <- covs[treat == focal, , drop = FALSE] * s.weights[subset][treat == focal] * sum(treat == focal)/sum(s.weights[subset][treat == focal])

      for (i in control.levels) {
        treat.in.i.focal <- treat %in% c(focal, i)
        treat_ <- ifelse(treat[treat.in.i.focal] == i, 0, 1)
        covs_ <- covs[treat.in.i.focal, , drop = FALSE]

        colinear.covs.to.remove <- colnames(covs_)[colnames(covs_) %nin% colnames(make.full.rank(covs_[treat_ == 0, , drop = FALSE]))]

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
        fit.list[[i]] <- ebal.out
      }
    }
    else if (estimand == "ATE") {
      w <- rep(1, length(treat))
      fit.list <- setNames(vector("list", nlevels(treat)), levels(treat))

      for (i in levels(treat)) {
        covs_i <- rbind(covs, covs[treat==i, , drop = FALSE])
        treat_i <- c(rep(1, nrow(covs)), rep(0, sum(treat==i)))

        colinear.covs.to.remove <- colnames(covs_i)[colnames(covs_i) %nin% colnames(make.full.rank(covs_i[treat_i == 0, , drop = FALSE]))]

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
        fit.list[[i]] <- ebal.out_i
      }
    }
  }

  obj <- list(w = w, fit.obj = fit.list)
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
      fit.list <- setNames(vector("list", length(control.levels)), control.levels)

      for (i in control.levels) {
        treat.in.i.focal <- treat %in% c(focal, i)
        treat_ <- ifelse(treat[treat.in.i.focal] == i, 0, 1)
        covs_ <- covs[treat.in.i.focal, , drop = FALSE]

        colinear.covs.to.remove <- colnames(covs_)[colnames(covs_) %nin% colnames(make.full.rank(covs_[treat_ == 0, , drop = FALSE]))]

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
        fit.list[[i]] <- ate.out
      }
    }
    else if (estimand == "ATE") {
      w <- rep(1, length(treat))
      fit.list <- setNames(vector("list", nlevels(treat)), levels(treat))

      for (i in levels(treat)) {
        covs_i <- rbind(covs, covs[treat==i, , drop = FALSE])
        treat_i <- c(rep(1, nrow(covs)), rep(0, sum(treat==i)))

        colinear.covs.to.remove <- colnames(covs_i)[colnames(covs_i) %nin% colnames(make.full.rank(covs_i[treat_i == 0, , drop = FALSE]))]

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
        fit.list[[i]] <- ate.out
      }
    }
    obj <- list(w = w, fit.obj = fit.list)
    return(obj)
  }
}

#PS weights using SuperLearner
weightit2super <- function(covs, treat, s.weights, subset, estimand, focal, stabilize, ...) {
  A <- list(...)

  check.package("SuperLearner")

  covs <- covs[subset, , drop = FALSE]
  t <- factor(treat[subset])

  if (any(vars.w.missing <- apply(covs, 2, function(x) any(is.na(x))))) {
    missing.ind <- apply(covs[, vars.w.missing, drop = FALSE], 2, function(x) as.numeric(is.na(x)))
    covs[is.na(covs)] <- 0
    covs <- cbind(covs, missing.ind)
  }
  covs <- data.frame(apply(covs, 2, make.closer.to.1))

  if (ncol(covs) > 1) {
    colinear.covs.to.remove <- colnames(covs)[colnames(covs) %nin% colnames(make.full.rank(covs))]
    covs <- covs[, colnames(covs) %nin% colinear.covs.to.remove, drop = FALSE]
  }

  for (f in names(formals(SuperLearner::SuperLearner))) {
    if (f == "method") {if (is_null(A[["SL.method"]])) A[["SL.method"]] <- formals(SuperLearner::SuperLearner)[["method"]]}
    else if (f == "env") {if (is_null(A[["env"]])) A[["env"]] <- environment(SuperLearner::SuperLearner)}
    else if (is_null(A[[f]])) A[[f]] <- formals(SuperLearner::SuperLearner)[[f]]
  }

  ps <- setNames(as.data.frame(matrix(NA_real_, ncol = nunique(t), nrow = length(t))),
                 levels(t))

  if (is_binary(t)) {

    fit.list <- do.call(SuperLearner::SuperLearner, list(Y = as.numeric(as.character(t)),
                                                         X = covs, newX = covs,
                                                         family = binomial(),
                                                         SL.library = A[["SL.library"]],
                                                         verbose = FALSE,
                                                         method = A[["SL.method"]],
                                                         id = NULL,
                                                         obsWeights = s.weights[subset],
                                                         control = A[["control"]],
                                                         cvControl = A[["cvControl"]],
                                                         env = A[["env"]]))
    ps[[2]] <- p.score <- fit.list$SL.predict
    ps[[1]] <- 1 - ps[[2]]

  }
  else {
    fit.list <- setNames(vector("list", nlevels(t)), levels(t))

    for (i in levels(t)) {
      t_i <- rep(0, length(t)); t_i[t == i] <- 1

      fit.list[[i]] <- do.call(SuperLearner::SuperLearner, list(Y = t_i,
                                                                X = covs, newX = covs,
                                                                family = binomial(),
                                                                SL.library = A[["SL.library"]],
                                                                verbose = FALSE,
                                                                method = A[["SL.method"]],
                                                                id = NULL,
                                                                obsWeights = s.weights[subset],
                                                                control = A[["control"]],
                                                                cvControl = A[["cvControl"]],
                                                                env = A[["env"]]))
      ps[[i]] <- fit.list[[i]]$SL.predict
    }

    p.score <- NULL
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
    w <- w/apply(ps, 1, function(x) sum(1/x)) #Li & Li (2018)
  }
  else if (toupper(estimand) == "ATM") {
    w <- w*apply(ps, 1, min)
  }
  else w <- NULL

  if (stabilize) {
    tab <- vapply(levels(t), function(x) mean(t == x), numeric(1L))
    w <- w * tab[t]
  }

  obj <- list(w = w, ps = p.score, fit.obj = fit.list)
  return(obj)
}
weightit2super.cont <- function(covs, treat, s.weights, subset, stabilize, ps, ...) {
  A <- B <- list(...)

  covs <- covs[subset, , drop = FALSE]
  t <- treat[subset]

  if (any(vars.w.missing <- apply(covs, 2, function(x) any(is.na(x))))) {
    missing.ind <- apply(covs[, vars.w.missing, drop = FALSE], 2, function(x) as.numeric(is.na(x)))
    covs[is.na(covs)] <- 0
    covs <- cbind(covs, missing.ind)
  }

  covs <- data.frame(apply(covs, 2, make.closer.to.1))

  if (ncol(covs) > 1) {
    colinear.covs.to.remove <- colnames(covs)[colnames(covs) %nin% colnames(make.full.rank(covs))]
    covs <- covs[, colnames(covs) %nin% colinear.covs.to.remove, drop = FALSE]
  }

  stabilize <- TRUE

  for (f in names(formals(SuperLearner::SuperLearner))) {
    if (f == "method") {if (is_null(B[["SL.method"]])) B[["SL.method"]] <- formals(SuperLearner::SuperLearner)[["method"]]}
    else if (f == "env") {if (is_null(B[["env"]])) B[["env"]] <- environment(SuperLearner::SuperLearner)}
    else if (is_null(B[[f]])) B[[f]] <- formals(SuperLearner::SuperLearner)[[f]]
  }

  fit <- do.call(SuperLearner::SuperLearner, list(Y = t,
                                                  X = covs, newX = covs,
                                                  family = gaussian(),
                                                  SL.library = B[["SL.library"]],
                                                  verbose = FALSE,
                                                  method = B[["SL.method"]],
                                                  id = NULL,
                                                  obsWeights = s.weights[subset],
                                                  control = B[["control"]],
                                                  cvControl = B[["cvControl"]],
                                                  env = B[["env"]]))
  p.denom <- t - fit$SL.predict

  if (isTRUE(A[["use.kernel"]])) {
    if (is_null(A[["bw"]])) A[["bw"]] <- "nrd0"
    if (is_null(A[["adjust"]])) A[["adjust"]] <- 1
    if (is_null(A[["kernel"]])) A[["kernel"]] <- "gaussian"
    if (is_null(A[["n"]])) A[["n"]] <- 10*length(t)

    d.d <- density(p.denom, n = A[["n"]],
                   weights = s.weights[subset]/sum(s.weights[subset]), give.Rkern = FALSE,
                   bw = A[["bw"]], adjust = A[["adjust"]], kernel = A[["kernel"]])
    if (isTRUE(A[["plot"]])) {
      par(mfrow=c(2,1))
      plot(d.d, main = "Denominator density")
    }
    dens.denom <- with(d.d, approxfun(x = x, y = y))(p.denom)
  }
  else {
    dens.denom <- dnorm(p.denom, 0, sd = sd(p.denom))
  }

  if (stabilize) {
    if (is_null(A$link)) A$link <- "identity"
    num.fit <- do.call("glm", c(list(t ~ 1,
                                     data = data.frame(t = t),
                                     weights = s.weights[subset],
                                     family = gaussian(link = A$link),
                                     control = list()),
                                A[names(A %nin% "family")]),
                       quote = TRUE)

    p.num <- t - num.fit$fitted.values

    if (isTRUE(A[["use.kernel"]])) {
      d.n <- density(p.num, n = A[["n"]],
                     weights = s.weights[subset]/sum(s.weights[subset]), give.Rkern = FALSE,
                     bw = A[["bw"]], adjust = A[["adjust"]], kernel = A[["kernel"]])
      if (isTRUE(A[["plot"]])) plot(d.n, main = "Numerator density")
      dens.num <- with(d.n, approxfun(x = x, y = y))(p.num)
    }
    else {
      dens.num <- dnorm(p.num, 0, sqrt(summary(num.fit)$dispersion))
    }
    w <- dens.num/dens.denom
  }
  else {
    w <- 1/dens.denom
  }

  obj <- list(w = w, fit.obj = fit)
  return(obj)
}

#Stable balancing weights with sbw
weightit2sbw <- function(...) {
  stop("method = \"sbw\" has been deprecated in place of method = \"optweight\". Please use that instead.", call. = FALSE)
  #stop("Stable balancing weights are not currently supported. Please choose another method.\n        The github version of WeightIt may allow stable balancing weights.\n        Install it with devtools::install_github(\"ngreifer/WeightIt\").", call. = FALSE)
}

#SBW--------
# weightit2sbw <- function(covs, treat, s.weights, subset, estimand, focal, moments, int, ...) {
#   A <- list(...)
#
#   if (check.package("sbw")) {
#     check.package("slam")
#     if (!"package:slam" %in% search()) {
#       need.to.detach.slam <- TRUE
#       attachNamespace("slam")
#     } else {
#       need.to.detach.slam <- FALSE
#     }
#
#     if (is_null(A$l_norm)) A$l_norm <- "l_2"
#     if (is_null(A$solver)) A$solver <- "quadprog"
#     if (is_null(A$max_iter)) A$max_iter <- 100000
#     if (is_null(A$rel_tol)) A$rel_tol <- 1e-4
#     if (is_null(A$abs_tol)) A$abs_tol <- 1e-4
#     if (is_null(A$gap_stop)) A$gap_stop <- TRUE
#     if (is_null(A$adaptive_rho)) A$adaptive_rho <- TRUE
#
#     solve.package = switch(A$solver, quadprog = "quadprog",
#                            cplex = "Rcplex",
#                            gurobi = "gurobi",
#                            pogs = "pogs")
#     check.package(solve.package)
#     if (!paste0("package:", solve.package) %in% search()) {
#       need.to.detach.solver <- TRUE
#       attachNamespace(solve.package)
#     } else {
#       need.to.detach.solver <- FALSE
#     }
#
#     covs <- covs[subset, , drop = FALSE]
#     treat <- factor(treat)[subset]
#
#     if (any(is.na(covs))) {
#       stop("Stable balancing weights are not compatible with missing values.", call. = FALSE)
#     }
#     covs <- cbind(covs, int.poly.f(covs, poly = moments, int = int))
#     covs <- apply(covs, 2, make.closer.to.1)
#
#     #new.data <- setNames(data.frame(treat, covs), as.character(seq_len(1+ncol(covs))))
#
#     binary.vars <- apply(covs, 2, is_binary)
#
#     if (is_null(A$bal_tols)) bal_tols <- .0001
#     else {
#       bal_tols <- A$bal_tols
#       if (length(bal_tols) != 1 && length(bal_tols) != ncol(covs)) {
#         stop(paste0("bal_tols needs to be of length 1 or equal to the number of covariates (", ncol(covs),
#                     ").\nThe covariates (in order) are:\n   ", paste0(colnames(covs), collapse = " ")), call.= FALSE)
#       }
#     }
#     if (is_null(A$bal_tols_sd)) bal_tols_sd <- TRUE
#     else bal_tols_sd <- A$bal_tols_sd
#
#     if (is_null(A$bal_tols) && is_null(A$bal_tols_sd)) {
#       message("Using bal_tols = 0.0001 and bal_tols_sd = TRUE.")
#     }
#     else if (is_null(A$bal_tols)) {
#       message("Using bal_tols = 0.0001.")
#     }
#     else if (is_null(A$bal_tols_sd)) {
#       message("Using bal_tols_sd = TRUE.")
#     }
#
#     if (estimand == "ATT") {
#       control.levels <- levels(treat)[levels(treat) != focal]
#
#       w <- rep(1, length(treat))
#
#       if (bal_tols_sd) {
#         bal_tols <- bal_tols * sapply(1:ncol(covs), function(x) {if (binary.vars[x]) 1 else sqrt(cov.wt(covs[treat == focal, x, drop = FALSE], s.weights[subset][treat == focal])$cov[1,1])})
#       }
#
#       covs[treat == focal,] <- covs[treat == focal,] * s.weights[subset][treat == focal] * sum(treat == focal)/sum(s.weights[subset][treat == focal])
#
#       for (i in control.levels) {
#
#         treat.in.i.focal <- treat %in% c(focal, i)
#         treat_ <- ifelse(treat[treat.in.i.focal] == i, 0, 1)
#         covs_ <- covs[treat.in.i.focal, , drop = FALSE]
#
#         new.data_ <- data.frame(treat_, covs_)
#         t_ind <- names(new.data_)[1]
#         bal_covs = names(new.data_)[-1]
#
#         sbw.fit <- sbw::sbw(new.data_,
#                             t_ind = t_ind,
#                             bal_covs = bal_covs,
#                             bal_tols = bal_tols,
#                             bal_tols_sd = FALSE,
#                             target = "treated",
#                             l_norm = A[["l_norm"]],
#                             w_min = 0,
#                             normalize = TRUE,
#                             solver = A[["solver"]],
#                             display = 1,
#                             max_iter = A[["max_iter"]],
#                             rel_tol = A[["rel_tol"]],
#                             abs_tol = A[["abs_tol"]],
#                             gap_stop = A[["gap_stop"]],
#                             adaptive_rho = A[["adaptive_rho"]])
#
#         w[treat==i] <- sbw.fit$data_frame_weights$weights[sbw.fit$data_frame_weights[[t_ind]] == 0]*sum(treat == i) / s.weights[subset][treat == i]
#       }
#       w[w < 0] <- 0
#     }
#     else if (estimand == "ATE") {
#       if (bal_tols_sd) {
#         bal_tols <- bal_tols * sapply(1:ncol(covs), function(x) {if (binary.vars[x]) 1 else sqrt(mean(sapply(unique(treat), function(t) cov.wt(covs[treat == t, x, drop = FALSE], s.weights[subset][treat == t])$cov[1,1])))})
#       }
#
#       bal_tols <- bal_tols/nunique(treat)
#
#       w <- rep(1, length(treat))
#
#       for (i in levels(treat)) {
#         covs_i <- rbind(covs, covs[treat==i, , drop = FALSE])
#         treat_i <- c(rep(1, nrow(covs)), rep(0, sum(treat==i)))
#
#         covs_i[treat_i == 1,] <- covs_i[treat_i == 1, , drop = FALSE] * s.weights[subset] * sum(treat_i == 1) / sum(s.weights[subset])
#
#         new.data_i <- data.frame(treat_i, covs_i)
#         t_ind <- names(new.data_i)[1]
#         bal_covs = names(new.data_i)[-1]
#
#         sbw.fit_i <- sbw::sbw(new.data_i, t_ind = t_ind,
#                               bal_covs = bal_covs,
#                               bal_tols = bal_tols,
#                               bal_tols_sd = FALSE,
#                               target = "treated",
#                               l_norm = A$l_norm,
#                               w_min = 0,
#                               normalize = TRUE,
#                               solver = A[["solver"]],
#                               display = 1,
#                               max_iter = A[["max_iter"]],
#                               rel_tol = A[["rel_tol"]],
#                               abs_tol = A[["abs_tol"]],
#                               gap_stop = A[["gap_stop"]],
#                               adaptive_rho = A[["adaptive_rho"]])
#
#         w[treat==i] <- sbw.fit_i$data_frame_weights$weights[sbw.fit_i$data_frame_weights[[t_ind]] == 0]*sum(treat == i) / s.weights[subset][treat == i]
#
#       }
#       w[w < 0] <- 0
#     }
#   }
#
#   if (need.to.detach.slam) detach("package:slam", character.only = TRUE)
#   if (need.to.detach.solver) detach(paste0("package:", solve.package), character.only = TRUE)
#
#   obj <- list(w = w)
#   return(obj)
#
# }

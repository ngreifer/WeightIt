#User-defined weighting function
weightit2user <- function(Fun, covs, treat, s.weights, subset, estimand, focal, stabilize, subclass, missing, ps, moments, int, ...) {
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

      if (length(Anames) > 1) warning(paste0("The following arguments were specified but are not suitable arguments to the provided function:\n\t", word_list(Anames)), call. = FALSE)
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
weightitMSM2user <- function(Fun, covs.list, treat.list, s.weights, subset, stabilize, missing, moments, int, ...) {
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

      if (length(Anames) > 1) warning(paste0("The following arguments were specified but are not suitable arguments to the provided function:\n\t", word_list(Anames)), call. = FALSE)
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
weightit2ps <- function(covs, treat, s.weights, subset, estimand, focal, stabilize, subclass, missing, ps, ...) {
  A <- list(...)

  fit.obj <- NULL

  if (is_null(ps)) {

    covs <- covs[subset, , drop = FALSE]
    treat_sub <- factor(treat[subset])
    bin.treat <- is_binary(treat_sub)
    ord.treat <- is.ordered(treat_sub)

    if (missing == "ind") {
      missing.ind <- apply(covs[, apply(covs, 2, anyNA), drop = FALSE], 2, function(x) as.numeric(is.na(x)))
      covs[is.na(covs)] <- 0
      covs <- cbind(covs, missing.ind)
    }

    covs <- apply(covs, 2, make.closer.to.1)

    if (ncol(covs) > 1) {
      if (missing == "saem") {
        covs0 <- covs
        for (i in colnames(covs)[apply(covs, 2, anyNA)]) covs0[is.na(covs0[,i]),i] <- covs0[!is.na(covs0[,i]),i][1]
        colinear.covs.to.remove <- colnames(covs)[colnames(covs) %nin% colnames(make_full_rank(covs0))]
      }
      else colinear.covs.to.remove <- colnames(covs)[colnames(covs) %nin% colnames(make_full_rank(covs))]
      covs <- covs[, colnames(covs) %nin% colinear.covs.to.remove, drop = FALSE]
    }

    if (is_null(A$link)) A$link <- "logit"
    else {
      if (bin.treat || isFALSE(A$use.mlogit)) {
        if (missing == "saem") acceptable.links <- "logit"
        else acceptable.links <- expand.grid_string(c("", "br."), c("logit", "probit", "cloglog", "identity", "log", "cauchit"))
      }
      else if (ord.treat) acceptable.links <- c("logit", "probit", "loglog", "cloglog", "cauchit")
      else acceptable.links <- c("logit", "probit", "bayes.probit", "br.logit")

      which.link <- acceptable.links[pmatch(A$link, acceptable.links, nomatch = 0)][1]
      if (is.na(which.link)) {
        A$link <- acceptable.links[1]
        warning(paste0("Only ", word_list(acceptable.links, quotes = TRUE, is.are = TRUE), " allowed as the link for ",
                       if (bin.treat) "binary" else if (ord.treat) "ordinal" else "multinomial",
                       " treatments", if (missing == "saem") " with missing = \"saem\"", ". Using link = ", word_list(acceptable.links[1], quotes = TRUE), "."),
                call. = FALSE, immediate. = TRUE)
      }
      else A$link <- which.link
    }

    if (startsWith(A$link, "br.")) {
      A$link <- substr(A$link, 4, nchar(A$link))
      use.br <- TRUE
    }
    else use.br <- FALSE

    if (bin.treat) {

      ps <- setNames(as.data.frame(matrix(NA_real_, ncol = nunique(treat_sub), nrow = length(treat_sub))),
                     levels(treat_sub))

      if (missing == "saem") {
        check.package("misaem")
        for (f in names(formals(misaem::miss.saem))) {
          if (f %in% c("X.obs", "y", "pos_var", "print_iter", "var_cal", "ll_obs_cal", "seed")) {A[[f]] <- NULL}
          else if (is_null(A[[f]])) A[[f]] <- formals(misaem::miss.saem)[[f]]
        }
        fit <- do.call(misaem::miss.saem, c(list(X.obs = covs,
                                                 y = binarize(treat_sub),
                                                 seed = NULL,
                                                 var_cal = FALSE),
                                            A[names(A) %in% names(formals(misaem::miss.saem))]))

        if (is_null(A[["saem_method"]])) A[["saem_method"]] <- formals(misaem::pred_saem)[["method"]]

        p.score <- misaem::pred_saem(covs, fit$beta, fit$mu, fit$sig2,
                                     seed = NULL, method = A[["saem_method"]])
        ps[[2]] <- p.score
        ps[[1]] <- 1 - ps[[2]]
      }
      else {
        data <- data.frame(binarize(treat_sub), covs)
        formula <- formula(data)

        if (use.br) {
          check.package("brglm2")
          control <- A[names(formals(brglm2::brglmControl))[pmatch(names(A), names(formals(brglm2::brglmControl)), 0)]]
          fit <- do.call("glm", c(list(formula, data = data,
                                       weights = s.weights[subset],
                                       family = binomial(link = A[["link"]]),
                                       method = brglm2::brglmFit),
                                  control), quote = TRUE)
        }
        else {
          if (is_null(A[["glm.method"]])) A[["glm.method"]] <- "glm.fit"
          control <- A[names(formals(glm.control))[pmatch(names(A), names(formals(glm.control)), 0)]]
          fit <- do.call("glm", c(list(formula, data = data,
                                       weights = s.weights[subset],
                                       family = quasibinomial(link = A[["link"]]),
                                       method = A[["glm.method"]]),
                                  control), quote = TRUE)
        }

        ps[[2]] <- p.score <- fit$fitted.values
        ps[[1]] <- 1 - ps[[2]]
      }

      fit.obj <- fit
    }
    else if (ord.treat) {
      if (A[["link"]] == "logit") A[["link"]] <- "logistic"
      check.package("MASS")
      # message(paste("Using ordinal", A$link, "regression."))
      data <- data.frame(treat_sub, covs)
      formula <- formula(data)
      # control <- A[names(A) %nin% c("link", names(formals(MASS::polr)))]
      tryCatch({fit <- do.call(MASS::polr,
                               c(list(formula,
                                      data = data,
                                      weights = s.weights,
                                      Hess = FALSE,
                                      model = FALSE,
                                      method = A[["link"]],
                                      contrasts = NULL)))},
               error = function(e) {stop(paste0("There was a problem fitting the ordinal ", A$link, " regressions with polr().\n       Try again with an un-ordered treatment.",
                                                "\n       Error message: ", conditionMessage(e)), call. = FALSE)})

      ps <- fit$fitted.values
      fit.obj <- fit
      p.score <- NULL
    }
    else {
      if (use.br) {
        check.package("brglm2")
        data <- data.frame(treat_sub, covs)
        formula <- formula(data)
        control <- A[names(formals(brglm2::brglmControl))[pmatch(names(A), names(formals(brglm2::brglmControl)), 0)]]
        tryCatch({fit <- do.call(brglm2::brmultinom,
                                 c(list(formula, data,
                                        weights = s.weights),
                                   control))},
                 error = function(e) stop("There was a problem with the bias-reduced logit regression. Try a different link.", call. = FALSE))

        ps <- fit$fitted.values
        fit.obj <- fit
      }
      else if (A$link %in% c("logit", "probit")) {
        if (!isFALSE(A$use.mlogit)) {
          if (check.package("mlogit")) {
            data <- data.frame(treat = treat_sub , s.weights = s.weights[subset], covs)
            covnames <- names(data)[-c(1,2)]
            mult <- mlogit::mlogit.data(data, varying = NULL, shape = "wide", sep = "", choice = "treat")
            tryCatch({fit <- mlogit::mlogit(as.formula(paste0("treat ~ 1 | ", paste(covnames, collapse = " + "),
                                                              " | 1")), data = mult, estimate = TRUE,
                                            probit = ifelse(A$link[1] == "probit", TRUE, FALSE),
                                            weights = s.weights, ...)},
                     error = function(e) {stop(paste0("There was a problem fitting the multinomial ", A$link, " regressions with mlogit().\n       Try again with use.mlogit = FALSE."), call. = FALSE)}
            )
            ps <- fitted(fit, outcome = FALSE)
            fit.obj <- fit
          }
        }
        else {
          ps <- setNames(as.data.frame(matrix(NA_real_, ncol = nunique(treat_sub), nrow = length(treat_sub))),
                         levels(treat_sub))

          control <- A[names(formals(glm.control))[pmatch(names(A), names(formals(glm.control)), 0)]]
          fit.list <- setNames(vector("list", nlevels(treat_sub)), levels(treat_sub))
          for (i in levels(treat_sub)) {
            t_i <- as.numeric(treat_sub == i)
            data_i <- data.frame(t_i, covs)
            fit.list[[i]] <- do.call(glm, c(list(formula(data_i), data = data_i,
                                                 family = quasibinomial(link = A$link),
                                                 weights = s.weights[subset]),
                                            control))
            ps[[i]] <- fit.list[[i]]$fitted.values
          }
          fit.obj <- fit.list
        }
      }
      else if (A$link == "bayes.probit") {
        check.package("MNP")
        data <- data.frame(treat_sub, covs)
        formula <- formula(data)
        tryCatch({fit <- MNP::mnp(formula, data, verbose = TRUE)},
                 error = function(e) stop("There was a problem with the Bayes probit regression. Try a different link.", call. = FALSE))
        ps <- MNP::predict.mnp(fit, type = "prob")$p
        fit.obj <- fit
      }
      else {
        stop('link must be "logit", "probit", "bayes.probit", or "br.logit.', call. = FALSE)
      }
      p.score <- NULL
    }
  }
  else {
    n <- length(treat)
    p.score <- NULL
    treat <- factor(treat)
    treat_sub <- factor(treat[subset])
    bin.treat <- is_binary(treat_sub)

    if (bin.treat) {
      bad.ps <- FALSE
      if (is.matrix(ps) || is.data.frame(ps)) {
        if (dim(ps) == c(n, 2)) {
          ps <- setNames(as.data.frame(ps), levels(treat))[subset, , drop = FALSE]
          p.score <- ps[[2]]
        }
        else if (dim(ps) == c(length(treat), 1)) {
          ps <- setNames(as.data.frame(matrix(c(1-ps[,1], ps[,1]), ncol = 2)),
                         levels(treat))[subset, , drop = FALSE]
          p.score <- ps[[2]]
        }
        else {
          bad.ps <- TRUE
        }
      }
      else if (is.numeric(ps)) {
        if (length(ps) == n) {
          ps <- setNames(as.data.frame(matrix(c(1-ps, ps), ncol = 2)),
                         levels(treat))[subset, , drop = FALSE]
          p.score <- ps[[2]]
        }
        else {
          bad.ps <- TRUE
        }
      }
      else bad.ps <- TRUE

      if (bad.ps) stop("ps must be a numeric vector with a propensity score for each unit.", call. = FALSE)

    }
    else {
      bad.ps <- FALSE
      if (is.matrix(ps) || is.data.frame(ps)) {
        if (dim(ps) == c(n, nunique(treat))) {
          ps <- setNames(as.data.frame(ps), levels(treat))[subset, , drop = FALSE]
        }
        else if (dim(ps) == c(n, 1)) {
          ps <- setNames(do.call("cbind", lapply(levels(treat), function(x) {
            p_ <- rep(1, length(treat))
            p_[treat == x] <- ps[treat == x, 1]
            return(p_)
          })), levels(treat))[subset, , drop = FALSE]
        }
        else {
          bad.ps <- TRUE
        }
      }
      else if (is.numeric(ps)) {
        if (length(ps) == n) {
          ps <- setNames(do.call("cbind", lapply(levels(treat), function(x) {
            p_ <- rep(1, length(treat))
            p_[treat == x] <- ps[treat == x]
            return(p_)
          })), levels(treat))[subset, , drop = FALSE]
        }
        else {
          bad.ps <- TRUE
        }
      }
      else bad.ps <- TRUE

      if (bad.ps) stop("ps must be a numeric vector with a propensity score for each unit or a matrix \n\twith the probability of being in each treatment for each unit.", call. = FALSE)

    }

  }

  #ps should be matrix of probs for each treat
  #Computing weights
  w <- get_w_from_ps(ps = ps, treat = treat_sub, estimand, focal, stabilize = stabilize, subclass = subclass)

  obj <- list(w = w, ps = p.score, fit.obj = fit.obj)
  return(obj)
}
weightit2ps.cont <- function(covs, treat, s.weights, subset, stabilize, missing, ps, ...) {
  A <- list(...)

  fit.obj <- NULL

  covs <- covs[subset, , drop = FALSE]
  treat <- treat[subset]

  if (missing == "ind") {
    missing.ind <- apply(covs[, apply(covs, 2, anyNA), drop = FALSE], 2, function(x) as.numeric(is.na(x)))
    covs[is.na(covs)] <- 0
    covs <- cbind(covs, missing.ind)
  }
  covs <- apply(covs, 2, make.closer.to.1)
  data <- data.frame(treat, covs)
  formula <- formula(data)

  stabilize <- TRUE

  if (is_null(ps)) {
    if (is_null(A$link)) A$link <- "identity"
    fit <- do.call("glm", c(list(formula, data = data,
                                 weights = s.weights[subset],
                                 family = gaussian(link = A$link),
                                 control = as.list(A$control))),
                   quote = TRUE)
    p.denom <- treat - fit$fitted.values

    fit.obj <- fit

    if (isTRUE(A[["use.kernel"]])) {
      if (is_null(A[["bw"]])) A[["bw"]] <- "nrd0"
      if (is_null(A[["adjust"]])) A[["adjust"]] <- 1
      if (is_null(A[["kernel"]])) A[["kernel"]] <- "gaussian"
      if (is_null(A[["n"]])) A[["n"]] <- 10*length(treat)

      d.d <- density(p.denom, n = A[["n"]],
                     weights = s.weights[subset]/sum(s.weights[subset]), give.Rkern = FALSE,
                     bw = A[["bw"]], adjust = A[["adjust"]], kernel = A[["kernel"]])
      dens.denom <- with(d.d, approxfun(x = x, y = y))(p.denom)
    }
    else {
      if (is_null(A[["density"]])) densfun <- dnorm
      else if (is.function(A[["density"]])) densfun <- A[["density"]]
      else if (is.character(A[["density"]]) && length(A[["density"]] == 1)) {
        splitdens <- strsplit(A[["density"]], "_", fixed = TRUE)[[1]]
        if (exists(splitdens[1], mode = "function", envir = parent.frame())) {
          if (length(splitdens) > 1 && !can_str2num(splitdens[-1])) {
            stop(paste(A[["density"]], "is not an appropriate argument to density because",
                       word_list(splitdens[-1], and.or = "or", quotes = TRUE), "cannot be coerced to numeric."), call. = FALSE)
          }
          densfun <- function(x) {
            tryCatch(do.call(get(splitdens[1]), c(list(x), as.list(str2num(splitdens[-1])))),
                     error = function(e) stop(paste0("Error in applying density:\n  ", conditionMessage(e)), call. = FALSE))
          }
        }
        else {
          stop(paste(A[["density"]], "is not an appropriate argument to density because",
                     splitdens[1], "is not an available function."), call. = FALSE)
        }
      }
      else stop("The argument to density cannot be evaluated as a density function.", call. = FALSE)
      dens.denom <- densfun(p.denom/sqrt(summary(fit)$dispersion))
      if (is_null(dens.denom) || !is.atomic(dens.denom) || anyNA(dens.denom)) {
        stop("There was a problem with the output of density. Try another density function or leave it blank to use the normal density.", call. = FALSE)
      }
      else if (any(dens.denom <= 0)) {
        warning("The input to density may not accept the full range of standardized residuals of the treatment model. Some weights may be NA or zero.", call. = FALSE)
      }
    }

    if (stabilize) {

      num.fit <- do.call("glm", c(list(treat ~ 1, data = data.frame(treat = treat),
                                       weights = s.weights[subset],
                                       family = gaussian(link = A$link),
                                       control = as.list(A$control))),
                         quote = TRUE)

      p.num <- treat - num.fit$fitted.values

      if (isTRUE(A[["use.kernel"]])) {
        d.n <- density(p.num, n = A[["n"]],
                       weights = s.weights[subset]/sum(s.weights[subset]), give.Rkern = FALSE,
                       bw = A[["bw"]], adjust = A[["adjust"]], kernel = A[["kernel"]])
        dens.num <- with(d.n, approxfun(x = x, y = y))(p.num)
      }
      else {
        dens.num <- densfun(p.num/sqrt(summary(num.fit)$dispersion))
        if (is_null(dens.num) || !is.atomic(dens.num) || anyNA(dens.num)) {
          stop("There was a problem with the output of density. Try another density function or leave it blank to use the normal density.", call. = FALSE)
        }
        else if (any(dens.num <= 0) && !any(dens.denom <= 0)) {
          warning("The input to density may not accept the full range of treatment values. Some weights may be NA or zero.", call. = FALSE)
        }
      }
      w <- dens.num/dens.denom
    }
    else {
      w <- 1/dens.denom
    }
  }

  if (isTRUE(A[["use.kernel"]]) && isTRUE(A[["plot"]])) {
    d.d_ <- cbind(as.data.frame(d.d[c("x", "y")]), dens = "Denominator Density", stringsAsfactors = FALSE)
    d.n_ <- cbind(as.data.frame(d.n[c("x", "y")]), dens = "Numerator Density", stringsAsfactors = FALSE)
    d.all <- rbind(d.d_, d.n_)
    d.all$dens <- factor(d.all$dens, levels = c("Numerator Density", "Denominator Density"))
    pl <- ggplot(d.all, aes(x=x,y=y)) + geom_line() + labs(title = "Weight Component Densities", x = "E[Treat|X]", y = "Density") +
      facet_grid(rows = vars(dens)) + theme(panel.background = element_rect(fill = "white"),
                                            panel.border = element_rect(fill = NA, color = "black"),
                                            axis.text.x = element_text(color = "black"),
                                            axis.text.y = element_text(color = "black"),
                                            panel.grid.major = element_blank(),
                                            panel.grid.minor = element_blank()
      )
    print(pl)
  }

  obj <- list(w = w, fit.obj = fit.obj)
  return(obj)
}

#MABW with optweight
weightit2optweight <- function(covs, treat, s.weights, subset, estimand, focal, missing, moments, int, ...) {
  A <- list(...)

  covs <- covs[subset, , drop = FALSE]
  treat <- factor(treat[subset])

  covs <- cbind(covs, int.poly.f(covs, poly = moments, int = int))
  covs <- apply(covs, 2, make.closer.to.1)

  if (missing == "ind") {
    missing.ind <- apply(covs[, apply(covs, 2, anyNA), drop = FALSE], 2, function(x) as.numeric(is.na(x)))
    covs[is.na(covs)] <- 0
    covs <- cbind(covs, missing.ind)
  }

  new.data <- data.frame(treat, covs)
  new.formula <- formula(new.data)

  if (check.package("optweight")) {

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


    out <- do.call(optweight::optweight, A, quote = TRUE)

    obj <- list(w = out[["weights"]], info = list(duals = out$duals), fit.obj = out)
    return(obj)
  }
}
weightit2optweight.cont <- function(covs, treat, s.weights, subset, missing, moments, int, ...) {
  A <- list(...)
  check.package("optweight")

  covs <- covs[subset, , drop = FALSE]
  treat <- treat[subset]

  covs <- cbind(covs, int.poly.f(covs, poly = moments, int = int))
  covs <- apply(covs, 2, make.closer.to.1)

  if (missing == "ind") {
    missing.ind <- apply(covs[, apply(covs, 2, anyNA), drop = FALSE], 2, function(x) as.numeric(is.na(x)))
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

  obj <- list(w = out[["weights"]], info = list(duals = out$duals), fit.obj = out)
  return(obj)

}
weightit2optweight.msm <- function(covs.list, treat.list, s.weights, subset, missing, moments, int, ...) {
  A <- list(...)
  check.package("optweight")
  if (is_not_null(covs.list)) {
    covs.list <- lapply(covs.list, function(c) {
      covs <- c[subset, , drop = FALSE]
      covs <- cbind(covs, int.poly.f(covs, poly = moments, int = int))
      covs <- apply(covs, 2, make.closer.to.1)

      if (missing == "ind") {
        missing.ind <- apply(covs[, apply(covs, 2, anyNA), drop = FALSE], 2, function(x) as.numeric(is.na(x)))
        covs[is.na(covs)] <- 0
        covs <- cbind(covs, missing.ind)
      }
      return(covs)
    })
  }
  if (is_not_null(treat.list)) {
    treat.list <- lapply(treat.list, function(t) {
      treat <- t[subset]
      if (get.treat.type(t) != "continuous") treat <- factor(treat)
      return(treat)
    })
  }
  if (is_not_null(s.weights)) {
    s.weights <- s.weights[subset]
  }

  baseline.data <- data.frame(treat.list[[1]], covs.list[[1]])
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
weightit2twang <- function(covs, treat, s.weights, estimand, focal, subset, stabilize, subclass, missing, ...) {
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
  if (is.na(s.m.matches) || s.m.matches == 0L) {stop(paste0("stop.method must be one of ", word_list(available.stop.methods, "or", quotes = TRUE), "."), call. = FALSE)}
  else A[["stop.method"]] <- available.stop.methods[s.m.matches]

  for (f in names(formals(twang::ps))) {
    if (is_null(A[[f]])) A[[f]] <- formals(twang::ps)[[f]]
  }

  covs <- covs[subset, , drop = FALSE]
  treat <- factor(treat[subset])

  covs <- apply(covs, 2, make.closer.to.1)

  if (missing == "ind") {
    missing.ind <- apply(covs[, apply(covs, 2, anyNA), drop = FALSE], 2, function(x) as.numeric(is.na(x)))
    colnames(missing.ind) <- paste0(colnames(missing.ind), ":<NA>")
    # covs[is.na(covs)] <- 0
    covs <- cbind(covs, missing.ind)
  }

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

        w[treat == i] <- get.w.ps(fit.list[[i]], stop.method = s)[treat_ == 0]

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
          w <- get.w.ps(fit.list[[i]], stop.method = s)
          ps <- fit.list[[i]][["ps"]][[1]]
          fit.list <- fit.list[[i]]
          break
        }
        else {
          w[treat == i] <- get.w.ps(fit.list[[i]], stop.method = s)[treat == i]
        }
      }
    }

    if (stabilize) {
      w <- stabilize_w(w, treat)
    }
  }
  obj <- list(w = w, ps = ps, fit.obj = fit.list)
  return(obj)
}
weightit2twang.cont <- function(covs, treat, s.weights, subset, stabilize, missing, ...) {
  A <- list(...)
  A[c("formula", "data", "sampw", "verbose")] <- NULL
  if (is_null(A[["stop.method"]])) {
    warning("No stop.method was entered. Using \"p.mean\", the mean of the absolute Pearson correlations.",
            call. = FALSE, immediate. = TRUE)
    A[["stop.method"]] <- "p.mean"
  }
  else if (length(A$stop.method) > 1) {
    warning("Only one stop.method is allowed at a time. Using just the first stop.method.",
            call. = FALSE, immediate. = TRUE)
    A$stop.method <- A$stop.method[1]
  }

  covs <- covs[subset, , drop = FALSE]
  treat <- treat[subset]

  covs <- apply(covs, 2, make.closer.to.1)

  if (missing == "ind") {
    missing.ind <- apply(covs[, apply(covs, 2, anyNA), drop = FALSE], 2, function(x) as.numeric(is.na(x)))
    colnames(missing.ind) <- paste0(colnames(missing.ind), ":<NA>")
    # covs[is.na(covs)] <- 0
    covs <- cbind(covs, missing.ind)
  }

  new.data <- data.frame(treat, covs)

  if (check.package("gbm")) {
    fit <- do.call(ps.cont, c(list(formula(new.data),
                                     data = new.data,
                                     sampw = s.weights[subset],
                                     verbose = TRUE), A))
    w <- get.w.ps.cont(fit, stop.method = A[["stop.method"]])
  }

  #ps <- fit[["ps"]][[A[["stop.method"]]]]

  obj <- list(w = w, fit.obj = fit)
  return(obj)
}

#Generalized boosted modeling with gbm and cobalt
weightit2gbm <- function(covs, treat, s.weights, estimand, focal, subset, stabilize, subclass, missing, ...) {

  check.package(c("gbm", "cobalt"))

  A <- list(...)

  covs <- covs[subset, , drop = FALSE]
  treat <- treat[subset]
  s.weights <- s.weights[subset]

  if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
  treat.type <- get.treat.type(treat)

  covs <- apply(covs, 2, make.closer.to.1)

  if (missing == "ind") {
    missing.ind <- apply(covs[, apply(covs, 2, anyNA), drop = FALSE], 2, function(x) as.numeric(is.na(x)))
    colnames(missing.ind) <- paste0(colnames(missing.ind), ":<NA>")
    covs <- cbind(covs, missing.ind)
  }

  bin.vars <- apply(covs, 2, is_binary)

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
  available.stop.methods <- c("ks.mean", "es.mean", "ks.max", "es.max", "ks.rms", "es.rms")
  s.m.matches <- charmatch(A[["stop.method"]], available.stop.methods)
  if (is.na(s.m.matches) || s.m.matches == 0L) {
    if (startsWith(A[["stop.method"]], "cv") && can_str2num(numcv <- substr(A[["stop.method"]], 3, nchar(A[["stop.method"]])))) {
      cv <- round(str2num(numcv))
      if (cv < 2) stop("At least 2 CV-folds must be specified in stop.method.", call. = FALSE)
    }
    else stop(paste0("stop.method must be one of ", word_list(c(available.stop.methods, "cv{#}"), "or", quotes = TRUE), "."), call. = FALSE)
  }
  else stop.method <- available.stop.methods[s.m.matches]

  if (cv == 0) {
    if (startsWith(stop.method, "es.")) stop.fun <- function(mat, treat, weights, s.d.denom, s.weights = NULL, bin.vars, subset = NULL) {
      cobalt::col_w_smd(mat, treat, weights, std = rep(TRUE, ncol(mat)), s.d.denom = s.d.denom, abs = TRUE,
                        s.weights = s.weights, bin.vars = bin.vars, subset = subset)
    }
    else if (startsWith(stop.method, "ks.")) stop.fun <- function(mat, treat, weights, s.d.denom, s.weights = NULL, bin.vars, subset = NULL) {
      cobalt::col_w_ks(mat, treat, weights, s.weights = s.weights, bin.vars = bin.vars, subset = subset)
    }

    if (endsWith(stop.method, ".mean")) stop.sum <- mean
    else if (endsWith(stop.method, ".max")) stop.sum <- max
    else if (endsWith(stop.method, ".rms")) stop.sum <- function(x, ...) sqrt(mean(x^2, ...))
  }

  if (is_null(A[["trim.at"]])) trim.at <- 0
  else trim.at <- A[["trim.at"]]

  for (f in names(formals(gbm::gbm.fit))) {
    if (is_null(A[[f]])) {
      if (f %in% c("x", "y", "misc", "w", "verbose", "var.names",
                   "response.name", "group", "distribution")) A[f] <- list(NULL)
      else A[f] <- list(switch(f, n.trees = 1e4,
                               interaction.depth = 3,
                               shrinkage = .01,
                               bag.fraction = 1,
                               keep.data = FALSE,
                               formals(gbm::gbm.fit)[[f]]))
    }
  }

  if (treat.type == "binary")  available.distributions <- c("bernoulli", "adaboost")
  else available.distributions <- "multinomial"

  if (cv == 0) {
    start.tree <- if_null_then(A[["start.tree"]], 1)
    if (is_null(A[["n.grid"]])) n.grid <- round(1+sqrt(2*(A[["n.trees"]]-start.tree+1)))
    else if (!is_(A[["n.grid"]], "numeric") || length(A[["n.grid"]]) > 1 ||
             !between(A[["n.grid"]], c(2, A[["n.trees"]]))) {
      stop("n.grid must be a numeric value between 2 and n.trees.", call. = FALSE)
    }
    else n.grid <- round(A[["n.grid"]])
  }

  A[["x"]] <- covs
  A[["y"]] <- treat
  A[["distribution"]] <- if (is_null(A[["distribution"]])) {
    available.distributions[1]} else {
      match_arg(A[["distribution"]], available.distributions)}
  A[["w"]] <- s.weights
  A[["verbose"]] <- FALSE

  if (treat.type == "binary") {
    fit <- do.call(gbm::gbm.fit, A[names(A) %in% names(formals(gbm::gbm.fit))])
    if (cv == 0) {
      n.trees <- fit[["n.trees"]]
      iters <- 1:n.trees
      iters.grid <- round(seq(start.tree, n.trees, length.out = n.grid))

      if (anyNA(iters.grid) || is_null(iters.grid) || any(iters.grid > n.trees)) stop("A problem has occurred")

      ps <- gbm::predict.gbm(fit, n.trees = iters.grid, type = "response", newdata = covs)
      w <- apply(ps, 2, get_w_from_ps, treat = treat, estimand = estimand, focal = focal, stabilize = stabilize, subclass = subclass)
      w <- suppressMessages(apply(w, 2, trim, at = trim.at, treat = treat))

      s.d.denom <- switch(estimand, ATT = "treated", ATC = "control", ATE = "all", "weighted")
      iter.grid.balance <- apply(w, 2, function(w_) {
        stop.sum(stop.fun(covs, treat, weights = w_, s.d.denom = s.d.denom,
                          s.weights = s.weights, bin.vars = bin.vars))
      })

      it <- which.min(iter.grid.balance) + c(-1, 1)
      it[1] <- iters.grid[max(1, it[1])]
      it[2] <- iters.grid[min(length(iters.grid), it[2])]
      iters.to.check <- iters[between(iters, iters[it])]

      if (anyNA(iters.to.check) || is_null(iters.to.check) || any(iters.to.check > n.trees)) stop("A problem has occurred")

      ps <- gbm::predict.gbm(fit, n.trees = iters.to.check, type = "response", newdata = covs)
      w <- apply(ps, 2, get_w_from_ps, treat = treat, estimand = estimand, focal = focal, stabilize = stabilize, subclass = subclass)
      w <- suppressMessages(apply(w, 2, trim, at = trim.at, treat = treat))
      iter.grid.balance.fine <- apply(w, 2, function(w_) {
        stop.sum(stop.fun(covs, treat, weights = w_, s.d.denom = s.d.denom,
                          s.weights = s.weights, bin.vars = bin.vars))
      })

      best.tree <- iters.to.check[which.min(iter.grid.balance.fine)]
      tree.val <- setNames(data.frame(c(iters.grid, iters.to.check),
                                      c(iter.grid.balance, iter.grid.balance.fine)),
                           c("tree", stop.method))
      tree.val <- unique(tree.val[order(tree.val$tree),])
      w <- w[,as.character(best.tree)]
      ps <- ps[,as.character(best.tree)]
    }
    else {
      A["data"] <- list(data.frame(treat, covs))
      A[["cv.folds"]] <- cv
      A["n.cores"] <- list(A[["n.cores"]])
      A["var.names"] <- list(A[["var.names"]])
      A["offset"] <- list(NULL)
      A[["nTrain"]] <- floor(nrow(covs))
      A[["class.stratify.cv"]] <- FALSE
      A[["y"]] <- treat
      A[["x"]] <- covs
      A[["distribution"]] <- list(name = A[["distribution"]])
      A[["w"]] <- s.weights
      cv.results <- do.call(gbm::gbmCrossVal,
                            A[names(A) %in% names(formals(gbm::gbmCrossVal))])
      best.tree <- which.min(cv.results$error)
      tree.val <- data.frame(tree = seq_along(cv.results$error),
                             error = cv.results$error)
      ps <- gbm::predict.gbm(fit, n.trees = best.tree, type = "response", newdata = covs)
      w <- get_w_from_ps(ps, treat = treat, estimand = estimand, focal = focal, stabilize = stabilize, subclass = subclass)
      w <- suppressMessages(trim(w, at = trim.at, treat = treat))
    }
  }
  else if (treat.type == "multinomial") {

    treat <- factor(treat)

    fit <- do.call(gbm::gbm.fit, A[names(A) %in% names(formals(gbm::gbm.fit))])

    if (cv == 0) {
      n.trees <- fit[["n.trees"]]
      iters <- 1:n.trees
      iters.grid <- round(seq(start.tree, n.trees, length.out = n.grid))

      if (anyNA(iters.grid) || is_null(iters.grid) || any(iters.grid > n.trees)) stop("A problem has occurred")

      ps <- gbm::predict.gbm(fit, n.trees = iters.grid, type = "response", newdata = covs)
      w <- apply(ps, 3, get_w_from_ps, treat = treat, estimand = estimand, focal = focal, stabilize = stabilize, subclass = subclass)
      w <- suppressMessages(apply(w, 2, trim, at = trim.at, treat = treat))

      iter.grid.balance <- {
        if (is_not_null(focal)) {
          bin.treat <- as.numeric(treat == focal)
          apply(w, 2, function(w_) {
            bal <- unlist(lapply(levels(treat)[levels(treat) != focal], function(t) {
              stop.fun(covs, bin.treat, weights = w_, s.d.denom = "treated", s.weights = s.weights,
                       bin.vars = bin.vars, subset = treat %in% c(t, focal))
            }))
            stop.sum(bal)
          })
        }
        else {
          covs_i_list <- setNames(lapply(levels(treat), function(i) {
            rbind(covs, covs[treat == i, , drop = FALSE])
          }), levels(treat))
          treat_i_list <- setNames(lapply(levels(treat), function(i) {
            c(rep(1, nrow(covs)), rep(0, sum(treat == i)))
          }), levels(treat))
          s.weights_i_list <- setNames(lapply(levels(treat), function(i) {
            if (is_not_null(s.weights)) c(s.weights, s.weights[treat == i])
            else NULL
          }), levels(treat))

          s.d.denom <- NULL
          if (estimand == "ATE" && startsWith(stop.method, "es.")) {
            s.d.denom <- cobalt::col_w_sd(covs, weights = NULL,
                                          s.weights = s.weights,
                                          bin.vars = bin.vars)
          }

          apply(w, 2, function(w_) {
            bal <- unlist(lapply(levels(treat), function(i) {
              if (estimand == "ATE") {
                w_i <- c(rep(1, nrow(covs)), w_[treat == i])
              }
              else {
                if (startsWith(stop.method, "es.")) {
                  s.d.denom <- cobalt::col_w_sd(covs, weights = w_,
                                                s.weights = s.weights,
                                                bin.vars = bin.vars)
                }
                w_i <- c(w_, w_[treat == i])
              }

              stop.fun(covs_i_list[[i]], treat_i_list[[i]], weights = w_i,
                       s.d.denom = s.d.denom, s.weights = s.weights_i_list[[i]],
                       bin.vars = bin.vars)
            }))
            stop.sum(bal)
          })
        }
      }

      it <- which.min(iter.grid.balance) + c(-1, 1)
      it[1] <- iters.grid[max(1, it[1])]
      it[2] <- iters.grid[min(length(iters.grid), it[2])]
      iters.to.check <- iters[between(iters, iters[it])]

      if (anyNA(iters.to.check) || is_null(iters.to.check) || any(iters.to.check > n.trees)) stop("A problem has occurred")

      ps <- gbm::predict.gbm(fit, n.trees = iters.to.check, type = "response", newdata = covs)
      w <- apply(ps, 3, get_w_from_ps, treat = treat, estimand = estimand, focal = focal, stabilize = stabilize, subclass = subclass)
      w <- suppressMessages(apply(w, 2, trim, at = trim.at, treat = treat))

      iter.grid.balance.fine <- {
        if (is_not_null(focal)) {
          bin.treat <- as.numeric(treat == focal)
          apply(w, 2, function(w_) {
            bal <- unlist(lapply(levels(treat)[levels(treat) != focal], function(t) {
              stop.fun(covs, bin.treat, weights = w_, s.d.denom = "treated", s.weights = s.weights,
                       bin.vars = bin.vars, subset = treat %in% c(t, focal))
            }))
            stop.sum(bal)
          })
        }
        else {
          covs_i_list <- setNames(lapply(levels(treat), function(i) {
            rbind(covs, covs[treat == i, , drop = FALSE])
          }), levels(treat))
          treat_i_list <- setNames(lapply(levels(treat), function(i) {
            c(rep(1, nrow(covs)), rep(0, sum(treat == i)))
          }), levels(treat))
          s.weights_i_list <- setNames(lapply(levels(treat), function(i) {
            if (is_not_null(s.weights)) c(s.weights, s.weights[treat == i])
            else NULL
          }), levels(treat))

          s.d.denom <- NULL
          if (estimand == "ATE" && startsWith(stop.method, "es.")) {
            s.d.denom <- cobalt::col_w_sd(covs, weights = NULL,
                                          s.weights = s.weights,
                                          bin.vars = bin.vars)
          }

          apply(w, 2, function(w_) {
            bal <- unlist(lapply(levels(treat), function(i) {
              if (estimand == "ATE") {
                w_i <- c(rep(1, nrow(covs)), w_[treat == i])
              }
              else {
                if (startsWith(stop.method, "es.")) {
                  s.d.denom <- cobalt::col_w_sd(covs, weights = w_,
                                                s.weights = s.weights,
                                                bin.vars = bin.vars)
                }
                w_i <- c(w_, w_[treat == i])
              }

              stop.fun(covs_i_list[[i]], treat_i_list[[i]], weights = w_i,
                       s.d.denom = s.d.denom, s.weights = s.weights_i_list[[i]],
                       bin.vars = bin.vars)
            }))
            stop.sum(bal)
          })
        }
      }

      best.tree <- iters.to.check[which.min(iter.grid.balance.fine)]
      tree.val <- setNames(data.frame(c(iters.grid, iters.to.check),
                                      c(iter.grid.balance, iter.grid.balance.fine)),
                           c("tree", stop.method))
      tree.val <- unique(tree.val[order(tree.val$tree),])

      ps <- ps[, , as.character(best.tree)]
      w <- get_w_from_ps(ps, treat, estimand, focal, stabilize = stabilize, subclass = subclass)
      ps <- NULL
    }
    else {
      A["data"] <- list(data.frame(treat, covs))
      A[["cv.folds"]] <- cv
      A["n.cores"] <- list(A[["n.cores"]])
      A["var.names"] <- list(A[["var.names"]])
      A["offset"] <- list(NULL)
      A[["nTrain"]] <- floor(nrow(covs))
      A[["class.stratify.cv"]] <- FALSE
      A[["y"]] <- treat
      A[["x"]] <- covs
      A[["distribution"]] <- list(name = A[["distribution"]])
      A[["w"]] <- s.weights
      cv.results <- do.call(gbm::gbmCrossVal,
                            A[names(A) %in% names(formals(gbm::gbmCrossVal))])
      best.tree <- which.min(cv.results$error)
      tree.val <- data.frame(tree = seq_along(cv.results$error),
                             error = cv.results$error)
      ps <- gbm::predict.gbm(fit, n.trees = best.tree, type = "response", newdata = covs)
      w <- get_w_from_ps(ps[, , 1], treat = treat, estimand = estimand, focal = focal, stabilize = stabilize, subclass = subclass)
      w <- suppressMessages(trim(w, at = trim.at, treat = treat))
      ps <- NULL
    }
  }

  obj <- list(w = w, ps = ps, info = list(best.tree = best.tree,
                                          tree.val = tree.val), fit.obj = fit)
  return(obj)
}
weightit2gbm.cont <- function(covs, treat, s.weights, subset, stabilize, missing, ...) {

  check.package(c("gbm", "cobalt"))

  A <- list(...)

  covs <- covs[subset, , drop = FALSE]
  treat <- treat[subset]

  covs <- apply(covs, 2, make.closer.to.1)

  if (missing == "ind") {
    missing.ind <- apply(covs[, apply(covs, 2, anyNA), drop = FALSE], 2, function(x) as.numeric(is.na(x)))
    colnames(missing.ind) <- paste0(colnames(missing.ind), ":<NA>")
    # covs[is.na(covs)] <- 0
    covs <- cbind(covs, missing.ind)
  }

  bin.vars <- apply(covs, 2, is_binary)

  if (is_null(A[["stop.method"]])) {
    warning("No stop.method was provided. Using \"p.mean\".",
            call. = FALSE, immediate. = TRUE)
    A[["stop.method"]] <- "p.mean"
  }
  else if (length(A[["stop.method"]]) > 1) {
    warning("Only one stop.method is allowed at a time. Using just the first stop.method.",
            call. = FALSE, immediate. = TRUE)
    A[["stop.method"]] <- A[["stop.method"]][1]
  }

  cv <- 0
  available.stop.methods <- c("p.mean", "s.mean", "p.max", "s.max", "p.rms", "s.rms")
  s.m.matches <- charmatch(A[["stop.method"]], available.stop.methods)
  if (is.na(s.m.matches) || s.m.matches == 0L) {
    if (startsWith(A[["stop.method"]], "cv") && can_str2num(numcv <- substr(A[["stop.method"]], 3, nchar(A[["stop.method"]])))) {
      cv <- round(str2num(numcv))
      if (cv < 2) stop("At least 2 CV-folds must be specified in stop.method.", call. = FALSE)
    }
    else stop(paste0("stop.method must be one of ", word_list(c(available.stop.methods, "cv{#}"), "or", quotes = TRUE), "."), call. = FALSE)
  }
  else stop.method <- available.stop.methods[s.m.matches]

  if (cv == 0) {
    if (startsWith(stop.method, "s.")) stop.fun <- function(mat, treat, weights, s.weights, bin.vars, subset = NULL) {
      cobalt::col_w_corr(mat, treat, weights, type = "spearman", abs = TRUE,
                         s.weights = s.weights, bin.vars = bin.vars, subset = subset)
    }
    else if (startsWith(stop.method, "p.")) stop.fun <- function(mat, treat, weights, s.weights, bin.vars, subset = NULL) {
      cobalt::col_w_corr(mat, treat, weights, type = "pearson", abs = TRUE,
                         s.weights = s.weights, bin.vars = bin.vars, subset = subset)
    }

    if (endsWith(stop.method, ".mean")) stop.sum <- mean
    else if (endsWith(stop.method, ".max")) stop.sum <- max
    else if (endsWith(stop.method, ".rms")) stop.sum <- function(x, ...) sqrt(mean(x^2, ...))
  }

  if (is_null(trim.at <- A[["trim.at"]])) trim.at <- 0

  for (f in names(formals(gbm::gbm.fit))) {
    if (is_null(A[[f]])) {
      if (f %in% c("x", "y", "misc", "w", "verbose", "var.names",
                   "response.name", "group", "distribution")) A[f] <- list(NULL)
      else A[f] <- list(switch(f, n.trees = 2e4,
                               interaction.depth = 4,
                               shrinkage = 0.0005,
                               bag.fraction = 1,
                               formals(gbm::gbm.fit)[[f]]))
    }
  }

  available.distributions <- c("gaussian", "laplace", "tdist", "poisson")

  if (cv == 0) {
    start.tree <- if_null_then(A[["start.tree"]], 1)
    if (is_null(A[["n.grid"]])) n.grid <- round(1+sqrt(2*(A[["n.trees"]]-start.tree+1)))
    else if (!is_(A[["n.grid"]], "numeric") || length(A[["n.grid"]]) > 1 ||
             !between(A[["n.grid"]], c(2, A[["n.trees"]]))) {
      stop("n.grid must be a numeric value between 2 and n.trees.", call. = FALSE)
    }
    else n.grid <- round(A[["n.grid"]])
  }

  A[["x"]] <- covs
  A[["y"]] <- treat
  A[["distribution"]] <- if (is_null(A[["distribution"]])) {
    available.distributions[1]} else {
      match_arg(A[["distribution"]], available.distributions)}
  A[["w"]] <- s.weights
  A[["verbose"]] <- FALSE

  if (isTRUE(A[["use.kernel"]])) {
    if (is_null(A[["bw"]])) A[["bw"]] <- "nrd0"
    if (is_null(A[["adjust"]])) A[["adjust"]] <- 1
    if (is_null(A[["kernel"]])) A[["kernel"]] <- "gaussian"
    if (is_null(A[["n"]])) A[["n"]] <- 10*length(treat)
  }
  else {
    if (is_null(A[["density"]])) densfun <- dnorm
    else if (is.function(A[["density"]])) densfun <- A[["density"]]
    else if (is.character(A[["density"]]) && length(A[["density"]] == 1)) {
      splitdens <- strsplit(A[["density"]], "_", fixed = TRUE)[[1]]
      if (exists(splitdens[1], mode = "function", envir = parent.frame())) {
        if (length(splitdens) > 1 && !can_str2num(splitdens[-1])) {
          stop(paste(A[["density"]], "is not an appropriate argument to density because",
                     word_list(splitdens[-1], and.or = "or", quotes = TRUE), "cannot be coerced to numeric."), call. = FALSE)
        }
        densfun <- function(x) {
          tryCatch(do.call(get(splitdens[1]), c(list(x), as.list(str2num(splitdens[-1])))),
                   error = function(e) stop(paste0("Error in applying density:\n  ", conditionMessage(e)), call. = FALSE))
        }
      }
      else {
        stop(paste(A[["density"]], "is not an appropriate argument to density because",
                   splitdens[1], "is not an available function."), call. = FALSE)
      }
    }
    else stop("The argument to density cannot be evaluated as a density function.", call. = FALSE)
  }

  #Stabilization
  p.num <- treat - mean(treat)

  if (isTRUE(A[["use.kernel"]])) {
    d.n <- density(p.num, n = A[["n"]],
                   weights = s.weights/sum(s.weights), give.Rkern = FALSE,
                   bw = A[["bw"]], adjust = A[["adjust"]], kernel = A[["kernel"]])
    dens.num <- with(d.n, approxfun(x = x, y = y))(p.num)
  }
  else {
    dens.num <- densfun(p.num/sd(treat))
    if (is_null(dens.num) || !is.atomic(dens.num) || anyNA(dens.num)) {
      stop("There was a problem with the output of density. Try another density function or leave it blank to use the normal density.", call. = FALSE)
    }
    else if (any(dens.num <= 0)) {
      stop("The input to density may not accept the full range of treatment values.", call. = FALSE)
    }
  }

  get_cont_weights <- function(ps, treat, s.weights, dens.num) {
    p.denom <- treat - ps

    if (isTRUE(A[["use.kernel"]])) {
      d.d <- density(p.denom, n = A[["n"]],
                     weights = s.weights/sum(s.weights), give.Rkern = FALSE,
                     bw = A[["bw"]], adjust = A[["adjust"]], kernel = A[["kernel"]])
      dens.denom <- with(d.d, approxfun(x = x, y = y))(p.denom)
    }
    else {
      dens.denom <- densfun(p.denom/sd(p.denom))
      if (is_null(dens.denom) || !is.atomic(dens.denom) || anyNA(dens.denom)) {
        stop("There was a problem with the output of density. Try another density function or leave it blank to use the normal density.", call. = FALSE)
      }
      else if (any(dens.denom <= 0)) {
        stop("The input to density may not accept the full range of standardized residuals of the treatment model.", call. = FALSE)
      }

    }

    w <- dens.num/dens.denom

    return(w)
  }

  fit <- do.call(gbm::gbm.fit, A[names(A) %in% names(formals(gbm::gbm.fit))])

  if (cv == 0) {
    n.trees <- fit[["n.trees"]]
    iters <- 1:n.trees
    iters.grid <- round(seq(start.tree, n.trees, length.out = n.grid))

    if (anyNA(iters.grid) || is_null(iters.grid) || any(iters.grid > n.trees)) stop("A problem has occurred")

    ps <- gbm::predict.gbm(fit, n.trees = iters.grid, newdata = covs)

    w <- apply(ps, 2, get_cont_weights, treat = treat, s.weights = s.weights, dens.num = dens.num)
    w <- suppressMessages(apply(w, 2, trim, at = trim.at, treat = treat))

    iter.grid.balance <- apply(w, 2, function(w_) stop.sum(stop.fun(covs, treat, weights = w_, s.weights, bin.vars)))

    it <- which.min(iter.grid.balance) + c(-1, 1)
    it[1] <- iters.grid[max(1, it[1])]
    it[2] <- iters.grid[min(length(iters.grid), it[2])]
    iters.to.check <- iters[between(iters, iters[it])]

    if (anyNA(iters.to.check) || is_null(iters.to.check) || any(iters.to.check > n.trees)) stop("A problem has occurred")

    ps <- gbm::predict.gbm(fit, n.trees = iters.to.check, newdata = covs)
    w <- apply(ps, 2, get_cont_weights, treat = treat, s.weights = s.weights, dens.num = dens.num)
    w <- suppressMessages(apply(w, 2, trim, at = trim.at, treat = treat))

    iter.grid.balance.fine <- apply(w, 2, function(w_) stop.sum(stop.fun(covs, treat, weights = w_, s.weights, bin.vars)))

    best.tree <- iters.to.check[which.min(iter.grid.balance.fine)]
    tree.val <- setNames(data.frame(c(iters.grid, iters.to.check),
                                    c(iter.grid.balance, iter.grid.balance.fine)),
                         c("tree", stop.method))
    tree.val <- unique(tree.val[order(tree.val$tree),])
    w <- w[,as.character(best.tree)]
  }
  else {
    A["data"] <- list(data.frame(treat, covs))
    A[["cv.folds"]] <- cv
    A["n.cores"] <- list(A[["n.cores"]])
    A["var.names"] <- list(A[["var.names"]])
    A["offset"] <- list(NULL)
    A[["nTrain"]] <- floor(nrow(covs))
    A[["class.stratify.cv"]] <- FALSE
    A[["distribution"]] <- list(name = A[["distribution"]])

    cv.results <- do.call(gbm::gbmCrossVal,
                          A[names(A) %in% names(formals(gbm::gbmCrossVal))])
    best.tree <- which.min(cv.results$error)
    tree.val <- data.frame(tree = seq_along(cv.results$error),
                           error = cv.results$error)
    ps <- gbm::predict.gbm(fit, n.trees = best.tree, newdata = covs)
    w <- get_cont_weights(ps, treat = treat, s.weights = s.weights, dens.num = dens.num)
    w <- suppressMessages(trim(w, at = trim.at, treat = treat))
  }
  obj <- list(w = w, info = list(best.tree = best.tree,
                                 tree.val = tree.val), fit.obj = fit)
  return(obj)
}

#CBPS
weightit2cbps <- function(covs, treat, s.weights, estimand, focal, subset, stabilize, subclass, missing, ...) {
  A <- list(...)

  covs <- covs[subset, , drop = FALSE]
  treat <- factor(treat[subset])

  if (missing == "ind") {
    missing.ind <- apply(covs[, apply(covs, 2, anyNA), drop = FALSE], 2, function(x) as.numeric(is.na(x)))
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


        if (length(control.levels) == 1) {
          ps <- fit.list[[i]][["fitted.values"]]
          fit.list <- fit.list[[1]]
          w <- get_w_from_ps(ps, treat_, estimand = "ATT",
                             treated = 1, subclass = subclass, stabilize = stabilize)
        }
        else {
          w[treat == i] <- get.w(fit.list[[i]], estimand = "ATT")[treat_ == 0] / s.weights[subset][treat.in.i.focal][treat_ == 0]
        }
      }
    }
    else {
      new.data <- data.frame(treat, covs)
      if (!nunique.gt(treat, 4)) {
        tryCatch({fit.list <- CBPS::CBPS(formula(new.data),
                                         data = new.data,
                                         method = if (isFALSE(A$over)) "exact" else "over",
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

        if (is_binary(treat)) {
          ps <- fit.list[["fitted.values"]]
          fit.list <- fit.list[[1]]
          w <- get_w_from_ps(ps, treat, estimand = "ATE",
                             subclass = subclass, stabilize = stabilize)
        }
        else {
          w <- get.w(fit.list, estimand = "ATE") / s.weights[subset]
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

          w[treat==i] <- get.w(fit.list[[i]], estimand = "ATE")[treat==i] / s.weights[subset][treat==i]
        }
      }
    }

  }

  if (stabilize && !is_binary(treat)) {
    w <- stabilize_w(w, treat)
  }

  obj <- list(w = w, ps = ps, fit.obj = fit.list)
  return(obj)
}
weightit2cbps.cont <- function(covs, treat, s.weights, subset, missing, ...) {
  A <- list(...)

  covs <- covs[subset, , drop = FALSE]
  treat <- treat[subset]

  if (missing == "ind") {
    missing.ind <- apply(covs[, apply(covs, 2, anyNA), drop = FALSE], 2, function(x) as.numeric(is.na(x)))
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
  w <- get.w(fit) #/ s.weights[subset]

  obj <- list(w = w, fit.obj = fit)
  return(obj)
}
weightit2cbps.msm <- function(covs.list, treat.list, s.weights, subset, missing, ...) {
  stop("CBMSM doesn't work yet.")
}
weightit2npcbps <- function(covs, treat, s.weights, subset, moments, int, missing, ...) {
  A <- list(...)
  if (!all_the_same(s.weights)) stop(paste0("Sampling weights cannot be used with method = \"npcbps\"."),
                                     call. = FALSE)

  covs <- covs[subset, , drop = FALSE]
  treat <- factor(treat[subset])

  if (missing == "ind") {
    missing.ind <- apply(covs[, apply(covs, 2, anyNA), drop = FALSE], 2, function(x) as.numeric(is.na(x)))
    covs[is.na(covs)] <- 0
    covs <- cbind(covs, missing.ind)
  }

  covs <- cbind(covs, int.poly.f(covs, poly = moments, int = int))
  covs <- apply(covs, 2, make.closer.to.1)
  colinear.covs.to.remove <- colnames(covs)[colnames(covs) %nin% colnames(make_full_rank(covs))]
  covs <- covs[, colnames(covs) %nin% colinear.covs.to.remove, drop = FALSE]

  new.data <- data.frame(treat = treat, covs)

  if (check.package("CBPS")) {
    fit <- do.call(CBPS::npCBPS, c(list(formula(new.data), data = new.data, print.level = 1), A),
                   quote = TRUE)
  }
  w <- get.w(fit)

  obj <- list(w = w, fit.obj = fit)

  return(obj)
}
weightit2npcbps.cont <- function(covs, treat, s.weights, subset, estimand, moments, int, missing, ...) {
  A <- list(...)

  if (!all_the_same(s.weights)) stop(paste0("Sampling weights cannot be used with method = \"npcbps\"."),
                                     call. = FALSE)

  covs <- covs[subset, , drop = FALSE]
  treat <- treat[subset]

  if (missing == "ind") {
    missing.ind <- apply(covs[, apply(covs, 2, anyNA), drop = FALSE], 2, function(x) as.numeric(is.na(x)))
    covs[is.na(covs)] <- 0
    covs <- cbind(covs, missing.ind)
  }

  covs <- cbind(covs, int.poly.f(covs, poly = moments, int = int))
  covs <- apply(covs, 2, make.closer.to.1)
  colinear.covs.to.remove <- colnames(covs)[colnames(covs) %nin% colnames(make_full_rank(covs))]
  covs <- covs[, colnames(covs) %nin% colinear.covs.to.remove, drop = FALSE]

  new.data <- data.frame(treat = treat, covs)

  if (check.package("CBPS")) {
    fit <- do.call(CBPS::npCBPS, c(list(formula(new.data), data = new.data, print.level = 1), A),
                   quote = TRUE)
  }
  w <- get.w(fit)

  obj <- list(w = w, fit.obj = fit)

  return(obj)
}

#Entropy balancing with ebal
weightit2ebal <- function(covs, treat, s.weights, subset, estimand, focal, stabilize, missing, moments, int, ...) {
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

        colinear.covs.to.remove <- colnames(covs_)[colnames(covs_) %nin% colnames(make_full_rank(covs_[treat_ == 0, , drop = FALSE]))]

        covs_ <- covs_[, colnames(covs_) %nin% colinear.covs.to.remove, drop = FALSE]

        if (is_not_null(A[["base.weight"]])) {
          A[["base.weight"]] <- A[["base.weight"]][treat == i]
        }

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

        colinear.covs.to.remove <- colnames(covs_i)[colnames(covs_i) %nin% colnames(make_full_rank(covs_i[treat_i == 0, , drop = FALSE]))]

        covs_i <- covs_i[, colnames(covs_i) %nin% colinear.covs.to.remove, drop = FALSE]

        covs_i[treat_i == 1,] <- covs_i[treat_i == 1,] * s.weights[subset] * sum(treat_i == 1) / sum(s.weights[subset])

        if (is_not_null(A[["base.weight"]])) {
          A[["base.weight"]] <- A[["base.weight"]][treat == i]
        }

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
weightit2ebcw <- function(covs, treat, s.weights, subset, estimand, focal, missing, moments, int, ...) {
  A <- list(...)

  covs <- covs[subset, , drop = FALSE]
  treat <- factor(treat[subset])

  if (missing == "ind") {
    missing.ind <- apply(covs[, apply(covs, 2, anyNA), drop = FALSE], 2, function(x) as.numeric(is.na(x)))
    covs[is.na(covs)] <- 0
    covs <- cbind(covs, missing.ind)
  }

  covs <- cbind(covs, int.poly.f(covs, poly = moments, int = int))
  covs <- apply(covs, 2, make.closer.to.1)

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

        colinear.covs.to.remove <- colnames(covs_)[colnames(covs_) %nin% colnames(make_full_rank(covs_[treat_ == 0, , drop = FALSE]))]

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

        colinear.covs.to.remove <- colnames(covs_i)[colnames(covs_i) %nin% colnames(make_full_rank(covs_i[treat_i == 0, , drop = FALSE]))]

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
    if (length(fit.list) == 1) fit.list <- fit.list[[1]]
    obj <- list(w = w, fit.obj = fit.list)
    return(obj)
  }
}

#PS weights using SuperLearner
weightit2super <- function(covs, treat, s.weights, subset, estimand, focal, stabilize, subclass, missing, ...) {
  A <- list(...)

  check.package("SuperLearner")

  covs <- covs[subset, , drop = FALSE]
  treat <- factor(treat[subset])

  if (missing == "ind") {
    missing.ind <- apply(covs[, apply(covs, 2, anyNA), drop = FALSE], 2, function(x) as.numeric(is.na(x)))
    covs[is.na(covs)] <- 0
    covs <- cbind(covs, missing.ind)
  }
  covs <- data.frame(apply(covs, 2, make.closer.to.1))

  if (ncol(covs) > 1) {
    colinear.covs.to.remove <- colnames(covs)[colnames(covs) %nin% colnames(make_full_rank(covs))]
    covs <- covs[, colnames(covs) %nin% colinear.covs.to.remove, drop = FALSE]
  }

  for (f in names(formals(SuperLearner::SuperLearner))) {
    if (f == "method") {if (is_null(A[["SL.method"]])) A[["SL.method"]] <- formals(SuperLearner::SuperLearner)[["method"]]}
    else if (f == "env") {if (is_null(A[["env"]])) A[["env"]] <- environment(SuperLearner::SuperLearner)}
    else if (is_null(A[[f]])) A[[f]] <- formals(SuperLearner::SuperLearner)[[f]]
  }

  discrete <- if_null_then(A[["discrete"]], FALSE)
  if (length(discrete) != 1 || !is_(discrete, "logical")) stop("discrete must be TRUE or FALSE.", call. = FALSE)

  ps <- setNames(as.data.frame(matrix(NA_real_, ncol = nunique(treat), nrow = length(treat))),
                 levels(treat))

  if (is_binary(treat)) {

    fit.list <- do.call(SuperLearner::SuperLearner, list(Y = binarize(treat),
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

    if (discrete) p.score <- fit.list$library.predict[,which.min(fit.list$cvRisk)]
    else p.score <- fit.list$SL.predict

    ps[[2]] <- p.score
    ps[[1]] <- 1 - ps[[2]]

    info <- list(coef = fit.list$coef,
                 cvRisk = fit.list$cvRisk)

  }
  else {
    fit.list <- info <- setNames(vector("list", nlevels(treat)), levels(treat))

    for (i in levels(treat)) {
      treat_i <- as.numeric(treat == i)

      fit.list[[i]] <- do.call(SuperLearner::SuperLearner, list(Y = treat_i,
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
      if (discrete) ps[[i]] <- fit.list[[i]]$library.predict[,which.min(fit.list[[i]]$cvRisk)]
      else ps[[i]] <- fit.list[[i]]$SL.predict

      info[[i]] <- list(coef = fit.list[[i]]$coef,
                        cvRisk = fit.list[[i]]$cvRisk)
    }

    p.score <- NULL
  }

  #ps should be matrix of probs for each treat
  #Computing weights
  w <- get_w_from_ps(ps = ps, treat = treat, estimand, focal, stabilize = stabilize, subclass = subclass)

  obj <- list(w = w, ps = p.score, info = info, fit.obj = fit.list)
  return(obj)
}
weightit2super.cont <- function(covs, treat, s.weights, subset, stabilize, missing, ps, ...) {
  A <- B <- list(...)

  covs <- covs[subset, , drop = FALSE]
  treat <- treat[subset]

  if (missing == "ind") {
    missing.ind <- apply(covs[, apply(covs, 2, anyNA), drop = FALSE], 2, function(x) as.numeric(is.na(x)))
    covs[is.na(covs)] <- 0
    covs <- cbind(covs, missing.ind)
  }

  covs <- data.frame(apply(covs, 2, make.closer.to.1))

  if (ncol(covs) > 1) {
    colinear.covs.to.remove <- colnames(covs)[colnames(covs) %nin% colnames(make_full_rank(covs))]
    covs <- covs[, colnames(covs) %nin% colinear.covs.to.remove, drop = FALSE]
  }

  stabilize <- TRUE

  for (f in names(formals(SuperLearner::SuperLearner))) {
    if (f == "method") {if (is_null(B[["SL.method"]])) B[["SL.method"]] <- formals(SuperLearner::SuperLearner)[["method"]]}
    else if (f == "env") {if (is_null(B[["env"]])) B[["env"]] <- environment(SuperLearner::SuperLearner)}
    else if (is_null(B[[f]])) B[[f]] <- formals(SuperLearner::SuperLearner)[[f]]
  }

  discrete <- if_null_then(A[["discrete"]], FALSE)
  if (length(discrete) != 1 || !is_(discrete, "logical")) stop("discrete must be TRUE or FALSE.", call. = FALSE)

  fit <- do.call(SuperLearner::SuperLearner, list(Y = treat,
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

  if (discrete) gp.score <- fit$library.predict[,which.min(fit$cvRisk)]
  else gp.score <- fit$SL.predict

  p.denom <- treat - gp.score

  info <- list(coef = fit$coef,
               cvRisk = fit$cvRisk)

  if (isTRUE(A[["use.kernel"]])) {
    if (is_null(A[["bw"]])) A[["bw"]] <- "nrd0"
    if (is_null(A[["adjust"]])) A[["adjust"]] <- 1
    if (is_null(A[["kernel"]])) A[["kernel"]] <- "gaussian"
    if (is_null(A[["n"]])) A[["n"]] <- 10*length(treat)

    d.d <- density(p.denom, n = A[["n"]],
                   weights = s.weights[subset]/sum(s.weights[subset]), give.Rkern = FALSE,
                   bw = A[["bw"]], adjust = A[["adjust"]], kernel = A[["kernel"]])
    dens.denom <- with(d.d, approxfun(x = x, y = y))(p.denom)
  }
  else {
    if (is_null(A[["density"]])) densfun <- dnorm
    else if (is.function(A[["density"]])) densfun <- A[["density"]]
    else if (is.character(A[["density"]]) && length(A[["density"]] == 1)) {
      splitdens <- strsplit(A[["density"]], "_", fixed = TRUE)[[1]]
      if (exists(splitdens[1], mode = "function", envir = parent.frame())) {
        if (length(splitdens) > 1 && !can_str2num(splitdens[-1])) {
          stop(paste(A[["density"]], "is not an appropriate argument to density because",
                     word_list(splitdens[-1], quotes = TRUE, and.or = "or"), "cannot be coerced to numeric."), call. = FALSE)
        }
        densfun <- function(x) {
          tryCatch(do.call(get(splitdens[1]), c(list(x), as.list(str2num(splitdens[-1])))),
                   error = function(e) stop(paste0("Error in applying density:\n  ", conditionMessage(e)), call. = FALSE))
        }
      }
      else {
        stop(paste(A[["density"]], "is not an appropriate argument to density because",
                   splitdens[1], "is not an available function."), call. = FALSE)
      }
    }
    else stop("The argument to density cannot be evaluated as a density function.", call. = FALSE)

    dens.denom <- densfun(p.denom/sd(p.denom))

    if (is_null(dens.denom) || !is.atomic(dens.denom) || anyNA(dens.denom)) {
      stop("There was a problem with the output of density. Try another density function or leave it blank to use the normal density.", call. = FALSE)
    }
    else if (any(dens.denom <= 0)) {
      warning("The input to density may not accept the full range of standardized residuals of the treatment model. Some weights may be NA or zero.", call. = FALSE)
    }
  }

  if (stabilize) {

    if (is_null(A[["link"]])) A[["link"]] <- "identity"
    num.fit <- do.call(glm, list(treat ~ 1,
                                 data = data.frame(treat = treat),
                                 weights = s.weights[subset],
                                 family = gaussian(link = A[["link"]]),
                                 control = as.list(A[["control"]])),
                       quote = TRUE)

    p.num <- treat - num.fit$fitted.values

    if (isTRUE(A[["use.kernel"]])) {
      d.n <- density(p.num, n = A[["n"]],
                     weights = s.weights[subset]/sum(s.weights[subset]), give.Rkern = FALSE,
                     bw = A[["bw"]], adjust = A[["adjust"]], kernel = A[["kernel"]])
      dens.num <- with(d.n, approxfun(x = x, y = y))(p.num)
    }
    else {
      dens.num <- densfun(p.num/sqrt(summary(num.fit)$dispersion))
      if (is_null(dens.num) || !is.atomic(dens.num) || anyNA(dens.num)) {
        stop("There was a problem with the output of density. Try another density function or leave it blank to use the normal density.", call. = FALSE)
      }
      else if (any(dens.num <= 0) && !any(dens.denom <= 0)) {
        warning("The input to density may not accept the full range of treatment values. Some weights may be NA or zero.", call. = FALSE)
      }
    }
    w <- dens.num/dens.denom
  }
  else {
    w <- 1/dens.denom
  }

  if (isTRUE(A[["use.kernel"]]) && isTRUE(A[["plot"]])) {
    d.d_ <- cbind(as.data.frame(d.d[c("x", "y")]), dens = "Denominator Density", stringsAsfactors = FALSE)
    d.n_ <- cbind(as.data.frame(d.n[c("x", "y")]), dens = "Numerator Density", stringsAsfactors = FALSE)
    d.all <- rbind(d.d_, d.n_)
    d.all$dens <- factor(d.all$dens, levels = c("Numerator Density", "Denominator Density"))
    pl <- ggplot(d.all, aes(x = x, y = y)) + geom_line() +
      labs(title = "Weight Component Densities", x = "E[Treat|X]", y = "Density") +
      facet_grid(rows = vars(dens)) + theme(panel.background = element_rect(fill = "white"),
                                            panel.border = element_rect(fill = NA, color = "black"),
                                            axis.text.x = element_text(color = "black"),
                                            axis.text.y = element_text(color = "black"),
                                            panel.grid.major = element_blank(),
                                            panel.grid.minor = element_blank()
      )
    print(pl)
  }

  obj <- list(w = w, info = info, fit.obj = fit)
  return(obj)
}




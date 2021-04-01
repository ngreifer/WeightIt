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
  if (has_dots <- ("..." %in% names(Fun_formal))) {
    Fun_formal[["..."]] <- NULL
  }

  fun_args <- Fun_formal
  for (i in names(fun_args)) {
    if (exists(i, inherits = FALSE)) fun_args[i] <- list(get0(i, inherits = FALSE))
    else if (i %in% names(A)) {
      fun_args[i] <- A[i]
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

      if (length(Anames) > 1) warning(paste0("The following arguments were specified but are not suitable arguments to the provided function:\n\t", word_list(Anames)), call. = FALSE, immediate. = TRUE)
      else warning(paste0("The following argument was specified but is not a suitable argument to the provided function:\n\t", Anames), call. = FALSE, immediate. = TRUE)
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
    if (exists(i, inherits = FALSE)) fun_args[i] <- list(get0(i, inherits = FALSE))
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

      if (length(Anames) > 1) warning(paste0("The following arguments were specified but are not suitable arguments to the provided function:\n\t", word_list(Anames)), call. = FALSE, immediate. = TRUE)
      else warning(paste0("The following argument was specified but is not a suitable argument to the provided function:\n\t", Anames), call. = FALSE, immediate. = TRUE)
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
weightit2ps <- function(covs, treat, s.weights, subset, estimand, focal, stabilize, subclass, missing, ps, .data, ...) {
  A <- list(...)

  fit.obj <- NULL

  if (is_null(ps)) {

    covs <- covs[subset, , drop = FALSE]
    treat_sub <- factor(treat[subset])
    s.weights <- s.weights[subset]
    bin.treat <- is_binary(treat_sub)
    ord.treat <- is.ordered(treat_sub)

    if (missing == "ind") {
      missing.ind <- apply(covs[, anyNA_col(covs), drop = FALSE], 2, function(x) as.numeric(is.na(x)))
      if (is_not_null(missing.ind)) {
        covs[is.na(covs)] <- 0
        covs <- cbind(covs, missing.ind)
      }
    }

    for (i in seq_col(covs)) covs[,i] <- make.closer.to.1(covs[,i])

    if (ncol(covs) > 1) {
      if (missing == "saem") {
        covs0 <- covs
        for (i in colnames(covs)[anyNA_col(covs)]) covs0[is.na(covs0[,i]),i] <- covs0[!is.na(covs0[,i]),i][1]
        colinear.covs.to.remove <- colnames(covs)[colnames(covs) %nin% colnames(make_full_rank(covs0))]
      }
      else colinear.covs.to.remove <- colnames(covs)[colnames(covs) %nin% colnames(make_full_rank(covs))]
      covs <- covs[, colnames(covs) %nin% colinear.covs.to.remove, drop = FALSE]
    }

    if (is_null(A$link)) A$link <- "logit"
    else {
      if (ord.treat) acceptable.links <- c("logit", "probit", "loglog", "cloglog", "cauchit")
      else if (bin.treat || isFALSE(A$use.mlogit)) {
        if (missing == "saem") acceptable.links <- "logit"
        else acceptable.links <- expand.grid_string(c("", "br."), c("logit", "probit", "cloglog", "identity", "log", "cauchit"))
      }
      else acceptable.links <- c("logit", "probit", "bayes.probit", "br.logit")

      which.link <- acceptable.links[pmatch(A$link, acceptable.links, nomatch = 0)][1]
      if (is.na(which.link)) {
        A$link <- acceptable.links[1]
        warning(paste0("Only ", word_list(acceptable.links, quotes = TRUE, is.are = TRUE), " allowed as the link for ",
                       if (bin.treat) "binary" else if (ord.treat) "ordinal" else "multinomial",
                       " treatments", if (missing == "saem") " with missing = \"saem\"",
                       ". Using link = ", word_list(acceptable.links[1], quotes = TRUE), "."),
                call. = FALSE, immediate. = TRUE)
      }
      else A$link <- which.link
    }

    use.br <- startsWith(A$link, "br.")
    # use.bayes <- startsWith(A$link, "bayes.")
    if (use.br) A$link <- substr(A$link, 4, nchar(A$link))
    # else if (use.bayes) A$link <- substr(A$link, 7, nchar(A$link))

    if (bin.treat) {

      t.lev <- get.treated.level(treat_sub)
      c.lev <- setdiff(levels(treat_sub), t.lev)

      family <- binomial(link = A[["link"]])

      ps <- make_df(levels(treat_sub), nrow = length(treat_sub))

      if (missing == "saem") {
        # check.package("misaem")
        #
        # newdata <- data.frame(binarize(treat_sub, one = t.lev), covs)
        # fit <- misaem::miss.glm(formula(newdata), newdata)
        #
        # if (is_null(A[["saem.method"]])) A[["saem.method"]] <- "map"
        #
        # ps[[t.lev]] <- p.score <- drop(predict(fit, newdata = covs, method = A[["saem.method"]]))
        # ps[[c.lev]] <- 1 - ps[[t.lev]]
      }
      else {
        if (use.br) {
          ctrl_fun <- brglm2::brglmControl
          glm_method <- brglm2::brglmFit
        }
        else {
          ctrl_fun <- stats::glm.control
          glm_method <- if_null_then(A[["glm.method"]], stats::glm.fit)
        }

        control <- A[names(formals(ctrl_fun))[pmatch(names(A), names(formals(ctrl_fun)), 0)]]

        treat <- binarize(treat_sub, one = t.lev)

        withCallingHandlers({
          if (isTRUE(A[["quick"]])) {
            fit <- do.call(glm_method, c(list(y = treat,
                                              x = cbind(1, covs),
                                              start = c(family$linkfun(mean(treat)), rep(0, ncol(covs))),
                                              weights = s.weights,
                                              family = family),
                                         control), quote = TRUE)
          }
          else {
            data <- data.frame(treat, covs)
            formula <- formula(data)
            fit <- do.call(stats::glm, c(list(formula, data = data,
                                              weights = s.weights,
                                              start = c(family$linkfun(mean(treat)), rep(0, ncol(covs))),
                                              family = family,
                                              method = glm_method),
                                         control), quote = TRUE)

          }
        },
        warning = function(w) {
          if (conditionMessage(w) != "non-integer #successes in a binomial glm!") warning(w)
          invokeRestart("muffleWarning")
        })

        ps[[t.lev]] <- p.score <- fit$fitted.values
        ps[[c.lev]] <- 1 - ps[[t.lev]]
      }

      fit[["call"]] <- NULL
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
                                      contrasts = NULL)), quote = TRUE)},
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
                                   control), quote = TRUE)},
                 error = function(e) stop("There was a problem with the bias-reduced multinomial logit regression. Try a different link.", call. = FALSE))

        ps <- fit$fitted.values
        fit.obj <- fit
      }
      else if (A$link %in% c("logit", "probit")) {
        if (isTRUE(A$use.mclogit)) {
          check.package("mclogit")

          if (is_not_null(A$random)) {
            random <- get.covs.and.treat.from.formula(A$random, data = .data)$reported.covs[subset,,drop = FALSE]
            data <- cbind(data.frame(random), data.frame(treat = treat_sub, .s.weights = s.weights, covs))
            covnames <- names(data)[-c(seq_col(random), ncol(random) + (1:2))]
            tname <- names(data)[ncol(random) + 1]
            ctrl_fun <- mclogit::mmclogit.control
          }
          else {
            data <- data.frame(treat = treat_sub, .s.weights = s.weights, covs)
            covnames <- names(data)[-c(1,2)]
            tname <- names(data)[1]
            ctrl_fun <- mclogit::mclogit.control
          }
          form <- reformulate(covnames, tname)

          control <- do.call(ctrl_fun, c(A[["control"]], A[names(formals(ctrl_fun))[pmatch(names(A), names(formals(ctrl_fun)), 0)]]))

          tryCatch({
            fit <- do.call(mclogit::mblogit, list(form,
                                    data = data,
                                    weights = quote(.s.weights),
                                    random = A[["random"]],
                                    method = A[["mclogit.method"]],
                                    estimator = if_null_then(A[["estimator"]], eval(formals(mclogit::mclogit)[["estimator"]])),
                                    dispersion = if_null_then(A[["dispersion"]], eval(formals(mclogit::mclogit)[["dispersion"]])),
                                    groups = A[["groups"]],
                                    control = control))

          },
          error = function(e) {stop(paste0("There was a problem fitting the multinomial ", A$link, " regression with mblogit().\n       Try again with use.mclogit = FALSE.\nError message (from mclogit):\n       ", conditionMessage(e)), call. = FALSE)}
          )

          ps <- fitted(fit)
          colnames(ps) <- levels(treat_sub)
          fit.obj <- fit
        }
        else if (!isFALSE(A$use.mlogit)) {
          check.package("mlogit")

          data <- data.frame(treat = treat_sub, .s.weights = s.weights, covs)
          covnames <- names(data)[-c(1,2)]
          tryCatch({
            fit <- mlogit::mlogit(as.formula(paste0("treat ~ 1 | ", paste(covnames, collapse = " + "))),
                                  data = data,
                                  estimate = TRUE,
                                  probit = A$link[1] == "probit",
                                  weights = .s.weights,
                                  varying = NULL,
                                  shape = "wide",
                                  sep = "",
                                  choice = "treat",
                                  ...)
          },
          error = function(e) {stop(paste0("There was a problem fitting the multinomial ", A$link, " regression with mlogit().\n       Try again with use.mlogit = FALSE.\nError message (from mlogit):\n       ", conditionMessage(e)), call. = FALSE)}
          )

          ps <- fitted(fit, outcome = FALSE)
          fit.obj <- fit
        }
        else {
          ps <- make_df(levels(treat_sub), nrow = length(treat_sub))

          control <- A[names(formals(stats::glm.control))[pmatch(names(A), names(formals(stats::glm.control)), 0)]]
          fit.list <- make_list(levels(treat_sub))

          for (i in levels(treat_sub)) {
            if (isTRUE(A[["test1"]])) {
              if (i == last(levels(treat_sub))) {
                ps[[i]] <- 1 - rowSums(ps[names(ps) != i])
                next
              }
            }
            t_i <- as.numeric(treat_sub == i)
            data_i <- data.frame(t_i, covs)
            fit.list[[i]] <- do.call(stats::glm, c(list(formula(data_i), data = data_i,
                                                        family = quasibinomial(link = A$link),
                                                        weights = s.weights),
                                                   control), quote = TRUE)
            ps[[i]] <- fit.list[[i]]$fitted.values
          }
          if (isTRUE(A[["test2"]])) ps <- ps/rowSums(ps)
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
      t.lev <- get.treated.level(treat)
      c.lev <- setdiff(levels(treat_sub), t.lev)

      if (is_(ps, c("matrix", "data.frame"))) {
        if (all(dim(ps) == c(n, 2))) {

          if (all(colnames(ps) %in% levels(treat_sub))) {
            ps <- as.data.frame(ps[subset, , drop = FALSE])
          }
          else {
            ps <- as.data.frame(ps[subset, , drop = FALSE])
            names(ps) <- levels(treat_sub)
          }

          p.score <- ps[[t.lev]]
        }
        else if (all(dim(ps) == c(length(treat), 1))) {

          ps <- data.frame(ps[subset,1], 1-ps[subset,1])

          names(ps) <- c(t.lev, c.lev)

          p.score <- ps[[t.lev]]
        }
      }
      else if (is.numeric(ps)) {
        if (length(ps) == n) {
          ps <- data.frame(ps[subset], 1-ps[subset])

          names(ps) <- c(t.lev, c.lev)

          p.score <- ps[[t.lev]]
        }
      }

      if (is_null(p.score)) stop("'ps' must be a numeric vector with a propensity score for each unit.", call. = FALSE)

    }
    else {
      bad.ps <- FALSE
      if (is_(ps, c("matrix", "data.frame"))) {
        if (all(dim(ps) == c(n, nunique(treat)))) {
          ps <- setNames(as.data.frame(ps), levels(treat))[subset, , drop = FALSE]
        }
        else if (all(dim(ps) == c(n, 1))) {
          ps <- setNames(list2DF(lapply(levels(treat), function(x) {
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
          ps <- setNames(list2DF(lapply(levels(treat), function(x) {
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
  w <- get_w_from_ps(ps = ps, treat = treat_sub, estimand, focal = focal,
                     stabilize = stabilize, subclass = subclass)

  obj <- list(w = w, ps = p.score, fit.obj = fit.obj)
  return(obj)
}
weightit2ps.cont <- function(covs, treat, s.weights, subset, stabilize, missing, ps, ...) {
  A <- list(...)

  fit.obj <- NULL

  covs <- covs[subset, , drop = FALSE]
  treat <- treat[subset]
  s.weights <- s.weights[subset]

  if (missing == "ind") {
    missing.ind <- apply(covs[, anyNA_col(covs), drop = FALSE], 2, function(x) as.numeric(is.na(x)))
    if (is_not_null(missing.ind)) {
      covs[is.na(covs)] <- 0
      covs <- cbind(covs, missing.ind)
    }
  }
  for (i in seq_col(covs)) covs[,i] <- make.closer.to.1(covs[,i])

  if (ncol(covs) > 1) {
    if (missing == "saem") {
      covs0 <- covs
      for (i in colnames(covs)[anyNA_col(covs)]) covs0[is.na(covs0[,i]),i] <- covs0[!is.na(covs0[,i]),i][1]
      colinear.covs.to.remove <- colnames(covs)[colnames(covs) %nin% colnames(make_full_rank(covs0))]
    }
    else colinear.covs.to.remove <- colnames(covs)[colnames(covs) %nin% colnames(make_full_rank(covs))]
    covs <- covs[, colnames(covs) %nin% colinear.covs.to.remove, drop = FALSE]
  }

  data <- data.frame(treat, covs)
  formula <- formula(data)

  #Process density params
  if (isTRUE(A[["use.kernel"]])) {
    if (is_null(A[["bw"]])) A[["bw"]] <- "nrd0"
    if (is_null(A[["adjust"]])) A[["adjust"]] <- 1
    if (is_null(A[["kernel"]])) A[["kernel"]] <- "gaussian"
    if (is_null(A[["n"]])) A[["n"]] <- 10*length(treat)
    use.kernel <- TRUE
    densfun <- NULL
  }
  else {
    if (is_null(A[["density"]])) densfun <- dnorm
    else if (is.function(A[["density"]])) densfun <- A[["density"]]
    else if (is.character(A[["density"]]) && length(A[["density"]] == 1)) {
      splitdens <- strsplit(A[["density"]], "_", fixed = TRUE)[[1]]
      if (is_not_null(splitdens1 <- get0(splitdens[1], mode = "function", envir = parent.frame()))) {
        if (length(splitdens) > 1 && !can_str2num(splitdens[-1])) {
          stop(paste(A[["density"]], "is not an appropriate argument to 'density' because",
                     word_list(splitdens[-1], and.or = "or", quotes = TRUE), "cannot be coerced to numeric."), call. = FALSE)
        }
        densfun <- function(x) {
          tryCatch(do.call(splitdens1, c(list(x), as.list(str2num(splitdens[-1])))),
                   error = function(e) stop(paste0("Error in applying density:\n  ", conditionMessage(e)), call. = FALSE))
        }
      }
      else {
        stop(paste(A[["density"]], "is not an appropriate argument to 'density' because",
                   splitdens[1], "is not an available function."), call. = FALSE)
      }
    }
    else stop("The argument to 'density' cannot be evaluated as a density function.", call. = FALSE)
    use.kernel <- FALSE
  }

  #Stabilization - get dens.num
  p.num <- treat - mean(treat)

  if (use.kernel) {
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

  #Estimate GPS
  if (is_null(ps)) {

    if (missing == "saem") {
      # check.package("misaem")
      #
      # tryCatch({fit <- misaem::miss.lm(formula, data)},
      #          error = function(e) {
      #            if (trimws(conditionMessage(e)) == "object 'p' not found") {
      #              stop("missing = \"saem\" cannot be used with continuous treatments due to a bug in the misaem package.", call. = FALSE)
      #            }
      #          })
      #
      # if (is_null(A[["saem.method"]])) A[["saem.method"]] <- "map"
      #
      # gp.score <- drop(predict(fit, newdata = covs, method = A[["saem.method"]]))

    }
    else {
      if (is_null(A$link)) A$link <- "identity"
      fit <- do.call("glm", c(list(formula, data = data,
                                   weights = s.weights,
                                   family = gaussian(link = A$link),
                                   control = as.list(A$control))),
                     quote = TRUE)

      gp.score <- fit$fitted.values
    }

    fit.obj <- fit
  }

  #Get weights
  w <- get_cont_weights(gp.score, treat = treat, s.weights = s.weights,
                        dens.num = dens.num, densfun = densfun,
                        use.kernel = use.kernel, densControl = A)

  if (use.kernel && isTRUE(A[["plot"]])) {
    d.d <- density(treat - gp.score, n = A[["n"]],
                   weights = s.weights/sum(s.weights), give.Rkern = FALSE,
                   bw = A[["bw"]], adjust = A[["adjust"]],
                   kernel = A[["kernel"]])
    plot_density(d.n, d.d)
  }

  obj <- list(w = w, fit.obj = fit.obj)
  return(obj)
}

#MABW with optweight
weightit2optweight <- function(covs, treat, s.weights, subset, estimand, focal, missing, moments, int, ...) {
  A <- list(...)

  check.package("optweight")

  covs <- covs[subset, , drop = FALSE]
  treat <- factor(treat[subset])
  s.weights <- s.weights[subset]

  covs <- cbind(covs, int.poly.f(covs, poly = moments, int = int))
  for (i in seq_col(covs)) covs[,i] <- make.closer.to.1(covs[,i])

  if (missing == "ind") {
    missing.ind <- apply(covs[, anyNA_col(covs), drop = FALSE], 2, function(x) as.numeric(is.na(x)))
    if (is_not_null(missing.ind)) {
      covs[is.na(covs)] <- 0
      covs <- cbind(covs, missing.ind)
    }
  }

  new.data <- data.frame(treat, covs)
  new.formula <- formula(new.data)

  for (f in names(formals(optweight::optweight))) {
    if (is_null(A[[f]])) A[[f]] <- formals(optweight::optweight)[[f]]
  }
  A[names(A) %in% names(formals(weightit2optweight))] <- NULL

  if ("tols" %in% names(A)) A[["tols"]] <- optweight::check.tols(new.formula, new.data, A[["tols"]], stop = TRUE)
  if ("targets" %in% names(A)) {
    warning("targets cannot be used through WeightIt and will be ignored.", call. = FALSE, immediate. = TRUE)
    A[["targets"]] <- NULL
  }

  A[["formula"]] <- new.formula
  A[["data"]] <- new.data
  A[["estimand"]] <- estimand
  A[["s.weights"]] <- s.weights
  A[["focal"]] <- focal
  A[["verbose"]] <- TRUE

  out <- do.call(optweight::optweight, A, quote = TRUE)

  obj <- list(w = out[["weights"]], info = list(duals = out$duals), fit.obj = out)
  return(obj)
}
weightit2optweight.cont <- function(covs, treat, s.weights, subset, missing, moments, int, ...) {
  A <- list(...)
  check.package("optweight")

  covs <- covs[subset, , drop = FALSE]
  treat <- treat[subset]
  s.weights <- s.weights[subset]

  covs <- cbind(covs, int.poly.f(covs, poly = moments, int = int))
  for (i in seq_col(covs)) covs[,i] <- make.closer.to.1(covs[,i])

  if (missing == "ind") {
    missing.ind <- apply(covs[, anyNA_col(covs), drop = FALSE], 2, function(x) as.numeric(is.na(x)))
    if (is_not_null(missing.ind)) {
      covs[is.na(covs)] <- 0
      covs <- cbind(covs, missing.ind)
    }
  }

  new.data <- data.frame(treat, covs)
  new.formula <- formula(new.data)

  for (f in names(formals(optweight::optweight))) {
    if (is_null(A[[f]])) A[[f]] <- formals(optweight::optweight)[[f]]
  }
  A[names(A) %in% names(formals(weightit2optweight.cont))] <- NULL

  if ("tols" %in% names(A)) A[["tols"]] <- optweight::check.tols(new.formula, new.data, A[["tols"]], stop = TRUE)
  if ("targets" %in% names(A)) {
    warning("targets cannot be used through WeightIt and will be ignored.", call. = FALSE, immediate. = TRUE)
    A[["targets"]] <- NULL
  }

  A[["formula"]] <- new.formula
  A[["data"]] <- new.data
  A[["s.weights"]] <- s.weights
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
      for (i in seq_col(covs)) covs[,i] <- make.closer.to.1(covs[,i])

      if (missing == "ind") {
        missing.ind <- apply(covs[, anyNA_col(covs), drop = FALSE], 2, function(x) as.numeric(is.na(x)))
        if (is_not_null(missing.ind)) {
          covs[is.na(covs)] <- 0
          covs <- cbind(covs, missing.ind)
        }
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
    warning("targets cannot be used through WeightIt and will be ignored.", call. = FALSE, immediate. = TRUE)
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

#Generalized boosted modeling with gbm and cobalt
weightit2gbm <- function(covs, treat, s.weights, estimand, focal, subset, stabilize, subclass, missing, ...) {

  check.package("gbm")

  A <- list(...)

  covs <- covs[subset, , drop = FALSE]
  treat <- treat[subset]
  s.weights <- s.weights[subset]

  if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
  treat.type <- get.treat.type(treat)

  for (i in seq_col(covs)) covs[,i] <- make.closer.to.1(covs[,i])

  if (missing == "ind") {
    missing.ind <- apply(covs[, anyNA_col(covs), drop = FALSE], 2, function(x) as.numeric(is.na(x)))
    if (is_not_null(missing.ind)) {
      colnames(missing.ind) <- paste0(colnames(missing.ind), ":<NA>")
      covs <- cbind(covs, missing.ind)
    }
  }

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
      if (cv < 2) stop("At least 2 CV-folds must be specified in stop.method.", call. = FALSE)
    }
    else stop(paste0("'stop.method' must be one of ", word_list(c(available.stop.methods, "cv{#}"), "or", quotes = TRUE), "."), call. = FALSE)
  }
  else stop.method <- available.stop.methods[s.m.matches]

  tunable <- c("interaction.depth", "shrinkage", "distribution")

  trim.at <- if_null_then(A[["trim.at"]], 0)

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

  if (treat.type == "binary")  {
    available.distributions <- c("bernoulli", "adaboost")
    t.lev <- get.treated.level(treat)
    treat <- binarize(treat, one = focal)
  }
  else {
    available.distributions <- "multinomial"
    treat <- factor(treat)
  }

  if (cv == 0) {
    start.tree <- if_null_then(A[["start.tree"]], 1)
    if (is_null(A[["n.grid"]])) n.grid <- round(1+sqrt(2*(A[["n.trees"]]-start.tree+1)))
    else if (!is_(A[["n.grid"]], "numeric") || length(A[["n.grid"]]) > 1 ||
             !between(A[["n.grid"]], c(2, A[["n.trees"]]))) {
      stop("'n.grid' must be a numeric value between 2 and n.trees.", call. = FALSE)
    }
    else n.grid <- round(A[["n.grid"]])

    crit <- bal_criterion(treat.type, stop.method)
    init <- crit$init(covs, treat, estimand = estimand, s.weights = s.weights, focal = focal, ...)
  }

  A[["x"]] <- covs
  A[["y"]] <- treat
  A[["distribution"]] <- if (is_null(distribution <- A[["distribution"]])) {
    available.distributions[1]} else {
      match_arg(distribution, available.distributions, several.ok = TRUE)}
  A[["w"]] <- s.weights
  A[["verbose"]] <- FALSE

  tune <- do.call("expand.grid", c(A[names(A) %in% tunable],
                                   list(stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)))

  current.best.loss <- Inf

  for (i in seq_row(tune)) {

    A[["distribution"]] <- list(name = tune[["distribution"]][i])
    tune_args <- as.list(tune[i, setdiff(tunable, "distribution")])

    fit <- do.call(gbm::gbm.fit, c(A[names(A) %in% setdiff(names(formals(gbm::gbm.fit)), names(tune_args))], tune_args), quote = TRUE)

    if (cv == 0) {

      n.trees <- fit[["n.trees"]]
      iters <- 1:n.trees
      iters.grid <- round(seq(start.tree, n.trees, length.out = n.grid))

      if (is_null(iters.grid) || anyNA(iters.grid) || any(iters.grid > n.trees)) stop("A problem has occurred")

      ps <- gbm::predict.gbm(fit, n.trees = iters.grid, type = "response", newdata = covs)
      w <- get.w.from.ps(ps, treat = treat, estimand = estimand, focal = focal, stabilize = stabilize, subclass = subclass)
      if (trim.at != 0) w <- suppressMessages(apply(w, 2, trim, at = trim.at, treat = treat))

      iter.grid.balance <- apply(w, 2, function(w_) {
        crit$fun(init = init, weights = w_)
      })

      if (n.grid == n.trees) {
        best.tree.index <- which.min(iter.grid.balance)
        best.loss <- iter.grid.balance[best.tree.index]
        best.tree <- iters.grid[best.tree.index]
        tree.val <- setNames(data.frame(iters.grid,
                                        iter.grid.balance),
                             c("tree", stop.method))
      }
      else {
        it <- which.min(iter.grid.balance) + c(-1, 1)
        it[1] <- iters.grid[max(1, it[1])]
        it[2] <- iters.grid[min(length(iters.grid), it[2])]
        iters.to.check <- iters[between(iters, iters[it])]

        if (is_null(iters.to.check) || anyNA(iters.to.check) || any(iters.to.check > n.trees)) stop("A problem has occurred")

        ps <- gbm::predict.gbm(fit, n.trees = iters.to.check, type = "response", newdata = covs)
        w <- get.w.from.ps(ps, treat = treat, estimand = estimand, focal = focal, stabilize = stabilize, subclass = subclass)
        if (trim.at != 0) w <- suppressMessages(apply(w, 2, trim, at = trim.at, treat = treat))

        iter.grid.balance.fine <- apply(w, 2, function(w_) {
          crit$fun(init = init, weights = w_)
        })

        best.tree.index <- which.min(iter.grid.balance.fine)
        best.loss <- iter.grid.balance.fine[best.tree.index]
        best.tree <- iters.to.check[best.tree.index]
        tree.val <- setNames(data.frame(c(iters.grid, iters.to.check),
                                        c(iter.grid.balance, iter.grid.balance.fine)),
                             c("tree", stop.method))
      }

      tree.val <- unique(tree.val[order(tree.val$tree),])
      w <- w[,best.tree.index]
      ps <- if (treat.type == "binary") ps[,best.tree.index] else NULL

      tune[[paste.("best", stop.method)]][i] <- best.loss
      tune[["best.tree"]][i] <- best.tree

      if (best.loss < current.best.loss) {
        best.fit <- fit
        best.w <- w
        best.ps <- ps
        current.best.loss <- best.loss
        best.tune.index <- i

        info <- list(best.tree = best.tree,
                     tree.val = tree.val)
      }
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
      A[["distribution"]] <- list(name = tune[["distribution"]][i])
      A[["w"]] <- s.weights

      tune_args <- as.list(tune[i, setdiff(tunable, "distribution")])
      cv.results <- do.call(gbm::gbmCrossVal,
                            c(A[names(A) %in% setdiff(names(formals(gbm::gbmCrossVal)), names(tune_args))],
                              tune_args), quote = TRUE)
      best.tree.index <- which.min(cv.results$error)
      best.loss <- cv.results$error[best.tree.index]
      best.tree <- best.tree.index

      tune[[paste.("best", names(fit$name))]][i] <- best.loss
      tune[["best.tree"]][i] <- best.tree

      if (best.loss < current.best.loss) {
        best.fit <- fit
        best.ps <- gbm::predict.gbm(fit, n.trees = best.tree, type = "response", newdata = covs)
        best.w <- drop(get.w.from.ps(best.ps, treat = treat, estimand = estimand, focal = focal, stabilize = stabilize, subclass = subclass))
        # if (trim.at != 0) best.w <- suppressMessages(trim(best.w, at = trim.at, treat = treat))
        current.best.loss <- best.loss
        best.tune.index <- i

        tree.val <- data.frame(tree = seq_along(cv.results$error),
                               error = cv.results$error)

        info <- list(best.tree = best.tree,
                     tree.val = tree.val)

        if (treat.type == "multinomial") best.ps <- NULL
      }
    }

    if (treat.type == "multinomial") ps <- NULL
  }

  tune[tunable[vapply(tune[tunable], all_the_same, logical(1L))]] <- NULL

  if (ncol(tune) > 2) {
    info[["tune"]] <- tune
    info[["best.tune"]] <- tune[best.tune.index,]
  }

  if (is_not_null(best.ps)) {
    if (is_not_null(focal) && focal != t.lev) best.ps <- 1 - best.ps
  }

  obj <- list(w = best.w, ps = best.ps, info = info, fit.obj = best.fit)
  return(obj)
}
weightit2gbm.cont <- function(covs, treat, s.weights, estimand, focal, subset, stabilize, subclass, missing, ...) {

  check.package("gbm")

  A <- list(...)

  covs <- covs[subset, , drop = FALSE]
  treat <- treat[subset]
  s.weights <- s.weights[subset]

  for (i in seq_col(covs)) covs[,i] <- make.closer.to.1(covs[,i])

  if (missing == "ind") {
    missing.ind <- apply(covs[, anyNA_col(covs), drop = FALSE], 2, function(x) as.numeric(is.na(x)))
    if (is_not_null(missing.ind)) {
      colnames(missing.ind) <- paste0(colnames(missing.ind), ":<NA>")
      covs <- cbind(covs, missing.ind)
    }
  }

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
  available.stop.methods <- bal_criterion("continuous", list = TRUE)
  s.m.matches <- charmatch(A[["stop.method"]], available.stop.methods)
  if (is.na(s.m.matches) || s.m.matches == 0L) {
    if (startsWith(A[["stop.method"]], "cv") && can_str2num(numcv <- substr(A[["stop.method"]], 3, nchar(A[["stop.method"]])))) {
      cv <- round(str2num(numcv))
      if (cv < 2) stop("At least 2 CV-folds must be specified in stop.method.", call. = FALSE)
    }
    else stop(paste0("'stop.method' must be one of ", word_list(c(available.stop.methods, "cv{#}"), "or", quotes = TRUE), "."), call. = FALSE)
  }
  else stop.method <- available.stop.methods[s.m.matches]

  tunable <- c("interaction.depth", "shrinkage", "distribution")

  trim.at <- if_null_then(A[["trim.at"]], 0)

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
      stop("'n.grid' must be a numeric value between 2 and n.trees.", call. = FALSE)
    }
    else n.grid <- round(A[["n.grid"]])

    crit <- bal_criterion("continuous", stop.method)
    init <- crit$init(covs, treat, s.weights = s.weights, ...)
  }

  A[["x"]] <- covs
  A[["y"]] <- treat
  A[["distribution"]] <- if (is_null(distribution <- A[["distribution"]])) {
    available.distributions[1]} else {
      match_arg(distribution, available.distributions, several.ok = TRUE)}
  A[["w"]] <- s.weights
  A[["verbose"]] <- FALSE

  if (!is.numeric(A[["n.trees"]]) || length(A[["n.trees"]]) > 1 || A[["n.trees"]] <= 1) {
    stop("'n.trees' must be a number greater than 1.", call. = FALSE)
  }

  tune <- do.call("expand.grid", c(A[names(A) %in% tunable],
                                   list(stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)))

  #Process density params
  if (isTRUE(A[["use.kernel"]])) {
    if (is_null(A[["bw"]])) A[["bw"]] <- "nrd0"
    if (is_null(A[["adjust"]])) A[["adjust"]] <- 1
    if (is_null(A[["kernel"]])) A[["kernel"]] <- "gaussian"
    if (is_null(A[["n"]])) A[["n"]] <- 10*length(treat)
    use.kernel <- TRUE
    densfun <- NULL
  }
  else {
    if (is_null(A[["density"]])) densfun <- dnorm
    else if (is.function(A[["density"]])) densfun <- A[["density"]]
    else if (is.character(A[["density"]]) && length(A[["density"]] == 1)) {
      splitdens <- strsplit(A[["density"]], "_", fixed = TRUE)[[1]]
      if (exists(splitdens[1], mode = "function", envir = parent.frame())) {
        if (length(splitdens) > 1 && !can_str2num(splitdens[-1])) {
          stop(paste(A[["density"]], "is not an appropriate argument to 'density' because",
                     word_list(splitdens[-1], and.or = "or", quotes = TRUE), "cannot be coerced to numeric."), call. = FALSE)
        }
        densfun <- function(x) {
          tryCatch(do.call(get(splitdens[1]), c(list(x), as.list(str2num(splitdens[-1])))),
                   error = function(e) stop(paste0("Error in applying density:\n  ", conditionMessage(e)), call. = FALSE))
        }
      }
      else {
        stop(paste(A[["density"]], "is not an appropriate argument to 'density' because",
                   splitdens[1], "is not an available function."), call. = FALSE)
      }
    }
    else stop("The argument to 'density' cannot be evaluated as a density function.", call. = FALSE)
    use.kernel <- FALSE
  }

  #Stabilization - get dens.num
  p.num <- treat - mean(treat)

  if (use.kernel) {
    d.n <- density(p.num, n = A[["n"]],
                   weights = s.weights/sum(s.weights), give.Rkern = FALSE,
                   bw = A[["bw"]], adjust = A[["adjust"]], kernel = A[["kernel"]])
    dens.num <- with(d.n, approxfun(x = x, y = y))(p.num)
  }
  else {
    dens.num <- densfun(p.num/sd(treat))
    if (is_null(dens.num) || !is.atomic(dens.num) || anyNA(dens.num)) {
      stop("There was a problem with the output of 'density'. Try another density function or leave it blank to use the normal density.", call. = FALSE)
    }
    else if (any(dens.num <= 0)) {
      stop("The input to 'density' may not accept the full range of treatment values.", call. = FALSE)
    }
  }

  current.best.loss <- Inf

  for (i in seq_row(tune)) {

    fit <- do.call(gbm::gbm.fit, c(A[names(A) %in% setdiff(names(formals(gbm::gbm.fit)), tunable)],
                                   tune[i, tunable[tunable %in% names(formals(gbm::gbm.fit))]]), quote = TRUE)

    if (cv == 0) {

      n.trees <- fit[["n.trees"]]
      iters <- 1:n.trees
      iters.grid <- round(seq(start.tree, n.trees, length.out = n.grid))

      if (is_null(iters.grid) || anyNA(iters.grid) || any(iters.grid > n.trees)) stop("A problem has occurred")

      gps <- gbm::predict.gbm(fit, n.trees = iters.grid, newdata = covs)
      w <- get_cont_weights(gps, treat = treat, s.weights = s.weights, dens.num = dens.num,
                            densfun = densfun, use.kernel = use.kernel, densControl = A)
      if (trim.at != 0) w <- suppressMessages(apply(w, 2, trim, at = trim.at, treat = treat))

      iter.grid.balance <- apply(w, 2, function(w_) {
        crit$fun(init = init, weights = w_)
      })

      if (n.grid == n.trees) {
        best.tree.index <- which.min(iter.grid.balance)
        best.loss <- iter.grid.balance[best.tree.index]
        best.tree <- iters.grid[best.tree.index]
        tree.val <- setNames(data.frame(iters.grid,
                                        iter.grid.balance),
                             c("tree", stop.method))
      }
      else {
        it <- which.min(iter.grid.balance) + c(-1, 1)
        it[1] <- iters.grid[max(1, it[1])]
        it[2] <- iters.grid[min(length(iters.grid), it[2])]
        iters.to.check <- iters[between(iters, iters[it])]

        if (is_null(iters.to.check) || anyNA(iters.to.check) || any(iters.to.check > n.trees)) stop("A problem has occurred")

        gps <- gbm::predict.gbm(fit, n.trees = iters.to.check, newdata = covs)
        w <- get_cont_weights(gps, treat = treat, s.weights = s.weights, dens.num = dens.num,
                              densfun = densfun, use.kernel = use.kernel, densControl = A)
        if (trim.at != 0) w <- suppressMessages(apply(w, 2, trim, at = trim.at, treat = treat))

        iter.grid.balance.fine <- apply(w, 2, function(w_) {
          crit$fun(init = init, weights = w_)
        })

        best.tree.index <- which.min(iter.grid.balance.fine)
        best.loss <- iter.grid.balance.fine[best.tree.index]
        best.tree <- iters.to.check[best.tree.index]
        tree.val <- setNames(data.frame(c(iters.grid, iters.to.check),
                                        c(iter.grid.balance, iter.grid.balance.fine)),
                             c("tree", stop.method))
      }

      tree.val <- unique(tree.val[order(tree.val$tree),])
      w <- w[,best.tree.index]
      gps <- gps[,as.character(best.tree)]

      tune[[paste.("best", stop.method)]][i] <- best.loss
      tune[["best.tree"]][i] <- best.tree

      if (best.loss < current.best.loss) {
        best.fit <- fit
        best.w <- w
        best.gps <- gps
        current.best.loss <- best.loss
        best.tune.index <- i

        info <- list(best.tree = best.tree,
                     tree.val = tree.val)
      }
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
                            c(A[names(A) %in% setdiff(names(formals(gbm::gbmCrossVal)), tunable)],
                              tune[i, tunable[tunable %in% names(formals(gbm::gbmCrossVal))]]), quote = TRUE)
      best.tree.index <- which.min(cv.results$error)
      best.loss <- cv.results$error[best.tree.index]
      best.tree <- best.tree.index

      tune[[paste.("best", "error")]][i] <- best.loss
      tune[["best.tree"]][i] <- best.tree

      if (best.loss < current.best.loss) {
        best.fit <- fit
        best.gps <- gbm::predict.gbm(fit, n.trees = best.tree, newdata = covs)
        best.w <- get_cont_weights(best.gps, treat = treat, s.weights = s.weights, dens.num = dens.num,
                                   densfun = densfun, use.kernel = use.kernel, densControl = A)
        # if (trim.at != 0) best.w <- suppressMessages(trim(best.w, at = trim.at, treat = treat))
        current.best.loss <- best.loss
        best.tune.index <- i

        tree.val <- data.frame(tree = seq_along(cv.results$error),
                               error = cv.results$error)

        info <- list(best.tree = best.tree,
                     tree.val = tree.val)

      }
    }

  }

  if (use.kernel && isTRUE(A[["plot"]])) {
    d.d <- density(treat - best.gps, n = A[["n"]],
                   weights = s.weights/sum(s.weights), give.Rkern = FALSE,
                   bw = A[["bw"]], adjust = A[["adjust"]],
                   kernel = A[["kernel"]])
    plot_density(d.n, d.d)
  }

  tune[tunable[vapply(tunable, function(x) length(A[[x]]) == 1, logical(1L))]] <- NULL

  if (ncol(tune) > 2) {
    info[["tune"]] <- tune
    info[["best.tune"]] <- tune[best.tune.index,]
  }

  obj <- list(w = best.w, info = info, fit.obj = best.fit)
  return(obj)
}

#CBPS
weightit2cbps <- function(covs, treat, s.weights, estimand, focal, subset, stabilize, subclass, missing, ...) {

  check.package("CBPS")

  A <- list(...)

  covs <- covs[subset, , drop = FALSE]
  treat <- factor(treat[subset])
  s.weights <- s.weights[subset]

  if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
  treat.type <- get.treat.type(treat)

  if (missing == "ind") {
    missing.ind <- apply(covs[, anyNA_col(covs), drop = FALSE], 2, function(x) as.numeric(is.na(x)))
    if (is_not_null(missing.ind)) {
      covs[is.na(covs)] <- 0
      covs <- cbind(covs, missing.ind)
    }
  }
  for (i in seq_col(covs)) covs[,i] <- make.closer.to.1(covs[,i])

  colinear.covs.to.remove <- colnames(covs)[colnames(covs) %nin% colnames(make_full_rank(covs))]
  covs <- covs[, colnames(covs) %nin% colinear.covs.to.remove, drop = FALSE]

  if (estimand == "ATT") {
    ps <- make_df(levels(treat), length(treat))

    control.levels <- levels(treat)[levels(treat) != focal]
    fit.list <- make_list(control.levels)

    for (i in control.levels) {
      treat.in.i.focal <- which(treat %in% c(focal, i))
      treat_ <- as.integer(treat[treat.in.i.focal] != i)
      covs_ <- covs[treat.in.i.focal, , drop = FALSE]
      new.data <- data.frame(treat_, covs_)

      tryCatch({fit.list[[i]] <- CBPS::CBPS(formula(new.data),
                                            data = new.data,
                                            method = if (is_not_null(A$over) && A$over == FALSE) "exact" else "over",
                                            standardize = FALSE,
                                            sample.weights = s.weights[treat.in.i.focal],
                                            ATT = 1,
                                            ...)},
               error = function(e) {
                 e. <- conditionMessage(e)
                 e. <- gsub("method = \"exact\"", "over = FALSE", e., fixed = TRUE)
                 stop(e., call. = FALSE)
               }

      )

      ps[[focal]][treat.in.i.focal] <- fit.list[[i]][["fitted.values"]]
      ps[[i]][treat.in.i.focal] <- 1 - ps[[focal]][treat.in.i.focal]

    }

  }
  else {
    new.data <- data.frame(treat, covs)
    if (treat.type == "binary" || !nunique.gt(treat, 4)) {
      tryCatch({fit.list <- CBPS::CBPS(formula(new.data),
                                       data = new.data,
                                       method = if (isFALSE(A$over)) "exact" else "over",
                                       standardize = FALSE,
                                       sample.weights = s.weights,
                                       ATT = 0,
                                       ...)},
               error = function(e) {
                 e. <- conditionMessage(e)
                 e. <- gsub("method = \"exact\"", "over = FALSE", e., fixed = TRUE)
                 stop(e., call. = FALSE)
               }
      )

      ps <- fit.list[["fitted.values"]]
    }
    else {
      ps <- rep(NA_real_, length(treat))

      fit.list <- make_list(levels(treat))

      for (i in levels(treat)) {
        new.data[[1]] <- as.integer(treat == i)
        fit.list[[i]] <- CBPS::CBPS(formula(new.data), data = new.data,
                                    method = if (isFALSE(A$over)) "exact" else "over",
                                    standardize = FALSE,
                                    sample.weights = s.weights,
                                    ATT = 0, ...)

        ps[treat==i] <- fit.list[[i]][["fitted.values"]][treat==i]
      }

    }
  }

  w <- get_w_from_ps(ps, treat, estimand = estimand, subclass = subclass,
                     focal = focal, stabilize = stabilize)

  if (treat.type != "binary") {
    p.score <- NULL
  }
  else if (is_not_null(dim(ps)) && length(dim(ps)) == 2) {
    p.score <- ps[[get.treated.level(treat)]]
  }
  else p.score <- ps

  obj <- list(w = w, ps = p.score, fit.obj = fit.list)
  return(obj)
}
weightit2cbps.cont <- function(covs, treat, s.weights, subset, missing, ...) {
  check.package("CBPS")

  A <- list(...)

  covs <- covs[subset, , drop = FALSE]
  treat <- treat[subset]
  s.weights <- s.weights[subset]

  if (missing == "ind") {
    missing.ind <- apply(covs[, anyNA_col(covs), drop = FALSE], 2, function(x) as.numeric(is.na(x)))
    if (is_not_null(missing.ind)) {
      covs[is.na(covs)] <- 0
      covs <- cbind(covs, missing.ind)
    }
  }
  for (i in seq_col(covs)) covs[,i] <- make.closer.to.1(covs[,i])

  colinear.covs.to.remove <- colnames(covs)[colnames(covs) %nin% colnames(make_full_rank(covs))]
  covs <- covs[, colnames(covs) %nin% colinear.covs.to.remove, drop = FALSE]

  new.data <- data.frame(treat = treat, covs)

  tryCatch({fit <- CBPS::CBPS(formula(new.data),
                              data = new.data,
                              method = if (isFALSE(A$over)) "exact" else "over",
                              standardize = FALSE,
                              sample.weights = s.weights,
                              ...)},
           error = function(e) {
             e. <- conditionMessage(e)
             e. <- gsub("method = \"exact\"", "over = FALSE", e., fixed = TRUE)
             stop(e., call. = FALSE)
           }
  )

  w <- fit$weights / s.weights

  obj <- list(w = w, fit.obj = fit)
  return(obj)
}
weightit2cbps.msm <- function(covs.list, treat.list, s.weights, subset, missing, ...) {
  stop("CBMSM doesn't work yet.")
}
weightit2npcbps <- function(covs, treat, s.weights, subset, moments, int, missing, ...) {
  check.package("CBPS")

  A <- list(...)
  if (!all_the_same(s.weights)) stop(paste0("Sampling weights cannot be used with method = \"npcbps\"."),
                                     call. = FALSE)

  covs <- covs[subset, , drop = FALSE]
  treat <- factor(treat[subset])

  if (missing == "ind") {
    missing.ind <- apply(covs[, anyNA_col(covs), drop = FALSE], 2, function(x) as.numeric(is.na(x)))
    if (is_not_null(missing.ind)) {
      covs[is.na(covs)] <- 0
      covs <- cbind(covs, missing.ind)
    }
  }

  covs <- cbind(covs, int.poly.f(covs, poly = moments, int = int))

  for (i in seq_col(covs)) covs[,i] <- make.closer.to.1(covs[,i])

  colinear.covs.to.remove <- colnames(covs)[colnames(covs) %nin% colnames(make_full_rank(covs))]
  covs <- covs[, colnames(covs) %nin% colinear.covs.to.remove, drop = FALSE]

  new.data <- data.frame(treat = treat, covs)

  fit <- do.call(CBPS::npCBPS, c(list(formula(new.data), data = new.data, print.level = 1), A),
                 quote = TRUE)

  w <- fit$weights

  for (i in levels(treat)) w[treat == i] <- w[treat == i]/mean(w[treat == i])

  obj <- list(w = w, fit.obj = fit)

  return(obj)
}
weightit2npcbps.cont <- function(covs, treat, s.weights, subset, moments, int, missing, ...) {
  check.package("CBPS")

  A <- list(...)

  if (!all_the_same(s.weights)) stop(paste0("Sampling weights cannot be used with method = \"npcbps\"."),
                                     call. = FALSE)

  covs <- covs[subset, , drop = FALSE]
  treat <- treat[subset]

  if (missing == "ind") {
    missing.ind <- apply(covs[, anyNA_col(covs), drop = FALSE], 2, function(x) as.numeric(is.na(x)))
    if (is_not_null(missing.ind)) {
      covs[is.na(covs)] <- 0
      covs <- cbind(covs, missing.ind)
    }
  }

  for (i in seq_col(covs)) covs[,i] <- make.closer.to.1(covs[,i])

  covs <- cbind(covs, int.poly.f(covs, poly = moments, int = int))

  colinear.covs.to.remove <- colnames(covs)[colnames(covs) %nin% colnames(make_full_rank(covs))]
  covs <- covs[, colnames(covs) %nin% colinear.covs.to.remove, drop = FALSE]

  new.data <- data.frame(treat = treat, covs)

  fit <- do.call(CBPS::npCBPS, c(list(formula(new.data), data = new.data, print.level = 1), A),
                 quote = TRUE)

  w <- fit$weights

  w <- w/mean(w)

  obj <- list(w = w, fit.obj = fit)

  return(obj)
}

#Entropy balancing
weightit2ebal <- function(covs, treat, s.weights, subset, estimand, focal, stabilize, missing, moments, int, ...) {
  A <- list(...)

  covs <- covs[subset, , drop = FALSE]
  treat <- factor(treat[subset])
  s.weights <- s.weights[subset]

  if (missing == "ind") {
    missing.ind <- apply(covs[, anyNA_col(covs), drop = FALSE], 2, function(x) as.numeric(is.na(x)))
    if (is_not_null(missing.ind)) {
      covs[is.na(covs)] <- 0
      covs <- cbind(covs, missing.ind)
    }
  }

  covs <- cbind(covs, int.poly.f(covs, poly = moments, int = int, center = TRUE))
  for (i in seq_col(covs)) covs[,i] <- make.closer.to.1(covs[,i])

  if (is_not_null(A[["base.weights"]])) A[["base.weight"]] <- A[["base.weights"]]
  if (is_null(A[["base.weight"]])) {
    bw <- rep(1, length(treat))
  }
  else {
    if (!is.numeric(A[["base.weight"]]) || length(A[["base.weight"]]) != length(treat)) {
      stop("The argument to base.weight must be a numeric vector with length equal to the number of units.", call. = FALSE)
    }
    else bw <- A[["base.weight"]]
  }

  eb <- function(C, M, s.weights_t, Q) {
    #X_t : covariates in control group;
    #Returns weights for control group

    n <- nrow(C)

    W <- function(Z) {
      drop(Q * exp(-C %*% Z))
    }

    objective.EB <- function(Z) {
      log(sum(W(Z))) + sum(M * Z)
    }

    gradient.EB <- function(Z) {
      w <- W(Z)
      drop(M - w %*% C/sum(w))
    }

    opt.out <- optim(par = rep(0, ncol(C)),
                     fn = objective.EB,
                     gr = gradient.EB,
                     method = "BFGS",
                     control = list(trace = 0,
                                    reltol = if_null_then(A[["reltol"]], sqrt(.Machine$double.eps)),
                                    maxit = if_null_then(A[["maxit"]], 200)))

    w <- W(opt.out$par)

    list(Z = setNames(opt.out$par, colnames(C)),
         w = w/(mean(w) * s.weights_t),
                opt.out = opt.out)
  }

  w <- rep(1, length(treat))

  if (estimand == "ATT") {
    groups_to_weight <- levels(treat)[levels(treat) != focal]
    targets <- cobalt::col_w_mean(covs, s.weights = s.weights, subset = treat == focal)
  }
  else if (estimand == "ATE") {
    groups_to_weight <- levels(treat)
    targets <- cobalt::col_w_mean(covs, s.weights = s.weights)
  }

  fit.list <- make_list(groups_to_weight)
  for (i in groups_to_weight) {

    fit.list[[i]] <- eb(covs[treat == i,,drop = FALSE], targets,
                        s.weights[treat == i], bw[treat == i])

    w[treat == i] <- fit.list[[i]]$w
  }

  obj <- list(w = w, fit.obj = lapply(fit.list, function(x) x[["opt.out"]]))
  return(obj)
}
weightit2ebal.cont <- function(covs, treat, s.weights, subset, moments, int, missing, ...) {
  A <- list(...)

  covs <- covs[subset, , drop = FALSE]
  treat <- treat[subset]
  s.weights <- s.weights[subset]

  if (missing == "ind") {
    missing.ind <- apply(covs[, anyNA_col(covs), drop = FALSE], 2, function(x) as.numeric(is.na(x)))
    if (is_not_null(missing.ind)) {
      covs[is.na(covs)] <- 0
      covs <- cbind(covs, missing.ind)
    }
  }

  d.moments <- max(if_null_then(A[["d.moments"]], 1), moments)
  k <- ncol(covs)

  poly.covs <- int.poly.f(covs, poly = d.moments)
  int.covs <- int.poly.f(covs, int = int)
  covs <- cbind(covs, poly.covs, int.covs)

  for (i in seq_col(covs)) covs[,i] <- make.closer.to.1(covs[,i])
  # colinear.covs.to.remove <- colnames(covs)[colnames(covs) %nin% colnames(make_full_rank(covs))]
  # covs <- covs[, colnames(covs) %nin% colinear.covs.to.remove, drop = FALSE]

  if (is_not_null(A[["base.weights"]])) A[["base.weight"]] <- A[["base.weights"]]
  if (is_null(A[["base.weight"]])) {
    q <- rep(1, length(treat))
  }
  else {
    if (!is.numeric(A[["base.weight"]]) || length(A[["base.weight"]]) != length(treat)) {
      stop("The argument to base.weight must be a numeric vector with length equal to the number of units.", call. = FALSE)
    }
    else q <- A[["base.weight"]]
  }

  t.mat <- poly(treat, degree = d.moments)

  treat_sc <- mat_div(center(t.mat, at = cobalt::col_w_mean(t.mat, s.weights)),
                     cobalt::col_w_sd(t.mat, s.weights))
  covs_sc <- mat_div(center(covs, at = cobalt::col_w_mean(covs, s.weights)),
                     cobalt::col_w_sd(covs, s.weights))

  kp <- ncol(poly.covs)/(d.moments-1)
  cov_include <- c(seq_len(k),
                   if (moments > 1) k + unlist(lapply(seq_len(moments - 1), function(i) i + (d.moments - 1)*(seq_len(kp)-1))),
                   if (int) seq_col(covs)[-seq_len(k + ncol(poly.covs))])

  gTX <- do.call("cbind", c(list(treat_sc, covs_sc, treat_sc[,1] * covs_sc[,cov_include])))

  #----Code written by Stefan Tubbicke---#
  #define objective function (Lagrange dual)
  objective.EBCT <- function(theta) {
    f <- log(mean(q*s.weights*exp(gTX %*% theta)))*nrow(gTX)
    return(f)
  }

  #define gradient function (LHS of equations 8 in Tubbicke (2020))
  gradient.EBCT<- function(theta) {
    g <- t(gTX) %*% (q*s.weights*exp(gTX %*% theta)/(mean(q*s.weights*exp(gTX %*% theta))))
    return(g)
  }

  opt.out <- optim(par = rep(0, ncol(gTX)),
                   fn = objective.EBCT,
                   gr = gradient.EBCT,
                   method = "BFGS",
                   control = list(trace = TRUE,
                                  reltol = if_null_then(A[["reltol"]], sqrt(.Machine$double.eps)),
                                  maxit = if_null_then(A[["maxit"]], 200)))

  w <-  q*exp(gTX %*% opt.out$par)/(mean(q*exp(gTX %*% opt.out$par)))
  #--------------------------------------#

  obj <- list(w = w, fit.obj = opt.out)

  return(obj)
}

#Empirical Balancing Calibration weights with ATE
weightit2ebcw <- function(covs, treat, s.weights, subset, estimand, focal, missing, moments, int, ...) {
  check.package("ATE")

  A <- list(...)

  covs <- covs[subset, , drop = FALSE]
  treat <- factor(treat[subset])
  s.weights <- s.weights[subset]

  if (missing == "ind") {
    missing.ind <- apply(covs[, anyNA_col(covs), drop = FALSE], 2, function(x) as.numeric(is.na(x)))
    if (is_not_null(missing.ind)) {
      covs[is.na(covs)] <- 0
      covs <- cbind(covs, missing.ind)
    }
  }

  covs <- cbind(covs, int.poly.f(covs, poly = moments, int = int))
  for (i in seq_col(covs)) covs[,i] <- make.closer.to.1(covs[,i])

  for (f in names(formals(ATE::ATE))) {
    if (is_null(A[[f]])) A[[f]] <- formals(ATE::ATE)[[f]]
  }

  if (estimand == "ATT") {
    w <- rep(1, length(treat))
    control.levels <- levels(treat)[levels(treat) != focal]
    fit.list <- make_list(control.levels)

    for (i in control.levels) {
      treat.in.i.focal <- treat %in% c(focal, i)
      treat_ <- as.integer(treat[treat.in.i.focal] != i)
      covs_ <- covs[treat.in.i.focal, , drop = FALSE]

      colinear.covs.to.remove <- colnames(covs_)[colnames(covs_) %nin% colnames(make_full_rank(covs_[treat_ == 0, , drop = FALSE]))]

      covs_ <- covs_[, colnames(covs_) %nin% colinear.covs.to.remove, drop = FALSE]

      covs_[treat_ == 1,] <- covs_[treat_ == 1,] * s.weights[treat == focal] * sum(treat == focal)/ sum(s.weights[treat == focal])

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
      w[treat == i] <- ate.out$weights.q[treat_ == 0] / s.weights[treat == i]
      fit.list[[i]] <- ate.out
    }
  }
  else if (estimand == "ATE") {
    w <- rep(1, length(treat))
    fit.list <- make_list(levels(treat))

    for (i in levels(treat)) {
      covs_i <- rbind(covs, covs[treat==i, , drop = FALSE])
      treat_i <- c(rep(1, nrow(covs)), rep(0, sum(treat==i)))

      colinear.covs.to.remove <- colnames(covs_i)[colnames(covs_i) %nin% colnames(make_full_rank(covs_i[treat_i == 0, , drop = FALSE]))]

      covs_i <- covs_i[, colnames(covs_i) %nin% colinear.covs.to.remove, drop = FALSE]

      covs_i[treat_i == 1,] <- covs_i[treat_i == 1,] * s.weights * sum(treat_i == 1) / sum(s.weights)

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
      w[treat == i] <- ate.out$weights.q[treat_i == 0] / s.weights[treat == i]
      fit.list[[i]] <- ate.out
    }
  }
  if (length(fit.list) == 1) fit.list <- fit.list[[1]]
  obj <- list(w = w, fit.obj = fit.list)
  return(obj)

}

#PS weights using SuperLearner
weightit2super <- function(covs, treat, s.weights, subset, estimand, focal, stabilize, subclass, missing, ...) {
  A <- list(...)

  check.package("SuperLearner")

  covs <- covs[subset, , drop = FALSE]
  treat <- factor(treat[subset])
  s.weights <- s.weights[subset]

  if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
  treat.type <- get.treat.type(treat)

  if (missing == "ind") {
    missing.ind <- apply(covs[, anyNA_col(covs), drop = FALSE], 2, function(x) as.numeric(is.na(x)))
    if (is_not_null(missing.ind)) {
      covs[is.na(covs)] <- 0
      covs <- cbind(covs, missing.ind)
    }
  }
  for (i in seq_col(covs)) covs[,i] <- make.closer.to.1(covs[,i])

  covs <- as.data.frame(covs)

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
  if (length(discrete) != 1 || !is_(discrete, "logical")) stop("'discrete' must be TRUE or FALSE.", call. = FALSE)

  if (identical(A[["SL.method"]], "method.balance")) {
    if (treat.type != "binary") stop("\"method.balance\" cannot be used with multi-category treatments.", call. = FALSE)

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

    available.stop.methods <- bal_criterion(treat.type, list = TRUE)
    s.m.matches <- charmatch(A[["stop.method"]], available.stop.methods)
    if (is.na(s.m.matches) || s.m.matches == 0L) {
      stop(paste0("'stop.method' must be one of ", word_list(available.stop.methods, "or", quotes = TRUE), "."), call. = FALSE)
    }
    else stop.method <- available.stop.methods[s.m.matches]

    crit <- bal_criterion("binary", stop.method)
    init <- crit$init(covs, treat, estimand = estimand, s.weights = s.weights, focal = focal, ...)
    bal_fun <- crit$fun

    sneaky <- 0
    attr(sneaky, "vals") <- list(init = init, bal_fun = bal_fun, estimand = estimand)
    A[["control"]] <- list(trimLogit = sneaky)

    A[["SL.method"]] <- method.balance(stop.method)
  }

  fit.list <- info <- make_list(levels(treat))
  ps <- make_df(levels(treat), nrow = length(treat))

  for (i in levels(treat)) {

    if (treat.type == "binary" && i == last(levels(treat))) {
      ps[[i]] <- 1 - ps[[1]]
      fit.list <- fit.list[[1]]
      info <- info[[i]]
      next
    }

    treat_i <- as.numeric(treat == i)

    fit.list[[i]] <- do.call(SuperLearner::SuperLearner, list(Y = treat_i,
                                                              X = as.data.frame(covs),
                                                              family = binomial(),
                                                              SL.library = A[["SL.library"]],
                                                              verbose = FALSE,
                                                              method = A[["SL.method"]],
                                                              id = NULL,
                                                              obsWeights = s.weights,
                                                              control = A[["control"]],
                                                              cvControl = A[["cvControl"]],
                                                              env = A[["env"]]))
    if (discrete) ps[[i]] <- fit.list[[i]]$library.predict[,which.min(fit.list[[i]]$cvRisk)]
    else ps[[i]] <- fit.list[[i]]$SL.predict

    info[[i]] <- list(coef = fit.list[[i]]$coef,
                      cvRisk = fit.list[[i]]$cvRisk)

  }

  #ps should be matrix of probs for each treat
  #Computing weights
  w <- get_w_from_ps(ps = ps, treat = treat, estimand, focal, stabilize = stabilize, subclass = subclass)

  p.score <- if (treat.type == "binary") ps[[get.treated.level(treat)]] else NULL

  obj <- list(w = w, ps = p.score, info = info, fit.obj = fit.list)
  return(obj)
}
weightit2super.cont <- function(covs, treat, s.weights, subset, stabilize, missing, ps, ...) {
  A <- B <- list(...)

  covs <- covs[subset, , drop = FALSE]
  treat <- treat[subset]
  s.weights <- s.weights[subset]

  if (missing == "ind") {
    missing.ind <- apply(covs[, anyNA_col(covs), drop = FALSE], 2, function(x) as.numeric(is.na(x)))
    if (is_not_null(missing.ind)) {
      covs[is.na(covs)] <- 0
      covs <- cbind(covs, missing.ind)
    }
  }

  for (i in seq_col(covs)) covs[,i] <- make.closer.to.1(covs[,i])

  if (ncol(covs) > 1) {
    colinear.covs.to.remove <- colnames(covs)[colnames(covs) %nin% colnames(make_full_rank(covs))]
    covs <- covs[, colnames(covs) %nin% colinear.covs.to.remove, drop = FALSE]
  }

  #Process density params
  if (isTRUE(A[["use.kernel"]])) {
    if (is_null(A[["bw"]])) A[["bw"]] <- "nrd0"
    if (is_null(A[["adjust"]])) A[["adjust"]] <- 1
    if (is_null(A[["kernel"]])) A[["kernel"]] <- "gaussian"
    if (is_null(A[["n"]])) A[["n"]] <- 10*length(treat)
    use.kernel <- TRUE
    densfun <- NULL
  }
  else {
    if (is_null(A[["density"]])) densfun <- dnorm
    else if (is.function(A[["density"]])) densfun <- A[["density"]]
    else if (is.character(A[["density"]]) && length(A[["density"]] == 1)) {
      splitdens <- strsplit(A[["density"]], "_", fixed = TRUE)[[1]]
      if (exists(splitdens[1], mode = "function", envir = parent.frame())) {
        if (length(splitdens) > 1 && !can_str2num(splitdens[-1])) {
          stop(paste(A[["density"]], "is not an appropriate argument to 'density' because",
                     word_list(splitdens[-1], and.or = "or", quotes = TRUE), "cannot be coerced to numeric."), call. = FALSE)
        }
        densfun <- function(x) {
          tryCatch(do.call(get(splitdens[1]), c(list(x), as.list(str2num(splitdens[-1])))),
                   error = function(e) stop(paste0("Error in applying density:\n  ", conditionMessage(e)), call. = FALSE))
        }
      }
      else {
        stop(paste(A[["density"]], "is not an appropriate argument to 'density' because",
                   splitdens[1], "is not an available function."), call. = FALSE)
      }
    }
    else stop("The argument to 'density' cannot be evaluated as a density function.", call. = FALSE)
    use.kernel <- FALSE
  }

  #Stabilization - get dens.num
  p.num <- treat - mean(treat)

  if (use.kernel) {
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

  #Estimate GPS
  for (f in names(formals(SuperLearner::SuperLearner))) {
    if (f == "method") {if (is_null(B[["SL.method"]])) B[["SL.method"]] <- formals(SuperLearner::SuperLearner)[["method"]]}
    else if (f == "env") {if (is_null(B[["env"]])) B[["env"]] <- environment(SuperLearner::SuperLearner)}
    else if (is_null(B[[f]])) B[[f]] <- formals(SuperLearner::SuperLearner)[[f]]
  }

  discrete <- if_null_then(A[["discrete"]], FALSE)
  if (length(discrete) != 1 || !is_(discrete, "logical")) stop("discrete must be TRUE or FALSE.", call. = FALSE)

  if (identical(B[["SL.method"]], "method.balance")) {

    if (is_null(B[["stop.method"]])) {
      warning("No stop.method was provided. Using \"p.mean\".",
              call. = FALSE, immediate. = TRUE)
      B[["stop.method"]] <- "p.mean"
    }
    else if (length(B[["stop.method"]]) > 1) {
      warning("Only one stop.method is allowed at a time. Using just the first stop.method.",
              call. = FALSE, immediate. = TRUE)
      B[["stop.method"]] <- B[["stop.method"]][1]
    }

    available.stop.methods <- bal_criterion("continuous", list = TRUE)
    s.m.matches <- charmatch(B[["stop.method"]], available.stop.methods)
    if (is.na(s.m.matches) || s.m.matches == 0L) {
      stop(paste0("'stop.method' must be one of ", word_list(available.stop.methods, "or", quotes = TRUE), "."), call. = FALSE)
    }
    else stop.method <- available.stop.methods[s.m.matches]

    crit <- bal_criterion("continuous", stop.method)
    init <- crit$init(covs, treat, s.weights = s.weights, ...)
    bal_fun <- crit$fun

    sneaky <- 0
    attr(sneaky, "vals") <- list(init = init,
                                 bal_fun = bal_fun,
                                 dens.num = dens.num,
                                 densfun = densfun,
                                 use.kernel = use.kernel,
                                 densControl = A)
    B[["control"]] <- list(trimLogit = sneaky)

    B[["SL.method"]] <- method.balance.cont(stop.method)
  }

  fit <- do.call(SuperLearner::SuperLearner, list(Y = treat,
                                                  X = as.data.frame(covs),
                                                  family = gaussian(),
                                                  SL.library = B[["SL.library"]],
                                                  verbose = FALSE,
                                                  method = B[["SL.method"]],
                                                  id = NULL,
                                                  obsWeights = s.weights,
                                                  control = B[["control"]],
                                                  cvControl = B[["cvControl"]],
                                                  env = B[["env"]]))

  if (discrete) gp.score <- fit$library.predict[,which.min(fit$cvRisk)]
  else gp.score <- fit$SL.predict

  #Get weights
  w <- get_cont_weights(gp.score, treat = treat, s.weights = s.weights,
                        dens.num = dens.num, densfun = densfun,
                        use.kernel = use.kernel, densControl = A)

  if (use.kernel && isTRUE(A[["plot"]])) {
    d.d <- density(treat - gp.score, n = A[["n"]],
                   weights = s.weights/sum(s.weights), give.Rkern = FALSE,
                   bw = A[["bw"]], adjust = A[["adjust"]],
                   kernel = A[["kernel"]])
    plot_density(d.n, d.d)
  }

  info <- list(coef = fit$coef,
               cvRisk = fit$cvRisk)

  obj <- list(w = w, info = info, fit.obj = fit)
  return(obj)
}

#PS weights using BART
weightit2bart <- function(covs, treat, s.weights, subset, estimand, focal, stabilize, subclass, missing, ...) {
  A <- list(...)

  check.package("dbarts")

  covs <- covs[subset, , drop = FALSE]
  treat <- factor(treat[subset])
  s.weights <- s.weights[subset]

  if (!all_the_same(s.weights)) stop("Sampling weights cannot be used with method = \"bart\".",
                                     call. = FALSE)

  if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
  treat.type <- get.treat.type(treat)

  if (missing == "ind") {
    missing.ind <- apply(covs[, anyNA_col(covs), drop = FALSE], 2, function(x) as.numeric(is.na(x)))
    if (is_not_null(missing.ind)) {
      covs[is.na(covs)] <- 0
      covs <- cbind(covs, missing.ind)
    }
  }

  if (ncol(covs) > 1) {
    colinear.covs.to.remove <- colnames(covs)[colnames(covs) %nin% colnames(make_full_rank(covs))]
    covs <- covs[, colnames(covs) %nin% colinear.covs.to.remove, drop = FALSE]
  }

  for (i in seq_col(covs)) covs[,i] <- make.closer.to.1(covs[,i])

  ps <- make_df(levels(treat), nrow = length(treat))

  A[["formula"]] <- covs
  A[["keepCall"]] <- FALSE
  A[["combineChains"]] <- TRUE
  A[["verbose"]] <- FALSE #necessary to prevent crash

  fit.list <- make_list(levels(treat))

  for (i in levels(treat)) {

    if (treat.type == "binary" && i == last(levels(treat))) {
      ps[[i]] <- 1 - ps[[1]]
      fit.list <- fit.list[[1]]
      next
    }

    A[["data"]] <- as.integer(treat == i)
    fit.list[[i]] <- do.call(dbarts::bart2, A[names(A) %in% setdiff(c(names(formals(dbarts::bart2)),
                                                                      names(formals(dbarts::dbartsControl))),
                                                                    c("offset.test", "weights", "subset", "test"))],
                             quote = TRUE)

    ps[[i]] <- fitted(fit.list[[i]])

  }

  info <- list()

  #ps should be matrix of probs for each treat
  #Computing weights
  w <- get_w_from_ps(ps = ps, treat = treat, estimand, focal, stabilize = stabilize, subclass = subclass)

  p.score <- if (treat.type == "binary") ps[[get.treated.level(treat)]] else NULL

  obj <- list(w = w, ps = p.score, info = info, fit.obj = fit.list)
  return(obj)
}
weightit2bart.cont <- function(covs, treat, s.weights, subset, stabilize, missing, ps, ...) {
  A <- list(...)

  check.package("dbarts")

  covs <- covs[subset, , drop = FALSE]
  treat <- treat[subset]
  s.weights <- s.weights[subset]

  if (!all_the_same(s.weights)) stop("Sampling weights cannot be used with method = \"bart\".",
                                     call. = FALSE)

  if (missing == "ind") {
    missing.ind <- apply(covs[, anyNA_col(covs), drop = FALSE], 2, function(x) as.numeric(is.na(x)))
    if (is_not_null(missing.ind)) {
      covs[is.na(covs)] <- 0
      covs <- cbind(covs, missing.ind)
    }
  }

  for (i in seq_col(covs)) covs[,i] <- make.closer.to.1(covs[,i])

  #Process density params
  if (isTRUE(A[["use.kernel"]])) {
    if (is_null(A[["bw"]])) A[["bw"]] <- "nrd0"
    if (is_null(A[["adjust"]])) A[["adjust"]] <- 1
    if (is_null(A[["kernel"]])) A[["kernel"]] <- "gaussian"
    if (is_null(A[["n"]])) A[["n"]] <- 10*length(treat)
    use.kernel <- TRUE
    densfun <- NULL
  }
  else {
    if (is_null(A[["density"]])) densfun <- dnorm
    else if (is.function(A[["density"]])) densfun <- A[["density"]]
    else if (is.character(A[["density"]]) && length(A[["density"]] == 1)) {
      splitdens <- strsplit(A[["density"]], "_", fixed = TRUE)[[1]]
      if (exists(splitdens[1], mode = "function", envir = parent.frame())) {
        if (length(splitdens) > 1 && !can_str2num(splitdens[-1])) {
          stop(paste(A[["density"]], "is not an appropriate argument to 'density' because",
                     word_list(splitdens[-1], and.or = "or", quotes = TRUE), "cannot be coerced to numeric."), call. = FALSE)
        }
        densfun <- function(x) {
          tryCatch(do.call(get(splitdens[1]), c(list(x), as.list(str2num(splitdens[-1])))),
                   error = function(e) stop(paste0("Error in applying density:\n  ", conditionMessage(e)), call. = FALSE))
        }
      }
      else {
        stop(paste(A[["density"]], "is not an appropriate argument to 'density' because",
                   splitdens[1], "is not an available function."), call. = FALSE)
      }
    }
    else stop("The argument to 'density' cannot be evaluated as a density function.", call. = FALSE)
    use.kernel <- FALSE
  }

  #Stabilization - get dens.num
  p.num <- treat - mean(treat)

  if (use.kernel) {
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

  A[["formula"]] <- covs
  A[["data"]] <- treat
  A[["keepCall"]] <- FALSE
  A[["combineChains"]] <- TRUE
  A[["verbose"]] <- FALSE #necessary to prevent crash

  #Estimate GPS

  fit <- do.call(dbarts::bart2, A[names(A) %in% setdiff(c(names(formals(dbarts::bart2)),
                                                          names(formals(dbarts::dbartsControl))),
                                                        c("offset.test", "weights", "subset", "test"))],
                 quote = TRUE)

  gp.score <- fitted(fit)

  #Get weights
  w <- get_cont_weights(gp.score, treat = treat, s.weights = s.weights,
                        dens.num = dens.num, densfun = densfun,
                        use.kernel = use.kernel, densControl = A)

  if (use.kernel && isTRUE(A[["plot"]])) {
    d.d <- density(treat - gp.score, n = A[["n"]],
                   weights = s.weights/sum(s.weights), give.Rkern = FALSE,
                   bw = A[["bw"]], adjust = A[["adjust"]],
                   kernel = A[["kernel"]])
    plot_density(d.n, d.d)
  }

  info <- list()

  obj <- list(w = w, info = info, fit.obj = fit)
  return(obj)
}

#Energy balancing
weightit2energy <- function(covs, treat, s.weights, subset, estimand, focal, missing, moments, int, ...) {
  check.package("osqp")

  A <- list(...)

  n <- length(treat)
  covs <- covs[subset, , drop = FALSE]
  treat <- factor(treat[subset])
  s.weights <- s.weights[subset]

  if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
  treat.type <- get.treat.type(treat)

  if (missing == "ind") {
    missing.ind <- apply(covs[, anyNA_col(covs), drop = FALSE], 2, function(x) as.numeric(is.na(x)))
    if (is_not_null(missing.ind)) {
      covs[is.na(covs)] <- 0
      covs <- cbind(covs, missing.ind)
    }
  }

  covs <- mat_div(center(covs, at = col.w.m(covs, s.weights)),
                  sqrt(col.w.v(covs, s.weights)))

  if (is_not_null(A[["dist.mat"]])) {
    if (inherits(A[["dist.mat"]], "dist")) A[["dist.mat"]] <- as.matrix(A[["dist.mat"]])

    if (is.matrix(A[["dist.mat"]]) && all(dim(A[["dist.mat"]]) == n) &&
        all(check_if_zero(diag(A[["dist.mat"]]))) && !any(A[["dist.mat"]] < 0) &&
        isSymmetric(unname(A[["dist.mat"]]))) {
      d <- unname(A[["dist.mat"]][subset, subset])
    }
    else stop("'dist.mat' must be a square, symmetric distance matrix with a value for all pairs of units.", call. = FALSE)
  }
  else d <- as.matrix(dist(covs))

  n <- length(treat)
  levels_treat <- levels(treat)
  diagn <- diag(n)

  min.w <- if_null_then(A[["min.w"]], 1e-8)
  if (!is.numeric(min.w) || length(min.w) != 1 || min.w < 0) {
    warning("'min.w' must be a nonnegative number. Setting min.w = 1e-8.", call. = FALSE, immediate. = TRUE)
    min.w <- 1e-8
  }

  for (t in levels_treat) s.weights[treat == t] <- s.weights[treat == t]/mean(s.weights[treat == t])

  tmat <- vapply(levels_treat, function(t) treat == t, logical(n))
  nt <- colSums(tmat)

  J <- setNames(lapply(levels_treat, function(t) s.weights*tmat[,t]/nt[t]), levels_treat)

  if (estimand == "ATE") {
    J0 <- as.matrix(s.weights/n)

    M2_array <- vapply(levels_treat, function(t) -2 * tcrossprod(J[[t]]) * d, diagn)
    M1_array <- vapply(levels_treat, function(t) 2 * J[[t]] * d %*% J0, J0)

    M2 <- rowSums(M2_array, dims = 2)
    M1 <- rowSums(M1_array)

    if (!isFALSE(A[["improved"]])) {
      all_pairs <- combn(levels_treat, 2, simplify = FALSE)
      M2_pairs_array <- vapply(all_pairs, function(p) -2 * tcrossprod(J[[p[1]]]-J[[p[2]]]) * d, diagn)
      M2 <- M2 + rowSums(M2_pairs_array, dims = 2)
    }

    #Constraints for positivity and sum of weights
    Amat <- rbind(diagn, t(s.weights * tmat))
    lvec <- c(rep(min.w, n), nt)
    uvec <- c(ifelse(check_if_zero(s.weights), min.w, Inf), nt)
  }
  else {

    J0_focal <- as.matrix(J[[focal]])
    clevs <- levels_treat[levels_treat != focal]

    M2_array <- vapply(clevs, function(t) -2 * tcrossprod(J[[t]]) * d, diagn)
    M1_array <- vapply(clevs, function(t) 2 * J[[t]] * d %*% J0_focal, J0_focal)

    M2 <- rowSums(M2_array, dims = 2)
    M1 <- rowSums(M1_array)

    #Constraints for positivity and sum of weights
    Amat <- rbind(diagn, t(s.weights*tmat))
    lvec <- c(ifelse_(check_if_zero(s.weights), min.w, treat == focal, 1, min.w), nt)
    uvec <- c(ifelse_(check_if_zero(s.weights), min.w, treat == focal, 1, Inf), nt)
  }

  #Add weight penalty
  if (is_not_null(A[["lambda"]])) diag(M2) <- diag(M2) + A[["lambda"]] / n

  if (moments != 0 || int) {
    #Exactly balance moments and/or interactions
    covs <- cbind(covs, int.poly.f(covs, poly = moments, int = int))

    if (estimand == "ATE") targets <- col.w.m(covs, s.weights)
    else targets <- col.w.m(covs[treat == focal, , drop = FALSE], s.weights[treat == focal])

    Amat <- do.call("rbind", c(list(Amat),
                               lapply(levels_treat, function(t) {
                                 if (is_null(focal) || t != focal) t(covs * J[[t]])
                               })))
    lvec <- do.call("c", c(list(lvec),
                           lapply(levels_treat, function(t) {
                             if (is_null(focal) || t != focal) targets
                           })))
    uvec <- do.call("c", c(list(uvec),
                           lapply(levels_treat, function(t) {
                             if (is_null(focal) || t != focal) targets
                           })))
  }

  if (is_not_null(A[["eps"]])) {
    if (is_null(A[["eps_abs"]])) A[["eps_abs"]] <- A[["eps"]]
    if (is_null(A[["eps_rel"]])) A[["eps_rel"]] <- A[["eps"]]
  }
  A[names(A) %nin% names(formals(osqp::osqpSettings))] <- NULL
  if (is_null(A[["max_iter"]])) A[["max_iter"]] <- 2E3L
  if (is_null(A[["eps_abs"]])) A[["eps_abs"]] <- 1E-8
  if (is_null(A[["eps_rel"]])) A[["eps_rel"]] <- 1E-8
  A[["verbose"]] <- TRUE

  options.list <- do.call(osqp::osqpSettings, A)

  opt.out <- do.call(osqp::solve_osqp, list(P = M2, q = M1, A = Amat, l = lvec, u = uvec,
                                            pars = options.list),
                     quote = TRUE)

  if (identical(opt.out$info$status, "maximum iterations reached")) {
    warning("The optimization failed to converge. See Notes section at ?method_energy for information.", call. = FALSE)
  }

  w <- opt.out$x

  if (estimand == "ATT") w[treat == focal] <- 1

  w[w <= min.w] <- min.w

  obj <- list(w = w, fit.obj = opt.out)
  return(obj)

}
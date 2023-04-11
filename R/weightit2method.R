#User-defined weighting function
weightit2user <- function(Fun, covs, treat, s.weights, subset, estimand, focal, stabilize, subclass, missing, ps, moments, int, verbose, ...) {
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

  obj <- do.call(Fun, fun_args)

  if (is.numeric(obj)) {
    obj <- list(w = obj)
  }
  else if (!is.list(obj) || !any(c("w", "weights") %nin% names(obj))) {
    .err("the output of the user-provided function must be a list with an entry named \"w\" or \"weights\" containing the estimated weights")
  }
  else {
    names(obj)[names(obj) == "weights"] <- "w"
  }

  if (is_null(obj[["w"]])) .err("no weights were estimated")
  if (!is.vector(obj[["w"]], mode = "numeric")) .err("the \"w\" or \"weights\" entry of the output of the user-provided function must be a numeric vector of weights")
  if (all(is.na(obj[["w"]]))) .err("all weights were generated as `NA`")
  if (length(obj[["w"]]) != length(treat)) {
    .err(sprintf("%s weights were estimated, but there are %s units",
                 length(obj[["w"]]), length(treat)))
  }

  obj
}
weightitMSM2user <- function(Fun, covs.list, treat.list, s.weights, subset, stabilize, missing, moments, int, verbose, ...) {
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

  obj <- do.call(Fun, fun_args)

  if (is.numeric(obj)) {
    obj <- list(w = obj)
  }
  else if (!is.list(obj) || !any(c("w", "weights") %nin% names(obj))) {
    .err("the output of the user-provided function must be a list with an entry named \"w\" or \"weights\" containing the estimated weights")
  }
  else {
    names(obj)[names(obj) == "weights"] <- "w"
  }
  if (is_null(obj[["w"]])) .err("no weights were estimated")
  if (!is.vector(obj[["w"]], mode = "numeric")) .err("the \"w\" or \"weights\" entry of the output of the user-provided function must be a numeric vector of weights")
  if (all(is.na(obj[["w"]]))) .err("All weights were generated as `NA`")
  if (length(obj[["w"]]) != length(treat.list[[1]])) {
    .err(sprintf("%s weights were estimated, but there are %s units",
                 length(obj[["w"]]), length(treat.list[[1]])))
  }

  obj
}

#Propensity score estimation with regression
weightit2glm <- function(covs, treat, s.weights, subset, estimand, focal, stabilize, subclass, missing, ps, .data, verbose, ...) {
  A <- list(...)

  fit.obj <- NULL

  if (is_null(ps)) {

    covs <- covs[subset, , drop = FALSE]
    treat_sub <- factor(treat[subset])
    s.weights <- s.weights[subset]

    bin.treat <- is_binary(treat_sub)
    ord.treat <- is.ordered(treat_sub)

    if (missing == "ind") {
      covs <- add_missing_indicators(covs)
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

    #Process link
    if (ord.treat) acceptable.links <- c("logit", "probit", "loglog", "cloglog", "cauchit", "br.logit")
    else if (bin.treat || isFALSE(A$use.mlogit)) {
      if (missing == "saem") acceptable.links <- "logit"
      else acceptable.links <- expand.grid_string(c("", "br."), c("logit", "probit", "cloglog", "identity", "log", "cauchit"))
    }
    else acceptable.links <- c("logit", "probit", "bayes.probit", "br.logit")

    if (is_null(A[["link"]])) A$link <- acceptable.links[1]
    else {
      which.link <- acceptable.links[pmatch(A[["link"]], acceptable.links, nomatch = 0)][1]
      if (is.na(which.link)) {
        .err(sprintf("Only %s allowed as the link for %s treatments%",
                     word_list(acceptable.links, quotes = TRUE, is.are = TRUE),
                     if (bin.treat) "binary" else if (ord.treat) "ordinal" else "multinomial",
                     if (missing == "saem") ' with `missing = "saem"`' else ""))
      }

      A[["link"]] <- which.link
    }

    use.br <- startsWith(A[["link"]], "br.")
    # use.bayes <- startsWith(A$link, "bayes.")
    if (use.br) A$link <- substr(A$link, 4, nchar(A$link))
    # else if (use.bayes) A$link <- substr(A$link, 7, nchar(A$link))

    if (bin.treat) {

      t.lev <- get.treated.level(treat_sub)
      c.lev <- setdiff(levels(treat_sub), t.lev)

      ps <- make_df(levels(treat_sub), nrow = length(treat_sub))

      treat <- binarize(treat_sub, one = t.lev)

      if (missing == "saem") {
        rlang::check_installed("misaem")

        data <- data.frame(treat, covs)

        withCallingHandlers({
          verbosely({
            fit <- misaem::miss.glm(formula(data), data = data, control = as.list(A[["control"]]))
          }, verbose = verbose)
        },
        warning = function(w) {
          if (conditionMessage(w) != "one argument not used by format '%i '") .wrn("(from misaem) ", w, tidy = FALSE)
          invokeRestart("muffleWarning")
        })

        if (is_null(A[["saem.method"]])) A[["saem.method"]] <- "map"

        ps[[t.lev]] <- p.score <- drop(predict(fit, newdata = covs, method = A[["saem.method"]]))
        ps[[c.lev]] <- 1 - ps[[t.lev]]
      }
      else {
        if (use.br) {
          rlang::check_installed("brglm2")
          ctrl_fun <- brglm2::brglmControl
          glm_method <- brglm2::brglmFit
          family <- binomial(link = A[["link"]])
        }
        else {
          ctrl_fun <- stats::glm.control
          glm_method <- if_null_then(A[["glm.method"]], stats::glm.fit)
          family <- quasibinomial(link = A[["link"]])
        }

        control <- do.call(ctrl_fun, c(A[["control"]],
                                       A[setdiff(names(formals(ctrl_fun))[pmatch(names(A), names(formals(ctrl_fun)), 0)],
                                                 names(A[["control"]]))]))

        start <- mustart <- NULL

        if (family$link %in% c("log", "identity")) {
          #Need starting values because links are unbounded
          start <- c(family$linkfun(w.m(treat, s.weights)), rep(0, ncol(covs)))
        }
        else {
          #Default starting values from glm.fit() without weights; these
          #work better with s.weights than usual default.
          mustart <- .25 + .5*treat
        }

        withCallingHandlers({verbosely({
          if (isTRUE(A[["quick"]])) {
            fit <- do.call(glm_method, list(y = treat,
                                            x = cbind(`(Intercept)` = 1, covs),
                                            mustart = mustart,
                                            start = start,
                                            weights = s.weights,
                                            family = family,
                                            control = control), quote = TRUE)
          }
          else {
            data <- data.frame(treat, covs)
            formula <- formula(data)
            fit <- do.call(stats::glm, list(formula, data = data,
                                            weights = s.weights,
                                            mustart = mustart,
                                            start = start,
                                            family = family,
                                            method = glm_method,
                                            control = control), quote = TRUE)
          }
        }, verbose = verbose)},
        warning = function(w) {
          if (conditionMessage(w) != "non-integer #successes in a binomial glm!") .wrn("(from `glm()`) ", w, tidy = FALSE)
          invokeRestart("muffleWarning")
        })

        ps[[t.lev]] <- p.score <- fit$fitted.values
        ps[[c.lev]] <- 1 - ps[[t.lev]]
      }

      fit[["call"]] <- NULL
      fit.obj <- fit
    }
    else if (ord.treat) {
      if (use.br) {
        rlang::check_installed("brglm2")

        ctrl_fun <- brglm2::brglmControl
        control <- do.call(ctrl_fun, c(A[["control"]],
                                       A[setdiff(names(formals(ctrl_fun))[pmatch(names(A), names(formals(ctrl_fun)), 0)],
                                                 names(A[["control"]]))]))

        data <- data.frame(treat_sub, covs)
        formula <- formula(data)

        tryCatch({verbosely({
          fit <- do.call(brglm2::bracl, list(formula,
                                             data = data,
                                             weights = s.weights,
                                             control = control,
                                             parallel = if_null_then(A[["parallel"]], FALSE)),
                         quote = TRUE)
        }, verbose = verbose)},
        error = function(e) {
          .err(sprintf("there was a problem with the bias-reduced ordinal logit regression.\n       Try a different link.\n       Error message: %s", conditionMessage(e)), tidy = FALSE)
        })
      }
      else {
        if (A[["link"]] == "logit") A[["link"]] <- "logistic"
        rlang::check_installed("MASS")
        # message(paste("Using ordinal", A$link, "regression."))
        data <- data.frame(treat_sub, covs)
        formula <- formula(data)

        tryCatch({verbosely({
          fit <- do.call(MASS::polr,
                         list(formula,
                              data = data,
                              weights = s.weights,
                              Hess = FALSE,
                              model = FALSE,
                              method = A[["link"]],
                              contrasts = NULL), quote = TRUE)
        }, verbose = verbose)},
        error = function(e) {
          .err(sprintf("There was a problem fitting the ordinal %s regressions with `polr()`.\n       Try again with an un-ordered treatment.\n       Error message: (from `MASS::polr()`) %s",
                       A$link, conditionMessage(e)), tidy = FALSE)})
      }

      ps <- fit$fitted.values
      fit.obj <- fit
      p.score <- NULL
    }
    else {
      if (use.br) {
        rlang::check_installed("brglm2")
        data <- data.frame(treat_sub, covs)
        formula <- formula(data)
        ctrl_fun <- brglm2::brglmControl
        control <- do.call(ctrl_fun, c(A[["control"]],
                                       A[setdiff(names(formals(ctrl_fun))[pmatch(names(A), names(formals(ctrl_fun)), 0)],
                                                 names(A[["control"]]))]))
        tryCatch({verbosely({
          fit <- do.call(brglm2::brmultinom,
                         list(formula, data,
                              weights = s.weights,
                              control = control), quote = TRUE)
        }, verbose = verbose)},
        error = function(e) {
          .err(sprintf("There was a problem with the bias-reduced multinomial logit regression. Try a different link.\n       Error message: (from brglm2) %s", conditionMessage(e)), tidy = FALSE)
        })

        ps <- fit$fitted.values
        fit.obj <- fit
      }
      else if (A$link %in% c("logit", "probit")) {
        if (isTRUE(A$use.mclogit)) {
          rlang::check_installed("mclogit")

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

          control <- do.call(ctrl_fun, c(A[["control"]],
                                         A[setdiff(names(formals(ctrl_fun))[pmatch(names(A), names(formals(ctrl_fun)), 0)],
                                                   names(A[["control"]]))]))
          tryCatch({verbose({
            fit <- do.call(mclogit::mblogit, list(form,
                                                  data = data,
                                                  weights = quote(.s.weights),
                                                  random = A[["random"]],
                                                  method = A[["mclogit.method"]],
                                                  estimator = if_null_then(A[["estimator"]], eval(formals(mclogit::mclogit)[["estimator"]])),
                                                  dispersion = if_null_then(A[["dispersion"]], eval(formals(mclogit::mclogit)[["dispersion"]])),
                                                  groups = A[["groups"]],
                                                  control = control))

          }, verbose = verbose)},
          error = function(e) {
            .err(sprintf("there was a problem fitting the multinomial %s regression with `mblogit()`.\n       Try again with `use.mclogit = FALSE`.\nError message: (from `mclogit::mblogit()`) %s",
                         A$link, conditionMessage(e)), tidy = FALSE)
          }
          )

          ps <- fitted(fit)
          colnames(ps) <- levels(treat_sub)
          fit.obj <- fit
        }
        else if (!isFALSE(A$use.mlogit)) {
          rlang::check_installed("mlogit")

          data <- data.frame(treat = treat_sub, .s.weights = s.weights, covs)
          covnames <- names(data)[-c(1,2)]
          tryCatch({verbosely({
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
          }, verbose = verbose)},
          error = function(e) {
            .err(sprintf("There was a problem fitting the multinomial %s regression with `mlogit()`.\n       Try again with `use.mlogit = FALSE`.\nError message: (from `mlogit::mlogit()`) %s",
                         A$link, conditionMessage(e)), tidy = FALSE)
          }
          )

          ps <- fitted(fit, outcome = FALSE)
          fit.obj <- fit
        }
        else {
          ps <- make_df(levels(treat_sub), nrow = length(treat_sub))

          ctrl_fun <- stats::glm.control
          control <- do.call(ctrl_fun, c(A[["control"]],
                                         A[setdiff(names(formals(ctrl_fun))[pmatch(names(A), names(formals(ctrl_fun)), 0)],
                                                   names(A[["control"]]))]))
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
            verbosely({
              fit.list[[i]] <- do.call(stats::glm, list(formula(data_i), data = data_i,
                                                        family = quasibinomial(link = A$link),
                                                        weights = s.weights,
                                                        control = control), quote = TRUE)
            }, verbose = verbose)
            ps[[i]] <- fit.list[[i]]$fitted.values
          }
          if (isTRUE(A[["test2"]])) ps <- ps/rowSums(ps)
          fit.obj <- fit.list
        }
      }
      else if (A$link == "bayes.probit") {
        rlang::check_installed("MNP")
        data <- data.frame(treat_sub, covs)
        formula <- formula(data)
        tryCatch({verbosely({
          fit <- MNP::mnp(formula, data, verbose = TRUE)
        }, verbose = verbose)},
        error = function(e) {
          .err(sprintf("There was a problem fitting the Bayes probit regression with `MNP()`.\n       Try a different link.\nError message: (from `MNP::MNP()`) %s",
                       A$link, conditionMessage(e)), tidy = FALSE)
        })
        ps <- predict(fit, type = "prob")$p
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

      if (is_null(p.score)) .err("`ps` must be a numeric vector with a propensity score for each unit")

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
            p_
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
            p_
          })), levels(treat))[subset, , drop = FALSE]
        }
        else {
          bad.ps <- TRUE
        }
      }
      else bad.ps <- TRUE

      if (bad.ps) .err("`ps` must be a numeric vector with a propensity score for each unit or a matrix \n\twith the probability of being in each treatment for each unit")

    }

  }

  #ps should be matrix of probs for each treat
  #Computing weights
  w <- get_w_from_ps(ps = ps, treat = treat_sub, estimand, focal = focal,
                     stabilize = stabilize, subclass = subclass)

  list(w = w, ps = p.score, fit.obj = fit.obj)
}
weightit2glm.cont <- function(covs, treat, s.weights, subset, stabilize, missing, ps, verbose, ...) {
  A <- list(...)

  fit.obj <- NULL

  covs <- covs[subset, , drop = FALSE]
  treat <- treat[subset]
  s.weights <- s.weights[subset]

  if (missing == "ind") {
    covs <- add_missing_indicators(covs)
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
          .err(sprintf("%s is not an appropriate argument to `density` because %s cannot be coerced to numeric",
                       A[["density"]], word_list(splitdens[-1], and.or = "or", quotes = TRUE)))
        }
        densfun <- function(x) {
          tryCatch(do.call(splitdens1, c(list(x), as.list(str2num(splitdens[-1])))),
                   error = function(e) .err(sprintf("Error in applying density:\n  %s", conditionMessage(e)), tidy = FALSE))
        }
      }
      else {
        .err(sprintf("%s is not an appropriate argument to `density` because %s is not an available function",
                     A[["density"]], splitdens[1]))
      }
    }
    else .err("the argument to `density` cannot be evaluated as a density function")
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
      .err("there was a problem with the output of `density`. Try another density function or leave it blank to use the normal density")
    }
    if (any(dens.num <= 0)) {
      .err("the input to density may not accept the full range of treatment values")
    }
  }

  #Estimate GPS
  if (is_null(ps)) {

    if (missing == "saem") {
      rlang::check_installed("misaem")

      withCallingHandlers({verbosely({
        fit <- misaem::miss.lm(formula, data, control = as.list(A[["control"]]))
      }, verbose = verbose)},
      warning = function(w) {
        if (conditionMessage(w) != "one argument not used by format '%i '") {
          .wrn("(from `misaem::miss.lm()`) ", w, tidy = FALSE)
        }
        invokeRestart("muffleWarning")
      })

      if (is_null(A[["saem.method"]])) A[["saem.method"]] <- "map"

      gp.score <- drop(predict(fit, newdata = covs, method = A[["saem.method"]]))
    }
    else {
      if (is_null(A[["link"]])) A[["link"]] <- "identity"
      else {
        if (missing == "saem") acceptable.links <- "identity"
        else acceptable.links <- c("identity", "log", "inverse")

        which.link <- acceptable.links[pmatch(A[["link"]], acceptable.links, nomatch = 0)][1]
        if (is.na(which.link)) {
          .err(sprintf("Only %s allowed as the link for continuous treatments%s",
                       word_list(acceptable.links, quotes = TRUE, is.are = TRUE),
                       if (missing == "saem") ' with missing = "saem"' else ""))
        }

        A[["link"]] <- which.link
      }

      verbosely({
        fit <- do.call("glm", c(list(formula, data = data,
                                     weights = s.weights,
                                     family = gaussian(link = A[["link"]]),
                                     control = as.list(A$control))),
                       quote = TRUE)
      }, verbose = verbose)

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

  list(w = w, fit.obj = fit.obj)
}

#MABW with optweight
weightit2optweight <- function(covs, treat, s.weights, subset, estimand, focal, missing, moments, int, verbose, ...) {
  A <- list(...)

  rlang::check_installed("optweight")

  covs <- covs[subset, , drop = FALSE]
  treat <- factor(treat[subset])
  s.weights <- s.weights[subset]

  covs <- cbind(covs, int.poly.f(covs, poly = moments, int = int))
  for (i in seq_col(covs)) covs[,i] <- make.closer.to.1(covs[,i])

  if (missing == "ind") {
    covs <- add_missing_indicators(covs)
  }

  new.data <- data.frame(treat, covs)
  new.formula <- formula(new.data)

  for (f in names(formals(optweight::optweight))) {
    if (is_null(A[[f]])) A[[f]] <- formals(optweight::optweight)[[f]]
  }
  A[names(A) %in% names(formals(weightit2optweight))] <- NULL

  if ("tols" %in% names(A)) A[["tols"]] <- optweight::check.tols(new.formula, new.data, A[["tols"]], stop = TRUE)
  if ("targets" %in% names(A)) {
    .wrn("`targets` cannot be used through WeightIt and will be ignored")
    A[["targets"]] <- NULL
  }

  A[["formula"]] <- new.formula
  A[["data"]] <- new.data
  A[["estimand"]] <- estimand
  A[["s.weights"]] <- s.weights
  A[["focal"]] <- focal
  A[["verbose"]] <- TRUE

  verbosely({
    out <- do.call(optweight::optweight, A, quote = TRUE)
  }, verbose = verbose)

  list(w = out[["weights"]], info = list(duals = out$duals), fit.obj = out)
}
weightit2optweight.cont <- function(covs, treat, s.weights, subset, missing, moments, int, verbose, ...) {
  A <- list(...)
  rlang::check_installed("optweight")

  covs <- covs[subset, , drop = FALSE]
  treat <- treat[subset]
  s.weights <- s.weights[subset]

  covs <- cbind(covs, int.poly.f(covs, poly = moments, int = int))
  for (i in seq_col(covs)) covs[,i] <- make.closer.to.1(covs[,i])

  if (missing == "ind") {
    covs <- add_missing_indicators(covs)
  }

  new.data <- data.frame(treat, covs)
  new.formula <- formula(new.data)

  for (f in names(formals(optweight::optweight))) {
    if (is_null(A[[f]])) A[[f]] <- formals(optweight::optweight)[[f]]
  }
  A[names(A) %in% names(formals(weightit2optweight.cont))] <- NULL

  if ("tols" %in% names(A)) A[["tols"]] <- optweight::check.tols(new.formula, new.data, A[["tols"]], stop = TRUE)
  if ("targets" %in% names(A)) {
    .wrn("`targets` cannot be used through WeightIt and will be ignored")
    A[["targets"]] <- NULL
  }

  A[["formula"]] <- new.formula
  A[["data"]] <- new.data
  A[["s.weights"]] <- s.weights
  A[["verbose"]] <- TRUE

  verbosely({
    out <- do.call(optweight::optweight, A, quote = TRUE)
  }, verbose = verbose)

  list(w = out[["weights"]], info = list(duals = out$duals), fit.obj = out)
}
weightit2optweight.msm <- function(covs.list, treat.list, s.weights, subset, missing, moments, int, verbose, ...) {
  A <- list(...)
  rlang::check_installed("optweight")
  if (is_not_null(covs.list)) {
    covs.list <- lapply(covs.list, function(c) {
      covs <- c[subset, , drop = FALSE]
      covs <- cbind(covs, int.poly.f(covs, poly = moments, int = int))
      for (i in seq_col(covs)) covs[,i] <- make.closer.to.1(covs[,i])

      if (missing == "ind") {
        covs <- add_missing_indicators(covs)
      }

      covs
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
    .wrn("`targets` cannot be used through WeightIt and will be ignored")
    A[["targets"]] <- NULL
  }

  verbosely({
    out <- do.call(optweight::optweight.fit, c(list(treat = treat.list,
                                                    covs = covs.list,
                                                    s.weights = s.weights,
                                                    verbose = TRUE),
                                               A), quote = TRUE)
  }, verbose = verbose)

  list(w = out$w, fit.obj = out)
}

#Generalized boosted modeling with gbm and cobalt
weightit2gbm <- function(covs, treat, s.weights, estimand, focal, subset, stabilize, subclass, missing, verbose, ...) {

  rlang::check_installed("gbm")

  A <- list(...)

  covs <- covs[subset, , drop = FALSE]
  treat <- treat[subset]
  s.weights <- s.weights[subset]

  if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
  treat.type <- get.treat.type(treat)

  for (i in seq_col(covs)) covs[,i] <- make.closer.to.1(covs[,i])

  if (missing == "ind") {
    covs <- add_missing_indicators(covs, replace_with = NA)
  }

  if (is_null(A[["criterion"]])) {
    A[["criterion"]] <- A[["stop.method"]]
  }

  if (is_null(A[["criterion"]])) {
    .wrn("no `criterion` was provided. Using \"smd.mean\"",
         call. = FALSE, immediate. = TRUE)
    A[["criterion"]] <- "smd.mean"
  }
  else if (length(A[["criterion"]]) > 1) {
    .wrn("only one `criterion` is allowed at a time. Using just the first `criterion`")
    A[["criterion"]] <- A[["criterion"]][1]
  }

  available.criteria <- available.stats(treat.type)

  if (is.character(A[["criterion"]]) &&
      startsWith(A[["criterion"]], "es.")) {
    subbed.crit <- sub("es.", "smd.", A[["criterion"]], fixed = TRUE)
    subbed.match <- charmatch(subbed.crit, available.criteria)
    if (!anyNA(subbed.match) && subbed.match != 0L) {
      A[["criterion"]] <- subbed.crit
    }
  }

  cv <- 0

  s.m.matches <- charmatch(A[["criterion"]], available.criteria)
  if (anyNA(s.m.matches) || s.m.matches == 0L) {
    if (startsWith(A[["criterion"]], "cv") &&
        can_str2num(numcv <- substr(A[["criterion"]], 3, nchar(A[["criterion"]])))) {
      cv <- round(str2num(numcv))
      if (cv < 2) .err("at least 2 CV-folds must be specified in `criterion`")
    }
    else {
      .err(sprintf("`criterion` must be one of %s.",
                   word_list(c(available.criteria, "cv{#}"), "or", quotes = TRUE)))
    }
  }
  else criterion <- available.criteria[s.m.matches]

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

  chk::chk_count(A[["n.trees"]], "`n.trees`")
  chk::chk_gt(A[["n.trees"]], 1, "`n.trees`")

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
    if (is_null(A[["n.grid"]])) n.grid <- round(1 + sqrt(2 * (A[["n.trees"]] - start.tree + 1)))
    else if (!is_(A[["n.grid"]], "numeric") || length(A[["n.grid"]]) > 1 ||
             !between(A[["n.grid"]], c(2, A[["n.trees"]]))) {
      .err("`n.grid` must be a numeric value between 2 and `n.trees`")
    }
    else n.grid <- round(A[["n.grid"]])

    init <- bal.init(
      if (!anyNA(covs)) covs
      else if (missing == "surr") add_missing_indicators(covs)
      else replace_na_with(covs),
      treat = treat, stat = criterion,
      estimand = estimand, s.weights = s.weights,
      focal = focal, ...)
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

    verbosely({
      fit <- do.call(gbm::gbm.fit, c(A[names(A) %in% setdiff(names(formals(gbm::gbm.fit)), names(tune_args))], tune_args),
                     quote = TRUE)
    }, verbose = verbose)

    if (cv == 0) {

      n.trees <- fit[["n.trees"]]
      iters <- 1:n.trees
      iters.grid <- round(seq(start.tree, n.trees, length.out = n.grid))

      if (is_null(iters.grid) || anyNA(iters.grid) || any(iters.grid > n.trees)) {
        .err("a problem has occurred")
      }

      ps <- gbm::predict.gbm(fit, n.trees = iters.grid, type = "response", newdata = covs)
      w <- get.w.from.ps(ps, treat = treat, estimand = estimand, focal = focal, stabilize = stabilize, subclass = subclass)
      if (trim.at != 0) w <- suppressMessages(apply(w, 2, trim, at = trim.at, treat = treat))

      iter.grid.balance <- apply(w, 2, bal.compute, x = init)

      if (n.grid == n.trees) {
        best.tree.index <- which.min(iter.grid.balance)
        best.loss <- iter.grid.balance[best.tree.index]
        best.tree <- iters.grid[best.tree.index]
        tree.val <- setNames(data.frame(iters.grid,
                                        iter.grid.balance),
                             c("tree", criterion))
      }
      else {
        it <- which.min(iter.grid.balance) + c(-1, 1)
        it[1] <- iters.grid[max(1, it[1])]
        it[2] <- iters.grid[min(length(iters.grid), it[2])]
        iters.to.check <- iters[between(iters, iters[it])]

        if (is_null(iters.to.check) || anyNA(iters.to.check) || any(iters.to.check > n.trees)) {
          .err("a problem has occurred")
        }

        ps <- gbm::predict.gbm(fit, n.trees = iters.to.check, type = "response", newdata = covs)
        w <- get.w.from.ps(ps, treat = treat, estimand = estimand, focal = focal, stabilize = stabilize, subclass = subclass)
        if (trim.at != 0) w <- suppressMessages(apply(w, 2, trim, at = trim.at, treat = treat))

        iter.grid.balance.fine <- apply(w, 2, bal.compute, x = init)

        best.tree.index <- which.min(iter.grid.balance.fine)
        best.loss <- iter.grid.balance.fine[best.tree.index]
        best.tree <- iters.to.check[best.tree.index]
        tree.val <- setNames(data.frame(c(iters.grid, iters.to.check),
                                        c(iter.grid.balance, iter.grid.balance.fine)),
                             c("tree", criterion))
      }

      tree.val <- unique(tree.val[order(tree.val$tree),])
      w <- w[,best.tree.index]
      ps <- if (treat.type == "binary") ps[,best.tree.index] else NULL

      tune[[paste.("best", criterion)]][i] <- best.loss
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
      verbosely({
        cv.results <- do.call(gbm::gbmCrossVal,
                              c(A[names(A) %in% setdiff(names(formals(gbm::gbmCrossVal)), names(tune_args))],
                                tune_args), quote = TRUE)
      }, verbose = verbose)

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

  tune[tunable[lengths(A[tunable]) == 1]] <- NULL

  if (ncol(tune) > 2) {
    info[["tune"]] <- tune
    info[["best.tune"]] <- tune[best.tune.index,]
  }

  if (is_not_null(best.ps) && is_not_null(focal) && focal != t.lev) {
    best.ps <- 1 - best.ps
  }

  list(w = best.w, ps = best.ps, info = info, fit.obj = best.fit)
}
weightit2gbm.cont <- function(covs, treat, s.weights, estimand, focal, subset, stabilize, subclass, missing, verbose, ...) {

  rlang::check_installed("gbm")

  A <- list(...)

  covs <- covs[subset, , drop = FALSE]
  treat <- treat[subset]
  s.weights <- s.weights[subset]

  for (i in seq_col(covs)) covs[,i] <- make.closer.to.1(covs[,i])

  if (missing == "ind") {
    covs <- add_missing_indicators(covs, replace_with = NA)
  }

  if (is_null(A[["criterion"]])) {
    A[["criterion"]] <- A[["stop.method"]]
  }
  if (is_null(A[["criterion"]])) {
    .err("no `criterion` was provided. Using \"p.mean\"")
    A[["criterion"]] <- "p.mean"
  }
  else if (length(A[["criterion"]]) > 1) {
    .wrn("only one `criterion` is allowed at a time. Using just the first `criterion`")
    A[["criterion"]] <- A[["criterion"]][1]
  }

  available.criteria <- available.stats("continuous")

  cv <- 0

  s.m.matches <- charmatch(A[["criterion"]], available.criteria)
  if (anyNA(s.m.matches) || s.m.matches == 0L) {
    if (startsWith(A[["criterion"]], "cv") &&
        can_str2num(numcv <- substr(A[["criterion"]], 3, nchar(A[["criterion"]])))) {
      cv <- round(str2num(numcv))
      if (cv < 2) .err("at least 2 CV-folds must be specified in `criterion`")
    }
    else .err(sprintf("`criterion` must be one of %s",
                      word_list(c(available.criteria, "cv{#}"), "or", quotes = TRUE)))
  }
  else criterion <- available.criteria[s.m.matches]

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

  chk::chk_count(A[["n.trees"]], "`n.trees`")
  chk::chk_gt(A[["n.trees"]], 1, "`n.trees`")

  available.distributions <- c("gaussian", "laplace", "tdist", "poisson")

  if (cv == 0) {
    start.tree <- if_null_then(A[["start.tree"]], 1)
    if (is_null(A[["n.grid"]])) n.grid <- round(1 + sqrt(2 * (A[["n.trees"]] - start.tree + 1)))
    else if (!is_(A[["n.grid"]], "numeric") || length(A[["n.grid"]]) > 1 ||
             !between(A[["n.grid"]], c(2, A[["n.trees"]]))) {
      .err("`n.grid` must be a numeric value between 2 and `n.trees`")
    }
    else n.grid <- round(A[["n.grid"]])

    init <- bal.init(
      if (!anyNA(covs)) covs
      else if (missing == "surr") add_missing_indicators(covs)
      else replace_na_with(covs),
      treat = treat, stat = criterion,
      s.weights = s.weights, ...)
  }

  A[["x"]] <- covs
  A[["y"]] <- treat
  A[["distribution"]] <- {
    if (is_null(distribution <- A[["distribution"]])) {
      available.distributions[1]
    }
    else {
      match_arg(distribution, available.distributions, several.ok = TRUE)
    }
  }
  A[["w"]] <- s.weights
  A[["verbose"]] <- FALSE

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
      if (is_not_null(splitdens1 <- get0(splitdens[1], mode = "function", envir = parent.frame()))) {
        if (length(splitdens) > 1 && !can_str2num(splitdens[-1])) {
          .err(sprintf("%s is not an appropriate argument to `density` because %s cannot be coerced to numeric",
                       A[["density"]], word_list(splitdens[-1], and.or = "or", quotes = TRUE)))
        }
        densfun <- function(x) {
          tryCatch(do.call(splitdens1, c(list(x), as.list(str2num(splitdens[-1])))),
                   error = function(e) .err(sprintf("Error in applying density:\n  %s", conditionMessage(e)), tidy = FALSE))
        }
      }
      else {
        .err(sprintf("%s is not an appropriate argument to `density` because %s is not an available function",
                     A[["density"]], splitdens[1]))
      }
    }
    else .err("the argument to `density` cannot be evaluated as a density function")
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
      .err("there was a problem with the output of `density`. Try another density function or leave it blank to use the normal density")
    }
    else if (any(dens.num <= 0)) {
      .err("the input to `density` may not accept the full range of treatment values")
    }
  }

  current.best.loss <- Inf

  for (i in seq_row(tune)) {

    verbosely({
      fit <- do.call(gbm::gbm.fit, c(A[names(A) %in% setdiff(names(formals(gbm::gbm.fit)), tunable)],
                                     tune[i, tunable[tunable %in% names(formals(gbm::gbm.fit))]]),
                     quote = TRUE)
    }, verbose = verbose)

    if (cv == 0) {

      n.trees <- fit[["n.trees"]]
      iters <- 1:n.trees
      iters.grid <- round(seq(start.tree, n.trees, length.out = n.grid))

      if (is_null(iters.grid) || anyNA(iters.grid) || any(iters.grid > n.trees)) {
        .err("a problem has occurred")
      }

      gps <- gbm::predict.gbm(fit, n.trees = iters.grid, newdata = covs)
      w <- get_cont_weights(gps, treat = treat, s.weights = s.weights, dens.num = dens.num,
                            densfun = densfun, use.kernel = use.kernel, densControl = A)
      if (trim.at != 0) w <- suppressMessages(apply(w, 2, trim, at = trim.at, treat = treat))

      iter.grid.balance <- apply(w, 2, bal.compute, x = init)

      if (n.grid == n.trees) {
        best.tree.index <- which.min(iter.grid.balance)
        best.loss <- iter.grid.balance[best.tree.index]
        best.tree <- iters.grid[best.tree.index]
        tree.val <- setNames(data.frame(iters.grid,
                                        iter.grid.balance),
                             c("tree", criterion))
      }
      else {
        it <- which.min(iter.grid.balance) + c(-1, 1)
        it[1] <- iters.grid[max(1, it[1])]
        it[2] <- iters.grid[min(length(iters.grid), it[2])]
        iters.to.check <- iters[between(iters, iters[it])]

        if (is_null(iters.to.check) || anyNA(iters.to.check) || any(iters.to.check > n.trees)) {
          .err("a problem has occurred")
        }

        gps <- gbm::predict.gbm(fit, n.trees = iters.to.check, newdata = covs)
        w <- get_cont_weights(gps, treat = treat, s.weights = s.weights, dens.num = dens.num,
                              densfun = densfun, use.kernel = use.kernel, densControl = A)
        if (trim.at != 0) w <- suppressMessages(apply(w, 2, trim, at = trim.at, treat = treat))

        iter.grid.balance.fine <- apply(w, 2, bal.compute, x = init)

        best.tree.index <- which.min(iter.grid.balance.fine)
        best.loss <- iter.grid.balance.fine[best.tree.index]
        best.tree <- iters.to.check[best.tree.index]
        tree.val <- setNames(data.frame(c(iters.grid, iters.to.check),
                                        c(iter.grid.balance, iter.grid.balance.fine)),
                             c("tree", criterion))
      }

      tree.val <- unique(tree.val[order(tree.val$tree),])
      w <- w[,best.tree.index]
      gps <- gps[,as.character(best.tree)]

      tune[[paste.("best", criterion)]][i] <- best.loss
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

      verbosely({
        cv.results <- do.call(gbm::gbmCrossVal,
                              c(A[names(A) %in% setdiff(names(formals(gbm::gbmCrossVal)), tunable)],
                                tune[i, tunable[tunable %in% names(formals(gbm::gbmCrossVal))]]), quote = TRUE)
      }, verbose = verbose)

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

  tune[tunable[lengths(A[tunable]) == 1]] <- NULL

  if (ncol(tune) > 2) {
    info[["tune"]] <- tune
    info[["best.tune"]] <- tune[best.tune.index,]
  }

  list(w = best.w, info = info, fit.obj = best.fit)
}

#CBPS
weightit2cbps <- function(covs, treat, s.weights, estimand, focal, subset, stabilize, subclass, missing, verbose, ...) {

  rlang::check_installed("CBPS")

  A <- list(...)

  covs <- covs[subset, , drop = FALSE]
  treat <- factor(treat[subset])
  s.weights <- s.weights[subset]

  if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
  treat.type <- get.treat.type(treat)

  if (missing == "ind") {
    covs <- add_missing_indicators(covs)
  }

  for (i in seq_col(covs)) covs[,i] <- make.closer.to.1(covs[,i])

  colinear.covs.to.remove <- colnames(covs)[colnames(covs) %nin% colnames(make_full_rank(covs))]
  covs <- covs[, colnames(covs) %nin% colinear.covs.to.remove, drop = FALSE]

  sw0 <- check_if_zero(s.weights)

  if (estimand == "ATT") {
    ps <- make_df(levels(treat), length(treat))

    control.levels <- levels(treat)[levels(treat) != focal]
    fit.list <- make_list(control.levels)

    for (i in control.levels) {
      treat.in.i.focal <- treat %in% c(focal, i)
      treat_ <- as.integer(treat[treat.in.i.focal] != i)
      covs_ <- covs[treat.in.i.focal, , drop = FALSE]
      new.data <- as.data.frame(cbind(treat_, covs_))

      tryCatch({verbosely({
        fit.list[[i]] <- CBPS::CBPS(formula(new.data),
                                    data = new.data[!sw0[treat.in.i.focal],],
                                    method = if (is_not_null(A$over) && A$over == FALSE) "exact" else "over",
                                    standardize = FALSE,
                                    sample.weights = s.weights[!sw0 & treat.in.i.focal],
                                    ATT = 1,
                                    ...)
      }, verbose = verbose)},
      error = function(e) {
        e. <- conditionMessage(e)
        e. <- gsub("method = \"exact\"", "`over = FALSE`", e., fixed = TRUE)
        .err("(from `CBPS::CBPS()`) ", e., tidy = FALSE)
      })

      if (!any(sw0[treat.in.i.focal])) {
        ps[[focal]][treat.in.i.focal] <- fit.list[[i]][["fitted.values"]]
        ps[[i]][treat.in.i.focal] <- 1 - ps[[focal]][treat.in.i.focal]
      }
      else {
        ps[[focal]][treat.in.i.focal] <- plogis(drop(cbind(1, covs_) %*% fit.list[[i]][["coefficients"]]))
        ps[[i]][treat.in.i.focal] <- 1 - ps[[focal]][treat.in.i.focal]
      }
    }
  }
  else {
    new.data <- cbind(treat, as.data.frame(covs))
    if (treat.type == "binary" || !nunique.gt(treat, 4)) {
      tryCatch({verbosely({
        fit.list <- CBPS::CBPS(formula(new.data),
                               data = new.data[!sw0,],
                               method = if (isFALSE(A$over)) "exact" else "over",
                               standardize = FALSE,
                               sample.weights = s.weights[!sw0],
                               ATT = 0,
                               ...)
      }, verbose = verbose)},
      error = function(e) {
        e. <- conditionMessage(e)
        e. <- gsub("method = \"exact\"", "`over = FALSE`", e., fixed = TRUE)
        .err("(from `CBPS::CBPS()`) ", e., tidy = FALSE)
      })

      if (!any(sw0)) {
        ps <- fit.list[["fitted.values"]]
      }
      else if (treat.type == "binary") {
        ps <- plogis(drop(cbind(1, covs) %*% fit.list[["coefficients"]]))
      }
      else {
        ps <- make_df(levels(treat), length(treat))
        base.lvl <- setdiff(levels(treat), colnames(fit.list$coefficients))

        lin.pred <- cbind(1, covs) %*% fit.list[["coefficients"]]
        ps[, base.lvl] <- 1/rowSums(exp(lin.pred))
        ps[, colnames(fit.list$coefficients)] <- exp(lin.pred[, colnames(fit.list$coefficients)]) * ps[, base.lvl]
        ps <- ps/rowSums(ps)
      }
    }
    else {
      ps <- make_df(levels(treat), length(treat))

      fit.list <- make_list(levels(treat))

      for (i in levels(treat)) {
        new.data[[1]] <- as.integer(treat == i)

        tryCatch({verbosely({
          fit.list[[i]] <- CBPS::CBPS(formula(new.data), data = new.data[!sw0,],
                                      method = if (isFALSE(A$over)) "exact" else "over",
                                      standardize = FALSE,
                                      sample.weights = s.weights[!sw0],
                                      ATT = 0, ...)
        }, verbose = verbose)},
        error = function(e) {
          e. <- conditionMessage(e)
          e. <- gsub("method = \"exact\"", "`over = FALSE`", e., fixed = TRUE)
          .err("(from `CBPS::CBPS()`) ", e., tidy = FALSE)
        })

        if (!any(sw0)) {
          ps[, i] <- fit.list[[i]][["fitted.values"]]
        }
        else {
          ps[, i] <- plogis(drop(cbind(1, covs) %*% fit.list[[i]][["coefficients"]]))
        }
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

  list(w = w, ps = p.score, fit.obj = fit.list)
}
weightit2cbps.cont <- function(covs, treat, s.weights, subset, missing, verbose, ...) {
  rlang::check_installed("CBPS")

  A <- list(...)

  covs <- covs[subset, , drop = FALSE]
  treat <- treat[subset]
  s.weights <- s.weights[subset]

  if (missing == "ind") {
    covs <- add_missing_indicators(covs)
  }

  for (i in seq_col(covs)) covs[,i] <- make.closer.to.1(covs[,i])

  colinear.covs.to.remove <- colnames(covs)[colnames(covs) %nin% colnames(make_full_rank(covs))]
  covs <- covs[, colnames(covs) %nin% colinear.covs.to.remove, drop = FALSE]

  sw0 <- check_if_zero(s.weights)

  new.data <- data.frame(treat = treat, covs)

  w <- rep(0, length(treat))

  tryCatch({verbosely({
    fit <- CBPS::CBPS(formula(new.data),
                      data = new.data[!sw0,],
                      method = if (isFALSE(A$over)) "exact" else "over",
                      standardize = FALSE,
                      sample.weights = s.weights[!sw0],
                      ...)
  }, verbose = verbose)},
  error = function(e) {
    e. <- conditionMessage(e)
    e. <- gsub("method = \"exact\"", "`over = FALSE`", e., fixed = TRUE)
    .err("(from `CBPS::CBPS()`) ", e., tidy = FALSE)
  })

  w[!sw0] <- fit$weights / s.weights[!sw0]

  list(w = w, fit.obj = fit)
}
weightit2cbps.msm <- function(covs.list, treat.list, s.weights, subset, missing, verbose, ...) {
  .err("CBMSM doesn't work yet")
}
weightit2npcbps <- function(covs, treat, s.weights, subset, missing, moments, int, verbose, ...) {
  rlang::check_installed("CBPS")

  A <- list(...)

  if (!all_the_same(s.weights)) {
    .err("sampling weights cannot be used with `method = \"npcbps\"`")
  }

  covs <- covs[subset, , drop = FALSE]
  treat <- factor(treat[subset])

  if (missing == "ind") {
    covs <- add_missing_indicators(covs)
  }

  covs <- cbind(covs, int.poly.f(covs, poly = moments, int = int))

  for (i in seq_col(covs)) covs[,i] <- make.closer.to.1(covs[,i])

  colinear.covs.to.remove <- colnames(covs)[colnames(covs) %nin% colnames(make_full_rank(covs))]
  covs <- covs[, colnames(covs) %nin% colinear.covs.to.remove, drop = FALSE]

  new.data <- data.frame(treat = treat, covs)

  tryCatch({verbosely({
    fit <- do.call(CBPS::npCBPS, c(list(formula(new.data), data = new.data, print.level = 1), A),
                   quote = TRUE)
  }, verbose = verbose)},
  error = function(e) {
    e. <- conditionMessage(e)
    .err("(from `CBPS::npCBPS()`) ", e., tidy = FALSE)
  })

  w <- fit$weights

  for (i in levels(treat)) w[treat == i] <- w[treat == i]/mean(w[treat == i])

  list(w = w, fit.obj = fit)
}
weightit2npcbps.cont <- function(covs, treat, s.weights, subset, missing, moments, int, verbose, ...) {
  rlang::check_installed("CBPS")

  A <- list(...)

  if (!all_the_same(s.weights)) {
    .err("sampling weights cannot be used with `method = \"npcbps\"`")
  }

  covs <- covs[subset, , drop = FALSE]
  treat <- treat[subset]

  if (missing == "ind") {
    covs <- add_missing_indicators(covs)
  }

  for (i in seq_col(covs)) covs[,i] <- make.closer.to.1(covs[,i])

  covs <- cbind(covs, int.poly.f(covs, poly = moments, int = int))

  colinear.covs.to.remove <- colnames(covs)[colnames(covs) %nin% colnames(make_full_rank(covs))]
  covs <- covs[, colnames(covs) %nin% colinear.covs.to.remove, drop = FALSE]

  new.data <- data.frame(treat = treat, covs)

  tryCatch({verbosely({
    fit <- do.call(CBPS::npCBPS, c(list(formula(new.data), data = new.data, print.level = 1), A),
                   quote = TRUE)
  }, verbose = verbose)},
  error = function(e) {
    e. <- conditionMessage(e)
    .err("(from `CBPS::npCBPS()`) ", e., tidy = FALSE)
  })

  w <- fit$weights

  w <- w/mean(w)

  list(w = w, fit.obj = fit)
}

#Entropy balancing
weightit2ebal <- function(covs, treat, s.weights, subset, estimand, focal, stabilize, missing, moments, int, verbose, ...) {
  A <- list(...)

  covs <- covs[subset, , drop = FALSE]
  treat <- factor(treat[subset])
  s.weights <- s.weights[subset]

  if (missing == "ind") {
    covs <- add_missing_indicators(covs)
  }

  covs <- cbind(covs, int.poly.f(covs, poly = moments, int = int, center = TRUE))
  for (i in seq_col(covs)) covs[,i] <- make.closer.to.1(covs[,i])

  if (is_not_null(A[["base.weights"]])) A[["base.weight"]] <- A[["base.weights"]]
  if (is_null(A[["base.weight"]])) {
    bw <- rep(1, length(treat))
  }
  else {
    if (!is.numeric(A[["base.weight"]]) || length(A[["base.weight"]]) != length(treat)) {
      .err("the argument to `base.weight` must be a numeric vector with length equal to the number of units")
    }
    else bw <- A[["base.weight"]]
  }

  reltol <- if_null_then(A[["reltol"]], sqrt(.Machine$double.eps))

  eb <- function(C, s.weights_t, Q) {
    n <- nrow(C)

    W <- function(Z) {
      drop(Q * exp(-C %*% Z))
    }

    objective.EB <- function(Z) {
      log(sum(W(Z)))
    }

    gradient.EB <- function(Z) {
      w <- W(Z)
      drop(-(w %*% C)/sum(w))
    }

    opt.out <- optim(par = rep(0, ncol(C)),
                     fn = objective.EB,
                     gr = gradient.EB,
                     method = "BFGS",
                     control = list(trace = 1,
                                    reltol = reltol,
                                    maxit = if_null_then(A[["maxit"]], 200)))

    w <- W(opt.out$par)
    opt.out$gradient <- gradient.EB(opt.out$par)

    if (opt.out$convergence == 1) {
      .wrn("the optimization failed to converge in the alotted number of iterations. Try increasing `maxit`")
    }
    else if (any(abs(opt.out$gradient) > 1e-3)) {
      .wrn("the estimated weights do not balance the covariates, indicating the optimization arrived at a degenerate solution. Try decreasing the number of variables supplied to the optimization")
    }

    if (sum(w) > n*.Machine$double.eps) w <- w*n/sum(w)

    list(Z = setNames(opt.out$par, colnames(C)),
         w = w/s.weights_t,
         opt.out = opt.out)
  }

  w <- rep(1, length(treat))
  sw0 <- check_if_zero(s.weights)

  if (estimand == "ATT") {
    groups_to_weight <- levels(treat)[levels(treat) != focal]
    targets <- cobalt::col_w_mean(covs, s.weights = s.weights, subset = treat == focal)
  }
  else if (estimand == "ATE") {
    groups_to_weight <- levels(treat)
    targets <- cobalt::col_w_mean(covs, s.weights = s.weights)
  }
  covs <- sweep(covs, 2, targets, check.margin = FALSE)

  fit.list <- make_list(groups_to_weight)
  for (i in groups_to_weight) {
    verbosely({
      fit.list[[i]] <- eb(covs[treat == i & !sw0,,drop = FALSE], s.weights[treat == i & !sw0],
                          bw[treat == i & !sw0])
    }, verbose = verbose)

    w[treat == i & !sw0] <- fit.list[[i]]$w
  }

  list(w = w, fit.obj = lapply(fit.list, function(x) x[["opt.out"]]))
}
weightit2ebal.cont <- function(covs, treat, s.weights, subset, missing, moments, int, verbose, ...) {
  A <- list(...)

  covs <- covs[subset, , drop = FALSE]
  treat <- treat[subset]
  s.weights <- s.weights[subset]

  if (missing == "ind") {
    covs <- add_missing_indicators(covs)
  }

  if (is_not_null(A[["base.weights"]])) A[["base.weight"]] <- A[["base.weights"]]
  if (is_null(A[["base.weight"]])) {
    bw <- rep(1, length(treat))
  }
  else {
    if (!is.numeric(A[["base.weight"]]) || length(A[["base.weight"]]) != length(treat)) {
      .err("the argument to `base.weight` must be a numeric vector with length equal to the number of units")
    }
    else bw <- A[["base.weight"]]
  }

  bw <- bw[subset]

  reltol <- if_null_then(A[["reltol"]], sqrt(.Machine$double.eps))

  d.moments <- max(if_null_then(A[["d.moments"]], 1), moments)
  k <- ncol(covs)

  poly.covs <- int.poly.f(covs, poly = moments)
  int.covs <- int.poly.f(covs, int = int)

  treat <- make.closer.to.1(treat)
  for (i in seq_col(poly.covs)) poly.covs[,i] <- make.closer.to.1(poly.covs[,i])
  for (i in seq_col(int.covs)) int.covs[,i] <- make.closer.to.1(int.covs[,i])
  if (d.moments == moments) {
    d.poly.covs <- poly.covs
  }
  else {
    d.poly.covs <- int.poly.f(covs, poly = d.moments)
    for (i in seq_col(d.poly.covs)) d.poly.covs[,i] <- make.closer.to.1(d.poly.covs[,i])
  }
  for (i in seq_col(covs)) covs[,i] <- make.closer.to.1(covs[,i])

  covs <- cbind(covs, poly.covs, int.covs, d.poly.covs)
  # colinear.covs.to.remove <- colnames(covs)[colnames(covs) %nin% colnames(make_full_rank(covs))]
  # covs <- covs[, colnames(covs) %nin% colinear.covs.to.remove, drop = FALSE]

  t.mat <- matrix(treat, ncol = 1, dimnames = list(NULL, "treat"))
  if (d.moments > 1) t.mat <- cbind(t.mat, int.poly.f(t.mat, poly = d.moments))

  treat_c <- sweep(t.mat, 2, cobalt::col_w_mean(t.mat, s.weights))
  covs_c <- sweep(covs, 2, cobalt::col_w_mean(covs, s.weights))

  covs.ind <- seq_len(k)
  poly.covs.ind <- k + seq_col(poly.covs)
  int.covs.ind <- k + ncol(poly.covs) + seq_col(int.covs)
  d.poly.covs.ind <- k + ncol(poly.covs) + ncol(int.covs) + seq_col(d.poly.covs)

  C <- cbind(treat_c, covs_c[, c(covs.ind, int.covs.ind, d.poly.covs.ind)],
             treat_c[,1] * covs_c[, c(covs.ind, int.covs.ind, poly.covs.ind)])

  colnames(C) <- c(paste(colnames(treat_c), "(mean)"),
                   paste(colnames(covs_c)[c(covs.ind, int.covs.ind, d.poly.covs.ind)], "(mean)"),
                   colnames(covs_c)[c(covs.ind, int.covs.ind, poly.covs.ind)])

  eb <- function(C, s.weights, Q) {
    n <- nrow(C)

    W <- function(Z) {
      drop(Q * exp(-C %*% Z))
    }

    objective.EB <- function(Z) {
      log(sum(W(Z)))
    }

    gradient.EB <- function(Z) {
      w <- W(Z)
      drop(-(w %*% C)/sum(w))
    }

    opt.out <- optim(par = rep(0, ncol(C)),
                     fn = objective.EB,
                     gr = gradient.EB,
                     method = "BFGS",
                     control = list(trace = 0,
                                    reltol = reltol,
                                    maxit = if_null_then(A[["maxit"]], 200)))

    w <- W(opt.out$par)
    opt.out$gradient <- gradient.EB(opt.out$par)

    if (opt.out$convergence == 1) {
      .wrn("the optimization failed to converge in the alotted number of iterations. Try increasing `maxit`")
    }
    else if (any(abs(opt.out$gradient) > 1e-3)) {
      .wrn("the estimated weights do not balance the covariates, indicating the optimization arrived at a degenerate solution. Try decreasing the number of variables supplied to the optimization")
    }

    if (sum(w) > n*.Machine$double.eps) w <- w*n/sum(w)

    list(Z = setNames(opt.out$par, colnames(C)),
         w = w/s.weights,
         opt.out = opt.out)
  }

  w <- rep(0, length(treat))
  sw0 <- check_if_zero(s.weights)

  verbosely({
    fit <- eb(C[!sw0,, drop = FALSE], s.weights[!sw0], bw[!sw0])
  }, verbose = verbose)

  w[!sw0] <- fit$w

  list(w = w, fit.obj = fit$opt.out)
}

#PS weights using SuperLearner
weightit2super <- function(covs, treat, s.weights, subset, estimand, focal, stabilize, subclass, missing, verbose, ...) {
  A <- list(...)

  rlang::check_installed("SuperLearner")

  covs <- covs[subset, , drop = FALSE]
  treat <- factor(treat[subset])
  s.weights <- s.weights[subset]

  if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
  treat.type <- get.treat.type(treat)

  if (missing == "ind") {
    covs <- add_missing_indicators(covs)
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
  chk::chk_flag(discrete)

  if (identical(A[["SL.method"]], "method.balance")) {
    if (treat.type != "binary") .err("\"method.balance\" cannot be used with multi-category treatments")

    if (is_null(A[["criterion"]])) {
      A[["criterion"]] <- A[["stop.method"]]
    }

    if (is_null(A[["criterion"]])) {
      .wrn("no `criterion` was provided. Using \"smd.mean\"")
      A[["criterion"]] <- "smd.mean"
    }
    else if (length(A[["criterion"]]) > 1) {
      .wrn("only one `criterion` is allowed at a time. Using just the first `criterion`")
      A[["criterion"]] <- A[["criterion"]][1]
    }

    available.criteria <- available.stats(treat.type)

    if (is.character(A[["criterion"]]) &&
        startsWith(A[["criterion"]], "es.")) {
      subbed.crit <- sub("es.", "smd.", A[["criterion"]], fixed = TRUE)
      subbed.match <- charmatch(subbed.crit, available.criteria)
      if (!anyNA(subbed.match) && subbed.match != 0L) {
        A[["criterion"]] <- subbed.crit
      }
    }

    criterion <- A[["criterion"]]
    criterion <- match_arg(criterion, available.criteria)

    init <- bal.init(covs, treat = treat, stat = criterion,
                     estimand = estimand, s.weights = s.weights,
                     focal = focal, ...)

    sneaky <- 0
    attr(sneaky, "vals") <- list(init = init, estimand = estimand)
    A[["control"]] <- list(trimLogit = sneaky)

    A[["SL.method"]] <- method.balance()
  }

  fit.list <- info <- make_list(levels(treat))
  ps <- make_df(levels(treat), nrow = length(treat))

  for (i in levels(treat)) {

    if (treat.type == "binary" && i == last(levels(treat))) {
      ps[[i]] <- 1 - ps[[1]]
      fit.list <- fit.list[[1]]
      info <- info[[1]]
      next
    }

    treat_i <- as.numeric(treat == i)
    tryCatch({verbosely({
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
    }, verbose = verbose)},
    error = function(e) {
      e. <- conditionMessage(e)
      .err("(from `SuperLearner::SuperLearner()`) ", e., tidy = FALSE)
    })

    if (discrete) ps[[i]] <- fit.list[[i]]$library.predict[,which.min(fit.list[[i]]$cvRisk)]
    else ps[[i]] <- fit.list[[i]]$SL.predict

    info[[i]] <- list(coef = fit.list[[i]]$coef,
                      cvRisk = fit.list[[i]]$cvRisk)
  }

  #ps should be matrix of probs for each treat
  #Computing weights
  w <- get_w_from_ps(ps = ps, treat = treat, estimand, focal, stabilize = stabilize, subclass = subclass)

  p.score <- if (treat.type == "binary") ps[[get.treated.level(treat)]] else NULL

  list(w = w, ps = p.score, info = info, fit.obj = fit.list)
}
weightit2super.cont <- function(covs, treat, s.weights, subset, stabilize, missing, ps, verbose, ...) {
  A <- B <- list(...)

  covs <- covs[subset, , drop = FALSE]
  treat <- treat[subset]
  s.weights <- s.weights[subset]

  if (missing == "ind") {
    covs <- add_missing_indicators(covs)
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
      if (is_not_null(splitdens1 <- get0(splitdens[1], mode = "function", envir = parent.frame()))) {
        if (length(splitdens) > 1 && !can_str2num(splitdens[-1])) {
          .err(sprintf("%s is not an appropriate argument to `density` because %s cannot be coerced to numeric",
                       A[["density"]], word_list(splitdens[-1], and.or = "or", quotes = TRUE)))
        }
        densfun <- function(x) {
          tryCatch(do.call(splitdens1, c(list(x), as.list(str2num(splitdens[-1])))),
                   error = function(e) .err(sprintf("Error in applying density:\n  %s", conditionMessage(e)), tidy = FALSE))
        }
      }
      else {
        .err(sprintf("%s is not an appropriate argument to `density` because %s is not an available function",
                     A[["density"]], splitdens[1]))
      }
    }
    else .err("the argument to `density` cannot be evaluated as a density function")
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
      .err("there was a problem with the output of `density`. Try another density function or leave it blank to use the normal density")
    }
    if (any(dens.num <= 0)) {
      .err("the input to density may not accept the full range of treatment values")
    }
  }

  #Estimate GPS
  for (f in names(formals(SuperLearner::SuperLearner))) {
    if (f == "method") {if (is_null(B[["SL.method"]])) B[["SL.method"]] <- formals(SuperLearner::SuperLearner)[["method"]]}
    else if (f == "env") {if (is_null(B[["env"]])) B[["env"]] <- environment(SuperLearner::SuperLearner)}
    else if (is_null(B[[f]])) B[[f]] <- formals(SuperLearner::SuperLearner)[[f]]
  }

  discrete <- if_null_then(A[["discrete"]], FALSE)
  chk::chk_flag(discrete)

  if (identical(B[["SL.method"]], "method.balance")) {

    if (is_null(A[["criterion"]])) {
      A[["criterion"]] <- A[["stop.method"]]
    }

    if (is_null(A[["criterion"]])) {
      .wrn("no `criterion` was provided. Using \"p.mean\"")
      A[["criterion"]] <- "p.mean"
    }
    else if (length(A[["criterion"]]) > 1) {
      .wrn("only one `criterion` is allowed at a time. Using just the first `criterion`")
      A[["criterion"]] <- A[["criterion"]][1]
    }

    available.criteria <- available.stats("continuous")

    criterion <- A[["criterion"]]
    criterion <- match_arg(criterion, available.criteria)

    init <- bal.init(covs, treat = treat, stat = criterion,
                     s.weights = s.weights, ...)

    sneaky <- 0
    attr(sneaky, "vals") <- list(init = init,
                                 dens.num = dens.num,
                                 densfun = densfun,
                                 use.kernel = use.kernel,
                                 densControl = A)
    B[["control"]] <- list(trimLogit = sneaky)

    B[["SL.method"]] <- method.balance.cont()
  }

  tryCatch({verbosely({
    fit <- do.call(SuperLearner::SuperLearner,
                   list(Y = treat,
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
  }, verbose = verbose)},
  error = function(e) {
    e. <- conditionMessage(e)
    .err("(from `SuperLearner::SuperLearner()`) ", e., tidy = FALSE)
  })


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

  list(w = w, info = info, fit.obj = fit)
}

#PS weights using BART
weightit2bart <- function(covs, treat, s.weights, subset, estimand, focal, stabilize, subclass, missing, verbose, ...) {
  A <- list(...)

  rlang::check_installed("dbarts")

  covs <- covs[subset, , drop = FALSE]
  treat <- factor(treat[subset])
  s.weights <- s.weights[subset]

  if (!all_the_same(s.weights)) {
    .err("sampling weights cannot be used with `method = \"bart\"`")
  }

  if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
  treat.type <- get.treat.type(treat)

  if (missing == "ind") {
    covs <- add_missing_indicators(covs)
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
    tryCatch({verbosely({
      fit.list[[i]] <- do.call(dbarts::bart2,
                               A[names(A) %in% setdiff(c(names(formals(dbarts::bart2)),
                                                         names(formals(dbarts::dbartsControl))),
                                                       c("offset.test", "weights", "subset", "test"))],
                               quote = TRUE)
    }, verbose = verbose)},
    error = function(e) {
      e. <- conditionMessage(e)
      .err("(from `dbarts::bart2()`) ", e., tidy = FALSE)
    })

    ps[[i]] <- fitted(fit.list[[i]])

  }

  info <- list()

  #ps should be matrix of probs for each treat
  #Computing weights
  w <- get_w_from_ps(ps = ps, treat = treat, estimand, focal, stabilize = stabilize, subclass = subclass)

  p.score <- if (treat.type == "binary") ps[[get.treated.level(treat)]] else NULL

  list(w = w, ps = p.score, info = info, fit.obj = fit.list)
}
weightit2bart.cont <- function(covs, treat, s.weights, subset, stabilize, missing, ps, verbose, ...) {
  A <- list(...)

  rlang::check_installed("dbarts")

  covs <- covs[subset, , drop = FALSE]
  treat <- treat[subset]
  s.weights <- s.weights[subset]

  if (!all_the_same(s.weights)) {
    .err("sampling weights cannot be used with `method = \"bart\"`")
  }

  if (missing == "ind") {
    covs <- add_missing_indicators(covs)
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
      if (is_not_null(splitdens1 <- get0(splitdens[1], mode = "function", envir = parent.frame()))) {
        if (length(splitdens) > 1 && !can_str2num(splitdens[-1])) {
          .err(sprintf("%s is not an appropriate argument to `density` because %s cannot be coerced to numeric",
                       A[["density"]], word_list(splitdens[-1], and.or = "or", quotes = TRUE)))
        }
        densfun <- function(x) {
          tryCatch(do.call(splitdens1, c(list(x), as.list(str2num(splitdens[-1])))),
                   error = function(e) .err(sprintf("Error in applying density:\n  %s", conditionMessage(e)), tidy = FALSE))
        }
      }
      else {
        .err(sprintf("%s is not an appropriate argument to `density` because %s is not an available function",
                     A[["density"]], splitdens[1]))
      }
    }
    else .err("the argument to `density` cannot be evaluated as a density function")
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
      .err("there was a problem with the output of `density`. Try another density function or leave it blank to use the normal density")
    }
    if (any(dens.num <= 0)) {
      .err("the input to density may not accept the full range of treatment values")
    }
  }

  A[["formula"]] <- covs
  A[["data"]] <- treat
  A[["keepCall"]] <- FALSE
  A[["combineChains"]] <- TRUE
  A[["verbose"]] <- FALSE #necessary to prevent crash

  #Estimate GPS

  tryCatch({verbosely({
    fit <- do.call(dbarts::bart2,
                   A[names(A) %in% setdiff(c(names(formals(dbarts::bart2)),
                                             names(formals(dbarts::dbartsControl))),
                                           c("offset.test", "weights", "subset", "test"))],
                   quote = TRUE)
  }, verbose = verbose)},
  error = function(e) {
    e. <- conditionMessage(e)
    .err("(from `dbarts::bart2()`) ", e., tidy = FALSE)
  })

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

  list(w = w, info = info, fit.obj = fit)
}

#Energy balancing
weightit2energy <- function(covs, treat, s.weights, subset, estimand, focal, missing, moments, int, verbose, ...) {
  rlang::check_installed("osqp")

  A <- list(...)

  if (missing == "ind") {
    covs <- add_missing_indicators(covs)
  }

  d <- if_null_then(A[["dist.mat"]], "scaled_euclidean")
  A[["dist.mat"]] <- NULL

  if (is.character(d) && length(d) == 1L) {
    dist.covs <- transform_covariates(data = covs, method = d,
                                      s.weights = s.weights, discarded = !subset)
    d <- unname(eucdist_internal(dist.covs))
  }
  else {
    if (inherits(d, "dist")) d <- as.matrix(d)

    if (!is.matrix(d) || !all(dim(d) == length(treat)) ||
        !all(check_if_zero(diag(d))) || any(d < 0) ||
        !isSymmetric(unname(d))) {
      .err(sprintf("`dist.mat` must be one of %s or a square, symmetric distance matrix with a value for all pairs of units",
                   word_list(weightit_distances(), "or", quotes = TRUE)))
    }

  }

  d <- unname(d[subset, subset])

  covs <- covs[subset, , drop = FALSE]
  treat <- factor(treat[subset])
  s.weights <- s.weights[subset]

  n <- length(treat)
  levels_treat <- levels(treat)
  diagn <- diag(n)

  min.w <- if_null_then(A[["min.w"]], 1e-8)
  if (!chk::vld_number(min.w)) {
    .wrn("`min.w` must be a single number. Setting `min.w = 1e-8`")
    min.w <- 1e-8
  }

  lambda <- if_null_then(A[["lambda"]], 1e-4)
  if (!chk::vld_number(lambda)) {
    .wrn("`lambda` must be a single number. Setting lambda = 1e-4.")
    lambda <- 1e-4
  }

  for (t in levels_treat) s.weights[treat == t] <- s.weights[treat == t]/mean(s.weights[treat == t])

  treat_t <- vapply(levels_treat, function(t) treat == t, logical(n))
  n_t <- colSums(treat_t)

  s.weights_n_t <- setNames(lapply(levels_treat, function(t) treat_t[,t] * s.weights / n_t[t]),
                            levels_treat)

  if (estimand == "ATE") {

    P <- -d * Reduce("+", lapply(s.weights_n_t, tcrossprod))

    q <- ((s.weights * 2 / n) %*% d) * Reduce("+", s.weights_n_t)

    if (!isFALSE(A[["improved"]])) {
      all_pairs <- combn(levels_treat, 2, simplify = FALSE)
      P <- P - d * Reduce("+", lapply(all_pairs, function(p) {
        tcrossprod(s.weights_n_t[[p[1]]] - s.weights_n_t[[p[2]]])
      }))
    }

    #Constraints for positivity and sum of weights
    Amat <- cbind(diagn, s.weights * treat_t)
    lvec <- c(rep(min.w, n), n_t)
    uvec <- c(ifelse(check_if_zero(s.weights), min.w, Inf), n_t)

    if (moments != 0 || int) {
      #Exactly balance moments and/or interactions
      covs <- cbind(covs, int.poly.f(covs, poly = moments, int = int))

      targets <- col.w.m(covs, s.weights)

      Amat <- cbind(Amat, do.call("cbind", lapply(s.weights_n_t, function(x) covs * x)))
      lvec <- c(lvec, rep(targets, length(levels_treat)))
      uvec <- c(uvec, rep(targets, length(levels_treat)))
    }
  }
  else {
    non_focal <- setdiff(levels_treat, focal)
    in_focal <- treat == focal

    P <- -d[!in_focal, !in_focal] *
      Reduce("+", lapply(s.weights_n_t[non_focal], function(x) tcrossprod(x[!in_focal])))

    q <- 2 * (s.weights_n_t[[focal]][in_focal] %*% d[in_focal, !in_focal]) *
      Reduce("+", lapply(s.weights_n_t[non_focal], function(x) x[!in_focal]))

    Amat <- cbind(diag(sum(!in_focal)), s.weights[!in_focal] * treat_t[!in_focal, non_focal])
    lvec <- c(rep(min.w, sum(!in_focal)), n_t[non_focal])
    uvec <- c(ifelse_(check_if_zero(s.weights[!in_focal]), min.w, Inf), n_t[non_focal])

    if (moments != 0 || int) {
      #Exactly balance moments and/or interactions
      covs <- cbind(covs, int.poly.f(covs, poly = moments, int = int))

      targets <- col.w.m(covs[in_focal,, drop = FALSE], s.weights[in_focal])

      Amat <- cbind(Amat, do.call("cbind", lapply(s.weights_n_t[non_focal],
                                                  function(x) covs[!in_focal,, drop = FALSE] * x[!in_focal])))
      lvec <- c(lvec, rep(targets, length(non_focal)))
      uvec <- c(uvec, rep(targets, length(non_focal)))
    }
  }

  #Add weight penalty
  if (lambda < 0) {
    #Find lambda to make P PSD
    e <- eigen(P, symmetric = TRUE, only.values = TRUE)
    e.min <- min(e$values)
    if (e.min < 0) {
      lambda <- -e.min*n^2
    }
  }

  diag(P) <- diag(P) + lambda / n^2

  if (is_not_null(A[["eps"]])) {
    if (is_null(A[["eps_abs"]])) A[["eps_abs"]] <- A[["eps"]]
    if (is_null(A[["eps_rel"]])) A[["eps_rel"]] <- A[["eps"]]
  }
  A[names(A) %nin% names(formals(osqp::osqpSettings))] <- NULL
  if (is_null(A[["max_iter"]])) A[["max_iter"]] <- 2e3L
  if (is_null(A[["eps_abs"]])) A[["eps_abs"]] <- 1e-8
  if (is_null(A[["eps_rel"]])) A[["eps_rel"]] <- 1e-6
  A[["verbose"]] <- TRUE

  options.list <- do.call(osqp::osqpSettings, A)

  verbosely({
    opt.out <- osqp::solve_osqp(P = 2 * P, q = q, A = t(Amat), l = lvec, u = uvec,
                                pars = options.list)
  }, verbose = verbose)

  if (identical(opt.out$info$status, "maximum iterations reached")) {
    .wrn("the optimization failed to converge. See Notes section at `?method_energy` for information")
  }

  if (estimand == "ATT") {
    w <- rep(1, length(treat))
    w[treat != focal] <- opt.out$x
  }
  else {
    w <- opt.out$x
  }

  w[w <= min.w] <- min.w

  opt.out$lambda <- lambda

  list(w = w, fit.obj = opt.out)
}
weightit2energy.cont <- function(covs, treat, s.weights, subset, missing, moments, int, verbose, ...) {
  rlang::check_installed("osqp")

  A <- list(...)

  if (missing == "ind") {
    covs <- add_missing_indicators(covs)
  }

  Xdist <- if_null_then(A[["dist.mat"]], "scaled_euclidean")
  A[["dist.mat"]] <- NULL

  if (is.character(Xdist) && length(Xdist) == 1L) {
    dist.covs <- transform_covariates(data = covs, method = Xdist,
                                      s.weights = s.weights, discarded = !subset)
    Xdist <- unname(eucdist_internal(X = dist.covs))
  }
  else {
    if (inherits(Xdist, "dist")) Xdist <- as.matrix(Xdist)

    if (!is.matrix(Xdist) || !all(dim(Xdist) == length(treat)) ||
        !all(check_if_zero(diag(Xdist))) || any(Xdist < 0) ||
        !isSymmetric(unname(Xdist))) {
      .err(sprintf("`dist.mat` must be one of %s or a square, symmetric distance matrix with a value for all pairs of units",
                   word_list(weightit_distances(), "or", quotes = TRUE)))
    }

  }

  Xdist <- unname(Xdist[subset, subset])

  covs <- covs[subset, , drop = FALSE]
  treat <- treat[subset]
  s.weights <- s.weights[subset]

  n <- length(treat)

  sw0 <- check_if_zero(s.weights)

  s.weights <- n * s.weights/sum(s.weights)

  min.w <- if_null_then(A[["min.w"]], 1e-8)
  if (!chk::vld_number(min.w)) {
    .wrn("`min.w` must be a single number. Setting `min.w = 1e-8`")
    min.w <- 1e-8
  }

  lambda <- if_null_then(A[["lambda"]], 1e-4)
  if (!chk::vld_number(lambda)) {
    .wrn("`lambda` must be a single number. Setting lambda = 1e-4")
    lambda <- 1e-4
  }

  d.moments <- max(if_null_then(A[["d.moments"]], 0), moments)
  if (!chk::vld_number(d.moments)) {
    .wrn(sprintf("`d.moments` must be a single number. Setting `lambda = %s`", moments))
    lambda <- 1e-4
  }

  dimension.adj <- if_null_then(A[["dimension.adj"]], TRUE)

  Adist <- eucdist_internal(X = treat/sqrt(col.w.v(treat, s.weights)))

  Xmeans <- colMeans(Xdist)
  Xgrand_mean <- mean(Xmeans)
  XA <- Xdist + Xgrand_mean - outer(Xmeans, Xmeans, "+")

  Ameans <- colMeans(Adist)
  Agrand_mean <- mean(Ameans)
  AA <- Adist + Agrand_mean - outer(Ameans, Ameans, "+")

  Pdcow <- XA * AA/n^2
  PebA <- -Adist/n^2
  PebX <- -Xdist/n^2

  qebA <- drop(s.weights %*% Adist)*2/n^2
  qebX <- drop(s.weights %*% Xdist)*2/n^2

  if (isFALSE(dimension.adj)) {
    Q_energy_A_adj <- 1 / 2
  }
  else {
    Q_energy_A_adj <- 1 / (1 + sqrt(ncol(covs)))
  }
  Q_energy_X_adj <- 1 - Q_energy_A_adj

  PebA <- PebA * Q_energy_A_adj
  PebX <- PebX * Q_energy_X_adj

  qebA <- qebA * Q_energy_A_adj
  qebX <- qebX * Q_energy_X_adj

  P <- Pdcow + PebA + PebX
  q <- qebA + qebX

  P <- P * outer(s.weights, s.weights)
  q <- q * s.weights

  Amat <- cbind(diag(n), s.weights)
  lvec <- c(rep(min.w, n), n)
  uvec <- c(ifelse(sw0, min.w, Inf), n)

  if (d.moments != 0) {
    d.covs <- covs
    d.treat <- cbind(poly(treat, degree = d.moments))

    if (d.moments > 1) {
      d.covs <- cbind(d.covs, int.poly.f(d.covs, poly = d.moments))
    }

    X.targets <- col.w.m(d.covs, s.weights)
    A.targets <- col.w.m(d.treat, s.weights)

    d.covs <- scale(d.covs, center = X.targets, scale = FALSE)
    d.treat <- scale(d.treat, center = A.targets, scale = FALSE)

    Amat <- cbind(Amat, d.covs * s.weights, d.treat * s.weights)
    lvec <- c(lvec, rep(0, ncol(d.covs)), rep(0, ncol(d.treat)))
    uvec <- c(uvec, rep(0, ncol(d.covs)), rep(0, ncol(d.treat)))
  }
  if (moments != 0 || int) {
    covs <- cbind(covs, int.poly.f(covs, poly = moments, int = int))

    X.means <- col.w.m(covs, s.weights)
    A.mean <- w.m(treat, s.weights)

    covs <- scale(covs, center = X.means, scale = FALSE)
    treat <- treat - A.mean

    Amat <- cbind(Amat, covs * treat * s.weights)

    lvec <- c(lvec, rep(0, ncol(covs)))
    uvec <- c(uvec, rep(0, ncol(covs)))
  }

  #Add weight penalty
  if (lambda < 0) {
    #Find lambda to make P PSD
    e <- eigen(P, symmetric = TRUE, only.values = TRUE)
    e.min <- min(e$values)

    lambda <- -e.min*n^2
  }

  diag(P) <- diag(P) + lambda / n^2

  if (is_not_null(A[["eps"]])) {
    if (is_null(A[["eps_abs"]])) A[["eps_abs"]] <- A[["eps"]]
    if (is_null(A[["eps_rel"]])) A[["eps_rel"]] <- A[["eps"]]
  }
  A[names(A) %nin% names(formals(osqp::osqpSettings))] <- NULL
  if (is_null(A[["max_iter"]])) A[["max_iter"]] <- 5e4L
  if (is_null(A[["eps_abs"]])) A[["eps_abs"]] <- 1e-8
  if (is_null(A[["eps_rel"]])) A[["eps_rel"]] <- 1e-6
  A[["verbose"]] <- TRUE
  options.list <- do.call(osqp::osqpSettings, A)

  verbosely({
    opt.out <- osqp::solve_osqp(P = 2 * P, q = q, A = t(Amat), l = lvec, u = uvec,
                                pars = options.list)
  }, verbose = verbose)

  if (identical(opt.out$info$status, "maximum iterations reached")) {
    .wrn("the optimization failed to converge. See Notes section at `?method_energy` for information")
  }

  w <- opt.out$x
  w[w <= min.w] <- min.w

  opt.out$lambda <- lambda

  list(w = w, fit.obj = opt.out)
}

#Versions that use formula and data instead of covs and treat
#
##Problems: not all functions process subset and na.action
##          CBPS cannot remove factors with 1 level
##          Using by yields error when by covariate is included in model formula
.weightit2glm <- function(covs, treat, s.weights, subset, estimand, focal,
                          stabilize, subclass, missing, .data, .formula, verbose, ...) {
  A <- list(...)

  fit.obj <- NULL

  treat_sub <- factor(treat[subset])
  na.action <- na.pass

  if (missing == "ind") {
    na.action <- .na.impute
    .formula <- .add_missing_indicators(.formula, .data, subset)
  }

  #Process link
  acceptable.links <- {
    if (missing == "saem") "logit"
    else expand_grid_string(c("", "br."), c("logit", "probit", "cloglog", "identity", "log", "cauchit"))
  }

  if (is_null(A[["link"]])) A$link <- acceptable.links[1]
  else {
    which.link <- acceptable.links[pmatch(A[["link"]], acceptable.links, nomatch = 0)][1]
    if (is.na(which.link)) {
      .err(sprintf("Only %s allowed as the link for binary treatments%",
                   word_list(acceptable.links, quotes = TRUE, is.are = TRUE),
                   if (missing == "saem") ' with `missing = "saem"`' else ""))
    }

    A[["link"]] <- which.link
  }

  use.br <- startsWith(A[["link"]], "br.")
  # use.bayes <- startsWith(A$link, "bayes.")
  if (use.br) A$link <- substr(A$link, 4, nchar(A$link))
  # else if (use.bayes) A$link <- substr(A$link, 7, nchar(A$link))

  t.lev <- get_treated_level(treat_sub)
  c.lev <- setdiff(levels(treat_sub), t.lev)

  ps <- make_df(levels(treat_sub), nrow = length(treat_sub))

  .formula <- update(.formula, sprintf("I(. == %s) ~ .", t.lev))

  treat <- model.response(model.frame(update(.formula, . ~ 1), data = .data))

  if (missing == "saem") {
    rlang::check_installed("misaem")
    if (is_null(A[["saem.method"]])) A[["saem.method"]] <- "map"
    if (is_null(A[["control"]])) A[["control"]] <- list()
    if (is_null(A[["control"]][["var_cal"]])) A[["control"]][["var_cal"]] <- FALSE
    if (is_null(A[["control"]][["ll_obs_cal"]])) A[["control"]][["ll_obs_cal"]] <- FALSE

    .data <- .data[subset,, drop = FALSE]
    withCallingHandlers({
      verbosely({
        fit <- misaem::miss.glm(.formula, data = .data, control = as.list(A[["control"]]))
      }, verbose = verbose)
    },
    warning = function(w) {
      if (conditionMessage(w) != "one argument not used by format '%i '") .wrn("(from misaem) ", w, tidy = FALSE)
      invokeRestart("muffleWarning")
    })

    mf <- model.frame(.formula, data = .data, na.action = na.pass,
                      drop.unused.levels = TRUE)

    .pred.data <- model.matrix(.formula, data = mf)

    if (attr(terms(mf), "intercept") == 1) {
      .pred.data <- .pred.data[, -1, drop = FALSE]
    }

    ps[[t.lev]] <- p.score <- drop(predict(fit, newdata = .pred.data, method = A[["saem.method"]]))
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
      mm <- model.matrix(.formula, data = .data)
      start <- setNames(rep(0, ncol(mm)), colnames(mm))
      if ("(Intercept)" %in% colnames(mm)) {
        start["(Intercept)"] <- family$linkfun(w.m(treat, s.weights))
      }
    }
    else {
      #Default starting values from glm.fit() without weights; these
      #work better with s.weights than usual default.
      mustart <- .25 + .5*treat
    }

    withCallingHandlers({verbosely({
      fit <- do.call(stats::glm, list(.formula, data = .data,
                                      weights = s.weights,
                                      mustart = mustart,
                                      start = start,
                                      family = family,
                                      method = glm_method,
                                      subset = subset,
                                      na.action = na.action,
                                      control = control), quote = TRUE)
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

  #ps should be matrix of probs for each treat
  #Computing weights
  w <- get_w_from_ps(ps = ps, treat = treat_sub, estimand, focal = focal,
                     stabilize = stabilize, subclass = subclass)

  list(w = w, ps = p.score, fit.obj = fit.obj)
}

.weightit2glm.multi <- function(covs, treat, s.weights, subset, estimand, focal,
                                stabilize, subclass, missing, .data, .formula, verbose, ...) {
  A <- list(...)

  fit.obj <- NULL

  treat_sub <- factor(treat[subset])
  s.weights <- s.weights[subset]

  ord.treat <- is.ordered(treat_sub)
  na.action <- na.pass

  if (missing == "ind") {
    na.action <- .na.impute
    .formula <- .add_missing_indicators(.formula, .data, subset)
  }

  #Process link
  acceptable.links <- {
    if (ord.treat) c("logit", "probit", "loglog", "cloglog", "cauchit", "br.logit")
    else if (!isFALSE(A$use.mlogit)) c("logit", "probit", "bayes.probit", "br.logit")
    else if (missing == "saem") "logit"
    else expand_grid_string(c("", "br."), c("logit", "probit", "cloglog", "identity", "log", "cauchit"))
  }

  if (is_null(A[["link"]])) A$link <- acceptable.links[1]
  else {
    which.link <- acceptable.links[pmatch(A[["link"]], acceptable.links, nomatch = 0)][1]
    if (is.na(which.link)) {
      .err(sprintf("Only %s allowed as the link for %s treatments%",
                   word_list(acceptable.links, quotes = TRUE, is.are = TRUE),
                   if (ord.treat) "ordinal" else "multinomial",
                   if (missing == "saem") ' with `missing = "saem"`' else ""))
    }

    A[["link"]] <- which.link
  }

  use.br <- startsWith(A[["link"]], "br.")
  # use.bayes <- startsWith(A$link, "bayes.")
  if (use.br) A$link <- substr(A$link, 4, nchar(A$link))
  # else if (use.bayes) A$link <- substr(A$link, 7, nchar(A$link))

  if (ord.treat) {
    if (use.br) {
      rlang::check_installed("brglm2")

      ctrl_fun <- brglm2::brglmControl
      control <- do.call(ctrl_fun, c(A[["control"]],
                                     A[setdiff(names(formals(ctrl_fun))[pmatch(names(A), names(formals(ctrl_fun)), 0)],
                                               names(A[["control"]]))]))

      tryCatch({verbosely({
        fit <- do.call(brglm2::bracl,
                       list(.formula,
                            data = .data,
                            weights = s.weights,
                            control = control,
                            na.action = na.action,
                            subset = subset,
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

      tryCatch({verbosely({
        fit <- do.call(MASS::polr,
                       list(.formula,
                            data = .data,
                            weights = s.weights,
                            subset = subset,
                            na.action = na.action,
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
  }
  else if (use.br) {
    rlang::check_installed("brglm2")

    ctrl_fun <- brglm2::brglmControl
    control <- do.call(ctrl_fun, c(A[["control"]],
                                   A[setdiff(names(formals(ctrl_fun))[pmatch(names(A), names(formals(ctrl_fun)), 0)],
                                             names(A[["control"]]))]))
    tryCatch({verbosely({
      fit <- do.call(brglm2::brmultinom,
                     list(.formula, .data,
                          weights = s.weights,
                          subset = subset,
                          na.action = na.action,
                          control = control), quote = TRUE)
    }, verbose = verbose)},
    error = function(e) {
      .err(sprintf("There was a problem with the bias-reduced multinomial logit regression. Try a different link.\n       Error message: (from brglm2) %s", conditionMessage(e)), tidy = FALSE)
    })

    ps <- fit$fitted.values
    fit.obj <- fit
  }
  else if (A$link == "bayes.probit") {
    rlang::check_installed("MNP")
    mf <- do.call("model.frame", list(.formula, .data, subset = subset,
                                      na.action = na.action))
    mm <- model.matrix(.formula, data = mf)
    if ("(Intercept)" %in% colnames(mm)) {
      mm <- mm[,colnames(mm) != "(Intercept)", drop = FALSE]
    }

    data <- cbind(model.response(mf), as.data.frame(mm))
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
  else if (missing == "saem") {
    rlang::check_installed("misaem")
    if (is_null(A[["saem.method"]])) A[["saem.method"]] <- "map"
    if (is_null(A[["control"]])) A[["control"]] <- list()
    if (is_null(A[["control"]][["var_cal"]])) A[["control"]][["var_cal"]] <- FALSE
    if (is_null(A[["control"]][["ll_obs_cal"]])) A[["control"]][["ll_obs_cal"]] <- FALSE

    ps <- make_df(levels(treat_sub), nrow = length(treat_sub))

    fit.list <- make_list(levels(treat_sub))

    .data <- .data[subset,, drop = FALSE]

    mf <- model.frame(.formula, data = .data, na.action = na.pass,
                      drop.unused.levels = TRUE)

    .pred.data <- model.matrix(.formula, data = mf)

    if (attr(terms(mf), "intercept") == 1) {
      .pred.data <- .pred.data[, -1, drop = FALSE]
    }

    for (i in levels(treat_sub)) {
      .formula_i <- update(.formula, sprintf("I(. == %s) ~ .", i))

      withCallingHandlers({
        verbosely({
          fit.list[[i]] <- do.call(misaem::miss.glm, list(.formula_i, data = .data, subset = subset,
                                                          control = A[["control"]]),
                                   quote = TRUE)
        }, verbose = verbose)
      },
      warning = function(w) {
        if (conditionMessage(w) != "one argument not used by format '%i '") .wrn("(from misaem) ", w, tidy = FALSE)
        invokeRestart("muffleWarning")
      })

      ps[[i]] <- drop(predict(fit, newdata = .data, method = A[["saem.method"]]))
    }
    fit.obj <- fit.list
  }
  else if (isTRUE(A$use.mclogit)) {
    rlang::check_installed("mclogit")
    stop("not ready for formula")
    if (is_not_null(A$random)) {
      random <- get_covs_and_treat_from_formula(A$random, data = .data)$reported.covs[subset,,drop = FALSE]
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

    tr_name <- deparse1(.formula[[2]])
    .env <- environment(.formula)
    .formula <- as.formula(paste(tr_name,
                                 "~ 1 |",
                                 deparse1(.formula[[3]])))
    environment(.formula) <- .env

    # browser()
    tryCatch({verbosely({
      fit <- do.call(mlogit::mlogit, list(.formula,
                                          data = .data,
                                          estimate = TRUE,
                                          probit = A$link[1] == "probit",
                                          weights = rep(s.weights, nlevels(treat_sub)),
                                          varying = NULL,
                                          shape = "wide",
                                          sep = "",
                                          choice = tr_name,
                                          subset = subset,
                                          na.action = na.action),
                     quote = FALSE)
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

      .formula_i <- update(.formula, sprintf("I(. == %s) ~ .", i))

      verbosely({
        fit.list[[i]] <- do.call(stats::glm, list(.formula_i, data = .data,
                                                  family = quasibinomial(link = A$link),
                                                  weights = s.weights,
                                                  subset = subset,
                                                  na.action = na.action,
                                                  control = control), quote = TRUE)
      }, verbose = verbose)

      ps[[i]] <- fit.list[[i]]$fitted.values

    }
    if (isTRUE(A[["test2"]])) ps <- ps/rowSums(ps)
    fit.obj <- fit.list
  }

  p.score <- NULL

  #ps should be matrix of probs for each treat
  #Computing weights
  w <- get_w_from_ps(ps = ps, treat = treat_sub, estimand, focal = focal,
                     stabilize = stabilize, subclass = subclass)

  list(w = w, ps = p.score, fit.obj = fit.obj)
}

.weightit2glm.cont <- function(covs, treat, s.weights, subset, stabilize, missing,
                               .data, .formula, verbose, ...) {
  A <- list(...)

  fit.obj <- NULL

  treat <- treat[subset]
  s.weights <- s.weights[subset]
  na.action <- na.pass

  if (missing == "ind") {
    na.action <- .na.impute
    .formula <- .add_missing_indicators(.formula, .data, subset)
  }

  #Process density params
  densfun <- get_dens_fun(use.kernel = isTRUE(A[["use.kernel"]]), bw = A[["bw"]],
                          adjust = A[["adjust"]], kernel = A[["kernel"]],
                          n = A[["n"]], treat = treat, density = A[["density"]])

  #Stabilization - get dens.num
  dens.num <- densfun(treat - mean(treat), s.weights)

  #Estimate GPS
  if (is_null(A[["link"]])) A[["link"]] <- "identity"

  if (missing == "saem") {
    rlang::check_installed("misaem")

    acceptable.links <- "identity"
    which.link <- acceptable.links[pmatch(A[["link"]], acceptable.links, nomatch = 0)][1]

    if (is.na(which.link)) {
      .err(sprintf("only %s allowed as the link for continuous treatments with missing = \"saem\"",
                   word_list(acceptable.links, quotes = TRUE, is.are = TRUE)))
    }

    .data <- .data[subset,, drop = FALSE]
    withCallingHandlers({verbosely({
      fit <- misaem::miss.lm(.formula, .data, control = as.list(A[["control"]]))
    }, verbose = verbose)},
    warning = function(w) {
      if (conditionMessage(w) != "one argument not used by format '%i '") {
        .wrn("(from `misaem::miss.lm()`) ", w, tidy = FALSE)
      }
      invokeRestart("muffleWarning")
    })

    if (is_null(A[["saem.method"]])) A[["saem.method"]] <- "map"

    mf <- model.frame(.formula, data = .data, na.action = na.pass,
                      drop.unused.levels = TRUE)

    .pred.data <- model.matrix(.formula, data = mf)

    if (attr(terms(mf), "intercept") == 1) {
      .pred.data <- .pred.data[, -1, drop = FALSE]
    }

    gp.score <- drop(predict(fit, newdata = .pred.data, method = A[["saem.method"]]))
  }
  else {
    acceptable.links <- c("identity", "log", "inverse")

    which.link <- acceptable.links[pmatch(A[["link"]], acceptable.links, nomatch = 0)][1]
    if (is.na(which.link)) {
      .err(sprintf("only %s allowed as the link for continuous treatments",
                   word_list(acceptable.links, quotes = TRUE, is.are = TRUE)))
    }

    A[["link"]] <- which.link

    verbosely({
      fit <- do.call("glm", c(list(.formula, data = .data,
                                   weights = s.weights,
                                   family = gaussian(link = A[["link"]]),
                                   subset = subset,
                                   na.action = na.action,
                                   control = as.list(A$control))),
                     quote = TRUE)
    }, verbose = verbose)

    gp.score <- fit$fitted.values
  }

  fit.obj <- fit

  #Get weights
  dens.denom <- densfun(treat - gp.score, s.weights)

  w <- dens.num/dens.denom

  if (isTRUE(A[["plot"]])) {
    d.n <- attr(dens.num, "density")
    d.d <- attr(dens.denom, "density")
    plot_density(d.n, d.d)
  }

  list(w = w, fit.obj = fit.obj)
}

### Note: bart can't handle na.action
.weightit2bart <- function(covs, treat, s.weights, subset, estimand, focal, stabilize, subclass, missing,
                           .data, .formula, verbose, ...) {
  A <- list(...)

  rlang::check_installed("dbarts", version = "0.9-23")

  treat <- factor(treat[subset])
  s.weights <- s.weights[subset]

  if (!all_the_same(s.weights)) {
    .err("sampling weights cannot be used with `method = \"bart\"`")
  }

  if (!has_treat_type(treat)) treat <- assign_treat_type(treat)
  treat.type <- get_treat_type(treat)

  na.action <- na.pass

  if (missing == "ind") {
    na.action <- .na.impute
    .formula <- .add_missing_indicators(.formula, .data, subset)
  }

  ps <- make_df(levels(treat), nrow = length(treat))

  A[["data"]] <- .data
  A[["subset"]] <- subset
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

    .formula_i <- update(.formula, sprintf("I(. == %s) ~ .", i))

    A[["formula"]] <- .formula_i

    tryCatch({verbosely({
      fit.list[[i]] <- do.call(dbarts::bart2,
                               A[names(A) %in% setdiff(c(names(formals(dbarts::bart2)),
                                                         names(formals(dbarts::dbartsControl))),
                                                       c("offset.test", "weights", "test"))],
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
  w <- get_w_from_ps(ps = ps, treat = treat, estimand, focal,
                     stabilize = stabilize, subclass = subclass)

  p.score <- if (treat.type == "binary") ps[[get_treated_level(treat)]] else NULL

  list(w = w, ps = p.score, info = info, fit.obj = fit.list)
}

# weightit2bart.multi <- weightit2bart

.weightit2bart.cont <- function(covs, treat, s.weights, subset, stabilize, missing, ps,
                                .data, .formula, verbose, ...) {
  A <- list(...)

  rlang::check_installed("dbarts")

  treat <- treat[subset]
  s.weights <- s.weights[subset]

  if (!all_the_same(s.weights)) {
    .err("sampling weights cannot be used with `method = \"bart\"`")
  }

  na.action <- na.pass

  if (missing == "ind") {
    na.action <- .na.impute
    .formula <- .add_missing_indicators(.formula, .data, subset)
  }

  #Process density params
  densfun <- get_dens_fun(use.kernel = isTRUE(A[["use.kernel"]]), bw = A[["bw"]],
                          adjust = A[["adjust"]], kernel = A[["kernel"]],
                          n = A[["n"]], treat = treat, density = A[["density"]])

  #Stabilization - get dens.num
  dens.num <- densfun(treat - mean(treat), s.weights)

  A[["formula"]] <- .formula
  A[["data"]] <- .data
  A[["subset"]] <- subset
  A[["keepCall"]] <- FALSE
  A[["combineChains"]] <- TRUE
  A[["verbose"]] <- FALSE #necessary to prevent crash

  #Estimate GPS

  tryCatch({verbosely({
    fit <- do.call(dbarts::bart2,
                   A[names(A) %in% setdiff(c(names(formals(dbarts::bart2)),
                                             names(formals(dbarts::dbartsControl))),
                                           c("offset.test", "weights", "test"))],
                   quote = TRUE)
  }, verbose = verbose)},
  error = function(e) {
    e. <- conditionMessage(e)
    .err("(from `dbarts::bart2()`) ", e., tidy = FALSE)
  })

  gp.score <- fitted(fit)

  #Get weights
  dens.denom <- densfun(treat - gp.score, s.weights)

  w <- dens.num/dens.denom

  if (isTRUE(A[["use.kernel"]]) && isTRUE(A[["plot"]])) {
    d.n <- attr(dens.num, "density")
    d.d <- attr(dens.denom, "density")
    plot_density(d.n, d.d)
  }

  info <- list()

  list(w = w, info = info, fit.obj = fit)
}

.weightit2cbps <- function(covs, treat, s.weights, estimand, focal, subset,
                          stabilize, subclass, missing, .data, .formula, verbose, ...) {

  rlang::check_installed("CBPS")

  A <- list(...)

  treat <- factor(treat[subset])
  na.action <- na.pass

  if (missing == "ind") {
    na.action <- .na.impute
    .formula <- .add_missing_indicators(.formula, .data, subset)
  }

  sw0 <- check_if_zero(s.weights)

  if (estimand == "ATT") {
    .formula <- update(.formula, sprintf("I(. == %s) ~ .", focal))
  }
  else {
    .formula <- update(.formula, sprintf("I(. == %s) ~ .", get_treated_level(treat)))
  }

  tryCatch({verbosely({
    # fit <- CBPS::CBPS(.formula,
    #                   data = .data[!sw0 & subset,],
    #                   method = if (is_not_null(A$over) && isFALSE(A$over)) "exact" else "over",
    #                   standardize = FALSE,
    #                   na.action = na.action,
    #                   sample.weights = s.weights[!sw0 & subset],
    #                   ATT = if (estimand = "ATT") 1 else 0,
    #                   ...)
    fit <- do.call(CBPS_wrapper, list(.formula,
                                      data = .data,
                                      method = if (is_not_null(A$over) && isFALSE(A$over)) "exact" else "over",
                                      standardize = FALSE,
                                      na.action = na.action,
                                      subset = subset,
                                      weights = s.weights,
                                      ATT = if (estimand == "ATT") 1 else 0,
                                      ...), quote = TRUE)
  }, verbose = verbose)},
  error = function(e) {
    e. <- conditionMessage(e)
    e. <- gsub("method = \"exact\"", "`over = FALSE`", e., fixed = TRUE)
    .err("(from `CBPS::CBPS()`) ", e., tidy = FALSE)
  })

  if (!any(sw0)) {
    ps <- fit[["fitted.values"]]
  }
  else {
    mf <- do.call("model.frame", list(.formula, data = .data, subset = subset),
                  quote = TRUE)

    ps <- plogis(drop(model.matrix(.formula, data = mf) %*% fit[["coefficients"]]))
  }

  w <- get_w_from_ps(ps, treat, estimand = estimand, subclass = subclass,
                     focal = focal, stabilize = stabilize)

  p.score <- {
    if (is_null(dim(ps)) || length(dim(ps)) != 2) ps
    else ps[[get_treated_level(treat)]]
  }

  list(w = w, ps = p.score, fit.obj = fit)
}

.weightit2cbps.multi <- function(covs, treat, s.weights, estimand, focal, subset,
                                 stabilize, subclass, missing, .data, .formula, verbose, ...) {

  rlang::check_installed("CBPS")

  A <- list(...)

  treat <- factor(treat[subset])
  s.weights <- s.weights[subset]
  na.action <- na.pass

  if (missing == "ind") {
    na.action <- .na.impute
    .formula <- .add_missing_indicators(.formula, .data, subset)
  }

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
                                    na.action = na.action,
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
    if (nlevels(treat) <= 4) {
      tryCatch({verbosely({
        fit.list <- CBPS::CBPS(formula(new.data),
                               data = new.data[!sw0,],
                               method = if (isFALSE(A$over)) "exact" else "over",
                               standardize = FALSE,
                               na.action = na.action,
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
          fit.list[[i]] <- CBPS::CBPS(formula(new.data),
                                      data = new.data[!sw0,],
                                      method = if (isFALSE(A$over)) "exact" else "over",
                                      standardize = FALSE,
                                      na.action = na.action,
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

  p.score <- NULL

  list(w = w, ps = p.score, fit.obj = fit.list)
}

.weightit2cbps.cont <- function(covs, treat, s.weights, subset, missing, .data, .formula, verbose, ...) {
  rlang::check_installed("CBPS")

  A <- list(...)

  treat <- treat[subset]
  s.weights <- s.weights[subset]
  na.action <- na.pass

  if (missing == "ind") {
    na.action <- .na.impute
    .formula <- .add_missing_indicators(.formula, .data, subset)
  }

  sw0 <- check_if_zero(s.weights)

  new.data <- data.frame(treat = treat, covs)

  w <- rep(0, length(treat))

  tryCatch({verbosely({
    fit <- CBPS::CBPS(formula(new.data),
                      data = new.data[!sw0,],
                      method = if (isFALSE(A$over)) "exact" else "over",
                      standardize = FALSE,
                      na.action = na.action,
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

CBPS_wrapper <- function(formula, data, na.action, subset, ATT = 1, iterations = 1000,
                         standardize = TRUE, method = "over", twostep = TRUE, weights = NULL,
                         baseline.formula = NULL, diff.formula = NULL, ...) {

  if (missing(data))
    data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  # mf$drop.unused.levels <- TRUE
  mf$drop.unused.levels <- FALSE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  Y <- model.response(mf, "any")
  if (length(dim(Y)) == 1L) {
    nm <- rownames(Y)
    dim(Y) <- NULL
    if (!is.null(nm))
      names(Y) <- nm
  }
  browser()
  X <- if (!is.empty.model(mt))
    model.matrix(mt, mf, contrasts = NULL)
  else matrix(, NROW(Y), 0L)
  X <- cbind(1, X[, apply(X, 2, sd) > 0])
  sample.weights <- as.vector(model.weights(mf))
  if (is.null(sample.weights))
    sample.weights <- rep(1, nrow(X))
  if (xor(is.null(baseline.formula), is.null(diff.formula))) {
    stop("Either baseline.formula or diff.formula not specified.  Both must be specified to use CBPSOptimal.  Otherwise, leave both NULL.")
  }
  if (!is.null(baseline.formula)) {
    baselineX <- model.matrix(terms(baseline.formula))
    baselineX <- baselineX[, apply(baselineX, 2, sd) > 0]
    diffX <- model.matrix(terms(diff.formula))
    diffX <- diffX[, apply(as.matrix(diffX), 2, sd) > 0]
  }
  else {
    baselineX <- NULL
    diffX <- NULL
  }
  fit <- CBPS::CBPS.fit(X = X, treat = Y, ATT = ATT,
                        intercept = attr(mt, "intercept") > 0L, method = method,
                        iterations = iterations, standardize = standardize, twostep = twostep,
                        baselineX = baselineX, diffX = diffX, sample.weights = sample.weights)

  fit$na.action <- attr(mf, "na.action")
  xlevels <- .getXlevels(mt, mf)
  fit$data <- data
  fit$call <- call
  fit$formula <- formula
  fit$terms <- mt
  fit
}

#########
.add_missing_indicators <- function(formula, data, subset) {
  tt <- terms(formula, data = data)
  mf <- model.frame(tt, data = data, na.action = na.pass)

  if (!anyNA(mf)) return(formula)

  for (i in which(vapply(mf, anyNA, logical(1L)))) {
    formula <- update(formula, sprintf(". ~ . + is.na(%s)", names(mf)[i]))
  }

  formula
}

.na.impute <- function(object, fun = median, ...) {

  fun <- match.fun(fun)
  vars <- names(object)

  for (j in vars[vapply(object, anyNA, logical(1L))]) {
    x <- object[[j]]
    if (!is.atomic(x)) next
    if (is.numeric(x)) {
      x[is.na(x)] <- fun(x, na.rm = TRUE)
    }
    else {
      x[is.na(x)] <- x[!is.na(x)][1]
    }
    object[[j]] <- x
  }
  object
}


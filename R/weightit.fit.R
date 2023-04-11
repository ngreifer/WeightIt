weightit.fit <- function(covs, treat, method = "glm", s.weights = NULL, by.factor = NULL,
                         estimand = "ATE", focal = NULL, stabilize = FALSE, ps = NULL, moments = NULL, int = FALSE,
                         subclass = NULL, is.MSM.method = FALSE, missing = NULL, verbose = FALSE, include.obj = FALSE, ...){

  #main function of weightit that dispatches to weightit2method and returns object containing weights and ps

  #Checks
  if (!check_if_call_from_fun(weightit) && !check_if_call_from_fun(weightitMSM)) {

    chk::chk_not_missing(covs, "`covs`")
    chk::chk_matrix(covs)
    chk::chk_numeric(covs)

    chk::chk_not_missing(treat, "`treat`")
    chk::chk_vector(treat)
    chk::chk_numeric(treat)
    chk::chk_not_any_na(treat)

    if (length(treat) != nrow(covs)) {
      .err("`treat` and `'covs` must contain the same number of units")
    }

    if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
    treat.type <- get.treat.type(treat)

    check.acceptable.method(method, msm = FALSE, force = FALSE)

    if (is_not_null(ps)) {
      chk::chk_vector(ps)
      chk::chk_numeric(ps)

      if (length(ps) != length(treat)) {
        .err("`ps` and `treat` must be the same length")
      }

      method <- "glm"
    }

    if (is_null(s.weights)) s.weights <- rep(1, length(treat))
    else {
      chk::chk_vector(s.weights)
      chk::chk_numeric(s.weights)

      if (length(s.weights) != length(treat)) {
        .err("`s.weights` and `treat` must be the same length")
      }
    }

    if (is_null(by.factor)) by.factor <- factor(rep(1, length(treat)), levels = 1)
    else {
      chk::chk_factor(by.factor)

      if (length(by.factor) != length(treat)) {
        .err("`by.factor` and `treat` must be the same length")
      }
    }

    #Process estimand and focal
    estimand <- process.estimand(estimand, method, treat.type)
    f.e.r <- process.focal.and.estimand(focal, estimand, treat)
    focal <- f.e.r[["focal"]]
    estimand <- f.e.r[["estimand"]]

    .chk_null_or(missing, chk::chk_string)

    if (is_null(missing)) {
      if (anyNA(covs)) missing <- "ind"
      else missing <- ""
    }
    else if (missing != "" && anyNA(covs)) {
      missing <- process.missing(missing, method, treat.type)
    }
    else missing <- ""

    #Check subclass
    if (is_not_null(subclass)) check.subclass(method, treat.type)

    #Process moments and int
    moments.int <- process.moments.int(moments, int, method)
    moments <- moments.int[["moments"]]; int <- moments.int[["int"]]
  }
  else {
    if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
    treat.type <- get.treat.type(treat)
  }

  out <- make_list(c("weights", "ps", "fit.obj", "info"))
  out$weights <- out$ps <- rep(NA_real_, length(treat))

  if (include.obj) fit.obj <- make_list(levels(by.factor))
  info <- make_list(levels(by.factor))

  obj <- NULL

  for (i in levels(by.factor)) {
    #Run method
    if (is.function(method)) {
      if (is.MSM.method) {
        obj <- weightitMSM2user(Fun = method,
                                covs.list = covs,
                                treat.list = treat,
                                s.weights = s.weights,
                                subset = by.factor == i,
                                stabilize = stabilize,
                                missing = missing,
                                verbose = verbose,
                                ...)
      }
      else {
        obj <- weightit2user(Fun = method,
                             covs = covs,
                             treat = treat,
                             s.weights = s.weights,
                             subset = by.factor == i,
                             estimand = estimand,
                             focal = focal,
                             stabilize = stabilize,
                             subclass = subclass,
                             ps = ps,
                             missing = missing,
                             moments = moments,
                             int = int,
                             verbose = verbose,
                             ...)
      }
    }
    else if (method == "glm") {
      if (treat.type %in% c("binary", "multinomial")) {
        obj <- weightit2glm(covs = covs,
                            treat = treat,
                            s.weights = s.weights,
                            subset = by.factor == i,
                            estimand = estimand,
                            focal = focal,
                            stabilize = stabilize,
                            subclass = subclass,
                            ps = ps,
                            missing = missing,
                            verbose = verbose,
                            ...)
      }
      else if (treat.type == "continuous") {
        obj <- weightit2glm.cont(covs = covs,
                                 treat = treat,
                                 s.weights = s.weights,
                                 subset = by.factor == i,
                                 ps = ps,
                                 missing = missing,
                                 verbose = verbose,
                                 ...)
      }
    }
    else if (method == "optweight") {
      if (is.MSM.method) {
        obj <- weightit2optweight.msm(covs.list = covs,
                                      treat.list = treat,
                                      s.weights = s.weights,
                                      subset = by.factor == i,
                                      moments = moments,
                                      int = int,
                                      missing = missing,
                                      verbose = verbose,
                                      ...)
      }
      else {
        if (treat.type %in% c("binary", "multinomial")) {
          obj <- weightit2optweight(covs = covs,
                                    treat = treat,
                                    s.weights = s.weights,
                                    subset = by.factor ==i,
                                    estimand = estimand,
                                    focal = focal,
                                    moments = moments,
                                    int = int,
                                    missing = missing,
                                    verbose = verbose,
                                    ...)
        }
        else if (treat.type == "continuous") {
          obj <- weightit2optweight.cont(covs = covs,
                                         treat = treat,
                                         subset = by.factor == i,
                                         s.weights = s.weights,
                                         moments = moments,
                                         int = int,
                                         missing = missing,
                                         verbose = verbose,
                                         ...)

        }
      }
    }
    else if (method == "gbm") {
      if (treat.type %in% c("binary", "multinomial")) {
        obj <- weightit2gbm(covs = covs,
                            treat = treat,
                            s.weights = s.weights,
                            estimand = estimand,
                            focal = focal,
                            subset = by.factor == i,
                            stabilize = stabilize,
                            subclass = subclass,
                            missing = missing,
                            verbose = verbose,
                            ...)
      }
      else {
        obj <- weightit2gbm.cont(covs = covs,
                                 treat = treat,
                                 s.weights = s.weights,
                                 subset = by.factor == i,
                                 missing = missing,
                                 verbose = verbose,
                                 ...)
      }

    }
    else if (method == "cbps") {
      if (is.MSM.method) {
        # obj <- weightit2cbps.msm()
      }
      else {
        if (treat.type %in% c("binary", "multinomial")) {
          obj <- weightit2cbps(covs = covs,
                               treat = treat,
                               subset = by.factor == i,
                               s.weights = s.weights,
                               stabilize = stabilize,
                               subclass = subclass,
                               estimand = estimand,
                               focal = focal,
                               missing = missing,
                               verbose = verbose,
                               ...)
        }
        else if (treat.type == "continuous") {
          obj <- weightit2cbps.cont(covs = covs,
                                    treat = treat,
                                    subset = by.factor == i,
                                    s.weights = s.weights,
                                    missing = missing,
                                    verbose = verbose,
                                    ...)

        }
      }

    }
    else if (method == "npcbps") {
      if (treat.type %in% c("binary", "multinomial")) {
        obj <- weightit2npcbps(covs = covs,
                               treat = treat,
                               subset = by.factor == i,
                               s.weights = s.weights,
                               moments = moments,
                               int = int,
                               missing = missing,
                               verbose = verbose,
                               ...)
      }
      else if (treat.type == "continuous") {
        obj <- weightit2npcbps.cont(covs = covs,
                                    treat = treat,
                                    subset = by.factor == i,
                                    s.weights = s.weights,
                                    moments = moments,
                                    int = int,
                                    missing = missing,
                                    verbose = verbose,
                                    ...)
      }

    }
    else if (method == "ebal") {
      if (treat.type %in% c("binary", "multinomial")) {
        obj <- weightit2ebal(covs = covs,
                             treat = treat,
                             s.weights = s.weights,
                             subset = by.factor == i,
                             estimand = estimand,
                             focal = focal,
                             stabilize = stabilize,
                             moments = moments,
                             int = int,
                             missing = missing,
                             verbose = verbose,
                             ...)
      }
      else {
        obj <- weightit2ebal.cont(covs = covs,
                                  treat = treat,
                                  subset = by.factor == i,
                                  s.weights = s.weights,
                                  moments = moments,
                                  int = int,
                                  missing = missing,
                                  verbose = verbose,
                                  ...)
      }
    }
    else if (method == "super") {
      if (treat.type %in% c("binary", "multinomial")) {
        obj <- weightit2super(covs = covs,
                              treat = treat,
                              s.weights = s.weights,
                              subset = by.factor == i,
                              estimand = estimand,
                              focal = focal,
                              stabilize = stabilize,
                              subclass = subclass,
                              missing = missing,
                              verbose = verbose,
                              ...)
      }
      else if (treat.type == "continuous") {
        obj <- weightit2super.cont(covs = covs,
                                   treat = treat,
                                   s.weights = s.weights,
                                   subset = by.factor == i,
                                   stabilize = stabilize,
                                   missing = missing,
                                   verbose = verbose,
                                   ...)
      }
    }
    # else if (method == "ebcw") {
    #   if (treat.type %in% c("binary", "multinomial")) {
    #     obj <- weightit2ebcw(covs = covs,
    #                          treat = treat,
    #                          s.weights = s.weights,
    #                          subset = by.factor == i,
    #                          estimand = estimand,
    #                          focal = focal,
    #                          #stabilize = stabilize,
    #                          moments = moments,
    #                          int = int,
    #                          missing = missing,
    #                          ...)
    #   }
    #   else {
    #     stop("Empirical balancing calibration weights are not compatible with continuous treatments")
    #   }
    # }
    else if (method == "energy") {
      if (treat.type %in% c("binary", "multinomial")) {
        obj <- weightit2energy(covs = covs,
                               treat = treat,
                               s.weights = s.weights,
                               subset = by.factor == i,
                               estimand = estimand,
                               focal = focal,
                               stabilize = stabilize,
                               missing = missing,
                               moments = moments,
                               int = int,
                               verbose = verbose,
                               ...)
      }
      else {
        obj <- weightit2energy.cont(covs = covs,
                                    treat = treat,
                                    subset = by.factor == i,
                                    s.weights = s.weights,
                                    missing = missing,
                                    moments = moments,
                                    int = int,
                                    verbose = verbose,
                                    ...)
      }
    }
    else if (method == "bart") {
      if (treat.type %in% c("binary", "multinomial")) {
        obj <- weightit2bart(covs = covs,
                             treat = treat,
                             s.weights = s.weights,
                             subset = by.factor == i,
                             estimand = estimand,
                             focal = focal,
                             stabilize = stabilize,
                             subclass = subclass,
                             missing = missing,
                             verbose = verbose,
                             ...)
      }
      else if (treat.type == "continuous") {
        obj <- weightit2bart.cont(covs = covs,
                                  treat = treat,
                                  s.weights = s.weights,
                                  subset = by.factor == i,
                                  stabilize = stabilize,
                                  missing = missing,
                                  verbose = verbose,
                                  ...)
      }
    }
    # else if (method == "kbal") {
    #   if (treat.type %in% c("binary", "multinomial")) {
    #     obj <- weightit2kbal(covs = covs,
    #                          treat = treat,
    #                          s.weights = s.weights,
    #                          subset = by.factor == i,
    #                          estimand = estimand,
    #                          focal = focal,
    #                          missing = missing,
    #                          ...)
    #   }
    #   else stop("Kernel balancing is not compatible with continuous treatments")
    # }
    else {
      .err("invalid argument to `method`")
    }

    #Extract weights
    if (is_null(obj)) {
      .err("no object was created. This is probably a bug,\n     and you should report it at https://github.com/ngreifer/WeightIt/issues")
    }
    if (is_null(obj$w) || all(is.na(obj$w))) {
      .wrn("No weights were estimated. This is probably a bug,\n     and you should report it at https://github.com/ngreifer/WeightIt/issues")
    }
    if (any(!is.finite(obj$w))) {
      .wrn("Some weights were estimated as `NA`, which means a value was impossible to compute (e.g., Inf). Check for extreme values of the treatment or covariates and try removing them. Non-finite weights will be set to 0")
      obj$w[!is.finite(obj$w)] <- 0
    }
    # else if (any(!is.finite(obj$w))) probably.a.bug()

    out$weights[by.factor == i] <- obj$w
    if (is_not_null(obj$ps)) out$ps[by.factor == i] <- obj$ps

    if (include.obj) fit.obj[[i]] <- obj$fit.obj
    info[[i]] <- obj$info

  }


  if (include.obj) {
    if (nlevels(by.factor) == 1) fit.obj <- fit.obj[[1]]
    out$fit.obj <- fit.obj
  }

  if (is_not_null(info) && nlevels(by.factor) == 1) info <- info[[1]]
  out$info <- info

  if (all(is.na(out$ps))) out$ps <- NULL

  class(out) <- "weightit.fit"

  out
}

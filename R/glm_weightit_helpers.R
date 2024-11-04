# Makes a phrase from vcov_type to be printed by print() and summary()
.vcov_to_phrase <- function(vcov_type, cluster = FALSE) {
  switch(vcov_type,
         "const" = "maximum likelihood",
         "asympt" = {
           if (!cluster)
             "HC0 robust (adjusted for estimation of weights)"
           else
             "HC0 cluster-robust (adjusted for estimation of weights)"
         },
         "HC0" = {
           if (!cluster)
             "HC0 robust"
           else
             "HC0 cluster-robust"
         },
         "BS" = {
           if (!cluster)
             "traditional bootstrap"
           else
             "traditional cluster bootstrap"
         },
         "FWB" = {
           if (!cluster)
             "fractional weighted bootstrap"
           else
             "fractional weighted cluster bootstrap"
         },
         "none" = "none",
         "user-supplied")
}

# Processes the vcov argument from glm_weightit(); returns vcov_type with attributes
.process_vcov <- function(vcov = NULL, weightit = NULL, R, fwb.args = list(),
                          m_est_supported = TRUE) {
  if (is_null(weightit)) {
    if (is_null(vcov)) {
      return("HC0")
    }

    allowable_vcov <- c("none", "const", "HC0", "BS", "FWB")
  }
  else {
    chk::chk_is(weightit, "weightit")

    if (!m_est_supported ||
        (is_null(attr(weightit, "Mparts", exact = TRUE)) &&
         is_null(attr(weightit, "Mparts.list", exact = TRUE)))) {

      if (is_null(vcov)) {
        return("HC0")
      }

      allowable_vcov <- c("none", "const", "HC0", "BS", "FWB")
    }
    else {
      if (is_null(vcov)) {
        vcov <- "asympt"
      }

      allowable_vcov <- c("none", "const", "asympt", "HC0", "BS", "FWB")
    }
  }

  chk::chk_string(vcov)
  vcov <- match_arg(vcov, allowable_vcov)

  if (is_not_null(weightit) && vcov == "const") {
    .wrn('`vcov = "const"` should not be used when `weightit` is supplied; the resulting standard errors are invalid and should not be interpreted')
  }

  bootstrap <- vcov %in% c("BS", "FWB")

  if (bootstrap) {
    chk::chk_count(R)
    chk::chk_gt(R, 0)
    attr(vcov, "R") <- R

    if (vcov == "FWB") {
      #Error for weighting methods that don't accept s.weights
      if (is_not_null(weightit)) {
        if (is.character(weightit$method) &&
            !.weightit_methods[[weightit$method]]$s.weights_ok) {
          .err(sprintf('`vcov = "FWB"` cannot be used with `method = %s`',
                       add_quotes(weightit$method)))
        }

        if (identical(weightit$missing, "saem")) {
          .err('`vcov = "FWB"` cannot be used with `missing = "saem"`')
        }
      }

      rlang::check_installed("fwb")

      chk::chk_list(fwb.args)
    }

    attr(vcov, "fwb.args") <- fwb.args
  }

  vcov
}

# Process a user-supplied vcov for use in summary.glm_weightit()
.process_vcov_summary <- function(object, vcov. = NULL, ...) {
  if (is_null(vcov.)) {
    V <- stats::vcov(object)
    vcov_type <- object$vcov_type
    cluster <- attr(object, "cluster")
  }
  else if (is.function(vcov.)) {
    V <- vcov.(object, ...)

    if (!is.cov_like(V)) {
      .err("the function supplied to `vcov` must return a variance matrix (i.e., a square, symmetric matrix with a non-negative diagonal")
    }

    vcov_type <- "custom"
    cluster <- NULL
  }
  else if (is.cov_like(vcov.)) {
    V <- vcov.

    vcov_type <- "custom"
    cluster <- NULL
  }
  else if (chk::vld_string(vcov.)) {
    V <- .vcov_glm_weightit.internal(object, vcov. = vcov., ...)

    vcov_type <- attr(V, "vcov_type")
    cluster <- attr(V, "cluster")
  }
  else {
    .err("`vcov` must be `NULL`, a string corresponding to a method of computing the parameter variance matrix, a function that computes a variance matrix, or a variance matrix itself")
  }

  .modify_vcov_info(V, vcov_type = vcov_type, cluster = cluster)
}

# Process a user-supplied vcov for use in anova.glm_weightit()
.process_vcov_anova <- function(object, object2 = NULL, vcov. = NULL, b1 = NULL, ...) {

  if (is_null(vcov.)) {
    if (identical(object[["vcov_type"]], "none") || is_null(stats::vcov(object))) {
      .err("no variance matrix found in `object`; please supply an argument to `vcov`")
    }

    if (!identical(object[["vcov_type"]], object2[["vcov_type"]]) &&
        !identical(object2[["vcov_type"]], "none")) {
      .wrn("different `vcov` types detected for each model; using the variance matrix from the larger model")
    }

    V <- stats::vcov(object)
    vcov_type <- object$vcov_type
    cluster <- attr(object, "cluster")
  }
  else if (is.function(vcov.)) {
    V <- vcov.(object, ...)

    if (!is.cov_like(V)) {
      .err("the function supplied to `vcov` must return a variance matrix (i.e., a square, symmetric matrix with a non-negative diagonal")
    }

    vcov_type <- "custom"
    cluster <- NULL
  }
  else if (is.cov_like(vcov.)) {
    V <- vcov.

    vcov_type <- "custom"
    cluster <- NULL
  }
  else if (chk::vld_string(vcov.)) {
    V <- .vcov_glm_weightit.internal(object, vcov. = vcov., ...)

    vcov_type <- attr(V, "vcov_type")
    cluster <- attr(V, "cluster")
  }
  else {
    .err("`vcov` must be `NULL`, a string corresponding to a method of computing the parameter variance matrix, a function that computes a variance matrix, or a variance matrix itself")
  }

  if (is_null(V)) {
    .err("no variance matrix was found. See the `vcov` argument at `help(\"anova.glm_weightit\")` for details")
  }

  if (is_not_null(b1)) {
    if (!all(names(b1) %in% rownames(V)) || !all(names(b1) %in% colnames(V))) {
      .err("all coefficients in `object` must have entries in the supplied variance matrix")
    }

    V <- V[names(b1), names(b1), drop = FALSE]
  }

  .modify_vcov_info(V, vcov_type = vcov_type, cluster = cluster)
}

# Sets the vcov, vcov_type, and cluster components when given a user-supplied vcov;
# for use in summary.glm_weightit()
.set_vcov <- function(object, vcov, vcov_type = NULL) {
  object$vcov <- vcov

  object$vcov_type <- {
    if (is_null(vcov)) "none"
    else if_null_then(vcov_type, attr(vcov, "vcov_type"))
  }

  attr(object, "cluster") <- attr(vcov, "cluster")

  object$vcov <- .modify_vcov_info(object[["vcov"]])

  object
}

# Dispatches computation of vcov; for using in vcov.glm_weightit() and summary.glm_weightit()
.vcov_glm_weightit.internal <- function(object, vcov. = NULL, ...) {
  dots <- list(...)

  if (is_null(vcov.)) {
    vcov. <- dots[["type"]]
    dots[["type"]] <- NULL
  }

  if (is_null(vcov.) && is_null(dots)) {
    if (is_not_null(object[["vcov"]])) {
      return(.modify_vcov_info(object[["vcov"]],
                               vcov_type = object[["vcov_type"]],
                               cluster = attr(object, "cluster")))
    }

    if (!identical(object[["vcov_type"]], "none")) {
      .err("no variance-covariance matrix was found in the supplied object; this is likely a bug")
    }

    .wrn('`vcov` was specified as `"none"` in the original fitting call, so no variance-covariance matrix will be returned')

    return(NULL)
  }

  R <- if_null_then(dots[["R"]], 500)
  fwb.args <- if_null_then(dots[["fwb.args"]], list())

  vcov. <- .process_vcov(vcov., object[["weightit"]], R, fwb.args)

  if (vcov. == "none") {
    return(NULL)
  }

  cluster <- {
    if ("cluster" %in% names(dots)) dots[["cluster"]]
    else attr(object, "cluster")
  }

  glm_call <- .build_model_call(object, vcov = vcov.)

  .compute_vcov(object, object[["weightit"]], vcov., cluster, glm_call)
}

.modify_vcov_info <- function(vcov, vcov_type = NULL, cluster = NULL) {
  attr(vcov, "vcov_type") <- vcov_type
  attr(vcov, "cluster") <- cluster

  vcov
}

# Custom printCoefmat for glm_weightit objects
.printCoefmat_glm_weightit <- function(x,
                                       digits = max(3L, getOption("digits") - 2L),
                                       signif.stars = TRUE,
                                       signif.legend = FALSE,
                                       dig.tst = max(1L, min(5L, digits - 1L)),
                                       cs.ind = NULL,
                                       tst.ind = NULL,
                                       p.ind = NULL,
                                       zap.ind = integer(),
                                       P.values = NULL,
                                       has.Pvalue = TRUE,
                                       eps.Pvalue = 1e-6,
                                       na.print = ".",
                                       quote = FALSE,
                                       right = TRUE,
                                       ...) {
  if (is_null(d <- dim(x)) || length(d) != 2L) {
    .err("'x' must be coefficient matrix/data frame")
  }

  nm <- colnames(x)

  chk::chk_flag(has.Pvalue)

  if (has.Pvalue) {
    if (is_null(p.ind)) {
      if (is_null(nm)) {
        .err("`has.Pvalue` set to `TRUE` but `p.ind` is `NULL` and no colnames present")
      }

      p.ind <- which(substr(nm, 1L, 3L) %in% c("Pr(", "p-v"))

      if (is_null(p.ind)) {
        .err("`has.Pvalue` set to `TRUE` but `p.ind` is `NULL` and no colnames match p-value strings")
      }
    }
    else {
      chk::chk_whole_number(p.ind)
      chk::chk_subset(p.ind, seq_col(x))
    }
  }
  else {
    if (is_not_null(p.ind)) {
      .err("`has.Pvalue` set to `FALSE` but `p.ind` is not `NULL`")
    }
  }

  if (is_null(P.values)) {
    scp <- getOption("show.coef.Pvalues")
    if (!is.logical(scp) || is.na(scp)) {
      .wrn("option \"show.coef.Pvalues\" is invalid: assuming `TRUE`")
      scp <- TRUE
    }
    P.values <- has.Pvalue && scp
  }
  else {
    chk::chk_flag(P.values)
  }

  if (P.values && !has.Pvalue) {
    .err("`P.values` is `TRUE`, but `has.Pvalue` is not")
  }

  if (is_null(cs.ind)) {
    cs.ind <- which(nm %in% c("Estimate", "Std. Error") | endsWith(nm, " %"))
  }
  else {
    chk::chk_whole_numeric(cs.ind)
    chk::chk_subset(cs.ind, seq_col(x))
  }

  cs.ind <- setdiff(cs.ind, p.ind)

  if (is_null(tst.ind)) {
    tst.ind <- which(endsWith(nm, "value"))
  }
  else {
    chk::chk_whole_numeric(tst.ind)
    chk::chk_subset(tst.ind, seq_col(x))
  }

  tst.ind <- setdiff(tst.ind, p.ind)

  if (any(tst.ind %in% cs.ind)) {
    .err("`tst.ind` must not overlap with `cs.ind`")
  }

  xm <- data.matrix(x)

  if (is_null(tst.ind)) {
    tst.ind <- setdiff(which(endsWith(nm, "value")), p.ind)
  }

  Cf <- array("", dim = d, dimnames = dimnames(xm))

  ok <- !(ina <- is.na(xm))

  for (i in zap.ind) {
    xm[, i] <- zapsmall(xm[, i], digits)
  }

  if (is_not_null(cs.ind)) {
    acs <- abs(coef.se <- xm[, cs.ind, drop = FALSE])
    if (any(ia <- is.finite(acs))) {
      digmin <- 1 + if (length(acs <- acs[ia & acs != 0]))
        floor(log10(range(acs[acs != 0], finite = TRUE)))
      else 0

      Cf[, cs.ind] <- format(round(coef.se, max(1L, digits - digmin)), digits = digits)
    }
  }

  if (is_not_null(tst.ind)) {
    Cf[, tst.ind] <- format(round(xm[, tst.ind], digits = dig.tst),
                            digits = digits)
  }

  if (is_not_null(r.ind <- setdiff(seq_col(xm), c(cs.ind, tst.ind, p.ind)))) {
    for (i in r.ind) {
      Cf[, i] <- format(xm[, i], digits = digits)
    }
  }

  ok[, tst.ind] <- FALSE
  okP <- if (has.Pvalue) ok[, -p.ind] else ok

  x1 <- Cf[okP]
  dec <- getOption("OutDec")

  if (dec != ".") {
    x1 <- chartr(dec, ".", x1)
  }

  x0 <- (xm[okP] == 0) != (as.numeric(x1) == 0)

  if (is_not_null(not.both.0 <- which(x0 & !is.na(x0)))) {
    Cf[okP][not.both.0] <- format(xm[okP][not.both.0], digits = max(1L, digits - 1L))
  }

  if (any(ina)) {
    Cf[ina] <- na.print
  }

  if (any(inan <- is.nan(xm))) {
    Cf[inan] <- "NaN"
  }

  if (P.values) {
    chk::chk_flag(signif.stars)

    if (any(okP <- ok[, p.ind])) {
      pv <- as.vector(xm[, p.ind])
      Cf[okP, p.ind] <- format.pval(pv[okP], digits = dig.tst,
                                    eps = eps.Pvalue)

      signif.stars <- signif.stars && any(pv[okP] < 0.1)

      if (signif.stars) {
        Signif <- symnum(pv, corr = FALSE, na = FALSE,
                         cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                         symbols = c("***", "**", "*", ".", " "))
        Cf <- cbind(Cf, format(Signif))
      }
    }
    else signif.stars <- FALSE
  }
  else {
    if (has.Pvalue) {
      Cf <- Cf[, -p.ind, drop = FALSE]
    }
    signif.stars <- FALSE
  }

  print.default(Cf, quote = quote, right = right, na.print = na.print, ...)

  if (signif.stars) {
    chk::chk_flag(signif.legend)

    if (signif.legend) {
      if ((w <- getOption("width")) < nchar(sleg <- attr(Signif, "legend"))) {
        sleg <- strwrap(sleg, width = w - 2, prefix = space(2))
      }

      cat0("---\nSignif. codes:  ", sleg, fill = w + 4 + max(nchar(sleg, "bytes") - nchar(sleg)))
    }
  }

  invisible(x)
}

# Constructs a model call, used in .compute_vcov() to compute variance
.build_model_call <- function(object = NULL, model = "glm", model_call,
                              weightit = NULL, vcov = NULL, br = FALSE) {

  if (is_not_null(object)) {
    model <- {
      if (inherits(object, "coxph_weightit")) "coxph"
      else if (inherits(object, "multinom_weightit")) "multinom"
      else if (inherits(object, "ordinal_weightit")) "ordinal"
      else if (inherits(object, "glm_weightit")) "glm"
      else .err("can't build model call. This is probably a bug")
    }

    if (missing(model_call) || is_null(model_call)) {
      model_call <- object[["call"]]
    }

    if (is_null(weightit)) {
      weightit <- object[["weightit"]]
    }

    if (is_null(vcov)) {
      vcov <- object[["vcov_type"]]
    }

    br <- object[["br"]]
  }

  if (is_not_null(weightit)) {
    if (!inherits(weightit, c("weightit", "weightitMSM"))) {
      .err("`weightit` must be a `weightit` or `weightitMSM` object")
    }

    if (is_null(weightit[["s.weights"]])) {
      weightit[["s.weights"]] <- rep(1, nobs(weightit))
    }
  }

  if (model == "glm") {
    chk::chk_flag(br)

    model_call[[1]] <- quote(stats::glm)

    if (is_not_null(weightit)) {
      model_call$weights <- weightit[["weights"]] * weightit[["s.weights"]]
    }
    model_call$x <- TRUE
    model_call$y <- TRUE
    model_call$model <- TRUE

    if (br) {
      rlang::check_installed("brglm2")
      model_call$method <- quote(brglm2::brglmFit)
      ctrl <- brglm2::brglmControl
    }
    else {
      model_call$method <- quote(stats::glm.fit)
      ctrl <- stats::glm.control
    }

    model_call[setdiff(names(model_call), c(names(formals(stats::glm)), names(formals(ctrl))))] <- NULL
  }
  else if (model == "ordinal") {
    model_call[[1]] <- .ordinal_weightit
    model_call[setdiff(names(model_call), names(formals(.ordinal_weightit)))] <- NULL

    if (is_not_null(weightit)) {
      model_call$weights <- weightit[["weights"]] * weightit[["s.weights"]]
    }
    model_call$x <- TRUE
    model_call$y <- TRUE
    model_call$model <- TRUE
    model_call$hess <- vcov %nin% c("none", "BS", "FWB")
  }
  else if (model == "multinom") {
    model_call[[1]] <- .multinom_weightit
    model_call[setdiff(names(model_call), names(formals(.multinom_weightit)))] <- NULL

    if (is_not_null(weightit)) {
      model_call$weights <- weightit[["weights"]] * weightit[["s.weights"]]
    }
    model_call$x <- TRUE
    model_call$y <- TRUE
    model_call$model <- TRUE
    model_call$hess <- vcov %nin% c("none", "BS", "FWB")
  }
  else if (model == "coxph") {
    model_call[[1]] <- quote(survival::coxph)

    if (is_not_null(weightit)) {
      model_call$weights <- weightit[["weights"]] * weightit[["s.weights"]]
    }
    model_call$x <- TRUE
    model_call$y <- TRUE
    model_call$model <- TRUE
    model_call$robust <- vcov == "HC0"

    model_call$cluster <- NULL
    model_call[setdiff(names(model_call), c(names(formals(survival::coxph)), names(formals(survival::coxph.control))))] <- NULL
  }

  model_call
}

# Computes variance from model fit and call; used in glm_weightit()
.compute_vcov <- function(fit, weightit = NULL, vcov, cluster = NULL, glm_call) {
  if (is_not_null(cluster) && vcov %in% c("none", "const")) {
    .wrn("`cluster` is not used when `vcov = %s`", add_quotes(vcov))
  }

  if (vcov == "none") {
    return(.modify_vcov_info(matrix(nrow = 0L, ncol = 0L), vcov_type = "none", cluster = cluster))
  }

  bout <- fit[["coefficients"]]
  aliased <- is.na(bout)

  if (vcov == "const") {
    if (inherits(fit, "ordinal_weightit")) {
      if (is_null(fit[["hessian"]])) {
        fit[["hessian"]] <- .get_hess_ordinal(fit)
      }

      V <- solve(-fit[["hessian"]])
    }
    else if (inherits(fit, "multinom_weightit")) {
      if (is_null(fit[["hessian"]])) {
        fit[["hessian"]] <- .get_hess_multinom(fit)
      }

      V <- solve(-fit[["hessian"]])
    }
    else if (inherits(fit, "coxph_weightit")) {
      V <- if_null_then(fit[["naive.var"]], fit[["var"]])
    }
    else {
      .declass <- function(obj) {
        class(obj) <- setdiff(class(obj), "glm_weightit")
        obj
      }

      V <- stats::vcov(.declass(fit))
    }

    colnames(V) <- rownames(V) <- names(aliased)[!aliased]

    return(.modify_vcov_info(V, vcov_type = "const", cluster = cluster))
  }

  Xout <- if_null_then(fit[["x"]], model.matrix(fit))
  Y <- if_null_then(fit[["y"]], model.response(model.frame(fit)))

  if (is_not_null(weightit)) {
    W <- weightit[["weights"]]
    SW <- weightit[["s.weights"]]
  }
  else {
    W <- SW <- NULL
  }

  if (is_null(SW)) {
    SW <- rep.int(1, length(Y))
  }

  offset <- if_null_then(fit[["offset"]], rep.int(0, length(Y)))

  if (any(aliased)) {
    if (!is_null(attr(fit[["qr"]][["qr"]], "aliased"))) {
      Xout <- Xout[, !attr(fit[["qr"]][["qr"]], "aliased"), drop = FALSE]
    }
    else {
      Xout <- make_full_rank(Xout, with.intercept = FALSE)
    }
    bout <- bout[!aliased]
  }

  pout <- sum(!aliased)

  if (is_not_null(cluster) && vcov %in% c("asympt", "HC0", "BS", "FWB")) {
    if (inherits(cluster, "formula")) {
      cluster_tmp <- expand.model.frame(fit, cluster, na.expand = FALSE)
      cluster <- model.frame(cluster, cluster_tmp, na.action = na.pass)
    }
    else {
      cluster <- as.data.frame(cluster)
    }

    if (nrow(cluster) != length(Y)) {
      .err("the number of observations in `cluster` must equal that in `data`")
    }

    chk::chk_not_any_na(cluster)

    p <- ncol(cluster)
    if (p > 1L) {
      clu <- lapply(seq_len(p), function(i) utils::combn(seq_len(p), i, simplify = FALSE))
      clu <- unlist(clu, recursive = FALSE)
      sign <- vapply(clu, function(i) (-1)^(length(i) + 1), numeric(1L))
      paste_ <- function(...) paste(..., sep = "_")
      for (i in (p + 1L):length(clu)) {
        cluster <- cbind(cluster, Reduce(paste_, unclass(cluster[, clu[[i]]])))
      }
    }
    else {
      clu <- list(1)
      sign <- 1
    }

    #Small sample adjustment (setting cadjust = TRUE in vcovCL)
    g <- vapply(seq_along(clu), function(i) {
      if (is.factor(cluster[[i]])) {
        length(levels(cluster[[i]]))
      }
      else {
        length(unique(cluster[[i]]))
      }
    }, numeric(1L))
  }

  if (vcov == "FWB") {
    R <- attr(vcov, "R")
    fwb.args <- attr(vcov, "fwb.args")

    glm_call$x <- FALSE
    glm_call$y <- FALSE
    glm_call$model <- FALSE

    if (is_not_null(weightit)) {
      wcall <- weightit[["call"]]
      wenv <- weightit[["env"]]
    }
    else {
      weightit_boot <- list(weights = 1)
    }

    genv <- environment(fit[["formula"]])

    fwbfun <- function(data, w) {
      if (is_not_null(weightit)) {
        wcall$s.weights <- SW * w

        withCallingHandlers({
          weightit_boot <- eval(wcall, wenv)
        },
        warning = function(w) {
          w <- conditionMessage(w)
          if (!startsWith(tolower(w), "some extreme weights were generated"))
            .wrn("(from `weightit()`) ", w, tidy = FALSE)
          invokeRestart("muffleWarning")
        })

        glm_call$weights <- weightit_boot[["weights"]] * SW * w
      }
      else {
        glm_call$weights <- SW * w
      }

      withCallingHandlers({
        fit_boot <- eval(glm_call, genv)
      },
      warning = function(w) {
        w <- conditionMessage(w)
        if (w != "non-integer #successes in a binomial glm!") .wrn("(from `glm()`) ", w, tidy = FALSE)
        invokeRestart("muffleWarning")
      })

      fit_boot[["coefficients"]]
    }

    #Sham data.frame to get weight length
    fwb.args$data <- data.frame(SW)
    fwb.args$statistic <- fwbfun
    fwb.args$R <- R
    if (is_null(fwb.args$verbose)) fwb.args$verbose <- FALSE
    fwb.args <- fwb.args[names(fwb.args) %in% names(formals(fwb::fwb))]

    if (is_null(cluster)) {
      fwb_out <- eval(as.call(c(list(quote(fwb::fwb)), fwb.args)))
      V <- cov(fwb_out$t)
    }
    else {
      V <- 0
      for (i in seq_along(clu)) {
        fwb.args$cluster <- cluster[[i]]
        fwb_out <- eval(as.call(c(list(quote(fwb::fwb)), fwb.args)))
        V <- V + sign[i] * cov(fwb_out[["t"]])
      }
    }
  }
  else if (vcov == "BS") {
    R <- attr(vcov, "R")

    glm_call$x <- FALSE
    glm_call$y <- FALSE
    glm_call$model <- FALSE

    genv <- environment(fit$formula)
    if (is_not_null(weightit)) {
      wcall <- weightit$call
      wenv <- weightit$env
      data <- eval(wcall$data, wenv)
      if (is_null(data)) {
        .err(sprintf('a dataset must have been supplied to `data` in the original call to `%s()` to use `vcov = "BS"`',
                     deparse1(wcall[[1]])))
      }
    }
    else {
      weightit_boot <- list(weights = 1)
      data <- eval(glm_call$data, genv)
      if (is_null(data)) {
        .err(sprintf('a dataset must have been supplied to `data` in the original call to `%s()` to use `vcov = "BS"`',
                     deparse1(glm_call[[1]])))
      }
    }

    bootfun <- function(data, ind) {
      if (is_not_null(weightit)) {
        wcall$data <- data[ind,]
        wcall$s.weights <- SW[ind]

        withCallingHandlers({
          weightit_boot <- eval(wcall, wenv)
        },
        warning = function(w) {
          w <- conditionMessage(w)
          if (!startsWith(tolower(w), "some extreme weights were generated"))
            .wrn("(from `weightit()`) ", w, tidy = FALSE)
          invokeRestart("muffleWarning")
        })
        glm_call$weights <- weightit_boot$weights * SW[ind]
      }
      else {
        glm_call$weights <- SW[ind]
      }

      glm_call$data <- data[ind,]

      withCallingHandlers({
        fit_boot <- eval(glm_call, genv)
      },
      warning = function(w) {
        w <- conditionMessage(w)
        if (w != "non-integer #successes in a binomial glm!") .wrn("(from `glm()`) ", w, tidy = FALSE)
        invokeRestart("muffleWarning")
      })

      fit_boot$coefficients
    }

    if (is_null(cluster)) {
      boot_out <- do.call("rbind", lapply(seq_len(R), function(i) {
        ind <- sample(seq_row(data), nrow(data), replace = TRUE)
        bootfun(data, ind)
      }))

      V <- cov(boot_out)
    }
    else {
      V <- 0
      for (i in seq_along(clu)) {
        cli <- split(seq_along(cluster[[i]]), cluster[[i]])

        boot_out <- do.call("rbind", lapply(seq_len(R), function(i) {
          ind <- unlist(cli[sample(seq_along(cli), length(cli), replace = TRUE)])
          bootfun(data, ind)
        }))

        V <- V + sign[i] * cov(boot_out)
      }
    }
  }
  else if (inherits(fit, "coxph")) {
    if (is_null(cluster)) {
      V <- fit$var
    }
    else {
      V <- 0

      for (i in seq_along(clu)) {
        adj <- g[i]/(g[i] - 1)
        temp <- residuals(fit, type = "dfbeta",
                          collapse = cluster[[i]],
                          weighted = TRUE)

        V <- V + sign[i] * adj * crossprod(temp)
      }
    }
  }
  else {
    hess <- NULL
    if (vcov == "asympt") {
      Mparts <- attr(weightit, "Mparts", exact = TRUE)
      Mparts.list <- attr(weightit, "Mparts.list", exact = TRUE)
    }

    psi_out <- function(Bout, w, Y, Xout, SW, offset) {
      fit$psi(Bout, Xout, Y, w * SW, offset = offset)
    }

    if (vcov == "HC0") {
      wfun <- {
        if (is_null(W)) function(Btreat, Xtreat, A) 1
        else function(Btreat, Xtreat, A) W
      }

      Xtreat <- NULL
      A <- NULL

      psi <- function(B, Xout, Y, Xtreat, A, SW, offset) {
        Bout <- B[seq_len(pout)]
        Btreat <- NULL

        w <- wfun(Btreat, Xtreat, A)
        psi_out(B, w, Y, Xout, SW, offset)
      }

      b <- bout

      psi_b <- if_null_then(fit[["gradient"]], psi(b, Xout, Y, Xtreat, A, SW, offset))

      if (is_not_null(fit[["hessian"]])) {
        hess <- fit$hessian
      }
      else if (inherits(fit, "ordinal_weightit")) {
        hess <- .get_hess_ordinal(fit)
      }
      else if (inherits(fit, "multinom_weightit")) {
        hess <- .get_hess_multinom(fit)
      }
    }
    else if (is_not_null(Mparts)) {
      # Mparts from weightit()
      psi_treat <- Mparts[["psi_treat"]]
      wfun <- Mparts[["wfun"]]
      Xtreat <- Mparts[["Xtreat"]]
      A <- Mparts[["A"]]
      btreat <- Mparts[["btreat"]]

      psi <- function(B, Xout, Y, Xtreat, A, SW, offset) {
        Bout <- B[seq_len(pout)]
        Btreat <- B[-seq_len(pout)]

        w <- wfun(Btreat, Xtreat, A)

        cbind(psi_out(Bout, w, Y, Xout, SW, offset),
              psi_treat(Btreat, A, Xtreat, SW))
      }

      b <- c(bout, btreat)

      psi_b <- psi(b, Xout, Y, Xtreat, A, SW, offset)
    }
    else {
      # Mparts.list from weightitMSM() or weightit()
      psi_treat.list <- grab(Mparts.list, "psi_treat")
      wfun.list <- grab(Mparts.list, "wfun")
      Xtreat.list <- grab(Mparts.list, "Xtreat")
      A.list <- grab(Mparts.list, "A")
      btreat.list <- grab(Mparts.list, "btreat")

      psi_treat <- function(Btreat, A, Xtreat, SW) {
        do.call("cbind", lapply(seq_along(Btreat), function(i) {
          psi_treat.list[[i]](Btreat[[i]], A[[i]], Xtreat[[i]], SW)
        }))
      }

      wfun <- function(Btreat, Xtreat, A) {
        Reduce("*", lapply(seq_along(Btreat), function(i) {
          wfun.list[[i]](Btreat[[i]], Xtreat[[i]], A[[i]])
        }), init = 1)
      }

      psi <- function(B, Xout, Y, Xtreat, A, SW, offset) {
        Bout <- B[seq_len(pout)]
        Btreat <- btreat.list
        k <- 0
        for (i in seq_along(btreat.list)) {
          Btreat[[i]] <- B[-seq_len(pout)][k + seq_along(btreat.list[[i]])]
          k <- k + length(btreat.list[[i]])
        }

        w <- wfun(Btreat, Xtreat, A)

        cbind(psi_out(Bout, w, Y, Xout, SW, offset),
              psi_treat(Btreat, A, Xtreat, SW))
      }

      b <- c(bout, unlist(btreat.list))

      psi_b <- psi(b, Xout, Y, Xtreat.list, A.list, SW, offset)
    }

    if (is_not_null(cluster)) {
      B <- 0

      for (i in seq_along(clu)) {
        adj <- g[i]/(g[i] - 1)

        B <- B + sign[i] * adj * crossprod(rowsum(psi_b, cluster[[i]], reorder = FALSE))
      }
    }
    else {
      B <- crossprod(psi_b)
    }

    # Gradient of gradfun -> hessian
    gradfun <- function(B, Xout, Y, Xtreat, A, SW, offset) {
      colSums(psi(B, Xout, Y, Xtreat, A, SW, offset))
    }

    if (is_null(hess)) {
      hess <- .gradient(gradfun,
                        .x = b,
                        Xout = Xout,
                        Y = Y,
                        Xtreat = {
                          if (is_not_null(Mparts.list)) Xtreat.list
                          else Xtreat
                        },
                        A = {
                          if (is_not_null(Mparts.list)) A.list
                          else A
                        },
                        SW = SW,
                        offset = offset)
    }

    A1 <- try(solve(hess, t(.chol2(B))), silent = TRUE)

    if (null_or_error(A1)) {
      .e <- conditionMessage(attr(A1, "condition"))
      if (startsWith(.e, "system is computationally singular") ||
          startsWith(.e, "Lapack routine dgesv: system is exactly singular")) {
        .err('the Hessian could not be inverted, which indicates an estimation failure, likely due to perfect separation. Estimates from this model should not be trusted. Investigate the problem by refitting with `vcov = "none"`. Simplifying the model can sometimes help')
      }

      .err(.e, tidy = FALSE)
    }

    #Subset to just the outcome model coefs
    A1 <- A1[seq_len(pout), , drop = FALSE]

    V <- tcrossprod(A1)
  }

  colnames(V) <- rownames(V) <- names(aliased)[!aliased]

  .modify_vcov_info(V, vcov_type = vcov, cluster = cluster)
}

# Processes fit for output; used in glm_weightit()
.process_fit <- function(fit, weightit = NULL, vcov, glm_weightit_call, x, y) {
  if (is_not_null(weightit) && is_not_null(fit[["model"]])) {
    fit$model[["(s.weights)"]] <- weightit[["s.weights"]]
    fit$model[["(weights)"]] <- weightit[["weights"]] * weightit[["s.weights"]]
  }

  fit$vcov_type <- attr(fit[["vcov"]], "vcov_type")

  fit$call <- glm_weightit_call

  if (is_not_null(weightit)) {
    fit$weightit <- weightit
  }

  if (isFALSE(x)) fit$x <- NULL
  if (isFALSE(y)) fit$y <- NULL

  attr(fit, "cluster") <- attr(fit[["vcov"]], "cluster")

  fit[["vcov"]] <- .modify_vcov_info(fit[["vcov"]])

  fit
}
#' Fitting Weighted Generalized Linear Models
#'
#' @description
#' `lm_weightit()` and `glm_weightit()` are used to fit (generalized) linear models with a variance matrix that accounts for estimation of weights, if supplied. By default, these functions use M-estimation to construct a robust covariance matrix using the estimation equations for the weighting model and the outcome model. `lm_weightit()` is a wrapper for `glm_weightit()` with the Gaussian family and identity link (i.e., a linear model). `coxph_weightit()` fits a Cox proportional hazards model accounting for the weights and is a wrapper for [survival::coxph()].
#'
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted. For `coxph_weightit()`, see [survival::coxph()] for how this should be specified.
#' @param data a data frame containing the variables in the model. If not found in data, the variables are taken from `environment(formula)`, typically the environment from which the function is called.
#' @param family a description of the error distribution and link function to be used in the model. This can be a character string naming a family function, a family function or the result of a call to a family function. See [family] for details of family functions. Can also be the string `"multinomial"` for multinomial logistic regression.
#' @param weightit a `weightit` or `weightitMSM` object; the output of a call to [weightit()] or [weightitMSM()]. If not supplied, an unweighted model will be fit.
#' @param vcov string; the method used to compute the variance of the estimated parameters. Allowable options include `"asympt"`, which uses the asymptotically correct M-estimation-based method that accounts for estimation of the weights when available; `"const"`, which uses the usual maximum likelihood estimates (only available when `weightit` is not supplied); `"HC0"`, which computes the robust sandwich variance treating weights (if supplied) as fixed; `"BS"`, which uses the traditional bootstrap (including re-estimation of the weights, if supplied); `"FWB"`, which uses the fractional weighted bootstrap as implemented in \pkgfun{fwb}{fwb} (including re-estimation of the weights, if supplied); and `"none"` to omit calculation of a variance matrix. If `NULL` (the default), will use `"asympt"` if `weightit` is supplied and M-estimation is available and `"HC0"` otherwise. See the `vcov_type` component of the outcome object to see which was used.
#' @param cluster optional; for computing a cluster-robust variance matrix, a variable indicating the clustering of observations, a list (or data frame) thereof, or a one-sided formula specifying which variable(s) from the fitted model should be used. Note the cluster-robust variance matrix uses a correction for small samples, as is done in `sandwich::vcovCL()` by default. Cluster-robust variance calculations are available only when `vcov` is `"asympt"`, `"HC0"`, `"BS"`, or `"FWB"`.
#' @param R the number of bootstrap replications when `vcov` is `"BS"` or `"FWB"`. Default is 500. Ignored otherwise.
#' @param offset optional; a numeric vector contain the model offset. See [offset()]. An offset can also be preset in the model formula.
#' @param start optional starting values for the coefficients.
#' @param etastart,mustart optional starting values for the linear predictor and vector of means when `family` is not `"multinomial"`. Passed to [glm()].
#' @param control a list of parameters for controlling the fitting process.
#' @param x,y logical values indicating whether the response vector and model matrix used in the fitting process should be returned as components of the returned value.
#' @param contrasts an optional list define contrasts for factor variables. See [model.matrix()].
#' @param fwb.args an optional list of further arguments to supply to \pkgfun{fwb}{fwb} when `vcov = "FWB"`.
#' @param object a `glm_weightit` object.
#' @param ci `logical` whether to display Wald confidence intervals for estimated coefficients. Default is `FALSE`.
#' @param level when `ci = TRUE`, the desired confidence level.
#' @param transform the function used to transform the coefficients, e.g., `exp` (which can also be supplied as a string, e.g., `"exp"`); passed to [match.fun()] before being used on the coefficients. When `ci = TRUE`, this is also applied to the confidence interval bounds. If specified, the standard error will be omitted from the output. Default is no transformation.
#' @param \dots for `glm_weightit()` and `lm_weightit()`, arguments to be used to form the default control argument if it is not supplied directly. Ignored otherwise.
#'
#' @returns
#' For `lm_weightit()` and `glm_weightit()`, a `glm_weightit` object, which inherits from `glm`. Unless `vcov = "none"`, the `vcov` component contains the covariance matrix adjusted for the estimation of the weights if requested and a compatible `weightit` object was supplied. The `vcov_type` component contains the type of variance matrix requested. If `cluster` is supplied, it will be stored in the `"cluster"` attribute of the output object, even if not used. For `coxph_weightit()`, a `coxph_weightit` object, which inherits from `glm_weightit` and `coxph`. See [survival::coxph()] for details.
#'
#' `print()`, `vcov()`, `predict()`, and `confint()` methods are also available; these generally follow the same pattern as the respect method for `glm` objects. `confint()` uses Wald confidence intervals (internally calling [confint.lm()]). When `family = "multinomial"`, predict() produces a matrix of predicted probabilities, one for each level of the outcome, and the `type` argument is ignored. `model.frame()` output (also the `model` component of the output object) will include two extra column when `weightit` is supplied: `(weights)` containing the weights used in the model (the product of the estimated weights and the sampling weights, if any) and `(s.weights)` containing the sampling weights, which will be 1 if `s.weights` is not supplied in the original `weightit()` call.
#'
#' @details
#' [glm_weightit()] is essentially a wrapper for [glm()] that optionally computes a coefficient variance matrix that can be adjusted to account for estimation of the weights if a `weightit` or `weightitMSM` object is supplied to the `weightit` argument. When no argument is supplied to `weightit` or there is no `"Mparts"` attribute in the supplied object, the default variance matrix returned will be the "HC0" sandwich variance matrix, which is robust to misspecification of the outcome family (including heteroscedasticity). Otherwise, the default variance matrix uses M-estimation to additionally adjust for estimation of the weights. When possible, this often yields smaller (and more accurate) standard errors. See the individual methods pages to see whether and when an `"Mparts"` attribute is included in the supplied object. To request that a variance matrix be computed that doesn't account for estimation of the weights even when a compatible `weightit` object is supplied, set `vcov = "HC0"`, which treats the weights as fixed.
#'
#' Bootstrapping can also be used to compute the coefficient variance matrix; when `vcov = "BS"` or `vcov = "FWB"`, which implement the traditional resampling-based and fractional weighted bootstrap, respectively, the entire process of estimating the weights and fitting the outcome model is repeated in bootstrap samples (if a `weightit` object is supplied). This accounts for estimation of the weights and can be used with any weighting method. It is important to set a seed using `set.seed()` to ensure replicability of the results. The fractional weighted bootstrap is more reliable but requires the weighting method to accept sampling weights (which most do, and you'll get an error if it doesn't). Setting `vcov = "FWB"` and supplying `fwb.args = list(wtype = "multinom")` also performs the resampling-based bootstrap but with the additional features \pkg{fwb} provides (e.g., a progress and parallelization) at the expense of needing to have \pkg{fwb} installed.
#'
#' When `family = "multinomial"`, multinomial logistic regression is fit using a custom function in \pkg{WeightIt} that uses M-estimation to estimate the model coefficients. This implementation is less robust to failures than other multinomial logistic regression solvers and should be used with caution. Estimation of coefficients should align with that from `mlogit::mlogit()` and `mclogit::mblogit()`.
#'
#' Functions in the \pkg{sandwich} package can be to compute standard errors after fitting, regardless of how `vcov` was specified, though these will ignore estimation of the weights, if any. When no adjustment is done for estimation of the weights (i.e., because no `weightit` argument was supplied or there was no `"Mparts"` component in the supplied object), the default variance matrix produced by `glm_weightit()` should align with that from `sandwich::vcovHC(. type = "HC0")` or `sandwich::vcovCL(., type = "HC0", cluster = cluster)` when `cluster` is supplied.
#'
#' `coxph_weightit()` is a wrapper for [survival::coxph()] to fit weighted survival models, optionally accounting for estimation of the weights. It differs from `coxph()` in a few ways: the `print()` and `summary()` methods are more like those for `glm` objects then for `coxph` objects, and the `cluster` argument should be specified as a one-sided formula (which can include multiple clustering variables) and uses a small sample correction for cluster variance estimates when specified. Currently, M-estimation is not supported, so bootstrapping (i.e., `vcov = "BS"` or `"FWB"`) is the only way to correctly adjust for estimation of the weights.
#'
#' @seealso
#' [lm()] and [glm()] for fitting generalized linear models without adjusting standard errors for estimation of the weights. [survival::coxph()] for fitting Cox proportional hazards models without adjust standard errors for estimation of the weights.
#'
#' @examples
#'
#' data("lalonde", package = "cobalt")
#'
#' # Logistic regression ATT weights
#' w.out <- weightit(treat ~ age + educ + married + re74,
#'                   data = lalonde, method = "glm",
#'                   estimand = "ATT")
#'
#' # Linear regression outcome model that adjusts
#' # for estimation of weights
#' fit1 <- lm_weightit(re78 ~ treat, data = lalonde,
#'                     weightit = w.out)
#'
#' summary(fit1)
#'
#' # Linear regression outcome model that treats weights
#' # as fixed
#' fit2 <- lm_weightit(re78 ~ treat, data = lalonde,
#'                     weightit = w.out, vcov = "HC0")
#'
#' summary(fit2)
#'
#' # Linear regression outcome model that bootstraps
#' # estimation of weights and outcome model fitting
#' # using fractional weighted bootstrap with "Mammen"
#' # weights
#' set.seed(123)
#' fit3 <- lm_weightit(re78 ~ treat, data = lalonde,
#'                     weightit = w.out,
#'                     vcov = "FWB",
#'                     R = 50,
#'                     fwb.args = list(wtype = "mammen"))
#'
#' summary(fit3)
#'
#' # Multinomial logistic regression outcome model
#' # that adjusts for estimation of weights
#' lalonde$re78_3 <- factor(findInterval(lalonde$re78,
#'                                       c(0, 5e3, 1e4)))
#'
#' fit4 <- glm_weightit(re78_3 ~ treat, data = lalonde,
#'                      weightit = w.out,
#'                      family = "multinomial")
#'
#' summary(fit4)

#' @export
#' @name glm_weightit
glm_weightit <- function(formula, data, family = gaussian, weightit,
                         vcov = NULL, cluster, R = 500,
                         offset, start = NULL, etastart, mustart,
                         control = list(...),
                         x = FALSE, y = TRUE,
                         contrasts = NULL, fwb.args = list(), ...) {

  args <- list(...)

  if (missing(weightit) || is.null(weightit)) {
    weightit <- list()

    allowable_vcov <- c("none", "const", "HC0", "BS", "FWB")

    if (is_null(vcov)) {
      vcov <- "HC0"
    }
  }
  else {
    chk::chk_is(weightit, "weightit")
    if (is_null(attr(weightit, "Mparts", exact = TRUE)) &&
        is_null(attr(weightit, "Mparts.list", exact = TRUE))) {
      allowable_vcov <- c("none", "const", "HC0", "BS", "FWB")

      if (is_null(vcov)) {
        vcov <- "HC0"
      }
    }
    else {
      allowable_vcov <- c("none", "const", "asympt", "HC0", "BS", "FWB")

      if (is_null(vcov)) {
        vcov <- "asympt"
      }
    }
  }

  chk::chk_string(vcov)
  vcov <- match_arg(vcov, allowable_vcov)

  if (inherits(weightit, "weightit") && vcov == "const") {
    .wrn("`vcov = \"const\"` should not be used when `weightit` is supplied; the resulting standard errors are invalid and should not be interpreted")
  }

  bootstrap <- vcov %in% c("BS", "FWB")
  if (bootstrap) {
    chk::chk_count(R)
    chk::chk_gt(R, 0)
  }

  if (missing(cluster)) {
    cluster <- NULL
  }

  if (is_not_null(cluster) && vcov %in% c("none", "const")) {
    .wrn("`cluster` is not used when `vcov = %s`", add_quotes(vcov))
  }

  glm_call <- glm_weightit_call <- match.call()

  if (is.character(family) && identical(family, "multinomial")) {
    glm_call[[1]] <- .mlogit_weightit

    if (is_not_null(weightit))
      glm_call$weights <- weightit$weights * weightit$s.weights
    glm_call$x <- TRUE
    glm_call$y <- TRUE
    glm_call$model <- TRUE
    glm_call[setdiff(names(glm_call), names(formals(.mlogit_weightit)))] <- NULL

    if (!bootstrap && vcov != "none") {
      glm_call$hess <- TRUE
    }

    fit <- eval.parent(glm_call)

    fit$family <- list(family = "multinomial",
                       link = "logit")

    psi_out <- function(Bout, w, Y, Xout, SW) {
      fit$psi(Bout, Xout, Y, w * SW)
    }
  }
  else {
    glm_call[[1]] <- quote(stats::glm)

    if (is_not_null(weightit))
      glm_call$weights <- weightit$weights * weightit$s.weights
    glm_call$x <- TRUE
    glm_call$y <- TRUE
    glm_call$model <- TRUE
    glm_call$family <- family

    glm_call[setdiff(names(glm_call), c(names(formals(stats::glm)), names(formals(stats::glm.control))))] <- NULL

    withCallingHandlers({
      fit <- eval.parent(glm_call)
    },
    warning = function(w) {
      w <- conditionMessage(w)
      if (w != "non-integer #successes in a binomial glm!") .wrn("(from `glm()`) ", w, tidy = FALSE)
      invokeRestart("muffleWarning")
    })

    fam <- fit$family

    psi_out <- function(Bout, w, Y, Xout, SW) {
      lin_pred <- drop(Xout %*% Bout)
      Xout * (SW * w * fam$mu.eta(lin_pred) * (Y - fam$linkinv(lin_pred)) / fam$variance(fam$linkinv(lin_pred)))
    }
  }

  if (is_not_null(weightit)) {
    fit$model[["(s.weights)"]] <- weightit$s.weights
    fit$model[["(weights)"]] <- weightit$weights * weightit$s.weights
  }

  Xout <- fit$x
  Y <- fit$y
  W <- weightit$weights
  SW <- weightit$s.weights
  if (is_null(SW)) SW <- rep(1, length(Y))
  bout <- fit$coefficients
  aliased <- is.na(bout)

  if (any(aliased)) {
    if (!is_null(attr(fit$qr$qr, "aliased"))) {
      Xout <- Xout[, !attr(fit$qr$qr, "aliased"), drop = FALSE]
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

  if (vcov == "none") {
    V <- NULL
  }
  else if (vcov == "const") {
    V <- vcov(fit)
  }
  else if (vcov == "FWB") {
    rlang::check_installed("fwb")

    chk::chk_list(fwb.args)

    if (is_not_null(weightit)) {
      wcall <- weightit$call
      # wenv <- environment(weightit$formula)
      wenv <- weightit$env
    }
    else {
      weightit_boot <- list(weights = 1)
    }
    genv <- environment(fit$formula)

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

        glm_call$weights <- weightit_boot$weights * SW * w
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

      fit_boot$coefficients
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
        V <- V + sign[i] * cov(fwb_out$t)
      }
    }
  }
  else if (vcov == "BS") {
    genv <- environment(fit$formula)
    if (is_not_null(weightit)) {
      wcall <- weightit$call
      wenv <- environment(weightit$formula)
      data <- eval(wcall$data, wenv)
      if (is_null(data)) {
        .err(sprintf("a dataset must have been supplied to `data` in the original call to `%s()` to use `vcov = \"BS\"`",
                     deparse1(wcall[[1]])))
      }
    }
    else {
      weightit_boot <- list(weights = 1)
      data <- eval(glm_call$data, genv)
      if (is_null(data)) {
        .err(sprintf("a dataset must have been supplied to `data` in the original call to `%s()` to use `vcov = \"BS\"`",
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
        ind <- sample(seq_len(nrow(data)), nrow(data), replace = TRUE)
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
  else {
    hess <- NULL
    if (vcov == "asympt") {
      Mparts <- attr(weightit, "Mparts", exact = TRUE)
      Mparts.list <- attr(weightit, "Mparts.list", exact = TRUE)
    }

    if (vcov == "HC0") {
      wfun <- {
        if (is.null(W)) function(Btreat, Xtreat, A) 1
        else function(Btreat, Xtreat, A) W
      }

      Xtreat <- NULL
      A <- NULL

      psi <- function(B, Xout, Y, Xtreat, A, SW) {
        Bout <- B[seq_len(pout)]
        Btreat <- B[-seq_len(pout)]

        w <- wfun(Btreat, Xtreat, A)
        psi_out(B, w, Y, Xout, SW)
      }

      b <- bout

      if (is_not_null(fit$gradient)) {
        psi_b <- fit$gradient
      }
      else {
        psi_b <- psi(b, Xout, Y, Xtreat, A, SW)
      }

      if (is_not_null(fit$hessian)) {
        hess <- fit$hessian
      }
    }
    else if (is_not_null(Mparts)) {
      # Mparts from weightit()
      psi_treat <- Mparts$psi_treat
      wfun <- Mparts$wfun
      Xtreat <- Mparts$Xtreat
      A <- Mparts$A
      btreat <- Mparts$btreat

      psi <- function(B, Xout, Y, Xtreat, A, SW) {
        Bout <- B[seq_len(pout)]
        Btreat <- B[-seq_len(pout)]

        w <- wfun(Btreat, Xtreat, A)

        cbind(psi_out(Bout, w, Y, Xout, SW),
              psi_treat(Btreat, A, Xtreat, SW))
      }

      b <- c(bout, btreat)

      psi_b <- psi(b, Xout, Y, Xtreat, A, SW)
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

      psi <- function(B, Xout, Y, Xtreat, A, SW) {
        Bout <- B[seq_len(pout)]
        Btreat <- btreat.list
        k <- 0
        for (i in seq_along(btreat.list)) {
          Btreat[[i]] <- B[-seq_len(pout)][k + seq_along(btreat.list[[i]])]
          k <- k + length(btreat.list[[i]])
        }

        w <- wfun(Btreat, Xtreat, A)

        cbind(psi_out(Bout, w, Y, Xout, SW),
              psi_treat(Btreat, A, Xtreat, SW))
      }

      b <- c(bout, unlist(btreat.list))

      psi_b <- psi(b, Xout, Y, Xtreat.list, A.list, SW)
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
    gradfun <- function(B, Xout, Y, Xtreat, A, SW) {
      colSums(psi(B, Xout, Y, Xtreat, A, SW))
    }

    if (is_null(hess)) {
      hess <- gradient(gradfun,
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
                       SW = SW)
    }

    A1 <- try(solve(-hess), silent = TRUE)

    if (null_or_error(A1)) {
      .e <- conditionMessage(attr(A1, "condition"))
      if (startsWith(.e, "system is computationally singular")) {
        .err("the negative Hessian could not be inverted, which indicates an estimation failure, likely due to perfect separation. Estimates from this model should not be trusted. Investigate the problem by refitting with `vcov = \"none\"`. Simplifying the model can sometimes help")
      }

      .err(.e, tidy = FALSE)
    }

    V <- A1 %*% tcrossprod(B, A1)
  }

  if (is_not_null(V)) {
    fit$vcov <- V[seq_len(pout), seq_len(pout), drop = FALSE]
    colnames(fit$vcov) <- rownames(fit$vcov) <- names(aliased)[!aliased]
  }

  fit$vcov_type <- vcov
  fit$call <- glm_weightit_call

  if (isFALSE(x)) fit$x <- NULL
  if (isFALSE(y)) fit$y <- NULL

  attr(fit, "cluster") <- cluster
  class(fit) <- c("glm_weightit", class(fit))

  fit
}

#' @export
#' @rdname glm_weightit
coxph_weightit <- function(formula, data, weightit,
                           vcov = NULL, cluster, R = 500,
                           x = FALSE, y = TRUE,
                           fwb.args = list(), ...) {

  rlang::check_installed("survival")

  if (missing(weightit) || is.null(weightit)) {
    weightit <- list()

    allowable_vcov <- c("none", "const", "HC0", "BS", "FWB")

    if (is_null(vcov)) {
      vcov <- "HC0"
    }
  }
  else {
    chk::chk_is(weightit, "weightit")

    allowable_vcov <- c("none", "const", "HC0", "BS", "FWB")

    if (is_null(vcov)) {
      vcov <- "HC0"
    }
  }

  chk::chk_string(vcov)
  vcov <- match_arg(vcov, allowable_vcov)

  if (inherits(weightit, "weightit") && vcov == "const") {
    .wrn("`vcov = \"const\"` should not be used when `weightit` is supplied; the resulting standard errors are invalid and should not be interpreted")
  }

  bootstrap <- vcov %in% c("BS", "FWB")
  if (bootstrap) {
    chk::chk_count(R)
    chk::chk_gt(R, 0)
  }

  if (missing(cluster)) {
    cluster <- NULL
  }

  if (is_not_null(cluster) && vcov %in% c("none", "const")) {
    .wrn("`cluster` is not used when `vcov = %s`", add_quotes(vcov))
  }

  coxph_call <- coxph_weightit_call <- match.call()

  coxph_call[[1]] <- quote(survival::coxph)

  if (is_not_null(weightit)) {
    coxph_call$weights <- weightit$weights * weightit$s.weights
  }
  coxph_call$x <- TRUE
  coxph_call$y <- TRUE
  coxph_call$model <- TRUE
  coxph_call$robust <- vcov == "HC0"

  coxph_call$cluster <- NULL
  coxph_call[setdiff(names(coxph_call), c(names(formals(survival::coxph)), names(formals(survival::coxph.control))))] <- NULL

  fit <- eval.parent(coxph_call)

  if (is_not_null(weightit)) {
    fit$model[["(s.weights)"]] <- weightit$s.weights
    fit$model[["(weights)"]] <- weightit$weights * weightit$s.weights
  }

  Xout <- fit$x
  Y <- fit$y
  W <- weightit$weights
  SW <- weightit$s.weights
  if (is_null(SW)) SW <- rep(1, length(Y))
  bout <- fit$coefficients
  aliased <- is.na(bout)

  if (any(aliased)) {
    if (!is_null(attr(fit$qr$qr, "aliased"))) {
      Xout <- Xout[, !attr(fit$qr$qr, "aliased"), drop = FALSE]
    }
    else {
      Xout <- make_full_rank(Xout, with.intercept = FALSE)
    }
    bout <- bout[!aliased]
  }

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

  if (vcov == "none") {
    V <- NULL
  }
  else if (vcov == "const") {
    V <- fit$naive.var
    if (is_null(V)) {
      V <- fit$var
    }
  }
  else if (vcov == "FWB") {
    rlang::check_installed("fwb")

    chk::chk_list(fwb.args)

    coxenv <- environment(fit$formula)

    coxph_call$robust <- FALSE
    coxph_call$x <- FALSE
    coxph_call$y <- FALSE
    coxph_call$model <- FALSE

    if (is_not_null(weightit)) {
      wcall <- weightit$call
      # wenv <- environment(weightit$formula)
      wenv <- weightit$env
    }
    else {
      weightit_boot <- list(weights = 1)
    }

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
        coxph_call$weights <- weightit_boot$weights * SW * w
      }
      else {
        coxph_call$weights <- SW * w
      }

      fit_boot <- eval(coxph_call, coxenv)

      fit_boot$coefficients
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
        V <- V + sign[i] * cov(fwb_out$t)
      }
    }
  }
  else if (vcov == "BS") {
    coxenv <- environment(fit$formula)

    coxph_call$robust <- FALSE
    coxph_call$x <- FALSE
    coxph_call$y <- FALSE
    coxph_call$model <- FALSE

    if (is_not_null(weightit)) {
      wcall <- weightit$call
      wenv <- environment(weightit$formula)
      data <- eval(wcall$data, wenv)

      if (is_null(data)) {
        .err(sprintf("a dataset must have been supplied to `data` in the original call to `%s()` to use `vcov = \"BS\"`",
                     deparse1(wcall[[1]])))
      }
    }
    else {
      weightit_boot <- list(weights = 1)

      data <- eval(coxph_call$data, coxenv)

      if (is_null(data)) {
        .err(sprintf("a dataset must have been supplied to `data` in the original call to `%s()` to use `vcov = \"BS\"`",
                     deparse1(coxph_call[[1]])))
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
        coxph_call$weights <- weightit_boot$weights * SW[ind]
      }
      else {
        coxph_call$weights <- SW[ind]
      }

      coxph_call$data <- data[ind,]

      fit_boot <- eval(coxph_call, coxenv)

      fit_boot$coefficients
    }

    if (is_null(cluster)) {
      boot_out <- do.call("rbind", lapply(seq_len(R), function(i) {
        ind <- sample(seq_len(nrow(data)), nrow(data), replace = TRUE)
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
  else if (vcov == "HC0") {

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

  if (is_not_null(V)) {
    fit$var <- V[seq_len(sum(!aliased)), seq_len(sum(!aliased)), drop = FALSE]
    colnames(fit$var) <- rownames(fit$var) <- names(aliased)[!aliased]
  }

  fit$vcov_type <- vcov
  fit$call <- coxph_weightit_call

  if (isFALSE(x)) fit$x <- NULL
  if (isFALSE(y)) fit$y <- NULL

  attr(fit, "cluster") <- cluster
  class(fit) <- c("coxph_weightit", "glm_weightit", class(fit))

  fit
}

#' @export
#' @rdname glm_weightit
lm_weightit <- function(formula, data, weightit,
                        vcov = NULL, cluster, R = 500,
                        offset, start = NULL, etastart, mustart,
                        control = list(...),
                        x = FALSE, y = TRUE,
                        contrasts = NULL, ...) {
  cal <- cal0 <- match.call()
  cal[[1]] <- quote(WeightIt::glm_weightit)
  cal$family = quote(stats::gaussian())
  fit <- eval.parent(cal)
  fit$call <- cal0
  fit
}

#' @exportS3Method print glm_weightit
print.glm_weightit <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {

  cat("\n", underline("Call:"), "\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n", sep = "")

  if (is_not_null(coef(x))) {
    cat("\n", underline(paste0("Coefficients",
                               if (is.character(co <- x$contrasts))
                                 paste("  [contrasts: ", apply(cbind(names(co), co), 1L, paste, collapse = "="), "]"), ":\n")), sep = "")

    print.default(format(x$coefficients, digits = digits),
                  print.gap = 2, quote = FALSE)
    cat(italic(sprintf("Standard error: %s\n",
                       .vcov_to_phrase(x$vcov_type,
                                           is_not_null(x$cluster)))))
  }
  else cat("No coefficients\n\n")

  invisible(x)
}

#' @exportS3Method stats::vcov glm_weightit
vcov.glm_weightit <- function(object, complete = TRUE, ...) {
  if (is_not_null(object[["vcov"]])) {
    return(.vcov.aliased(is.na(object$coefficients), object[["vcov"]],
                         complete = complete))
  }

  NextMethod("vcov")
}

#' @exportS3Method summary glm_weightit
#' @name glm_weightit
summary.glm_weightit <- function(object, ci = FALSE, level = .95, transform = NULL, ...) {
  chk::chk_flag(ci)

  df.r <- object$df.residual

  fam <- object$family
  dispersion <- {
    if (is.null(fam)) NaN
    else if (!is.null(fam$dispersion) && !is.na(fam$dispersion)) fam$dispersion
    else if (fam$family %in% c("poisson", "binomial", "multinomial")) 1
    else if (df.r > 0) {
      if (any(object$weights == 0))
        .wrn("observations with zero weight not used for calculating dispersion")
      sum((object$weights * object$residuals^2)[object$weights > 0])/df.r
    }
    else {
      NaN
    }
  }

  if (is.null(transform)) transform <- base::identity
  else transform <- match.fun(transform)

  coef.p <- coef(object)
  aliased <- is.na(coef.p)
  p <- sum(!aliased)
  if (p > 0) {
    p1 <- seq_len(p)
    Qr <- object$qr
    coef.p <- object$coefficients[!aliased]

    if (inherits(object, "coxph_weightit")) {
      covmat.unscaled <- NULL
      df.f <- NULL
    }
    else if (identical(fam$family, "multinomial")) {
      covmat.unscaled <- NULL
      df.f <- NCOL(Qr$qr)
    }
    else {
      covmat.unscaled <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
      dimnames(covmat.unscaled) <- list(names(coef.p), names(coef.p))
      df.f <- NCOL(Qr$qr)
    }

    covmat <- vcov(object, complete = TRUE)
    if (is_null(covmat)) {
      ci <- FALSE
      s.err <- rep(NA_real_, length(coef.p))
      pvalue <- tvalue <- s.err
    }
    else {
      s.err <- sqrt(diag(covmat)[!aliased])
      tvalue <- coef.p/s.err
      pvalue <- 2 * pnorm(-abs(tvalue))
    }
    dn <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
    dimnames(coef.table) <- list(names(coef.p), dn)
  }
  else {
    coef.table <- matrix(NA_real_, 0L, 4L)
    dimnames(coef.table) <- list(NULL, c("Estimate", "Std. Error",
                                         "z value", "Pr(>|z|)"))
    covmat <- matrix(NA_real_, 0L, 0L)
  }

  transformed_coefs <- transform(coef.table[,"Estimate"])
  if (!is.numeric(transformed_coefs) || length(transformed_coefs) != length(coef.table[,"Estimate"])) {
    .err("`transform` must return a numeric vector")
  }

  identity_transform <- all(transformed_coefs == coef.table[,"Estimate"])
  if (!identity_transform) {
    coef.table[,"Estimate"] <- transformed_coefs
    coef.table <- coef.table[, -2, drop = FALSE]
  }

  if (ci) {
    conf <- confint(object, parm = which(!aliased), level = level)
    if (!identity_transform) conf[] <- transform(conf)
    coef.table <- cbind(coef.table, conf)
  }

  keep <- match(c("call", "terms", "family", "deviance", "aic", "contrasts", "df.residual",
                  "null.deviance", "df.null", "iter", "na.action", "vcov_type"), names(object), 0L)
  ans <- c(object[keep], list(deviance.resid = residuals(object, type = "deviance"),
                              coefficients = coef.table, aliased = aliased,
                              dispersion = dispersion, df = c(object$rank, df.r, df.f),
                              cov.unscaled = covmat.unscaled, cov.scaled = covmat,
                              cluster = attr(object, "cluster"),
                              transformed = !identity_transform))

  class(ans) <- c("summary.glm_weightit", "summary.glm")

  ans
}

#' @exportS3Method print summary.glm_weightit
print.summary.glm_weightit <- function(x, digits = max(3L, getOption("digits") - 3L),
                                       signif.stars = getOption("show.signif.stars"),
                                       ...) {
  cat("\n", underline("Call:"), "\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n", sep = "")

  if (length(x$aliased) == 0L) {
    cat("\nNo Coefficients\n\n")
    return(invisible(x))
  }

  cat("\n", underline(paste0("Coefficients", if (x$transformed) " (transformed):" else ":")),
      "\n", sep = "")

  coefs <- x$coefficients

  if (!is.null(aliased <- x$aliased) && any(aliased)) {
    cn <- names(aliased)
    coefs <- matrix(NA, length(aliased), ncol(coefs), dimnames = list(cn, colnames(coefs)))
    coefs[!aliased, ] <- x$coefficients
  }
  printCoefmat(coefs, digits = digits, signif.legend = FALSE,
               na.print = ".", cs.ind = -(3:4), tst.ind = 3,
               ...)

  cat(italic(sprintf("Standard error: %s\n",
                     .vcov_to_phrase(x$vcov_type,
                                         is_not_null(x$cluster)))))

  invisible(x)
}

#' @exportS3Method stats::confint glm_weightit
confint.glm_weightit <- function(object, parm, level = 0.95, ...) {
  chk::chk_number(level)
  chk::chk_gt(level, .5)
  chk::chk_lt(level, 1)
  object$df.residual <- Inf
  confint.lm(object, parm = parm, level = level, ...)
}

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
         })
}
#' Propensity Score Weighting Using SuperLearner
#' @name method_super
#' @aliases method_super
#' @usage NULL
#'
#' @description
#' This page explains the details of estimating weights from SuperLearner-based propensity scores by setting `method = "super"` in the call to [weightit()] or [weightitMSM()]. This method can be used with binary, multi-category, and continuous treatments.
#'
#' In general, this method relies on estimating propensity scores using the SuperLearner algorithm for stacking predictions and then converting those propensity scores into weights using a formula that depends on the desired estimand. For binary and multi-category treatments, one or more binary classification algorithms are used to estimate the propensity scores as the predicted probability of being in each treatment given the covariates. For continuous treatments, regression algorithms are used to estimate generalized propensity scores as the conditional density of treatment given the covariates. This method relies on \pkgfun{SuperLearner}{SuperLearner} from the \CRANpkg{SuperLearner} package.
#'
#' ## Binary Treatments
#'
#' For binary treatments, this method estimates the propensity scores using \pkgfun{SuperLearner}{SuperLearner}. The following estimands are allowed: ATE, ATT, ATC, ATO, ATM, and ATOS. Weights can also be computed using marginal mean weighting through stratification for the ATE, ATT, and ATC. See [get_w_from_ps()] for details.
#'
#' ## Multi-Category Treatments
#'
#' For multi-category treatments, the propensity scores are estimated using several calls to \pkgfun{SuperLearner}{SuperLearner}, one for each treatment group; the treatment probabilities are not normalized to sum to 1. The following estimands are allowed: ATE, ATT, ATC, ATO, and ATM. The weights for each estimand are computed using the standard formulas or those mentioned above. Weights can also be computed using marginal mean weighting through stratification for the ATE, ATT, and ATC. See [get_w_from_ps()] for details.
#'
#' ## Continuous Treatments
#'
#' For continuous treatments, the generalized propensity score is estimated using \pkgfun{SuperLearner}{SuperLearner}. In addition, kernel density estimation can be used instead of assuming a normal density for the numerator and denominator of the generalized propensity score by setting `use.kernel = TRUE`. Other arguments to [density()] can be specified to refine the density estimation parameters. `plot = TRUE` can be specified to plot the density for the numerator and denominator, which can be helpful in diagnosing extreme weights.
#'
#' ## Longitudinal Treatments
#'
#' For longitudinal treatments, the weights are the product of the weights estimated at each time point.
#'
#' ## Sampling Weights
#'
#' Sampling weights are supported through `s.weights` in all scenarios.
#'
#' ## Missing Data
#'
#' In the presence of missing data, the following value(s) for `missing` are allowed:
#'     \describe{
#'       \item{`"ind"` (default)}{First, for each variable with missingness, a new missingness indicator variable is created which takes the value 1 if the original covariate is `NA` and 0 otherwise. The missingness indicators are added to the model formula as main effects. The missing values in the covariates are then replaced with the covariate medians (this value is arbitrary and does not affect estimation). The weight estimation then proceeds with this new formula and set of covariates. The covariates output in the resulting `weightit` object will be the original covariates with the `NA`s.
#'       }
#'     }
#'
#' ## M-estimation
#'
#' M-estimation is not supported.
#'
#' @section Additional Arguments:
#' \describe{
#'   \item{`discrete`}{if `TRUE`, uses discrete SuperLearner, which simply selects the best performing method. Default `FALSE`, which finds the optimal combination of predictions for the libraries using `SL.method`.}
#' }
#'
#' An argument to `SL.library` **must** be supplied. To see a list of available entries, use \pkgfun{SuperLearner}{listWrappers}.
#'
#' All arguments to \pkgfun{SuperLearner}{SuperLearner} can be passed through `weightit()` or `weightitMSM()`, with the following exceptions:
#'
#'   * `obsWeights` is ignored because sampling weights are passed using `s.weights`.
#'   * `method` in `SuperLearner()` is replaced with the argument `SL.method` in `weightit()`.
#'
#' For continuous treatments only, the following arguments may be supplied:
#'   \describe{
#'     \item{`density`}{A function corresponding to the conditional density of the treatment. The standardized residuals of the treatment model will be fed through this function to produce the numerator and denominator of the generalized propensity score weights. If blank, [dnorm()] is used as recommended by Robins et al. (2000). This can also be supplied as a string containing the name of the function to be called. If the string contains underscores, the call will be split by the underscores and the latter splits will be supplied as arguments to the second argument and beyond. For example, if `density = "dt_2"` is specified, the density used will be that of a t-distribution with 2 degrees of freedom. Using a t-distribution can be useful when extreme outcome values are observed (Naimi et al., 2014). Ignored if `use.kernel = TRUE` (described below).
#'     }
#'     \item{`use.kernel`}{If `TRUE`, uses kernel density estimation through the [density()] function to estimate the numerator and denominator densities for the weights. If `FALSE`, the argument to the `density` parameter is used instead.
#'     }
#'     \item{`bw`, `adjust`, `kernel`, `n`}{If `use.kernel = TRUE`, the arguments to the [density()] function. The defaults are the same as those in `density` except that `n` is 10 times the number of units in the sample.
#'     }
#'     \item{`plot`}{If `use.kernel = TRUE`, whether to plot the estimated density.
#'     }
#'   }
#'
#' ## Balance SuperLearner
#'
#' In addition to the methods allowed by `SuperLearner()`, one can specify `SL.method = "method.balance"` to use "Balance SuperLearner" as described by Pirracchio and Carone (2018), wherein covariate balance is used to choose the optimal combination of the predictions from the methods specified with `SL.library`. Coefficients are chosen (one for each prediction method) so that the weights generated from the weighted combination of the predictions optimize a balance criterion, which must be set with the `criterion` argument, described below.
#'   \describe{
#'     \item{`criterion`}{A string describing the balance criterion used to select the best weights. See \pkgfun{cobalt}{bal.compute} for allowable options for each treatment type. For binary and multi-category treatments, the default is `"smd.mean"`, which minimizes the average absolute standard mean difference among the covariates between treatment groups. For continuous treatments, the default is `"p.mean"`, which minimizes the average absolute Pearson correlation between the treatment and covariates.
#'     }
#'   }
#'   Note that this implementation differs from that of Pirracchio and Carone (2018) in that here, balance is measured only on the terms included in the model formula (i.e., and not their interactions unless specifically included), and balance results from a sample weighted using the estimated predicted values as propensity scores, not a sample matched using propensity score matching on the predicted values. Binary and continuous treatments are supported, but currently multi-category treatments are not.
#'
#' @section Additional Outputs:
#' \describe{
#'   \item{`info`}{
#'     For binary and continuous treatments, a list with two entries, `coef` and `cvRisk`. For multi-category treatments, a list of lists with these two entries, one for each treatment level.
#'     \describe{
#'       \item{`coef`}{
#'         The coefficients in the linear combination of the predictions from each method in `SL.library`. Higher values indicate that the corresponding method plays a larger role in determining the resulting predicted value, and values close to zero indicate that the method plays little role in determining the predicted value. When `discrete = TRUE`, these correspond to the coefficients that would have been estimated had `discrete` been `FALSE`.
#'       }
#'       \item{`cvRisk`}{
#'         The cross-validation risk for each method in `SL.library`. Higher values indicate that the method has worse cross-validation accuracy. When `SL.method = "method.balance"`, the sample weighted balance statistic requested with `criterion`. Higher values indicate worse balance.
#'       }
#'     }
#'   }
#'   \item{`obj`}{
#'     When `include.obj = TRUE`, the SuperLearner fit(s) used to generate the predicted values. For binary and continuous treatments, the output of the call to \pkgfun{SuperLearner}{SuperLearner}. For multi-category treatments, a list of outputs to calls to `SuperLearner::SuperLearner()`.
#'   }
#' }
#'
#' @details
#' SuperLearner works by fitting several machine learning models to the treatment and covariates and then taking a weighted combination of the generated predicted values to use as the propensity scores, which are then used to construct weights. The machine learning models used are supplied using the `SL.library` argument; the more models are supplied, the higher the chance of correctly modeling the propensity score. The predicted values are combined using the method supplied in the `SL.method` argument (which is nonnegative least squares by default). A benefit of SuperLearner is that, asymptotically, it is guaranteed to perform as well as or better than the best-performing method included in the library. Using Balance SuperLearner by setting `SL.method = "method.balance"` works by selecting the combination of predicted values that minimizes an imbalance measure.
#'
#' @note
#' Some methods formerly available in \pkg{SuperLearner} are now in \pkg{SuperLearnerExtra}, which can be found on GitHub at \url{https://github.com/ecpolley/SuperLearnerExtra}.
#'
#' The `criterion` argument used to be called `stop.method`, which is its name in \pkg{twang}. `stop.method` still works for backward compatibility. Additionally, the criteria formerly named as `es.mean`, `es.max`, and `es.rms` have been renamed to `smd.mean`, `smd.max`, and `smd.rms`. The former are used in \pkg{twang} and will still work with `weightit()` for backward compatibility.
#'
#' @seealso
#' [weightit()], [weightitMSM()], [get_w_from_ps()]
#'
#' @references
#' ## Binary treatments
#'
#' Pirracchio, R., Petersen, M. L., & van der Laan, M. (2015). Improving Propensity Score Estimators’ Robustness to Model Misspecification Using Super Learner. American Journal of Epidemiology, 181(2), 108–119. \doi{10.1093/aje/kwu253}
#'
#' ## Multi-Category Treatments
#'
#' Imai, K., & Ratkovic, M. (2014). Covariate balancing propensity score. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 76(1), 243–263.
#'
#' ## Continuous treatments
#'
#' Kreif, N., Grieve, R., Díaz, I., & Harrison, D. (2015). Evaluation of the Effect of a Continuous Treatment: A Machine Learning Approach with an Application to Treatment for Traumatic Brain Injury. Health Economics, 24(9), 1213–1228. \doi{10.1002/hec.3189}
#'
#' ## Balance SuperLearner (`SL.method = "method.balance"`)
#'
#' Pirracchio, R., & Carone, M. (2018). The Balance Super Learner: A robust adaptation of the Super Learner to improve estimation of the average treatment effect in the treated based on propensity score matching. Statistical Methods in Medical Research, 27(8), 2504–2518. \doi{10.1177/0962280216682055}
#'
#' See [`method_glm`] for additional references.
#'
#' @examplesIf all(sapply(c("SuperLearner", "MASS"), requireNamespace, quietly = TRUE))
#' library("cobalt")
#' data("lalonde", package = "cobalt")
#'
#' #Balancing covariates between treatment groups (binary)
#' (W1 <- weightit(treat ~ age + educ + married +
#'                   nodegree + re74, data = lalonde,
#'                 method = "super", estimand = "ATT",
#'                 SL.library = c("SL.glm", "SL.stepAIC",
#'                                "SL.glm.interaction")))
#' summary(W1)
#' bal.tab(W1)
#' \donttest{
#'   #Balancing covariates with respect to race (multi-category)
#'   (W2 <- weightit(race ~ age + educ + married +
#'                     nodegree + re74, data = lalonde,
#'                   method = "super", estimand = "ATE",
#'                   SL.library = c("SL.glm", "SL.stepAIC",
#'                                  "SL.glm.interaction")))
#'   summary(W2)
#'   bal.tab(W2)
#'
#'   #Balancing covariates with respect to re75 (continuous)
#'   #assuming t(8) conditional density for treatment
#'   (W3 <- weightit(re75 ~ age + educ + married +
#'                     nodegree + re74, data = lalonde,
#'                   method = "super", density = "dt_8",
#'                   SL.library = c("SL.glm", "SL.ridge",
#'                                  "SL.glm.interaction")))
#'   summary(W3)
#'   bal.tab(W3)
#' }
#' #Balancing covariates between treatment groups (binary)
#' # using balance SuperLearner to minimize the maximum
#' # KS statistic
#' (W4 <- weightit(treat ~ age + educ + married +
#'                   nodegree + re74, data = lalonde,
#'                 method = "super", estimand = "ATT",
#'                 SL.library = c("SL.glm", "SL.stepAIC",
#'                                "SL.lda"),
#'                 SL.method = "method.balance",
#'                 criterion = "ks.max"))
#' summary(W4)
#' bal.tab(W4, stats = c("m", "ks"))
NULL

weightit2super <- function(covs, treat, s.weights, subset, estimand, focal,
                           stabilize, subclass, missing, verbose, ...) {
  rlang::check_installed("SuperLearner")

  A <- list(...)

  covs <- covs[subset, , drop = FALSE]
  treat <- factor(treat[subset])
  s.weights <- s.weights[subset]

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

    available.criteria <- cobalt::available.stats("binary")

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

    init <- cobalt::bal.init(covs,
                             treat = treat,
                             stat = criterion,
                             estimand = estimand,
                             s.weights = s.weights,
                             focal = focal,
                             ...)

    sneaky <- 0
    attr(sneaky, "vals") <- list(init = init, estimand = estimand)
    A[["control"]] <- list(trimLogit = sneaky)

    A[["SL.method"]] <- .method.balance()
  }

  t.lev <- get_treated_level(treat)
  treat <- binarize(treat, one = t.lev)

  tryCatch({verbosely({
    fit <- do.call(SuperLearner::SuperLearner, list(Y = treat,
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

  ps <- {
    if (discrete) fit$library.predict[,which.min(fit$cvRisk)]
    else fit$SL.predict
  }

  info <- list(coef = fit$coef,
               cvRisk = fit$cvRisk)

  #ps should be matrix of probs for each treat
  #Computing weights
  w <- .get_w_from_ps_internal_bin(ps = ps, treat = treat, estimand,
                                   stabilize = stabilize, subclass = subclass)

  list(w = w, ps = ps, info = info, fit.obj = fit)
}

weightit2super.multi <- function(covs, treat, s.weights, subset, estimand, focal,
                                 stabilize, subclass, missing, verbose, ...) {
  rlang::check_installed("SuperLearner")

  A <- list(...)

  covs <- covs[subset, , drop = FALSE]
  treat <- factor(treat[subset])
  s.weights <- s.weights[subset]

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
    .err("\"method.balance\" cannot be used with multi-category treatments")
  }

  fit.list <- info <- make_list(levels(treat))
  ps <- make_df(levels(treat), nrow = length(treat))

  for (i in levels(treat)) {

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

    ps[[i]] <- {
      if (discrete) fit.list[[i]]$library.predict[,which.min(fit.list[[i]]$cvRisk)]
      else fit.list[[i]]$SL.predict
    }

    info[[i]] <- list(coef = fit.list[[i]]$coef,
                      cvRisk = fit.list[[i]]$cvRisk)
  }

  #ps should be matrix of probs for each treat
  #Computing weights
  w <- get_w_from_ps(ps = ps, treat = treat, estimand, focal, stabilize = stabilize, subclass = subclass)

  p.score <- NULL

  list(w = w, ps = p.score, info = info, fit.obj = fit.list)
}

weightit2super.cont <- function(covs, treat, s.weights, subset, stabilize, missing, ps, verbose, ...) {
  rlang::check_installed("SuperLearner")

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
  densfun <- get_dens_fun(use.kernel = isTRUE(A[["use.kernel"]]), bw = A[["bw"]],
                          adjust = A[["adjust"]], kernel = A[["kernel"]],
                          n = A[["n"]], treat = treat, density = A[["density"]])

  #Stabilization - get dens.num
  dens.num <- densfun(treat - mean(treat), s.weights)

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

    available.criteria <- cobalt::available.stats("continuous")

    criterion <- A[["criterion"]]
    criterion <- match_arg(criterion, available.criteria)

    init <- cobalt::bal.init(covs,
                             treat = treat,
                             stat = criterion,
                             s.weights = s.weights,
                             ...)

    sneaky <- 0
    attr(sneaky, "vals") <- list(init = init,
                                 dens.num = dens.num,
                                 densfun = densfun)
    B[["control"]] <- list(trimLogit = sneaky)

    B[["SL.method"]] <- .method.balance.cont()
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
  dens.denom <- densfun(treat - gp.score, s.weights)

  w <- dens.num/dens.denom

  if (isTRUE(A[["use.kernel"]]) && isTRUE(A[["plot"]])) {
    d.n <- attr(dens.num, "density")
    d.d <- attr(dens.denom, "density")
    plot_density(d.n, d.d)
  }

  info <- list(coef = fit$coef,
               cvRisk = fit$cvRisk)

  list(w = w, info = info, fit.obj = fit)
}

#For balance SuperLearner
.method.balance <- function() {

  out <- list(
    # require allows you to pass a character vector with required packages
    # use NULL if no required packages
    require = "cobalt",

    # computeCoef is a function that returns a list with two elements:
    # 1) coef: the weights (coefficients) for each algorithm
    # 2) cvRisk: the V-fold CV risk for each algorithm
    computeCoef = function(Z, Y, libraryNames, obsWeights, control, verbose, ...) {
      estimand <- attr(control$trimLogit, "vals")$estimand
      init <- attr(control$trimLogit, "vals")$init

      tol <- .001
      for (i in seq_col(Z)) {
        Z[Z[,i] < tol, i] <- tol
        Z[Z[,i] > 1-tol, i] <- 1-tol
      }
      w_mat <- .get_w_from_ps_internal_array(Z, treat = Y, estimand = estimand)
      cvRisk <- apply(w_mat, 2, cobalt::bal.compute, x = init)

      names(cvRisk) <- libraryNames

      loss <- function(coefs) {
        ps <- crossprod(t(Z), coefs/sum(coefs))
        w <- get_w_from_ps(ps, Y, estimand)
        cobalt::bal.compute(init, weights = w)
      }
      fit <- optim(rep(1/ncol(Z), ncol(Z)), loss, method = "L-BFGS-B", lower = 0, upper = 1)
      coef <- fit$par

      list(cvRisk = cvRisk, coef = coef/sum(coef))
    },

    # computePred is a function that takes the weights and the predicted values
    # from each algorithm in the library and combines them based on the model to
    # output the super learner predicted values
    computePred = function(predY, coef, control, ...) {
      crossprod(t(predY), coef/sum(coef))
    }
  )

  out
}

.method.balance.cont <- function() {

  out <- list(
    # require allows you to pass a character vector with required packages
    # use NULL if no required packages
    require = "cobalt",

    # computeCoef is a function that returns a list with two elements:
    # 1) coef: the weights (coefficients) for each algorithm
    # 2) cvRisk: the V-fold CV risk for each algorithm
    computeCoef = function(Z, Y, libraryNames, obsWeights, control, verbose, ...) {
      dens.num <- attr(control$trimLogit, "vals")$dens.num
      densfun <- attr(control$trimLogit, "vals")$densfun
      init <- attr(control$trimLogit, "vals")$init

      w_mat <- apply(Z, 2, function(gp.score) {
        dens.num/densfun(Y - gp.score, obsWeights)
      })

      cvRisk <- apply(w_mat, 2, cobalt::bal.compute, x = init)
      names(cvRisk) <- libraryNames

      loss <- function(coefs) {
        gp.score <- crossprod(t(Z), coefs/sum(coefs))
        w <- dens.num/densfun(Y - gp.score, obsWeights)
        cobalt::bal.compute(init, weights = w)
      }

      fit <- optim(rep(1/ncol(Z), ncol(Z)), loss, method = "L-BFGS-B", lower = 0, upper = 1)
      coef <- fit$par

      list(cvRisk = cvRisk, coef = coef/sum(coef))
    },

    # computePred is a function that takes the weights and the predicted values
    # from each algorithm in the library and combines them based on the model to
    # output the super learner predicted values
    computePred = function(predY, coef, control, ...) {
      crossprod(t(predY), coef/sum(coef))
    }
  )

  out
}

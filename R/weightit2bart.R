#' Propensity Score Weighting Using BART
#' @name method_bart
#' @aliases method_bart
#' @usage NULL
#'
#' @description This page explains the details of estimating weights from
#' Bayesian additive regression trees (BART)-based propensity scores by setting
#' `method = "bart"` in the call to [weightit()] or [weightitMSM()]. This method
#' can be used with binary, multi-category, and continuous treatments.
#'
#' In general, this method relies on estimating propensity scores using BART and
#' then converting those propensity scores into weights using a formula that
#' depends on the desired estimand. This method relies on \pkgfun{dbarts}{bart2}
#' from the \CRANpkg{dbarts} package.
#'
#' ## Binary Treatments
#'
#' For binary treatments, this method estimates the propensity scores using
#' \pkgfun{dbarts}{bart2}. The following estimands are allowed: ATE, ATT, ATC,
#' ATO, ATM, and ATOS. Weights can also be computed using marginal mean
#' weighting through stratification for the ATE, ATT, and ATC. See
#' [get_w_from_ps()] for details.
#'
#' ## Multi-Category Treatments
#'
#' For multi-category treatments, the propensity scores are estimated using
#' several calls to \pkgfun{dbarts}{bart2}, one for each treatment group; the
#' treatment probabilities are not normalized to sum to 1. The following
#' estimands are allowed: ATE, ATT, ATC, ATO, and ATM. The weights for each
#' estimand are computed using the standard formulas or those mentioned above.
#' Weights can also be computed using marginal mean weighting through
#' stratification for the ATE, ATT, and ATC. See [get_w_from_ps()] for details.
#'
#' ## Continuous Treatments For continuous treatments, weights are estimated as
#' \eqn{w_i = f_A(a_i) / f_{A|X}(a_i)}, where \eqn{f_A(a_i)} (known as the
#' stabilization factor) is the unconditional density of treatment evaluated the
#' observed treatment value and \eqn{f_{A|X}(a_i)} (known as the generalized
#' propensity score) is the conditional density of treatment given the
#' covariates evaluated at the observed value of treatment. The shape of
#' \eqn{f_A(.)} and \eqn{f_{A|X}(.)} is controlled by the `density` argument
#' described below (normal distributions by default), and the predicted values
#' used for the mean of the conditional density are estimated using BART as
#' implemented in \pkgfun{dbarts}{bart2}. Kernel density estimation can be used
#' instead of assuming a specific density for the numerator and denominator by
#' setting `density = "kernel"`. Other arguments to [density()] can be specified
#' to refine the density estimation parameters.
#'
#' ## Longitudinal Treatments
#'
#' For longitudinal treatments, the weights are the product of the weights
#' estimated at each time point.
#'
#' ## Sampling Weights
#'
#' Sampling weights are not supported.
#'
#' ## Missing Data
#'
#' In the presence of missing data, the following value(s) for `missing` are
#' allowed:
#'     \describe{
#'       \item{`"ind"` (default)}{
#'         First, for each variable with missingness, a new missingness indicator variable is created which takes the value 1 if the original covariate is `NA` and 0 otherwise. The missingness indicators are added to the model formula as main effects. The missing values in the covariates are then replaced with the covariate medians. The weight estimation then proceeds with this new formula and set of covariates. The covariates output in the resulting `weightit` object will be the original covariates with the `NA`s.
#'       }
#'     }
#' ## M-estimation
#'
#' M-estimation is not supported.
#'
#' @section Additional Arguments:
#'
#'   All arguments to \pkgfun{dbarts}{bart2} can be passed through `weightit()`
#'   or `weightitMSM()`, with the following exceptions:
#'
#'  * `test`, `weights`,`subset`, `offset.test` are ignored
#'  * `combine.chains` is always set to `TRUE`
#'  * `sampleronly` is always set to `FALSE`
#'
#'   For continuous treatments only, the following arguments may be supplied:
#'   \describe{
#'     \item{`density`}{A function corresponding to the conditional density of the treatment. The standardized residuals of the treatment model will be fed through this function to produce the numerator and denominator of the generalized propensity score weights. If blank, [dnorm()] is used as recommended by Robins et al. (2000). This can also be supplied as a string containing the name of the function to be called. If the string contains underscores, the call will be split by the underscores and the latter splits will be supplied as arguments to the second argument and beyond. For example, if `density = "dt_2"` is specified, the density used will be that of a t-distribution with 2 degrees of freedom. Using a t-distribution can be useful when extreme outcome values are observed (Naimi et al., 2014).
#'
#' Can also be `"kernel"` to use kernel density estimation, which calls [density()] to estimate the numerator and denominator densities for the weights. (This used to be requested by setting `use.kernel = TRUE`, which is now deprecated.)}
#'     \item{`bw`, `adjust`, `kernel`, `n`}{If `density = "kernel"`, the arguments to [density()]. The defaults are the same as those in `density()` except that `n` is 10 times the number of units in the sample.}
#'     \item{`plot`}{If `density = "kernel"`, whether to plot the estimated densities.}
#'   }
#'
#' @section Additional Outputs:
#' \describe{
#'   \item{`obj`}{
#'     When `include.obj = TRUE`, the `bart2` fit(s) used to generate the predicted values. With multi-category treatments, this will be a list of the fits; otherwise, it will be a single fit. The predicted probabilities used to compute the propensity scores can be extracted using \pkgfun2{dbarts}{bart}{fitted}.
#'   }
#' }
#'
#' @details BART works by fitting a sum-of-trees model for the treatment or
#' probability of treatment. The number of trees is determined by the `n.trees`
#' argument. Bayesian priors are used for the hyperparameters, so the result is
#' a posterior distribution of predicted values for each unit. The mean of these
#' for each unit is taken for use in computing the (generalized) propensity
#' score. Although the hyperparameters governing the priors can be modified by
#' supplying arguments to `weightit()` that are passed to the BART fitting
#' function, the default values tend to work well and require little
#' modification (though the defaults differ for continuous and categorical
#' treatments; see the \pkgfun{dbarts}{bart2} documentation for details). Unlike
#' many other machine learning methods, no loss function is optimized and the
#' hyperparameters do not need to be tuned (e.g., using cross-validation),
#' though performance can benefit from tuning. BART tends to balance sparseness
#' with flexibility by using very weak learners as the trees, which makes it
#' suitable for capturing complex functions without specifying a particular
#' functional form and without overfitting.
#'
#' ## Reproducibility
#'
#' BART has a random component, so some work must be done to ensure
#' reproducibility across runs. See the *Reproducibility* section at
#' \pkgfun{dbarts}{bart2} for more details. To ensure reproducibility, one can
#' do one of two things: 1) supply an argument to `seed`, which is passed to
#' `dbarts::bart2()` and sets the seed for single- and multi-threaded uses, or
#' 2) call [set.seed()], though this only ensures reproducibility when using
#' single-threading, which can be requested by setting `n.threads = 1`. Note
#' that to ensure reproducibility on any machine, regardless of the number of
#' cores available, one should use single-threading and either supply `seed` or
#' call `set.seed()`.
#'
#' @seealso [weightit()], [weightitMSM()], [get_w_from_ps()]
#'
#' [`method_super`] for stacking predictions from several machine learning
#' methods, including BART.
#'
#' @references Hill, J., Weiss, C., & Zhai, F. (2011). Challenges With
#' Propensity Score Strategies in a High-Dimensional Setting and a Potential
#' Alternative. *Multivariate Behavioral Research*, 46(3), 477–513.
#' \doi{10.1080/00273171.2011.570161}
#'
#' Chipman, H. A., George, E. I., & McCulloch, R. E. (2010). BART: Bayesian
#' additive regression trees. *The Annals of Applied Statistics*, 4(1), 266–298.
#' \doi{10.1214/09-AOAS285}
#'
#' Note that many references that deal with BART for causal inference focus on
#' estimating potential outcomes with BART, not the propensity scores, and so
#' are not directly relevant when using BART to estimate propensity scores for
#' weights.
#'
#' See [`method_glm`] for additional references on propensity score weighting
#' more generally.
#'
#' @examplesIf requireNamespace("dbarts", quietly = TRUE)
#' \donttest{library("cobalt")
#' data("lalonde", package = "cobalt")
#'
#' #Balancing covariates between treatment groups (binary)
#' (W1 <- weightit(treat ~ age + educ + married +
#'                   nodegree + re74, data = lalonde,
#'                 method = "bart", estimand = "ATT"))
#' summary(W1)
#' bal.tab(W1)
#'
#' #Balancing covariates with respect to race (multi-category)
#' (W2 <- weightit(race ~ age + educ + married +
#'                 nodegree + re74, data = lalonde,
#'                 method = "bart", estimand = "ATE"))
#' summary(W2)
#' bal.tab(W2)
#'
#' #Balancing covariates with respect to re75 (continuous)
#' #assuming t(3) conditional density for treatment
#' (W3 <- weightit(re75 ~ age + educ + married +
#'                   nodegree + re74, data = lalonde,
#'                 method = "bart", density = "dt_3"))
#'  summary(W3)
#'  bal.tab(W3)
#' }
NULL

weightit2bart <- function(covs, treat, s.weights, subset, estimand, focal, stabilize,
                          subclass, missing, verbose, ...) {

  rlang::check_installed("dbarts", version = "0.9-23")

  covs <- covs[subset, , drop = FALSE]
  treat <- treat[subset]
  s.weights <- s.weights[subset]

  missing <- .process_missing2(missing, covs)

  if (missing == "ind") {
    covs <- add_missing_indicators(covs)
  }

  if (ncol(covs) > 1L) {
    colinear.covs.to.remove <- colnames(covs)[colnames(covs) %nin% colnames(make_full_rank(covs))]
    covs <- covs[, colnames(covs) %nin% colinear.covs.to.remove, drop = FALSE]
  }

  for (i in seq_col(covs)) {
    covs[, i] <- .make_closer_to_1(covs[, i])
  }

  t.lev <- get_treated_level(treat, estimand, focal)
  treat <- binarize(treat, one = t.lev)

  A <- ...mget(setdiff(c(names(formals(dbarts::bart2)),
                         names(formals(dbarts::dbartsControl))),
                       c("offset.test", "weights", "subset", "test")))

  A[["data"]] <- treat
  A[["formula"]] <- covs
  A[["keepCall"]] <- FALSE
  A[["combineChains"]] <- TRUE
  A[["verbose"]] <- FALSE #necessary to prevent crash

  bart.call <- as.call(c(list(quote(dbarts::bart2)), A))

  tryCatch({verbosely({
    fit <- eval(bart.call)
  }, verbose = verbose)},
  error = function(e) {
    e. <- conditionMessage(e)
    .err("(from `dbarts::bart2()`) ", e., tidy = FALSE)
  })

  p.score <- fitted(fit)

  w <- .get_w_from_ps_internal_bin(ps = p.score, treat = treat, estimand,
                                   stabilize = stabilize, subclass = subclass)

  list(w = w, ps = p.score, fit.obj = fit)
}

weightit2bart.multi <-  function(covs, treat, s.weights, subset, estimand, focal, stabilize,
                                 subclass, missing, verbose, ...) {
  rlang::check_installed("dbarts", version = "0.9-23")

  covs <- covs[subset, , drop = FALSE]
  treat <- factor(treat[subset])
  s.weights <- s.weights[subset]

  missing <- .process_missing2(missing, covs)

  if (missing == "ind") {
    covs <- add_missing_indicators(covs)
  }

  if (ncol(covs) > 1L) {
    colinear.covs.to.remove <- colnames(covs)[colnames(covs) %nin% colnames(make_full_rank(covs))]
    covs <- covs[, colnames(covs) %nin% colinear.covs.to.remove, drop = FALSE]
  }

  for (i in seq_col(covs)) {
    covs[, i] <- .make_closer_to_1(covs[, i])
  }

  ps <- make_df(levels(treat), nrow = length(treat))

  A <- ...mget(setdiff(c(names(formals(dbarts::bart2)),
                         names(formals(dbarts::dbartsControl))),
                       c("offset.test", "weights", "subset", "test")))

  A[["formula"]] <- covs
  A[["keepCall"]] <- FALSE
  A[["combineChains"]] <- TRUE
  A[["verbose"]] <- FALSE #necessary to prevent crash

  fit.list <- make_list(levels(treat))

  for (i in levels(treat)) {
    A[["data"]] <- as.integer(treat == i)
    bart.call <- as.call(c(list(quote(dbarts::bart2)), A))

    tryCatch({verbosely({
      fit.list[[i]] <- eval(bart.call)
    }, verbose = verbose)},
    error = function(e) {
      e. <- conditionMessage(e)
      .err("(from `dbarts::bart2()`) ", e., tidy = FALSE)
    })

    ps[[i]] <- fitted(fit.list[[i]])
  }

  #ps should be matrix of probs for each treat
  #Computing weights
  w <- .get_w_from_ps_internal_multi(ps = ps, treat = treat, estimand, focal,
                                     stabilize = stabilize, subclass = subclass)

  list(w = w, fit.obj = fit.list)
}

weightit2bart.cont <- function(covs, treat, s.weights, subset, stabilize, missing, ps, verbose, ...) {

  covs <- covs[subset, , drop = FALSE]
  treat <- treat[subset]
  s.weights <- s.weights[subset]

  missing <- .process_missing2(missing, covs)

  if (missing == "ind") {
    covs <- add_missing_indicators(covs)
  }

  for (i in seq_col(covs)) {
    covs[, i] <- .make_closer_to_1(covs[, i])
  }

  #Process density params
  densfun <- .get_dens_fun(use.kernel = isTRUE(...get("use.kernel")),
                           bw = ...get("bw"),
                           adjust = ...get("adjust"),
                           kernel = ...get("kernel"),
                           n = ...get("n"),
                           treat = treat,
                           density = ...get("density"),
                           weights = s.weights)

  A <- ...mget(setdiff(c(names(formals(dbarts::bart2)),
                         names(formals(dbarts::dbartsControl))),
                       c("offset.test", "weights", "subset", "test")))

  A[["formula"]] <- covs
  A[["data"]] <- treat
  A[["keepCall"]] <- FALSE
  A[["combineChains"]] <- TRUE
  A[["verbose"]] <- FALSE #necessary to prevent crash

  bart.call <- as.call(c(list(quote(dbarts::bart2)), A))

  #Estimate GPS
  tryCatch({verbosely({
    fit <- eval(bart.call)
  }, verbose = verbose)},
  error = function(e) {
    e. <- conditionMessage(e)
    .err("(from `dbarts::bart2()`) ", e., tidy = FALSE)
  })

  r <- residuals(fit)

  #Get weights
  log.dens.num <- densfun(scale_w(treat, s.weights), log = TRUE)
  log.dens.denom <- densfun(r / sqrt(col.w.v(r, s.weights)), log = TRUE)

  w <- exp(log.dens.num - log.dens.denom)

  if (isTRUE(...get("plot"))) {
    d.n <- attr(log.dens.num, "density")
    d.d <- attr(log.dens.denom, "density")
    plot_density(d.n, d.d, log = TRUE)
  }

  list(w = w, fit.obj = fit)
}

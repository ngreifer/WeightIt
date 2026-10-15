#' Propensity Score Weighting Using BART
#' @name method_bart
#' @usage NULL
#'
#' @description
#' This page explains the details of estimating weights from
#' Bayesian additive regression trees (BART)-based propensity scores by setting
#' `method = "bart"` in the call to [weightit()] or [weightitMSM()]. This method
#' can be used with binary, multi-category, and continuous treatments, as well as for estimating censoring weights.
#'
#' In general, this method relies on estimating propensity scores using BART and
#' then converting those propensity scores into weights using a formula that
#' depends on the desired estimand. This method relies on \pkgfun{bartisan}{bartisan}
#' from the \CRANpkg{bartisan} package, which fits BART models with soft decision
#' rules and with likelihoods beyond the Gaussian using the sampler of Linero
#' (2025).
#'
#' ## Binary Treatments
#'
#' For binary treatments, this method estimates the propensity scores using
#' \pkgfun{bartisan}{bartisan} with a binomial likelihood (logistic by default;
#' see `link` below). The following estimands are allowed: ATE, ATT, ATC,
#' ATO, ATM, and ATOS. Weights can also be computed using marginal mean
#' weighting through stratification for the ATE, ATT, and ATC. See
#' [get_w_from_ps()] for details.
#'
#' ## Multi-Category Treatments
#'
#' For multi-category treatments, the propensity scores are estimated using a
#' single call to \pkgfun{bartisan}{bartisan} with a multinomial likelihood
#' (\pkgfun{bartisan}{multinomial}), which fits one forest per treatment level,
#' so the predicted probabilities for each unit sum to 1. The following
#' estimands are allowed: ATE, ATT, ATC, ATO, and ATM. The weights for each
#' estimand are computed using the standard formulas. Weights can also be
#' computed using marginal mean weighting through stratification for the ATE,
#' ATT, and ATC. See [get_w_from_ps()] for details.
#'
#' ## Continuous Treatments
#'
#' For continuous treatments, weights are estimated as
#' \eqn{w_i = f_A(a_i) / f_{A|X}(a_i)}, where \eqn{f_A(a_i)} (known as the
#' stabilization factor) is the unconditional density of treatment evaluated at the
#' observed treatment value and \eqn{f_{A|X}(a_i)} (known as the generalized
#' propensity score) is the conditional density of treatment given the
#' covariates evaluated at the observed treatment value. The conditional density
#' is estimated directly by a BART model whose mean is a sum of trees and whose
#' error distribution is a Dirichlet process mixture of normals estimated from
#' the data (\pkgfun{bartisan}{dpm}; George et al., 2019), so no parametric form
#' for the density needs to be assumed. \eqn{f_A(.)} is estimated by
#' marginalizing over \eqn{f_{A|X}(.)}. A different family can be requested with
#' the `family` argument described below; the `density` argument used by other
#' methods does not apply.
#'
#' ## Multilevel Treatment Models
#'
#' When the model `formula` contains \CRANpkg{lme4}-style random intercept terms
#' (e.g., `treat ~ x1 + x2 + (1 | school)`), the random intercepts are included
#' in the BART model directly; see *Random Intercepts* at
#' \pkgfun{bartisan}{bartisan}. This can improve balance when units are
#' clustered (e.g., patients within hospitals or students within schools).
#' Multiple grouping factors (e.g., `(1 | school) + (1 | district)`) and nesting
#' (`(1 | district/school)`) are allowed, but random slopes are not; a variable
#' whose effect varies by group belongs in the fixed part of the formula, where
#' a tree can split on it and on the group together. The grouping variables
#' are taken from `data`. The estimated propensity scores (or conditional
#' densities for continuous treatments) include the estimated random effects;
#' i.e., they are cluster-specific. As with the single-level case, M-estimation
#' is not supported.
#'
#' ## Censoring Weights
#'
#' For censoring weights, requested by wrapping the censoring indicator in [.cens()], a single BART model of the probability of being censored is fit, and the weights are `1/P(C = 0 | X)` for the units still under observation and 0 for the censored units. The returned propensity score is the probability of *being censored*.
#'
#' ## Longitudinal Treatments
#'
#' For longitudinal treatments, the weights are the product of the weights
#' estimated at each time point.
#'
#' ## Sampling Weights
#'
#' Sampling weights are supported through `s.weights` for binary and
#' multi-category treatments and for censoring. Before being supplied to
#' \pkgfun{bartisan}{bartisan} as prior weights, they are scaled to have a mean
#' of 1 among the units in the model, so that the fit does not depend on their
#' scale. For continuous treatments, the default `family = bartisan::dpm()`
#' cannot take fractional weights, so the sampling weights are used only in
#' computing the marginal density of the treatment in the numerator of the
#' weights and not in fitting the conditional density; to include them in the
#' fit, supply `family = gaussian()` or `family = bartisan::gaussian_ls()`.
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
#'
#' ## M-estimation
#'
#' M-estimation is not supported.
#'
#' @section Additional Arguments:
#'
#' All arguments to \pkgfun{bartisan}{bartisan_control} can be passed through
#' `weightit()` or `weightitMSM()` to control the sampler and the priors, including
#' `num_trees`, `num_draws`, `num_burn`, `chains`, `gate`, and `sparsity`. Their
#' defaults are those of `bartisan_control()`. For compatibility with earlier
#' versions of this method, which used \pkg{dbarts}, the following arguments are
#' translated to their `bartisan_control()` equivalents when supplied: `n.trees`
#' to `num_trees`, `n.samples` to `num_draws`, `n.burn` to `num_burn`, and
#' `n.chains` to `chains`. A natively named argument takes precedence over its
#' translation. `n.threads` is ignored (see *Parallel Processing* below), and
#' `seed` is ignored with a warning (see *Reproducibility* below).
#'
#' The following additional argument may be supplied:
#' \describe{
#'   \item{`link`}{for binary treatments and censoring weights, the link function of the binomial model. Allowable options include `"logit"` (the default), `"probit"`, and `"cloglog"`. For multi-category treatments, the link function of the multinomial model, which can be `"logit"` (the default) or `"probit"`. Ignored for continuous treatments.}
#'   \item{`use.offset`}{`logical`; whether to use the linear predictor resulting from a generalized linear model as an offset to the BART model. If `TRUE`, this fits a binomial regression model with the link specified by `link` (for binary treatments) or a linear regression model (for continuous treatments) and supplies the linear predictor to the `offset` argument of `bartisan()`, overriding any supplied argument to `offset`. Default is `FALSE` to omit this offset. Only allowed for binary and continuous treatments.
#'   }
#' }
#'
#' For multi-category treatments, the following arguments may be supplied:
#' \describe{
#'   \item{`reference`}{the treatment level to hold as the reference category in the multinomial model. Default is `NULL`, which with the logit link fits one forest per level and leaves the prior symmetric in the levels; see \pkgfun{bartisan}{multinomial}.}
#'   \item{`replicates`}{for `link = "probit"`, the number of simulation draws used to compute the category probabilities, which have no closed form. Default is 200. See \pkgfun{bartisan}{multinomial}.}
#' }
#'
#' For binary and multi-category treatments, the following arguments may be supplied:
#' \describe{
#'   \item{`subclass`}{`integer`; the number of subclasses to use for computing
#'   weights using marginal mean weighting through stratification (MMWS). If `NULL`,
#'   standard inverse probability weights (and their extensions) will be
#'   computed; if a number greater than 1, subclasses will be formed and weights
#'   will be computed based on subclass membership. See [get_w_from_ps()] for
#'   details and references.}
#' }
#'
#' For continuous treatments, the following arguments may be supplied:
#'   \describe{
#'     \item{`family`}{the response family used to model the conditional density of the treatment, given as a family object or the name of one; see [`bartisan-families`][bartisan::bartisan-families]. Default is \pkgfun{bartisan}{dpm} for a Dirichlet process mixture of normals for the error distribution, which estimates the shape of the density from the data. Other options include `gaussian()` for a normal density with a single variance and \pkgfun{bartisan}{gaussian_ls} for a normal density whose variance depends on the covariates; these two accept sampling weights, whereas `dpm()` does not. Arguments to the family (e.g., `bartisan::dpm(nu = 5)`) are supplied by constructing the family object.}
#'   }
#'
#' The `density`, `bw`, `adjust`, `kernel`, and `n` arguments that other methods
#' use to specify the conditional density are not used by this method and are
#' ignored with a warning, since the density is estimated by the model itself.
#'
#' @section Additional Outputs:
#' \describe{
#'   \item{`obj`}{
#'     When `include.obj = TRUE`, the `bartisan_fit` object used to generate the predicted values. For multi-category treatments, this is a single fit, since one multinomial model is fit for all treatment levels. The predicted probabilities used to compute the propensity scores can be extracted using \pkgfun{bartisan}{predict.bartisan_fit} with `type = "prob"`, and the conditional density of a continuous treatment with `type = "density"`.
#'   }
#' }
#'
#' @details
#' BART works by fitting a sum-of-trees model for the treatment or
#' probability of treatment. The number of trees is determined by the `num_trees`
#' argument. Bayesian priors are used for the hyperparameters, so the result is
#' a posterior distribution of predicted values for each unit. The mean of these
#' for each unit is taken for use in computing the (generalized) propensity
#' score. Although the hyperparameters governing the priors can be modified by
#' supplying arguments to [weightit()] that are passed to the BART fitting
#' function, the default values tend to work well and require little
#' modification. Unlike many other machine learning methods, no loss function is
#' optimized and the hyperparameters do not need to be tuned (e.g., using
#' cross-validation), though performance can benefit from tuning. BART tends to
#' balance sparseness with flexibility by using very weak learners as the trees,
#' which makes it suitable for capturing complex functions without specifying a
#' particular functional form and without overfitting.
#'
#' By default, the decision rules in each tree are soft, as in the SoftBart
#' model of Linero and Yang (2018): an observation reaches every leaf with some
#' weight, so the fitted propensity score is a smooth function of the covariates
#' rather than a step function. Smoother propensity scores tend to be better
#' calibrated and to yield less extreme weights than those from standard BART.
#' Soft rules are slower than hard rules; the standard BART sampler can be
#' requested with `gate = "hard"`. See \pkgfun{bartisan}{bartisan_control} for
#' details.
#'
#' ## Reproducibility
#'
#' BART has a random component. \pkg{bartisan} uses R's random number generator,
#' so calling [set.seed()] before `weightit()` makes the results reproducible,
#' including when the chains run in parallel, since each chain is given its own
#' random number stream. The `seed` argument accepted by earlier versions of this
#' method, which used \pkg{dbarts}, is ignored with a warning.
#'
#' ## Parallel Processing
#'
#' When `chains` is greater than 1 and the \CRANpkg{future.apply} package is
#' installed, the chains run in parallel under whatever \pkgfun{future}{plan} has
#' been set. See *Parallelization with future* at \pkgfun{bartisan}{bartisan_control}.
#'
#' @seealso [weightit()], [weightitMSM()], [get_w_from_ps()]
#'
#' [`method_super`] for stacking predictions from several machine learning
#' methods, including BART.
#'
#' @references
#' Chipman, H. A., George, E. I., & McCulloch, R. E. (2010). BART: Bayesian
#' additive regression trees. *The Annals of Applied Statistics*, 4(1), 266–298.
#' \doi{10.1214/09-AOAS285}
#'
#' Hill, J., Weiss, C., & Zhai, F. (2011). Challenges With
#' Propensity Score Strategies in a High-Dimensional Setting and a Potential
#' Alternative. *Multivariate Behavioral Research*, 46(3), 477–513.
#' \doi{10.1080/00273171.2011.570161}
#'
#' Linero, A. R. (2025). Generalized Bayesian additive regression trees models:
#' beyond conditional conjugacy. *Journal of the American Statistical
#' Association*, 120(549), 356–369. \doi{10.1080/01621459.2024.2337156}
#'
#' Linero, A. R., & Yang, Y. (2018). Bayesian regression tree ensembles that
#' adapt to smoothness and sparsity. *Journal of the Royal Statistical Society
#' Series B*, 80(5), 1087–1110. \doi{10.1111/rssb.12293}
#'
#' - Continuous treatments
#'
#' George, E., Laud, P., Logan, B., McCulloch, R., & Sparapani, R. (2019). Fully
#' Nonparametric Bayesian Additive Regression Trees. In *Topics in
#' Identification, Limited Dependent Variables, Partial Observability,
#' Experimentation, and Flexible Modeling: Part B* (Advances in Econometrics,
#' Vol. 40B, pp. 89–110). Emerald Publishing Limited.
#' \doi{10.1108/S0731-90532019000040B006}
#'
#' Note that many references that deal with BART for causal inference focus on
#' estimating potential outcomes with BART, not the propensity scores, and so
#' are not directly relevant when using BART to estimate propensity scores for
#' weights.
#'
#' See [`method_glm`] for additional references on propensity score weighting
#' more generally.
#'
#' @examplesIf rlang::is_installed("bartisan")
#' \donttest{data("lalonde", package = "cobalt")
#'
#' set.seed(123)
#'
#' (W1 <- weightit(treat ~ age + educ + married +
#'                   nodegree + re74, data = lalonde,
#'                 method = "bart", estimand = "ATT"))
#'
#' summary(W1)
#'
#' cobalt::bal.tab(W1)
#'
#' (W2 <- weightit(race ~ age + educ + married +
#'                 nodegree + re74, data = lalonde,
#'                 method = "bart", estimand = "ATE"))
#'
#' summary(W2)
#'
#' cobalt::bal.tab(W2)
#'
#' (W3 <- weightit(re75 ~ age + educ + married +
#'                   nodegree + re74, data = lalonde,
#'                 method = "bart"))
#'
#' summary(W3)
#'
#' cobalt::bal.tab(W3)}
NULL

#Assemble the control arguments for bartisan::bartisan() from the arguments
#passed through `...`. Natively named `bartisan_control()` arguments are collected
#as they are; the `dbarts::bart2()`-style names accepted by earlier versions of
#this method are translated to their equivalents, with a native name taking
#precedence. `seed` is dropped with a warning, since bartisan uses R's own random
#number generator; `n.threads` is dropped silently, since parallelism is set with
#a future plan. The returned list excludes `verbose`, which the caller sets.
.make_bartisan_args <- function(...) {
  native <- setdiff(rlang::fn_fmls_names(bartisan::bartisan_control), "verbose")

  A <- ...mget(intersect(...names(), native))

  translations <- c(n.trees = "num_trees", n.samples = "num_draws",
                    n.burn = "num_burn", n.chains = "chains")

  for (i in intersect(...names(), names(translations))) {
    if (is_null(A[[translations[i]]])) {
      A[[translations[i]]] <- ...get(i)
    }
  }

  if (is_not_null(...get("seed"))) {
    arg::wrn("{.arg seed} is ignored with {.code method = \"bart\"}; call {.fn set.seed} before {.fn weightit} to make the results reproducible")
  }

  A
}

#Scale the sampling weights to have a mean of 1 so the fit does not depend on
#their scale.
.scale_s_weights <- function(s.weights) {
  s.weights / mean(s.weights)
}

#Assemble the data frame and formula for bartisan::bartisan(). The covariates are
#renamed `.x1`, `.x2`, ... so that any covariate name is a valid formula term,
#and the treatment, sampling weights, and offset are stored as columns, so that
#the call bartisan() records does not carry the vectors themselves. With random
#effects terms, `.make_re_data_formula()` adds the grouping variables and the
#terms to the formula.
.make_bartisan_data <- function(covs, treat, s.weights, offset, re.bars = NULL,
                                .data = NULL, subset = NULL) {
  if (is_not_null(re.bars)) {
    df <- .make_re_data_formula(covs, treat, s.weights, re.bars, .data,
                                subset, bart = FALSE)
  }
  else {
    fe.string <- "1"

    if (ncol(covs) > 0L) {
      fe.data <- as.data.frame(covs)
      names(fe.data) <- paste0(".x", seq_col(covs))
      data <- cbind(data.frame(.treat = treat, .s.weights = s.weights),
                    fe.data)
      fe.string <- paste(names(fe.data), collapse = " + ")
    }
    else {
      data <- data.frame(.treat = treat, .s.weights = s.weights)
    }

    df <- list(data = data,
               formula = stats::as.formula(sprintf(".treat ~ %s", fe.string)))
  }

  df[["data"]][[".offset"]] <- offset

  df
}

#Fit the bartisan model, referring to the weights and offset by the names of
#their columns in the data so that the call bartisan() stores stays small.
.fit_bartisan <- function(df, family, A, verbose, weights = TRUE) {
  A[["formula"]] <- df[["formula"]]
  A[["data"]] <- df[["data"]]
  A[["family"]] <- family
  A[["offset"]] <- quote(.offset)
  A[["verbose"]] <- isTRUE(verbose)

  if (weights) {
    A[["weights"]] <- quote(.s.weights)
  }

  rlang::try_fetch({
    fit <- eval(as.call(c(list(quote(bartisan::bartisan)), A)))
  },
  error = function(e) {
    arg::err("(from {.fun bartisan::bartisan}): {conditionMessage(e)}")
  })

  fit
}

#Whether a bartisan family accepts prior weights. The Dirichlet process mixture
#families refuse them, since a weight would have to be a multiplicity.
.bartisan_family_takes_weights <- function(family) {
  name <- {
    if (is.character(family)) family
    else family[["family"]]
  }

  !identical(name, "dpm") && !identical(name, "dpm_aft")
}

weightit2bart <- function(covs, treat, s.weights, subset, estimand, focal, stabilize,
                          missing, .data, verbose, ...) {
  #Read before `subset` narrows `treat`; see `.recorded_treated_level()`.
  t.lev <- .recorded_treated_level(treat)

  covs <- covs[subset, , drop = FALSE]
  treat <- treat[subset]
  s.weights <- s.weights[subset]

  missing <- .process_missing2(missing, covs)

  if (missing == "ind") {
    covs <- add_missing_indicators(covs)
  }

  covs <- covs |>
    .make_covs_closer_to_1() |>
    .make_covs_full_rank()

  t.lev <- t.lev %or% .get_treated_level(treat, estimand, focal)
  treat <- binarize(treat, one = t.lev)

  link <- ...get("link", "logit")
  link <- arg::match_arg(link, c("logit", "probit", "cloglog"))

  s.weights <- .scale_s_weights(s.weights)

  use.offset <- ...get("use.offset", FALSE)
  arg::arg_flag(use.offset)

  if (use.offset) {
    rlang::try_fetch({
      fit <- glm.fit(x = as.matrix(cbind(1, covs)), y = treat,
                     weights = s.weights, family = quasibinomial(link))
    },
    warning = function(w) {
      w <- conditionMessage(w)
      if (w != "non-integer #successes in a binomial glm!") {
        arg::wrn("(from {.fun stats::glm.fit}): {w}")
      }
      invokeRestart("muffleWarning")
    })

    offset <- fit$linear.predictors
  }
  else {
    offset <- ...get("offset", rep.int(0, length(treat)))
  }

  df <- .make_bartisan_data(covs, as.numeric(treat), s.weights, offset,
                            re.bars = ...get(".random"), .data = .data,
                            subset = subset)

  A <- .make_bartisan_args(...)

  fit <- .fit_bartisan(df, family = binomial(link), A = A, verbose = verbose)

  p.score <- unname(fitted(fit))

  w <- .get_w_from_ps_internal_bin(ps = p.score, treat = treat, estimand,
                                   stabilize = stabilize, subclass = ...get("subclass"))

  list(w = w, ps = p.score, fit.obj = fit)
}

weightit2bart.cens <- function(covs, treat, s.weights, subset, missing, verbose,
                               estimand = NULL, focal = NULL, stabilize = FALSE, ...) {

  C <- .make_cens_treat(treat)

  out <- .cens_degenerate_out(C[subset])

  if (is_not_null(out)) {
    return(out)
  }

  #The censoring model P(C = 1 | X) is the binary propensity score model with the
  #censored units as the focal group; see `weightit2glm.cens()`.
  out <- weightit2bart(covs = covs, treat = C, s.weights = s.weights,
                       subset = subset, estimand = "ATT", focal = 1,
                       stabilize = FALSE, missing = missing,
                       verbose = verbose, ...)

  .att_out_to_cens(out, C[subset])
}

weightit2bart.multi <- function(covs, treat, s.weights, subset, estimand, focal, stabilize,
                                missing, .data, verbose, ...) {
  covs <- covs[subset, , drop = FALSE]
  treat <- factor(treat[subset])
  s.weights <- s.weights[subset]

  missing <- .process_missing2(missing, covs)

  if (missing == "ind") {
    covs <- add_missing_indicators(covs)
  }

  covs <- covs |>
    .make_covs_closer_to_1() |>
    .make_covs_full_rank()

  if (!isFALSE(...get("use.offset", FALSE))) {
    arg::err("{.arg use.offset} cannot be used with multi-category treatments")
  }

  link <- ...get("link", "logit")
  link <- arg::match_arg(link, c("logit", "probit"))

  family <- bartisan::multinomial(link = link,
                                  reference = ...get("reference"),
                                  replicates = ...get("replicates", 200L))

  s.weights <- .scale_s_weights(s.weights)

  #A plain factor, so that bartisan() reads the response as multi-category
  treat.plain <- factor(as.character(treat), levels = levels(treat))

  df <- .make_bartisan_data(covs, treat.plain, s.weights,
                            offset = rep.int(0, length(treat)),
                            re.bars = ...get(".random"), .data = .data,
                            subset = subset)

  A <- .make_bartisan_args(...)

  fit <- .fit_bartisan(df, family = family, A = A, verbose = verbose)

  #ps is a matrix of probabilities for each treatment level, one column per level
  ps <- predict(fit, type = "prob")[, levels(treat), drop = FALSE] |>
    as.data.frame()

  #Computing weights
  w <- .get_w_from_ps_internal_multi(ps = ps, treat = treat, estimand, focal,
                                     stabilize = stabilize, subclass = ...get("subclass"))

  list(w = w, fit.obj = fit)
}

weightit2bart.cont <- function(covs, treat, s.weights, subset, stabilize, missing, ps, .data, verbose, ...) {

  covs <- covs[subset, , drop = FALSE]
  treat <- treat[subset]
  s.weights <- s.weights[subset]

  missing <- .process_missing2(missing, covs)

  if (missing == "ind") {
    covs <- add_missing_indicators(covs)
  }

  covs <- covs |>
    .make_covs_closer_to_1()

  #The conditional density comes from the model, so the density arguments used
  #by the other methods have no role here
  dens.args <- intersect(...names(), c("density", "bw", "adjust", "kernel", "n", "use.kernel"))
  if (is_not_null(dens.args)) {
    arg::wrn("{.arg {dens.args}} {?is/are} ignored with {.code method = \"bart\"} for continuous treatments; the conditional density of the treatment is estimated by the model specified by {.arg family}")
  }

  family <- ...get("family") %or% bartisan::dpm()

  if (!is.character(family) && !inherits(family, "family")) {
    arg::err("{.arg family} must be a family object or the name of one; see {.help bartisan::bartisan-families}")
  }

  s.weights <- .scale_s_weights(s.weights)

  use.offset <- ...get("use.offset", FALSE)
  arg::arg_flag(use.offset)

  if (use.offset) {
    rlang::try_fetch({
      fit <- glm.fit(x = as.matrix(cbind(1, covs)), y = treat,
                     weights = s.weights, family = gaussian())
    },
    warning = function(w) {
      invokeRestart("muffleWarning")
    })

    offset <- fit$linear.predictors
  }
  else {
    offset <- ...get("offset", rep.int(0, length(treat)))
  }

  df <- .make_bartisan_data(covs, as.numeric(treat), s.weights, offset,
                            re.bars = ...get(".random"), .data = .data,
                            subset = subset)

  A <- .make_bartisan_args(...)

  #Estimate the conditional density of the treatment
  fit <- .fit_bartisan(df, family = family, A = A, verbose = verbose,
                       weights = .bartisan_family_takes_weights(family))

  #Generalized propensity score: the conditional density at each unit's own
  #treatment value. Computed with `newdata` so that it is evaluated exactly as the
  #marginal density below is, with a prior weight of 1 for every unit.
  dens_denom <- .bartisan_density(fit, df[["data"]])

  #Marginal density: the conditional density at each of a grid of treatment
  #values, averaged over the units (with the sampling weights) and interpolated
  #to the observed treatment values, as for the other methods
  n_grid <- 50L
  t_grid <- seq(min(treat), max(treat), length.out = n_grid)

  dens_grid <- .bartisan_grid_density(fit, df[["data"]], t_grid, s.weights)

  dens_num <- spline(x = t_grid, y = dens_grid, xout = treat)[["y"]]

  if (any(!is.finite(dens_denom)) || any(dens_denom <= 0) ||
      any(!is.finite(dens_num)) || any(dens_num <= 0)) {
    arg::err("the estimated density of the treatment is not positive and finite for all units. Try another {.arg family} or more draws")
  }

  w <- dens_num / dens_denom

  list(w = w, fit.obj = fit)
}

#The marginal density of the treatment at each value in `t_grid`: the conditional
#density at that value, evaluated for every unit's covariates and averaged over
#the units with the sampling weights. Evaluating the density costs draws by rows
#by trees, so two approximations keep it affordable: the posterior draws are
#thinned to at most `max_draws`, since the average over units and draws is smooth,
#and at most `max_units` units (a random subset) are averaged over, since the
#marginal density is a stabilization factor rather than something that has to be
#exact. With the defaults, the average at each grid value is over up to 50,000
#density evaluations, and the cost no longer grows with the sample size beyond
#`max_units`. Each predict() call computes a matrix of draws by rows, so the grid is
#split into blocks whose stacked copies of the data keep that matrix small.
.bartisan_grid_density <- function(fit, data, t_grid, s.weights, max_draws = 25L,
                                   max_units = 2000L, max_rows = 1e5) {
  n <- nrow(data)

  if (n > max_units) {
    units <- sort(sample.int(n, max_units))
    data <- data[units, , drop = FALSE]
    s.weights <- s.weights[units]
    n <- max_units
  }

  num_draws <- length(fit[["loglik"]])
  iterations <- {
    if (num_draws > max_draws) unique(round(seq(1, num_draws, length.out = max_draws)))
    else NULL
  }

  per_call <- max(1L, min(length(t_grid), floor(max_rows / n)))
  blocks <- split(t_grid, ceiling(seq_along(t_grid) / per_call))

  out <- lapply(blocks, function(block) {
    nd <- data[rep(seq_len(n), times = length(block)), , drop = FALSE]
    nd[[".treat"]] <- rep(block, each = n)
    d <- matrix(.bartisan_density(fit, nd, iterations = iterations), nrow = n)
    apply(d, 2L, w.m, w = s.weights)
  })

  unlist(out, use.names = FALSE)
}

#The conditional density of the treatment at the values in `newdata`, averaged
#over the posterior draws. The offset column is passed back, since a model fit
#with an offset requires it at prediction time.
.bartisan_density <- function(fit, newdata, iterations = NULL) {
  rlang::try_fetch({
    d <- predict(fit, newdata = newdata, type = "density",
                 offset = newdata[[".offset"]], iterations = iterations)
  },
  error = function(e) {
    arg::err("(from {.fun bartisan::predict.bartisan_fit}): {conditionMessage(e)}")
  })

  unname(d)
}

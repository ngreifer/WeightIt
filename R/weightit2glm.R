#' Propensity Score Weighting Using Generalized Linear Models
#' @name method_glm
#' @aliases method_glm
#' @usage NULL
#'
#' @description
#' This page explains the details of estimating weights from generalized linear model-based propensity scores by setting `method = "glm"` in the call to [weightit()] or [weightitMSM()]. This method can be used with binary, multi-category, and continuous treatments.
#'
#' In general, this method relies on estimating propensity scores with a parametric generalized linear model and then converting those propensity scores into weights using a formula that depends on the desired estimand. For binary and multi-category treatments, a binomial or multinomial regression model is used to estimate the propensity scores as the predicted probability of being in each treatment given the covariates. For ordinal treatments, an ordinal regression model is used to estimate generalized propensity scores. For continuous treatments, a generalized linear model is used to estimate generalized propensity scores as the conditional density of treatment given the covariates.
#'
#' ## Binary Treatments
#'
#' For binary treatments, this method estimates the propensity scores using [glm()]. An additional argument is `link`, which uses the same options as `link` in [family()]. The default link is "logit", but others, including "probit", are allowed. The following estimands are allowed: ATE, ATT, ATC, ATO, ATM, and ATOS. Weights can also be computed using marginal mean weighting through stratification for the ATE, ATT, and ATC. See [get_w_from_ps()] for details.
#'
#' ## Multi-Category Treatments
#'
#' For multi-category treatments, the propensity scores are estimated using multinomial regression from one of a few functions depending on the argument supplied to `multi.method` (see Additional Arguments below). The following estimands are allowed: ATE, ATT, ATC, ATO, and ATM. The weights for each estimand are computed using the standard formulas or those mentioned above. Weights can also be computed using marginal mean weighting through stratification for the ATE, ATT, and ATC. See [get_w_from_ps()] for details. Ordinal treatments are treated exactly the same as non-order multi-category treatments except that additional models are available to estimate the generalized propensity score (e.g., ordinal logistic regression).
#'
#' ## Continuous Treatments
#'
#' For continuous treatments, the generalized propensity score is estimated using linear regression. The conditional density can be specified as normal or another distribution. In addition, kernel density estimation can be used instead of assuming a specific density for the numerator and denominator of the generalized propensity score by setting `density = "kernel"`. Other arguments to [density()] can be specified to refine the density estimation parameters. `plot = TRUE` can be specified to plot the density for the numerator and denominator, which can be helpful in diagnosing extreme weights.
#'
#' ## Longitudinal Treatments
#'
#' For longitudinal treatments, the weights are the product of the weights estimated at each time point.
#'
#' ## Sampling Weights
#'
#' Sampling weights are supported through `s.weights` in all scenarios except for multi-category treatments with `link = "bayes.probit"` and for binary and continuous treatments with `missing = "saem"` (see below). Warning messages may appear otherwise about non-integer successes, and these can be ignored.
#'
#' ## Missing Data
#'
#' In the presence of missing data, the following value(s) for `missing` are allowed:
#'     \describe{
#'       \item{`"ind"` (default)}{First, for each variable with missingness, a new missingness indicator variable is created which takes the value 1 if the original covariate is `NA` and 0 otherwise. The missingness indicators are added to the model formula as main effects. The missing values in the covariates are then replaced with the covariate medians (this value is arbitrary and does not affect estimation). The weight estimation then proceeds with this new formula and set of covariates. The covariates output in the resulting `weightit` object will be the original covariates with the `NA`s.
#'       }
#'       \item{`"saem"`}{For binary treatments with `link = "logit"` or continuous treatments, a stochastic approximation version of the EM algorithm (SAEM) is used via the \CRANpkg{misaem} package. No additional covariates are created. See Jiang et al. (2019) for information on this method. In some cases, this is a suitable alternative to multiple imputation.
#'       }
#'    }
#'
#' ## M-estimation
#'
#' For binary treatments, M-estimation is supported when `link` is neither `"flic"` nor `"flac"` (see below). For multi-category treatments, M-estimation is supported when `multi.method` is `"weightit"` (the default for non-ordered treatments) or `"glm"`. For continuous treatments, M-estimation is supported when `density` is not `"kernel"`. The conditional treatment variance and unconditional treatment mean and variance are included as parameters to estimate, as these all go into calculation of the weights. For all treatment type, M-estimation is not supported when `missing = "saem"`. See [glm_weightit()] and `vignette("estimating-effects")` for details. For longitudinal treatments, M-estimation is supported whenever the underlying methods are.
#'
#' @section Additional Arguments:
#' For binary treatments, the following additional argument can be specified:
#' \describe{
#'   \item{`link`}{the link used in the generalized linear model for the propensity scores. `link` can be any of those allowed by [binomial()]. A `br.` prefix can be added (e.g., `"br.logit"`); this changes the fitting method to the bias-corrected generalized linear models implemented in the \CRANpkg{brglm2} package. `link` can also be either `"flic"` or `"flac"` to fit the corresponding Firth corrected logistic regression models implemented in the \CRANpkg{logistf} package.}
#' }
#'
#' For multi-category treatments, the following additional arguments can be specified:
#' \describe{
#'   \item{`multi.method`}{the method used to estimate the generalized propensity scores. Allowable options include `"weightit"` to use an M-estimation-based method of multinomial logistic regression implemented in \pkg{WeightIt}, `"glm"` to use a series of binomial models using [glm()], `"mclogit"` to use multinomial logistic regression as implemented in \pkgfun{mclogit}{mblogit}, `"mnp"` to use Bayesian multinomial probit regression as implemented in \pkgfun{MNP}{MNP}, and `"brmultinom"` to use bias-reduced multinomial logistic regression as implemented in \pkgfun{brglm2}{brmultinom}. For ordered treatments, `"polr"` can be supplied to use ordinal regression as implemented in \pkgfun{MASS}{polr} unless `link` is `"br.logit"`, in which case bias-reduce ordinal logistic regression as implemented in \pkgfun{brglm2}{bracl} is used. `"weightit"` and `"mclogit"` should give near-identical results, the main difference being increased robustness and customizability when using `"mclogit"` at the expense of not being able to use M-estimation to compute standard errors after weighting. The default is `"weightit"` for un-ordered treatments and `"polr"` for ordered treatments. Ignored when `missing = "saem"`.}
#'   \item{`link`}{The link used in the multinomial, binomial, or ordered regression model for the generalized propensity scores depending on the argument supplied to `multi.method`. When `multi.method = "glm"`, `link` can be any of those allowed by [binomial()]. When treatment is ordered and `multi.method = "polr"`, `link` can be any of those allowed by `MASS::polr()` or `"br.logit"`. Otherwise, `link` should be `"logit"` or not specified.}
#' }
#'
#' For continuous treatments, the following additional arguments may be supplied:
#'   \describe{
#'     \item{`density`}{A function corresponding the conditional density of the treatment. The standardized residuals of the treatment model will be fed through this function to produce the numerator and denominator of the generalized propensity score weights. If blank, [dnorm()] is used as recommended by Robins et al. (2000). This can also be supplied as a string containing the name of the function to be called. If the string contains underscores, the call will be split by the underscores and the latter splits will be supplied as arguments to the second argument and beyond. For example, if `density = "dt_2"` is specified, the density used will be that of a t-distribution with 2 degrees of freedom. Using a t-distribution can be useful when extreme outcome values are observed (Naimi et al., 2014).
#'
#' Can also be `"kernel"` to use kernel density estimation, which calls [density()] to estimate the numerator and denominator densities for the weights. (This used to be requested by setting `use.kernel = TRUE`, which is now deprecated.)}
#'     \item{`bw`, `adjust`, `kernel`, `n`}{If `density = "kernel"`, the arguments to [density()]. The defaults are the same as those in `density()` except that `n` is 10 times the number of units in the sample.}
#'     \item{`plot`}{If `density = "kernel"`, whether to plot the estimated densities.}
#'     \item{`link`}{The link used to fit the linear model for the generalized propensity score. Can be any allowed by [gaussian()].
#'     }
#'   }
#'
#' Additional arguments to `glm()` can be specified as well when it is used for fitting. The `method` argument in `glm()` is renamed to `glm.method`. This can be used to supply alternative fitting functions, such as those implemented in the \CRANpkg{glm2} package. Other arguments to `weightit()` are passed to `...` in `glm()`. In the presence of missing data with `link = "logit"` and `missing = "saem"`, additional arguments are passed to \pkgfun{misaem}{miss.glm} and \pkgfun{misaem}{predict.miss.glm}, except the `method` argument in \pkgfun{misaem}{predict.miss.glm} is replaced with `saem.method`.
#'
#' For continuous treatments in the presence of missing data with `missing = "saem"`, additional arguments are passed to \pkgfun{misaem}{miss.lm} and \pkgfun{misaem}{predict.miss.lm}.
#'
#' @section Additional Outputs:
#' \describe{
#'   \item{`obj`}{When `include.obj = TRUE`, the (generalized) propensity score model fit. For binary treatments, the output of the call to [glm()] or the requested fitting function. For multi-category treatments, the output of the call to the fitting function (or a list thereof if `multi.method = "glm"`). For continuous treatments, the output of the call to `glm()` for the predicted values in the denominator density.
#'   }
#' }
#'
#' @details NULL
#'
#' @seealso
#' [weightit()], [weightitMSM()], [get_w_from_ps()]
#'
#' @references
#' ## Binary treatments
#'
#' - `estimand = "ATO"`
#'
#' Li, F., Morgan, K. L., & Zaslavsky, A. M. (2018). Balancing covariates via propensity score weighting. *Journal of the American Statistical Association*, 113(521), 390–400. \doi{10.1080/01621459.2016.1260466}
#'
#' - `estimand = "ATM"`
#'
#' Li, L., & Greene, T. (2013). A Weighting Analogue to Pair Matching in Propensity Score Analysis. *The International Journal of Biostatistics*, 9(2). \doi{10.1515/ijb-2012-0030}
#'
#' - `estimand = "ATOS"`
#'
#' Crump, R. K., Hotz, V. J., Imbens, G. W., & Mitnik, O. A. (2009). Dealing with limited overlap in estimation of average treatment effects. *Biometrika*, 96(1), 187–199. \doi{10.1093/biomet/asn055}
#'
#' - Other estimands
#'
#' Austin, P. C. (2011). An Introduction to Propensity Score Methods for Reducing the Effects of Confounding in Observational Studies. *Multivariate Behavioral Research*, 46(3), 399–424. \doi{10.1080/00273171.2011.568786}
#'
#' - Marginal mean weighting through stratification
#'
#' Hong, G. (2010). Marginal mean weighting through stratification: Adjustment for selection bias in multilevel data. *Journal of Educational and Behavioral Statistics*, 35(5), 499–531. \doi{10.3102/1076998609359785}
#'
#' - Bias-reduced logistic regression
#'
#' See references for the \pkg{brglm2} \pkgfun2{brglm2}{brglm2}{package}.
#'
#' - Firth corrected logistic regression
#'
#' Puhr, R., Heinze, G., Nold, M., Lusa, L., & Geroldinger, A. (2017). Firth’s logistic regression with rare events: Accurate effect estimates and predictions? *Statistics in Medicine*, 36(14), 2302–2317. \doi{10.1002/sim.7273}
#'
#' - SAEM logistic regression for missing data
#'
#' Jiang, W., Josse, J., & Lavielle, M. (2019). Logistic regression with missing covariates — Parameter estimation, model selection and prediction within a joint-modeling framework. *Computational Statistics & Data Analysis*, 106907. \doi{10.1016/j.csda.2019.106907}
#'
#' ## Multi-Category Treatments
#'
#' - `estimand = "ATO"`
#'
#' Li, F., & Li, F. (2019). Propensity score weighting for causal inference with multiple treatments. *The Annals of Applied Statistics*, 13(4), 2389–2415. \doi{10.1214/19-AOAS1282}
#'
#' - `estimand = "ATM"`
#'
#' Yoshida, K., Hernández-Díaz, S., Solomon, D. H., Jackson, J. W., Gagne, J. J., Glynn, R. J., & Franklin, J. M. (2017). Matching weights to simultaneously compare three treatment groups: Comparison to three-way matching. *Epidemiology* (Cambridge, Mass.), 28(3), 387–395. \doi{10.1097/EDE.0000000000000627}
#'
#' - Other estimands
#'
#' McCaffrey, D. F., Griffin, B. A., Almirall, D., Slaughter, M. E., Ramchand, R., & Burgette, L. F. (2013). A Tutorial on Propensity Score Estimation for Multiple Treatments Using Generalized Boosted Models. *Statistics in Medicine*, 32(19), 3388–3414. \doi{10.1002/sim.5753}
#'
#' - Marginal mean weighting through stratification
#'
#' Hong, G. (2012). Marginal mean weighting through stratification: A generalized method for evaluating multivalued and multiple treatments with nonexperimental data. *Psychological Methods*, 17(1), 44–60. \doi{10.1037/a0024918}
#'
#' ## Continuous treatments
#'
#' Robins, J. M., Hernán, M. Á., & Brumback, B. (2000). Marginal Structural Models and Causal Inference in Epidemiology. *Epidemiology*, 11(5), 550–560.
#'
#' - Using non-normal conditional densities
#'
#' Naimi, A. I., Moodie, E. E. M., Auger, N., & Kaufman, J. S. (2014). Constructing Inverse Probability Weights for Continuous Exposures: A Comparison of Methods. *Epidemiology*, 25(2), 292–299. \doi{10.1097/EDE.0000000000000053}
#'
#' - SAEM linear regression for missing data
#'
#' Jiang, W., Josse, J., & Lavielle, M. (2019). Logistic regression with missing covariates — Parameter estimation, model selection and prediction within a joint-modeling framework. *Computational Statistics & Data Analysis*, 106907. \doi{10.1016/j.csda.2019.106907}
#'
#' @examples
#' library("cobalt")
#' data("lalonde", package = "cobalt")
#'
#' #Balancing covariates between treatment groups (binary)
#' (W1 <- weightit(treat ~ age + educ + married +
#'                   nodegree + re74, data = lalonde,
#'                 method = "glm", estimand = "ATT",
#'                 link = "probit"))
#' summary(W1)
#' bal.tab(W1)
#'
#' #Balancing covariates with respect to race (multi-category)
#' (W2 <- weightit(race ~ age + educ + married +
#'                   nodegree + re74, data = lalonde,
#'                 method = "glm", estimand = "ATE"))
#' summary(W2)
#' bal.tab(W2)
#'
#' #Balancing covariates with respect to re75 (continuous)
#' #with kernel density estimate
#' (W3 <- weightit(re75 ~ age + educ + married +
#'                   nodegree + re74, data = lalonde,
#'                 method = "glm", density = "kernel"))
#' summary(W3)
#' bal.tab(W3)
NULL

weightit2glm <- function(covs, treat, s.weights, subset, estimand, focal,
                         stabilize, subclass, missing, .data, verbose, ...) {
  A <- list(...)

  fit.obj <- NULL

  covs <- covs[subset, , drop = FALSE]
  treat <- treat[subset]
  s.weights <- s.weights[subset]

  if (missing == "ind") {
    covs <- add_missing_indicators(covs)
  }

  for (i in seq_col(covs)) covs[,i] <- .make_closer_to_1(covs[,i])

  if (ncol(covs) > 1) {
    if (missing == "saem") {
      covs0 <- covs
      for (i in colnames(covs)[anyNA_col(covs)]) covs0[is.na(covs0[,i]),i] <- covs0[!is.na(covs0[,i]),i][1]
      colinear.covs.to.remove <- colnames(covs)[colnames(covs) %nin% colnames(make_full_rank(covs0))]
    }
    else colinear.covs.to.remove <- colnames(covs)[colnames(covs) %nin% colnames(make_full_rank(covs))]
    covs <- covs[, colnames(covs) %nin% colinear.covs.to.remove, drop = FALSE]
  }

  t.lev <- get_treated_level(treat)
  treat <- binarize(treat, one = t.lev)

  #Process link
  acceptable.links <- {
    if (missing == "saem") "logit"
    else c(expand_grid_string(c("", "br."), c("logit", "probit", "cloglog", "identity", "log", "cauchit")),
           "flic", "flac")
  }

  if (is_null(A[["link"]])) {
    A[["link"]] <- acceptable.links[1]
  }
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
  if (use.br) A[["link"]] <- substr(A[["link"]], 4, nchar(A[["link"]]))
  use.logistf <- A[["link"]] %in% c("flic", "flac")

  if (missing == "saem") {
    rlang::check_installed("misaem")
    if (is_null(A[["saem.method"]])) A[["saem.method"]] <- "map"
    if (is_null(A[["control"]])) A[["control"]] <- list()
    if (is_null(A[["control"]][["var_cal"]])) A[["control"]][["var_cal"]] <- FALSE
    if (is_null(A[["control"]][["ll_obs_cal"]])) A[["control"]][["ll_obs_cal"]] <- FALSE

    data <- data.frame(treat, covs)

    withCallingHandlers({
      verbosely({
        fit <- misaem::miss.glm(formula(data), data = data, control = as.list(A[["control"]]))
      }, verbose = verbose)
    },
    warning = function(w) {
      w <- conditionMessage(w)
      if (w != "one argument not used by format '%i '") .wrn("(from `misaem::miss.glm()`) ", w, tidy = FALSE)
      invokeRestart("muffleWarning")
    })

    p.score <- drop(predict(fit, newdata = covs, method = A[["saem.method"]]))
  }
  else if (use.logistf) {
    rlang::check_installed("logistf")
    fit_fun <- switch(A[["link"]],
                      "flic" = logistf::flic,
                      "flac" = logistf::flac)

    ctrl_fun <- logistf::logistf.control
    control <- do.call(ctrl_fun, c(A[["control"]],
                                   A[setdiff(names(formals(ctrl_fun))[pmatch(names(A), names(formals(ctrl_fun)), 0)],
                                             names(A[["control"]]))]))

    modctrl_fun <- logistf::logistf.mod.control
    modcontrol <- do.call(modctrl_fun, c(A[["modcontrol"]],
                                         A[setdiff(names(formals(modctrl_fun))[pmatch(names(A), names(formals(modctrl_fun)), 0)],
                                                   names(A[["modcontrol"]]))]))

    withCallingHandlers({verbosely({
      data <- data.frame(treat, covs)
      formula <- if (ncol(covs) > 0) formula(data) else treat ~ 1

      fit <- do.call(fit_fun, list(formula, data = data,
                                   weights = s.weights,
                                   control = control,
                                   modcontrol = modcontrol,
                                   pl = FALSE), quote = TRUE)
    }, verbose = verbose)},
    warning = function(w) {
      w <- conditionMessage(w)
      if (w != "non-integer #successes in a binomial glm!")
        .wrn(sprintf("(from `logistf::%s()`) ", A[["link"]]),
             w, tidy = FALSE)
      invokeRestart("muffleWarning")
    })

    p.score <- fit$predict
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
      start <- c(family$linkfun(w.m(treat, s.weights)), rep.int(0, ncol(covs)))
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
        formula <- if (ncol(covs) > 0) formula(data) else treat ~ 1

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
      w <- conditionMessage(w)
      if (w != "non-integer #successes in a binomial glm!")
        .wrn("(from `glm()`) ", w, tidy = FALSE)
      invokeRestart("muffleWarning")
    })

    p.score <- fit$fitted.values

    if (any(p.score <= 1e-14) || any(p.score >= 1 - 1e-14)) {
      .wrn('propensity scores numerically equal to 0 or 1 were estimated, indicating perfect separation and infinite parameter estimates. These may yield problems with inference. Consider trying a different `link`. See `help("method_glm", package = "WeightIt")` for details')
    }

    .psi <- .get_glm_psi(fit)
  }

  fit[["call"]] <- NULL
  fit.obj <- fit

  #Computing weights
  w <- .get_w_from_ps_internal_bin(ps = p.score, treat = treat, estimand,
                                   stabilize = stabilize, subclass = subclass)

  Mparts <- NULL
  if (missing != "saem" && !use.logistf &&
      !(use.br && identical(fit$type, "correction"))) {
    Mparts <- list(
      psi_treat = function(Btreat, A, Xtreat, SW) {
        .psi(B = Btreat, X = Xtreat, y = A, weights = SW)
      },
      wfun = function(Btreat, Xtreat, A) {
        XB <- drop(Xtreat %*% Btreat)
        p <- family$linkinv(XB)
        .get_w_from_ps_internal_bin(ps = p, treat = A, estimand,
                                    stabilize = stabilize, subclass = subclass)
      },
      Xtreat = cbind(`(Intercept)` = 1, covs),
      A = treat,
      btreat = fit$coefficients
    )
  }

  list(w = w, ps = p.score, fit.obj = fit.obj,
       Mparts = Mparts)
}

weightit2glm.multi <- function(covs, treat, s.weights, subset, estimand, focal,
                               stabilize, subclass, missing, .data, verbose, ...) {
  A <- list(...)

  fit.obj <- NULL

  covs <- covs[subset, , drop = FALSE]
  treat <- factor(treat[subset])
  s.weights <- s.weights[subset]

  ord.treat <- is.ordered(treat)

  if (missing == "ind") {
    covs <- add_missing_indicators(covs)
  }

  for (i in seq_col(covs)) covs[,i] <- .make_closer_to_1(covs[,i])

  if (ncol(covs) > 1) {
    if (missing == "saem") {
      covs0 <- covs
      for (i in colnames(covs)[anyNA_col(covs)]) covs0[is.na(covs0[,i]),i] <- covs0[!is.na(covs0[,i]),i][1]
      colinear.covs.to.remove <- colnames(covs)[colnames(covs) %nin% colnames(make_full_rank(covs0))]
    }
    else colinear.covs.to.remove <- colnames(covs)[colnames(covs) %nin% colnames(make_full_rank(covs))]
    covs <- covs[, colnames(covs) %nin% colinear.covs.to.remove, drop = FALSE]
  }

  if (is_not_null(A$use.mlogit)) {
    .wrn("`use.mlogit` is no longer accepted and will be ignored; use `multi.method` instead. See `help(\"method_glm\")` for details")
    A$use.mlogit <- NULL
  }
  if (is_not_null(A$use.mclogit)) {
    .wrn("`use.mclogit` is no longer accepted and will be ignored; use `multi.method` instead. See `help(\"method_glm\")` for details")
    A$use.mclogit <- NULL
  }

  if (missing == "saem") {
    if (is_not_null(A$multi.method) && !identical(A$multi.method, "saem")) {
      .wrn("`multi.method` is ignored when `missing = \"saem\"`")
    }
    multi.method <- "saem"
    ord.treat <- FALSE
  }
  else if (ord.treat) {
    if (is_null(A$multi.method)) {
      multi.method <- "weightit"
    }
    else {
      multi.method <- A$multi.method
      chk::chk_string(multi.method)
      multi.method <- tolower(multi.method)
      if (multi.method == "mblogit") multi.method <- "mclogit"
      allowable.multi.methods <- c("polr", "weightit", "glm", "mclogit", "mnp", "brmultinom")
      multi.method <- match_arg(multi.method, allowable.multi.methods)
      if (!multi.method %in% c("weightit", "polr")) {
        ord.treat <- FALSE
      }
    }
  }
  else if (is_not_null(A$multi.method)) {
    multi.method <- A$multi.method
    chk::chk_string(multi.method)
    multi.method <- tolower(multi.method)
    if (multi.method == "polr") {
      .err("`multi.method = \"polr\"` can only be used with an ordered treatment")
    }
    if (multi.method == "mblogit") multi.method <- "mclogit"
    allowable.multi.methods <- c("weightit", "glm", "mclogit", "mnp", "brmultinom")
    multi.method <- match_arg(multi.method, allowable.multi.methods)
  }
  else {
    multi.method <- NULL
  }

  # Process link
  link <- NULL
  use.br <- FALSE
  if (is_null(multi.method)) {
    if (is_null(A[["link"]]) || identical(A[["link"]], "logit")) {
      multi.method <- "weightit"
    }
    else if (identical(A[["link"]], "bayes.probit")) {
      multi.method <- "mnp"
    }
    else if (identical(A[["link"]], "br.logit")) {
      multi.method <- "brmultinom"
    }
    else {
      .wrn("`link` is ignored when `multi.method` is not specified")
    }
  }
  else if (ord.treat) {
    acceptable.links <- {
      if (multi.method == "polr") c("logit", "probit", "loglog", "cloglog", "cauchit", "br.logit")
      else c("logit", "probit", "cloglog", "identity", "log", "cauchit")
    }

    if (is_null(A[["link"]])) {
      link <- acceptable.links[1]
    }
    else {
      link <- acceptable.links[pmatch(A[["link"]], acceptable.links, nomatch = 0)][1]
      if (is.na(link)) {
        .err(sprintf("only %s allowed as the link for ordinal multi-category treatments with `multi.method = \"%s\"",
                     word_list(acceptable.links, quotes = TRUE, is.are = TRUE),
                     multi.method))
      }
    }

    use.br <- startsWith(link, "br.")
    if (use.br) link <- substr(link, 4, nchar(link))
  }
  else if (multi.method == "saem") {
    if (is_not_null(A[["link"]]) && !identical(A[["link"]], "logit")) {
      .wrn("`link` is ignored when `missing = \"saem\"`")
    }
  }
  else if (multi.method == "weightit") {
    acceptable.links <- c("logit")
    if (is_null(A[["link"]])) {
      link <- acceptable.links[1]
    }
    else {
      link <- acceptable.links[pmatch(A[["link"]], acceptable.links, nomatch = 0)][1]
      if (is.na(link)) {
        .err(sprintf("only %s allowed as the link for non-ordinal multi-category treatments with `multi.method = \"weightit\"`",
                     word_list(acceptable.links, quotes = TRUE, is.are = TRUE)))
      }
    }
  }
  else if (multi.method %in% c("mclogit", "mnp", "brmultinom")) {
    if (is_not_null(A[["link"]]) && !identical(A[["link"]], "logit")) {
      .wrn(sprintf("`link` is ignored when `multi.method` is %s",
                   add_quotes(multi.method)))
    }
  }
  else if (multi.method == "glm") {
    acceptable.links <- c("logit", "probit", "cloglog", "identity", "log", "cauchit")
    if (is_null(A[["link"]])) {
      link <- acceptable.links[1]
    }
    else {
      link <- acceptable.links[pmatch(A[["link"]], acceptable.links, nomatch = 0)][1]
      if (is.na(link)) {
        .err(sprintf("only %s allowed as the link for multi-category treatments with `multi.method = \"glm\"`",
                     word_list(acceptable.links, quotes = TRUE, is.are = TRUE)))
      }
    }
  }

  # Fit model
  if (multi.method == "weightit") {
    fit.fun <- if (ord.treat) ".ordinal_weightit.fit" else ".multinom_weightit.fit"

    verbosely({
      fit.obj <- do.call(fit.fun,
                         list(x = cbind(`(Intercept)` = 1, covs),
                              y = treat,
                              weights = s.weights,
                              hess = FALSE,
                              link = link),
                         quote = TRUE)
    }, verbose = verbose)

    ps <- fit.obj$fitted.values
  }
  else if (multi.method == "polr") {
    if (use.br) {
      rlang::check_installed("brglm2")

      ctrl_fun <- brglm2::brglmControl
      control <- do.call(ctrl_fun, c(A[["control"]],
                                     A[setdiff(names(formals(ctrl_fun))[pmatch(names(A), names(formals(ctrl_fun)), 0)],
                                               names(A[["control"]]))]))

      data <- data.frame(treat, covs)
      formula <- formula(data)

      tryCatch({verbosely({
        fit <- do.call(brglm2::bracl,
                       list(formula,
                            data = data,
                            weights = s.weights,
                            control = control,
                            parallel = if_null_then(A[["parallel"]], FALSE)),
                       quote = TRUE)
      }, verbose = verbose)},
      error = function(e) {
        .err(sprintf("there was a problem with the bias-reduced ordinal logit regression.\n       Try a different link.\n       Error message: (from `brglm2::bracl()`) %s", conditionMessage(e)), tidy = FALSE)
      })
    }
    else {
      rlang::check_installed("MASS")
      if (link == "logit") link <- "logistic"

      data <- data.frame(treat, covs)
      formula <- formula(data)

      tryCatch({verbosely({
        fit <- do.call(MASS::polr,
                       list(formula,
                            data = data,
                            weights = s.weights,
                            Hess = FALSE,
                            model = FALSE,
                            method = link,
                            contrasts = NULL), quote = TRUE)
      }, verbose = verbose)},
      error = function(e) {
        .err(sprintf("There was a problem fitting the ordinal %s regressions with `polr()`.\n       Try again with an un-ordered treatment.\n       Error message: (from `MASS::polr()`) %s",
                     link, conditionMessage(e)), tidy = FALSE)})
    }

    ps <- fit$fitted.values
    fit.obj <- fit
  }
  else if (multi.method == "saem") {
    rlang::check_installed("misaem")
    if (is_null(A[["saem.method"]])) A[["saem.method"]] <- "map"
    if (is_null(A[["control"]])) A[["control"]] <- list()
    if (is_null(A[["control"]][["var_cal"]])) A[["control"]][["var_cal"]] <- FALSE
    if (is_null(A[["control"]][["ll_obs_cal"]])) A[["control"]][["ll_obs_cal"]] <- FALSE

    ps <- make_df(levels(treat), nrow = length(treat))

    fit.list <- make_list(levels(treat))

    for (i in levels(treat)) {
      t_i <- as.numeric(treat == i)
      data_i <- data.frame(t_i, covs)

      withCallingHandlers({
        verbosely({
          fit.list[[i]] <- misaem::miss.glm(formula(data_i), data = data_i,
                                            control = A[["control"]])
        }, verbose = verbose)
      },
      warning = function(w) {
        w <- conditionMessage(w)
        if (w != "one argument not used by format '%i '") .wrn("(from `misaem::miss.glm()`) ", w, tidy = FALSE)
        invokeRestart("muffleWarning")
      })

      ps[[i]] <- drop(predict(fit, newdata = covs, method = A[["saem.method"]]))
    }
    fit.obj <- fit.list
  }
  else if (multi.method == "mnp") {
    rlang::check_installed("MNP")

    data <- data.frame(treat, covs)
    formula <- formula(data)
    tryCatch({verbosely({
      fit.obj <- MNP::mnp(formula, data, verbose = TRUE)
    }, verbose = verbose)},
    error = function(e) {
      .err(sprintf("There was a problem fitting the multinomial regression with `MNP()`.\n       Try a different `multi.method`.\nError message: (from `MNP::mnp()`) %s",
                   conditionMessage(e)), tidy = FALSE)
    })
    ps <- predict(fit.obj, type = "prob")$p
  }
  else if (multi.method == "mclogit") {
    rlang::check_installed("mclogit")

    if (is_not_null(A$random)) {
      random <- get_covs_and_treat_from_formula(A$random, data = .data)$reported.covs[subset,,drop = FALSE]
      data <- cbind(data.frame(random), data.frame(treat = treat, .s.weights = s.weights, covs))
      covnames <- names(data)[-c(seq_col(random), ncol(random) + (1:2))]
      tname <- names(data)[ncol(random) + 1]
      ctrl_fun <- mclogit::mmclogit.control
    }
    else {
      data <- data.frame(treat = treat, .s.weights = s.weights, covs)
      covnames <- names(data)[-c(1,2)]
      tname <- names(data)[1]
      ctrl_fun <- mclogit::mclogit.control
    }
    form <- reformulate(covnames, tname)

    control <- do.call(ctrl_fun, c(A[["control"]],
                                   A[setdiff(names(formals(ctrl_fun))[pmatch(names(A), names(formals(ctrl_fun)), 0)],
                                             names(A[["control"]]))]))
    tryCatch({verbosely({
      fit.obj <- do.call(mclogit::mblogit,
                         list(form,
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
      .err(sprintf("there was a problem fitting the multinomial %s regression with `mblogit()`.\n       Try a different `multi.method`.\nError message: (from `mclogit::mblogit()`) %s",
                   A[["link"]], conditionMessage(e)), tidy = FALSE)
    }
    )

    ps <- fitted(fit.obj)
    colnames(ps) <- levels(treat)
  }
  else if (multi.method == "brmultinom") {
    rlang::check_installed("brglm2")
    data <- data.frame(treat, covs)
    formula <- formula(data)
    ctrl_fun <- brglm2::brglmControl
    control <- do.call(ctrl_fun, c(A[["control"]],
                                   A[setdiff(names(formals(ctrl_fun))[pmatch(names(A), names(formals(ctrl_fun)), 0)],
                                             names(A[["control"]]))]))
    tryCatch({verbosely({
      fit.obj <- do.call(brglm2::brmultinom,
                         list(formula, data,
                              weights = s.weights,
                              control = control), quote = TRUE)
    }, verbose = verbose)},
    error = function(e) {
      .err(sprintf("There was a problem with the bias-reduced multinomial logit regression. Try a different `multi.method`.\n       Error message: (from `brglm2::brmultinom()`) %s", conditionMessage(e)), tidy = FALSE)
    })

    ps <- fit.obj$fitted.values
  }
  else if (multi.method == "glm") {
    ps <- make_df(levels(treat), nrow = length(treat))

    ctrl_fun <- stats::glm.control
    control <- do.call(ctrl_fun, c(A[["control"]],
                                   A[setdiff(names(formals(ctrl_fun))[pmatch(names(A), names(formals(ctrl_fun)), 0)],
                                             names(A[["control"]]))]))
    family <- quasibinomial(link = link)

    fit.list <- .psi.list <- make_list(levels(treat))

    for (i in levels(treat)) {
      t_i <- as.numeric(treat == i)
      data_i <- data.frame(t_i, covs)

      withCallingHandlers({verbosely({
        fit.list[[i]] <- do.call(stats::glm, list(formula(data_i), data = data_i,
                                                  family = family,
                                                  weights = s.weights,
                                                  control = control), quote = TRUE)
      }, verbose = verbose)},
      warning = function(w) {
        w <- conditionMessage(w)
        if (w != "non-integer #successes in a binomial glm!") .wrn("(from `glm()`) ", w, tidy = FALSE)
        invokeRestart("muffleWarning")
      })

      ps[[i]] <- fit.list[[i]]$fitted.values
      .psi.list[[i]] <- .get_glm_psi(fit.list[[i]])
    }

    fit.obj <- fit.list
  }

  #ps should be matrix of probs for each treat
  #Computing weights
  w <- .get_w_from_ps_internal_multi(ps = ps, treat = treat, estimand, focal = focal,
                                     stabilize = stabilize, subclass = subclass)

  #Get Mparts
  Mparts <- NULL
  if (multi.method == "weightit") {
    Mparts <- list(
      psi_treat = function(Btreat, A, Xtreat, SW) {
        fit.obj$psi(Btreat, Xtreat, A, SW)
      },
      wfun = function(Btreat, Xtreat, A) {
        ps <- fit.obj$get_p(Btreat, Xtreat)
        w <- .get_w_from_ps_internal_multi(ps = ps, treat = A, estimand, focal = focal,
                                           stabilize = stabilize, subclass = subclass)
      },
      Xtreat = fit.obj$x,
      A = treat,
      btreat = fit.obj$coefficients
    )
  }
  else if (multi.method == "glm") {
    Mparts <- list(
      psi_treat = function(Btreat, A, Xtreat, SW) {
        Btreat <- matrix(Btreat, nrow = ncol(Xtreat))

        do.call("cbind", lapply(seq_along(levels(A)), function(i) {
          .psi.list[[i]](Btreat[,i], Xtreat, A, SW)
        }))
      },
      wfun = function(Btreat, Xtreat, A) {
        ps <- family$linkinv(Xtreat %*% matrix(Btreat, nrow = ncol(Xtreat)))
        colnames(ps) <- levels(A)

        .get_w_from_ps_internal_multi(ps, A, estimand = estimand, focal = focal,
                                      subclass = subclass, stabilize = stabilize)
      },
      Xtreat = cbind(`(Intercept)` = 1, covs),
      A = treat,
      btreat = unlist(grab(fit.list, "coefficients"))
    )
  }

  list(w = w, fit.obj = fit.obj,
       Mparts = Mparts)
}

weightit2glm.cont <- function(covs, treat, s.weights, subset, stabilize, missing, verbose, ...) {
  A <- list(...)

  covs <- covs[subset, , drop = FALSE]
  treat <- treat[subset]
  s.weights <- s.weights[subset]

  if (missing == "ind") {
    covs <- add_missing_indicators(covs)
  }

  for (i in seq_col(covs)) covs[,i] <- .make_closer_to_1(covs[,i])

  if (ncol(covs) > 1) {
    if (missing == "saem") {
      covs0 <- covs
      for (i in colnames(covs)[anyNA_col(covs)]) {
        covs0[is.na(covs0[,i]),i] <- covs0[!is.na(covs0[,i]),i][1]
      }
      colinear.covs.to.remove <- colnames(covs)[colnames(covs) %nin% colnames(make_full_rank(covs0))]
    }
    else {
      colinear.covs.to.remove <- colnames(covs)[colnames(covs) %nin% colnames(make_full_rank(covs))]
    }
    covs <- covs[, colnames(covs) %nin% colinear.covs.to.remove, drop = FALSE]
  }

  s.weights <- s.weights / mean_fast(s.weights)

  #Process density params
  densfun <- .get_dens_fun(use.kernel = isTRUE(A[["use.kernel"]]), bw = A[["bw"]],
                           adjust = A[["adjust"]], kernel = A[["kernel"]],
                           n = A[["n"]], treat = treat, density = A[["density"]],
                           weights = s.weights)

  #Stabilization - get dens.num
  un_p <- mean_fast(s.weights * treat)
  un_s2 <- mean_fast(s.weights * (treat - un_p)^2)

  dens.num <- densfun((treat - un_p) / sqrt(un_s2))

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

    data <- data.frame(treat, covs)
    formula <- formula(data)

    withCallingHandlers({verbosely({
      fit <- misaem::miss.lm(formula, data, control = as.list(A[["control"]]))
    }, verbose = verbose)},
    warning = function(w) {
      w <- conditionMessage(w)
      if (w != "one argument not used by format '%i '") {
        .wrn("(from `misaem::miss.lm()`) ", w, tidy = FALSE)
      }
      invokeRestart("muffleWarning")
    })

    if (is_null(A[["saem.method"]])) A[["saem.method"]] <- "map"

    p <- drop(predict(fit, newdata = covs, method = A[["saem.method"]]))
  }
  else {
    acceptable.links <- c("identity", "log", "inverse")

    link <- acceptable.links[pmatch(A[["link"]], acceptable.links, nomatch = 0)][1]
    if (is.na(link)) {
      .err(sprintf("only %s allowed as the link for continuous treatments",
                   word_list(acceptable.links, quotes = TRUE, is.are = TRUE)))
    }

    family <- gaussian(link = link)

    verbosely({
      if (isTRUE(A[["quick"]])) {
        fit <- do.call(stats::glm.fit, list(y = treat,
                                            x = cbind(`(Intercept)` = 1, covs),
                                            weights = s.weights,
                                            family = family,
                                            control = as.list(A$control)), quote = TRUE)
      }
      else {
        data <- data.frame(treat, covs)
        formula <- if (ncol(covs) > 0) formula(data) else treat ~ 1

        fit <- do.call(stats::glm, list(formula, data = data,
                                        weights = s.weights,
                                        family = family,
                                        control = as.list(A$control)),
                       quote = TRUE)
      }
    }, verbose = verbose)

    p <- fit$fitted.values
  }

  s2 <- mean_fast(s.weights * (treat - p)^2)

  #Get weights
  dens.denom <- densfun((treat - p) / sqrt(s2))

  w <- dens.num / dens.denom

  if (isTRUE(A[["plot"]])) {
    d.n <- attr(dens.num, "density")
    d.d <- attr(dens.denom, "density")
    plot_density(d.n, d.d)
  }

  Mparts <- NULL
  if (missing != "saem" && !identical(A[["density"]], "kernel")) {
    Mparts <- list(
      psi_treat = function(Btreat, A, Xtreat, SW) {
        un_s2 <- exp(Btreat[1])
        un_p <- Btreat[2]

        s2 <- exp(Btreat[3])
        lin_pred <- drop(Xtreat %*% Btreat[-(1:3)])
        p <- family$linkinv(lin_pred)

        SW <- SW / mean_fast(SW)

        cbind(SW * (A - un_p)^2 - un_s2, #unconditional variance
              SW * (A - un_p), #unconditional mean
              SW * (A - p)^2 - s2, #conditional variance
              Xtreat * (SW * family$mu.eta(lin_pred) * (A - p) / family$variance(p))) #conditional mean
      },
      wfun = function(Btreat, Xtreat, A) {
        un_s2 <- exp(Btreat[1])
        un_p <- Btreat[2]
        dens.num <- densfun((A - un_p)/sqrt(un_s2))

        s2 <- exp(Btreat[3])
        lin_pred <- drop(Xtreat %*% Btreat[-(1:3)])
        p <- family$linkinv(lin_pred)
        dens.denom <- densfun((A - p) / sqrt(s2))

        dens.num / dens.denom
      },
      Xtreat = cbind(`(Intercept)` = 1, covs),
      A = treat,
      btreat = c("log(s^2)" = log(un_s2), "E[A]" = un_p, "log(s_r^2)" = log(s2),
                 fit$coefficients)
    )
  }

  list(w = w, fit.obj = fit,
       Mparts = Mparts)
}

#' Propensity Score Weighting Using Generalized Linear Models
#' @name method_glm
#' @aliases method_glm
#' @usage NULL
#'
#' @description This page explains the details of estimating weights from
#' generalized linear model-based propensity scores by setting `method = "glm"`
#' in the call to [weightit()] or [weightitMSM()]. This method can be used with
#' binary, multi-category, and continuous treatments.
#'
#' In general, this method relies on estimating propensity scores with a
#' parametric generalized linear model and then converting those propensity
#' scores into weights using a formula that depends on the desired estimand. For
#' binary and multi-category treatments, a binomial or multinomial regression
#' model is used to estimate the propensity scores as the predicted probability
#' of being in each treatment given the covariates. For ordinal treatments, an
#' ordinal regression model is used to estimate generalized propensity scores.
#' For continuous treatments, a generalized linear model is used to estimate
#' generalized propensity scores as the conditional density of treatment given
#' the covariates.
#'
#' ## Binary Treatments
#'
#' For binary treatments, this method estimates the propensity scores using
#' [glm()]. An additional argument is `link`, which uses the same options as
#' `link` in [family()]. The default link is `"logit"`, but others, including
#' `"probit"`, are allowed. The following estimands are allowed: ATE, ATT, ATC,
#' ATO, ATM, and ATOS. Weights can also be computed using marginal mean
#' weighting through stratification for the ATE, ATT, and ATC. See
#' [get_w_from_ps()] for details.
#'
#' ## Multi-Category Treatments
#'
#' For multi-category treatments, the propensity scores are estimated using
#' multinomial regression from one of a few functions depending on the argument
#' supplied to `multi.method` (see Additional Arguments below). The following
#' estimands are allowed: ATE, ATT, ATC, ATO, and ATM. The weights for each
#' estimand are computed using the standard formulas or those mentioned above.
#' Weights can also be computed using marginal mean weighting through
#' stratification for the ATE, ATT, and ATC. See [get_w_from_ps()] for details.
#' Ordinal treatments are treated exactly the same as non-order multi-category
#' treatments except that additional models are available to estimate the
#' generalized propensity score (e.g., ordinal logistic regression).
#'
#' ## Continuous Treatments
#'
#' For continuous treatments, weights are estimated as
#' \eqn{w_i = f_A(a_i) / f_{A|X}(a_i)}, where \eqn{f_A(a_i)} (known as the
#' stabilization factor) is
#' the unconditional density of treatment evaluated the observed treatment value
#' and \eqn{f_{A|X}(a_i)} (known as the generalized propensity score) is the
#' conditional density of treatment given the covariates evaluated at the
#' observed value of treatment. The shape of \eqn{f_A(.)} and \eqn{f_{A|X}(.)}
#' is controlled by the `density` argument described below (normal distributions
#' by default), and the predicted values used for the mean of the conditional
#' density are estimated using linear regression. Kernel density estimation can
#' be used instead of assuming a specific density for the numerator and
#' denominator by setting `density = "kernel"`. Other arguments to [density()]
#' can be specified to refine the density estimation parameters.
#'
#' ## Longitudinal Treatments
#'
#' For longitudinal treatments, the weights are the product of the weights
#' estimated at each time point.
#'
#' ## Sampling Weights
#'
#' Sampling weights are supported through `s.weights` in all scenarios except
#' for multi-category treatments with `link = "bayes.probit"` and for binary and
#' continuous treatments with `missing = "saem"` (see below). Warning messages
#' may appear otherwise about non-integer successes, and these can be ignored.
#'
#' ## Missing Data
#'
#' In the presence of missing data, the following value(s) for `missing` are
#' allowed:
#'     \describe{
#'       \item{`"ind"` (default)}{First, for each variable with missingness, a new missingness indicator variable is created which takes the value 1 if the original covariate is `NA` and 0 otherwise. The missingness indicators are added to the model formula as main effects. The missing values in the covariates are then replaced with the covariate medians (this value is arbitrary and does not affect estimation). The weight estimation then proceeds with this new formula and set of covariates. The covariates output in the resulting `weightit` object will be the original covariates with the `NA`s.
#'       }
#'       \item{`"saem"`}{For binary treatments with `link = "logit"` or continuous treatments, a stochastic approximation version of the EM algorithm (SAEM) is used via the \CRANpkg{misaem} package. No additional covariates are created. See Jiang et al. (2019) for information on this method. In some cases, this is a suitable alternative to multiple imputation.
#'       }
#'    }
#'
#' ## M-estimation
#'
#' For binary treatments, M-estimation is supported when `link` is neither
#' `"flic"` nor `"flac"` (see below). For multi-category treatments,
#' M-estimation is supported when `multi.method` is `"weightit"` (the default)
#' or `"glm"`. For continuous treatments, M-estimation is supported when
#' `density` is not `"kernel"`. The conditional treatment variance and
#' unconditional treatment mean and variance are included as parameters to
#' estimate, as these all go into calculation of the weights. For all treatment
#' types, M-estimation is not supported when `missing = "saem"`. See
#' [glm_weightit()] and `vignette("estimating-effects")` for details. For
#' longitudinal treatments, M-estimation is supported whenever the underlying
#' methods are.
#'
#' @section Additional Arguments: For binary treatments, the following
#'   additional argument can be specified:
#' \describe{
#'   \item{`link`}{the link used in the generalized linear model for the propensity scores. `link` can be any of those allowed by [binomial()] as well as `"loglog"` and `"clog"`. A `br.` prefix can be added (e.g., `"br.logit"`); this changes the fitting method to the bias-corrected generalized linear models implemented in the \CRANpkg{brglm2} package. `link` can also be either `"flic"` or `"flac"` to fit the corresponding Firth corrected logistic regression models implemented in the \CRANpkg{logistf} package.}
#' }
#'
#'   For multi-category treatments, the following additional arguments can be
#'   specified:
#' \describe{
#'   \item{`multi.method`}{the method used to estimate the generalized propensity scores. Allowable options include `"weightit"` (the default) to use multinomial logistic regression implemented in \pkg{WeightIt}, `"glm"` to use a series of binomial models using [glm()], `"mclogit"` to use multinomial logistic regression as implemented in \pkgfun{mclogit}{mblogit}, `"mnp"` to use Bayesian multinomial probit regression as implemented in \pkgfun{MNP}{MNP}, and `"brmultinom"` to use bias-reduced multinomial logistic regression as implemented in \pkgfun{brglm2}{brmultinom}. `"weightit"` and `"mclogit"` should give near-identical results, the main difference being increased robustness and customizability when using `"mclogit"` at the expense of not being able to use M-estimation to compute standard errors after weighting. For ordered treatments, allowable options include `"weightit"` (the default) to use ordinal regression implemented in \pkg{WeightIt} or `"polr"` to use ordinal regression implemented in \pkgfun{MASS}{polr}, unless `link` is `"br.logit"`, in which case bias-reduce ordinal logistic regression as implemented in \pkgfun{brglm2}{bracl} is used. Ignored when `missing = "saem"`. Using the defaults allows for the use of M-estimation and requires no additional dependencies, but other packages may provide benefits such as speed and flexibility.}
#'   \item{`link`}{The link used in the multinomial, binomial, or ordered regression model for the generalized propensity scores depending on the argument supplied to `multi.method`. When `multi.method = "glm"`, `link` can be any of those allowed by [binomial()]. When treatment is ordered and `multi.method` is `"weightit"` or `"polr"`, `link` can be any of those allowed by `MASS::polr()` or `"br.logit"`. Otherwise, `link` should be `"logit"` or not specified.}
#' }
#'
#'   For continuous treatments, the following additional arguments may be
#'   supplied:
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
#'   Additional arguments to `glm()` can be specified as well when it is used
#'   for fitting. The `method` argument in `glm()` is renamed to `glm.method`.
#'   This can be used to supply alternative fitting functions, such as those
#'   implemented in the \CRANpkg{glm2} package. Other arguments to `weightit()`
#'   are passed to `...` in `glm()`. In the presence of missing data with `link
#'   = "logit"` and `missing = "saem"`, additional arguments are passed to
#'   \pkgfun{misaem}{miss.glm} and \pkgfun{misaem}{predict.miss.glm}, except the
#'   `method` argument in \pkgfun{misaem}{predict.miss.glm} is replaced with
#'   `saem.method`.
#'
#'   For continuous treatments in the presence of missing data with `missing =
#'   "saem"`, additional arguments are passed to \pkgfun{misaem}{miss.lm} and
#'   \pkgfun{misaem}{predict.miss.lm}.
#'
#' @section Additional Outputs:
#' \describe{
#'   \item{`obj`}{When `include.obj = TRUE`, the (generalized) propensity score model fit. For binary treatments, the output of the call to [glm()] or the requested fitting function. For multi-category treatments, the output of the call to the fitting function (or a list thereof if `multi.method = "glm"`). For continuous treatments, the output of the call to `glm()` for the predicted values in the denominator density.
#'   }
#' }
#'
#' @details NULL
#'
#' @seealso [weightit()], [weightitMSM()], [get_w_from_ps()]
#'
#' @references ## Binary treatments
#'
#' - `estimand = "ATO"`
#'
#' Li, F., Morgan, K. L., & Zaslavsky, A. M. (2018). Balancing covariates via
#' propensity score weighting. *Journal of the American Statistical
#' Association*, 113(521), 390–400. \doi{10.1080/01621459.2016.1260466}
#'
#' - `estimand = "ATM"`
#'
#' Li, L., & Greene, T. (2013). A Weighting Analogue to Pair Matching in
#' Propensity Score Analysis. *The International Journal of Biostatistics*,
#' 9(2). \doi{10.1515/ijb-2012-0030}
#'
#' - `estimand = "ATOS"`
#'
#' Crump, R. K., Hotz, V. J., Imbens, G. W., & Mitnik, O. A. (2009). Dealing
#' with limited overlap in estimation of average treatment effects.
#' *Biometrika*, 96(1), 187–199. \doi{10.1093/biomet/asn055}
#'
#' - Other estimands
#'
#' Austin, P. C. (2011). An Introduction to Propensity Score Methods for
#' Reducing the Effects of Confounding in Observational Studies. *Multivariate
#' Behavioral Research*, 46(3), 399–424. \doi{10.1080/00273171.2011.568786}
#'
#' - Marginal mean weighting through stratification
#'
#' Hong, G. (2010). Marginal mean weighting through stratification: Adjustment
#' for selection bias in multilevel data. *Journal of Educational and Behavioral
#' Statistics*, 35(5), 499–531. \doi{10.3102/1076998609359785}
#'
#' - Bias-reduced logistic regression
#'
#' See references for the \CRANpkg{brglm2} package.
#'
#' - Firth corrected logistic regression
#'
#' Puhr, R., Heinze, G., Nold, M., Lusa, L., & Geroldinger, A. (2017). Firth’s
#' logistic regression with rare events: Accurate effect estimates and
#' predictions? *Statistics in Medicine*, 36(14), 2302–2317.
#' \doi{10.1002/sim.7273}
#'
#' - SAEM logistic regression for missing data
#'
#' Jiang, W., Josse, J., & Lavielle, M. (2019). Logistic regression with missing
#' covariates — Parameter estimation, model selection and prediction within a
#' joint-modeling framework. *Computational Statistics & Data Analysis*, 106907.
#' \doi{10.1016/j.csda.2019.106907}
#'
#' ## Multi-Category Treatments
#'
#' - `estimand = "ATO"`
#'
#' Li, F., & Li, F. (2019). Propensity score weighting for causal inference with
#' multiple treatments. *The Annals of Applied Statistics*, 13(4), 2389–2415.
#' \doi{10.1214/19-AOAS1282}
#'
#' - `estimand = "ATM"`
#'
#' Yoshida, K., Hernández-Díaz, S., Solomon, D. H., Jackson, J. W., Gagne, J.
#' J., Glynn, R. J., & Franklin, J. M. (2017). Matching weights to
#' simultaneously compare three treatment groups: Comparison to three-way
#' matching. *Epidemiology* (Cambridge, Mass.), 28(3), 387–395.
#' \doi{10.1097/EDE.0000000000000627}
#'
#' - Other estimands
#'
#' McCaffrey, D. F., Griffin, B. A., Almirall, D., Slaughter, M. E., Ramchand,
#' R., & Burgette, L. F. (2013). A Tutorial on Propensity Score Estimation for
#' Multiple Treatments Using Generalized Boosted Models. *Statistics in
#' Medicine*, 32(19), 3388–3414. \doi{10.1002/sim.5753}
#'
#' - Marginal mean weighting through stratification
#'
#' Hong, G. (2012). Marginal mean weighting through stratification: A
#' generalized method for evaluating multivalued and multiple treatments with
#' nonexperimental data. *Psychological Methods*, 17(1), 44–60.
#' \doi{10.1037/a0024918}
#'
#' ## Continuous treatments
#'
#' Robins, J. M., Hernán, M. Á., & Brumback, B. (2000). Marginal Structural
#' Models and Causal Inference in Epidemiology. *Epidemiology*, 11(5), 550–560.
#'
#' - Using non-normal conditional densities
#'
#' Naimi, A. I., Moodie, E. E. M., Auger, N., & Kaufman, J. S. (2014).
#' Constructing Inverse Probability Weights for Continuous Exposures: A
#' Comparison of Methods. *Epidemiology*, 25(2), 292–299.
#' \doi{10.1097/EDE.0000000000000053}
#'
#' - SAEM linear regression for missing data
#'
#' Jiang, W., Josse, J., & Lavielle, M. (2019). Logistic regression with missing
#' covariates — Parameter estimation, model selection and prediction within a
#' joint-modeling framework. *Computational Statistics & Data Analysis*, 106907.
#' \doi{10.1016/j.csda.2019.106907}
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
  fit.obj <- NULL

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

  if (ncol(covs) > 1L) {
    if (missing == "saem") {
      covs0 <- covs
      for (i in colnames(covs)[anyNA_col(covs)]) covs0[is.na(covs0[, i]), i] <- covs0[!is.na(covs0[, i]), i][1L]
      colinear.covs.to.remove <- colnames(covs)[colnames(covs) %nin% colnames(make_full_rank(covs0))]
    }
    else {
      colinear.covs.to.remove <- colnames(covs)[colnames(covs) %nin% colnames(make_full_rank(covs))]
    }

    covs <- covs[, colnames(covs) %nin% colinear.covs.to.remove, drop = FALSE]
  }

  t.lev <- get_treated_level(treat, estimand, focal)
  treat <- binarize(treat, one = t.lev)

  #Process link
  acceptable.links <- {
    if (missing == "saem") "logit"
    else c(expand_grid_string(c("", "br."), c("logit", "probit", "cloglog", "loglog", "identity", "log", "clog", "cauchit")),
           "flic", "flac")
  }

  link <- ...get("link")

  if (is_null(link)) {
    link <- acceptable.links[1L]
  }
  else {
    link <- acceptable.links[pmatch(link, acceptable.links, nomatch = 0L)][1L]

    if (is.na(link)) {
      .err(sprintf("only %s allowed as the link for binary treatments%s",
                   word_list(acceptable.links, quotes = TRUE, is.are = TRUE),
                   if (missing == "saem") ' with `missing = "saem"`' else ""))
    }
  }

  use.br <- startsWith(link, "br.")
  if (use.br) {
    link <- substr(link, 4L, nchar(link))
  }

  use.logistf <- link %in% c("flic", "flac")

  if (missing == "saem") {

    if (!all_the_same(s.weights)) {
      .err('sampling weights cannot be used with `missing = "saem"`')
    }

    rlang::check_installed("misaem")

    saem.method <- ...get("saem.method", "map")
    control <- ...get("control", list())

    if (is_null(control[["var_cal"]])) control[["var_cal"]] <- TRUE #need TRUE to bypass error in miss.glm.fit()
    if (is_null(control[["ll_obs_cal"]])) control[["ll_obs_cal"]] <- FALSE

    data <- data.frame(treat, covs)

    withCallingHandlers({
      verbosely({
        fit <- misaem::miss.glm(formula(data), data = data, control = as.list(control))
      }, verbose = verbose)
    },
    warning = function(w) {
      w <- conditionMessage(w)
      if (w != "one argument not used by format '%i '") {
        .wrn("(from `misaem::miss.glm()`) ", w, tidy = FALSE)
      }
      invokeRestart("muffleWarning")
    })

    p.score <- drop(predict(fit, newdata = covs, method = saem.method))
  }
  else if (use.logistf) {
    rlang::check_installed("logistf")

    fit_fun <- switch(link,
                      "flic" = logistf::flic,
                      "flac" = logistf::flac)

    ctrl_fun <- logistf::logistf.control

    control <- do.call(ctrl_fun, c(as.list(...get("control")),
                                   ...mget(setdiff(names(formals(ctrl_fun))[pmatch(...names(), names(formals(ctrl_fun)), 0L)],
                                                   names(...get("control"))))))

    modctrl_fun <- logistf::logistf.mod.control
    modcontrol <- do.call(modctrl_fun, c(as.list(...get("modcontrol")),
                                         ...mget(setdiff(names(formals(modctrl_fun))[pmatch(...names(), names(formals(modctrl_fun)), 0L)],
                                                         names(...get("modctrl_fun"))))))

    withCallingHandlers({verbosely({
      data <- data.frame(treat, covs)
      formula <- if (ncol(covs) > 0L) formula(data) else treat ~ 1

      fit <- do.call(fit_fun, list(formula, data = data,
                                   weights = s.weights,
                                   control = control,
                                   modcontrol = modcontrol,
                                   pl = FALSE),
                     quote = TRUE)
    }, verbose = verbose)},
    warning = function(w) {
      w <- conditionMessage(w)
      if (w != "non-integer #successes in a binomial glm!")
        .wrn(sprintf("(from `logistf::%s()`) ", link),
             w, tidy = FALSE)
      invokeRestart("muffleWarning")
    })

    p.score <- fit$predict
  }
  else {
    link <- .make_link(link)

    if (use.br) {
      rlang::check_installed("brglm2")

      ctrl_fun <- brglm2::brglmControl
      glm_method <- brglm2::brglmFit
      family <- binomial(link = link)
    }
    else {
      ctrl_fun <- stats::glm.control
      glm_method <- ...get("glm.method", stats::glm.fit)
      family <- quasibinomial(link = link)
    }

    control <- do.call(ctrl_fun, c(as.list(...get("control")),
                                   ...mget(setdiff(names(formals(ctrl_fun))[pmatch(...names(), names(formals(ctrl_fun)), 0L)],
                                                   names(...get("control"))))))

    start <- mustart <- NULL

    if (family$link %in% c("log", "clog", "identity")) {
      #Need starting values because links are unbounded
      start <- c(family$linkfun(w.m(treat, s.weights)), rep.int(0, ncol(covs)))
    }
    else {
      #Default starting values from glm.fit() without weights; these
      #work better with s.weights than usual default.
      mustart <- .25 + .5 * treat
    }

    withCallingHandlers({verbosely({
      if (isTRUE(...get("quick"))) {
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
        formula <- if (ncol(covs) > 0L) formula(data) else treat ~ 1

        fit <- do.call(stats::glm, list(formula, data = data,
                                        weights = s.weights,
                                        mustart = mustart,
                                        start = start,
                                        family = family,
                                        method = glm_method,
                                        control = control),
                       quote = TRUE)
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
      psi_treat = function(Btreat, Xtreat, A, SW) {
        .psi(B = Btreat, X = Xtreat, y = A, weights = SW)
      },
      wfun = function(Btreat, Xtreat, A) {
        XB <- drop(Xtreat %*% Btreat)
        p <- family$linkinv(XB)
        .get_w_from_ps_internal_bin(ps = p, treat = A, estimand,
                                    stabilize = stabilize, subclass = subclass)
      },
      dw_dBtreat = switch(estimand,
                          "ATOS" = NULL,
                          function(Btreat, Xtreat, A, SW) {
                            XB <- drop(Xtreat %*% Btreat)
                            ps <- family$linkinv(XB)
                            .dw_dp_bin(ps, A, estimand = estimand) * family$mu.eta(XB) * Xtreat
                          }),
      # hess_treat = function(Btreat, Xtreat, A, SW) {
      #   #Using expected information because observed is too complicated
      #   XB <- drop(Xtreat %*% Btreat)
      #   ps <- family$linkinv(XB)
      #
      #   d1mus <- family$mu.eta(XB)
      #   varmus <- family$variance(ps)
      #
      #   crossprod(Xtreat, Xtreat * (d1mus^2 * SW / -varmus))
      # },
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
  fit.obj <- NULL

  covs <- covs[subset, , drop = FALSE]
  treat <- factor(treat[subset])
  s.weights <- s.weights[subset]


  missing <- .process_missing2(missing, covs)

  if (missing == "ind") {
    covs <- add_missing_indicators(covs)
  }

  for (i in seq_col(covs)) {
    covs[, i] <- .make_closer_to_1(covs[, i])
  }

  if (ncol(covs) > 1L) {
    if (missing == "saem") {
      covs0 <- covs
      for (i in colnames(covs)[anyNA_col(covs)]) covs0[is.na(covs0[, i]), i] <- covs0[!is.na(covs0[, i]), i][1L]
      colinear.covs.to.remove <- colnames(covs)[colnames(covs) %nin% colnames(make_full_rank(covs0))]
    }
    else {
      colinear.covs.to.remove <- colnames(covs)[colnames(covs) %nin% colnames(make_full_rank(covs))]
    }

    covs <- covs[, colnames(covs) %nin% colinear.covs.to.remove, drop = FALSE]
  }

  for (i in c("use.mlogit", "use.mclogit")) {
    if (is_not_null(...get(i))) {
      .wrn(sprintf('`%s` is no longer accepted and will be ignored; use `multi.method` instead. See `help("method_glm")` for details', i))
    }
  }

  ord.treat <- is.ordered(treat) && missing != "saem"

  multi.method <- ...get("multi.method")
  link <- ...get("link")
  link.ignored <- FALSE

  #Process multi.method
  if (is_null(multi.method)) {
    if (ord.treat) {
      multi.method <- {
        if (identical(link, "br.logit")) "bracl"
        else "weightit"
      }
    }
    else {
      multi.method <- {
        if (missing == "saem") "saem"
        else if (identical(link, "bayes.probit")) "mnp"
        else if (identical(link, "br.logit")) "brmultinom"
        else "weightit"
      }
    }
  }
  else {
    if (missing == "saem") {
      if (!identical(multi.method, "saem") &&
          !identical(multi.method, "glm")) {
        .wrn('`multi.method` is ignored when `missing = "saem"`')
      }

      multi.method <- "saem"
    }
    else {
      chk::chk_string(multi.method)
      multi.method <- tolower(multi.method)
      if (multi.method == "mblogit") multi.method <- "mclogit"

      if (ord.treat) {
        allowable.multi.methods <- c("weightit", "polr", "glm", "mclogit", "mnp", "brmultinom", "bracl")
        multi.method <- match_arg(multi.method, allowable.multi.methods)

        ord.treat <- (multi.method %in% c("weightit", "polr", "bracl"))
      }
      else {
        allowable.multi.methods <- c("weightit", "glm", "mclogit", "mnp", "brmultinom")
        multi.method <- match_arg(multi.method, allowable.multi.methods)
      }
    }
  }

  link.ignored <- multi.method %in% c("mclogit", "mnp", "brmultinom", "bracl")

  #Process link
  acceptable.links <- switch(
    multi.method,
    "weightit" = {
      if (ord.treat) c("logit", "probit", "cloglog", "loglog",
                       "identity", "log", "clog", "cauchit")
      else "logit"
    },
    "glm" = c("logit", "probit", "cloglog", "loglog", "identity",
              "log", "clog", "cauchit"),
    "polr" = c("logit", "probit", "loglog", "cloglog", "cauchit"),
    "brmultinom" = c("br.logit", "logit"),
    "bracl" = c("br.logit", "logit"),
    "mnp" = c("bayes.probit", "probit"),
    "mclogit" = "logit",
    "saem" = "logit"
  )

  if (is_null(link)) {
    link <- acceptable.links[1L]
  }
  else if (link.ignored) {
    if (!any_apply(acceptable.links, identical, link)) {
      .wrn(sprintf("`link` is ignored when `%s = %s`",
                   if (missing == "saem") "missing" else "multi.method",
                   add_quotes(multi.method)))
    }

    link <- acceptable.links[1L]
  }
  else {
    link <- acceptable.links[pmatch(link, acceptable.links, nomatch = 0L)][1L]
    if (is.na(link)) {
      .err(sprintf("only %s allowed as the link for %smulti-category treatments with `multi.method = %s`",
                   if (ord.treat) "ordinal " else "",
                   word_list(acceptable.links, quotes = TRUE, is.are = TRUE),
                   add_quotes(multi.method)))
    }
  }

  # Fit model
  if (multi.method == "weightit") {
    if (ord.treat) {
      verbosely({
        fit.obj <- .ordinal_weightit.fit(x = cbind(`(Intercept)` = 1, covs),
                                         y = treat,
                                         weights = s.weights,
                                         hess = FALSE,
                                         link = link)
      }, verbose = verbose)
    }
    else {
      verbosely({
        fit.obj <- .multinom_weightit.fit(x = cbind(`(Intercept)` = 1, covs),
                                          y = treat,
                                          weights = s.weights,
                                          hess = FALSE)
      }, verbose = verbose)
    }

    ps <- fit.obj$fitted.values
  }
  else if (multi.method == "glm") {
    ps <- make_df(levels(treat), nrow = length(treat))

    ctrl_fun <- stats::glm.control

    control <- do.call(ctrl_fun, c(as.list(...get("control")),
                                   ...mget(setdiff(names(formals(ctrl_fun))[pmatch(...names(), names(formals(ctrl_fun)), 0L)],
                                                   names(...get("control"))))))

    link <- .make_link(link)

    family <- quasibinomial(link = link)

    fit.obj <- .psi.list <- make_list(levels(treat))

    for (i in levels(treat)) {
      t_i <- as.numeric(treat == i)
      data_i <- data.frame(t_i, covs)

      start <- mustart <- NULL

      if (family$link %in% c("log", "clog", "identity")) {
        #Need starting values because links are unbounded
        start <- c(family$linkfun(w.m(t_i, s.weights)), rep.int(0, ncol(covs)))
      }
      else {
        #Default starting values from glm.fit() without weights; these
        #work better with s.weights than usual default.
        mustart <- .25 + .5 * t_i
      }

      withCallingHandlers({verbosely({
        fit.obj[[i]] <- do.call(stats::glm,
                                list(formula(data_i),
                                     data = data_i,
                                     family = family,
                                     weights = s.weights,
                                     control = control,
                                     mustart = mustart,
                                     start = start),
                                quote = TRUE)
      }, verbose = verbose)},
      warning = function(w) {
        w <- conditionMessage(w)
        if (w != "non-integer #successes in a binomial glm!") {
          .wrn("(from `glm()`) ", w, tidy = FALSE)
        }
        invokeRestart("muffleWarning")
      })

      ps[[i]] <- fit.obj[[i]]$fitted.values
      .psi.list[[i]] <- .get_glm_psi(fit.obj[[i]])
    }
  }
  else if (multi.method == "polr") {
    rlang::check_installed("MASS")
    if (link == "logit") link <- "logistic"

    data <- data.frame(treat, covs)
    formula <- formula(data)

    tryCatch({verbosely({
      fit.obj <- do.call(MASS::polr,
                         list(formula,
                              data = data,
                              weights = s.weights,
                              Hess = FALSE,
                              model = FALSE,
                              method = link,
                              contrasts = NULL),
                         quote = TRUE)
    }, verbose = verbose)},
    error = function(e) {
      .err(sprintf("There was a problem fitting the ordinal %s regressions with `polr()`.\n       Try again with an un-ordered treatment.\n       Error message: (from `MASS::polr()`) %s",
                   link, conditionMessage(e)), tidy = FALSE)})

    ps <- fit.obj$fitted.values
  }
  else if (multi.method == "bracl") {
    rlang::check_installed("brglm2")

    ctrl_fun <- brglm2::brglmControl

    control <- do.call(ctrl_fun, c(as.list(...get("control")),
                                   ...mget(setdiff(names(formals(ctrl_fun))[pmatch(...names(), names(formals(ctrl_fun)), 0L)],
                                                   names(...get("control"))))))

    data <- data.frame(treat, covs)
    formula <- formula(data)

    tryCatch({verbosely({
      fit.obj <- do.call(brglm2::bracl,
                         list(formula,
                              data = data,
                              weights = s.weights,
                              control = control,
                              parallel = ...get("parallel", FALSE)),
                         quote = TRUE)
    }, verbose = verbose)},
    error = function(e) {
      .err(sprintf("there was a problem with the bias-reduced ordinal logit regression.\n       Try a different link.\n       Error message: (from `brglm2::bracl()`) %s", conditionMessage(e)), tidy = FALSE)
    })

    ps <- fit.obj$fitted.values
  }
  else if (multi.method == "saem") {

    if (!all_the_same(s.weights)) {
      .err('sampling weights cannot be used with `missing = "saem"`')
    }

    rlang::check_installed("misaem")

    saem.method <- ...get("saem.method", "map")
    control <- ...get("control", list())

    if (is_null(control[["var_cal"]])) control[["var_cal"]] <- TRUE #need TRUE to bypass error in miss.glm.fit()
    if (is_null(control[["ll_obs_cal"]])) control[["ll_obs_cal"]] <- FALSE

    ps <- make_df(levels(treat), nrow = length(treat))

    fit.obj <- make_list(levels(treat))

    for (i in levels(treat)) {
      t_i <- as.numeric(treat == i)
      data_i <- data.frame(t_i, covs)

      withCallingHandlers({
        verbosely({
          fit.obj[[i]] <- misaem::miss.glm(formula(data_i), data = data_i,
                                           control = as.list(control))
        }, verbose = verbose)
      },
      warning = function(w) {
        w <- conditionMessage(w)
        if (w != "one argument not used by format '%i '") {
          .wrn("(from `misaem::miss.glm()`) ", w, tidy = FALSE)
        }
        invokeRestart("muffleWarning")
      })

      ps[[i]] <- drop(predict(fit.obj[[i]], newdata = covs, method = saem.method))
    }
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

    if (is_not_null(...get("random"))) {
      random <- get_covs_and_treat_from_formula(...get("random"), data = .data)$reported.covs[subset, , drop = FALSE]
      data <- cbind(data.frame(random), data.frame(treat = treat, .s.weights = s.weights, covs))
      covnames <- names(data)[-c(seq_col(random), ncol(random) + (1:2))]
      tname <- names(data)[ncol(random) + 1L]
      ctrl_fun <- mclogit::mmclogit.control
    }
    else {
      data <- data.frame(treat = treat, .s.weights = s.weights, covs)
      covnames <- names(data)[-(1:2)]
      tname <- names(data)[1L]
      ctrl_fun <- mclogit::mclogit.control
    }
    form <- reformulate(covnames, tname)

    control <- c(as.list(...get("control")),
                 ...mget(setdiff(names(formals(ctrl_fun))[pmatch(...names(), names(formals(ctrl_fun)), 0L)],
                                 names(...get("control")))))
    tryCatch({verbosely({
      fit.obj <- do.call(mclogit::mblogit,
                         list(form,
                              data = data,
                              weights = quote(.s.weights),
                              random = ...get("random"),
                              method = ...get("mclogit.method"),
                              estimator = ...get("estimator", eval(formals(mclogit::mclogit)[["estimator"]])),
                              dispersion = ...get("dispersion", eval(formals(mclogit::mclogit)[["dispersion"]])),
                              groups = ...get("groups"),
                              control = control))

    }, verbose = verbose)},
    error = function(e) {
      .err(sprintf("there was a problem fitting the multinomial %s regression with `mblogit()`.\n       Try a different `multi.method`.\nError message: (from `mclogit::mblogit()`) %s",
                   link, conditionMessage(e)), tidy = FALSE)
    })

    ps <- fitted(fit.obj)
    colnames(ps) <- levels(treat)
  }
  else if (multi.method == "brmultinom") {
    rlang::check_installed("brglm2")
    data <- data.frame(treat, covs)
    formula <- formula(data)
    ctrl_fun <- brglm2::brglmControl

    control <- do.call(ctrl_fun, c(as.list(...get("control")),
                                   ...mget(setdiff(names(formals(ctrl_fun))[pmatch(...names(), names(formals(ctrl_fun)), 0L)],
                                                   names(...get("control"))))))

    tryCatch({verbosely({
      fit.obj <- do.call(brglm2::brmultinom,
                         list(formula, data,
                              weights = s.weights,
                              control = control),
                         quote = TRUE)
    }, verbose = verbose)},
    error = function(e) {
      .err(sprintf("There was a problem with the bias-reduced multinomial logit regression. Try a different `multi.method`.\n       Error message: (from `brglm2::brmultinom()`) %s", conditionMessage(e)), tidy = FALSE)
    })

    ps <- fit.obj$fitted.values
  }

  #ps should be matrix of probs for each treat
  #Computing weights
  w <- .get_w_from_ps_internal_multi(ps = ps, treat = treat, estimand, focal = focal,
                                     stabilize = stabilize, subclass = subclass)

  #Get Mparts
  Mparts <- NULL
  if (multi.method == "weightit") {
    Mparts <- list(
      psi_treat = function(Btreat, Xtreat, A, SW) {
        fit.obj$psi(Btreat, Xtreat, A, SW)
      },
      wfun = function(Btreat, Xtreat, A) {
        ps <- fit.obj$get_p(Btreat, Xtreat)
        .get_w_from_ps_internal_multi(ps = ps, treat = A, estimand, focal = focal,
                                      stabilize = stabilize, subclass = subclass)
      },
      Xtreat = fit.obj$x,
      A = treat,
      btreat = fit.obj$coefficients
    )
  }
  else if (multi.method == "glm") {
    Mparts <- list(
      psi_treat = function(Btreat, Xtreat, A, SW) {
        Btreat <- matrix(Btreat, nrow = ncol(Xtreat))

        do.call("cbind", lapply(seq_len(nlevels(A)), function(i) {
          .psi.list[[i]](Btreat[, i], Xtreat, A, SW)
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
      btreat = unlist(grab(fit.obj, "coefficients"))
    )
  }

  list(w = w, fit.obj = fit.obj,
       Mparts = Mparts)
}

weightit2glm.cont <- function(covs, treat, s.weights, subset, stabilize, missing, verbose, ...) {

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

  if (ncol(covs) > 1L) {
    if (missing == "saem") {
      covs0 <- covs
      for (i in colnames(covs)[anyNA_col(covs)]) {
        covs0[is.na(covs0[, i]), i] <- covs0[!is.na(covs0[, i]), i][1L]
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
  densfun <- .get_dens_fun(use.kernel = isTRUE(...get("use.kernel")), bw = ...get("bw"),
                           adjust = ...get("adjust"), kernel = ...get("kernel"),
                           n = ...get("n"), treat = treat, density = ...get("density"),
                           weights = s.weights)

  #Stabilization - get dens.num
  un_p <- mean_fast(s.weights * treat)
  un_s2 <- mean_fast(s.weights * (treat - un_p)^2)

  log.dens.num <- densfun((treat - un_p) / sqrt(un_s2), log = TRUE)

  #Estimate GPS
  link <- ...get("link", "identity")

  if (missing == "saem") {

    if (!all_the_same(s.weights)) {
      .err('sampling weights cannot be used with `missing = "saem"`')
    }

    rlang::check_installed("misaem")

    acceptable.links <- "identity"
    which.link <- acceptable.links[pmatch(link, acceptable.links, nomatch = 0L)][1L]

    if (is.na(which.link)) {
      .err(sprintf('only %s allowed as the link for continuous treatments with `missing = "saem"`',
                   word_list(acceptable.links, quotes = TRUE, is.are = TRUE)))
    }

    data <- data.frame(treat, covs)
    formula <- formula(data)

    withCallingHandlers({verbosely({
      fit <- misaem::miss.lm(formula, data = data, control = as.list(...get("control")))
    }, verbose = verbose)},
    warning = function(w) {
      w <- conditionMessage(w)
      if (w != "one argument not used by format '%i '") {
        .wrn("(from `misaem::miss.lm()`) ", w, tidy = FALSE)
      }
      invokeRestart("muffleWarning")
    })

    saem.method <- ...get("saem.method", "map")

    p <- drop(predict(fit, newdata = covs, method = saem.method))
  }
  else {
    acceptable.links <- c("identity", "log", "inverse")

    link <- acceptable.links[pmatch(link, acceptable.links, nomatch = 0L)][1L]
    if (is.na(link)) {
      .err(sprintf("only %s allowed as the link for continuous treatments",
                   word_list(acceptable.links, quotes = TRUE, is.are = TRUE)))
    }

    family <- gaussian(link = link)

    verbosely({
      if (isTRUE(...get("quick"))) {
        fit <- do.call(stats::glm.fit, list(y = treat,
                                            x = cbind(`(Intercept)` = rep.int(1, length(treat)), covs),
                                            weights = s.weights,
                                            family = family,
                                            control = as.list(...get("control"))),
                       quote = TRUE)
      }
      else {
        data <- data.frame(treat, covs)
        formula <- if (ncol(covs) > 0L) formula(data) else treat ~ 1

        fit <- do.call(stats::glm, list(formula, data = data,
                                        weights = s.weights,
                                        family = family,
                                        control = as.list(...get("control"))),
                       quote = TRUE)
      }
    }, verbose = verbose)

    p <- fit$fitted.values
  }

  s2 <- mean_fast(s.weights * (treat - p)^2)

  #Get weights
  log.dens.denom <- densfun((treat - p) / sqrt(s2), log = TRUE)

  w <- exp(log.dens.num - log.dens.denom)

  if (isTRUE(...get("plot"))) {
    d.n <- attr(log.dens.num, "density")
    d.d <- attr(log.dens.denom, "density")
    plot_density(d.n, d.d, log = TRUE)
  }

  Mparts <- NULL
  if (missing != "saem" && !identical(...get("density"), "kernel")) {
    Mparts <- list(
      psi_treat = function(Btreat, Xtreat, A, SW) {
        un_s2 <- exp(Btreat[1L])
        un_p <- Btreat[2L]

        s2 <- exp(Btreat[3L])
        lin_pred <- drop(Xtreat %*% Btreat[-(1:3)])
        p <- family$linkinv(lin_pred)

        SW <- SW / mean_fast(SW)

        cbind(SW * (A - un_p)^2 - un_s2, #unconditional variance
              SW * (A - un_p), #unconditional mean
              SW * (A - p)^2 - s2, #conditional variance
              Xtreat * (SW * family$mu.eta(lin_pred) * (A - p) / family$variance(p))) #conditional mean
      },
      wfun = function(Btreat, Xtreat, A) {
        un_s2 <- exp(Btreat[1L])
        un_p <- Btreat[2L]
        log.dens.num <- densfun((A - un_p) / sqrt(un_s2), log = TRUE)

        s2 <- exp(Btreat[3L])
        lin_pred <- drop(Xtreat %*% Btreat[-(1:3)])
        p <- family$linkinv(lin_pred)
        log.dens.denom <- densfun((A - p) / sqrt(s2), log = TRUE)

        exp(log.dens.num - log.dens.denom)
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

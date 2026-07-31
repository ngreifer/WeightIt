#' Propensity Score Weighting Using Generalized Linear Models
#' @name method_glm
#' @usage NULL
#'
#' @description
#' This page explains the details of estimating weights from
#' generalized linear model-based propensity scores by setting `method = "glm"`
#' in the call to [weightit()] or [weightitMSM()]. This method can be used with
#' binary, multi-category, and continuous treatments, as well as for estimating censoring weights.
#'
#' In general, this method relies on estimating propensity scores with a
#' parametric generalized linear model and then converting those propensity
#' scores into weights using a formula that depends on the desired estimand. For
#' binary and multi-category treatments, a binomial or multinomial regression
#' model (logistic, by default) is used to estimate the propensity scores as the predicted probability
#' of being in each treatment given the covariates. For ordinal treatments, an
#' ordinal regression model is used to estimate generalized propensity scores.
#' For continuous treatments, a generalized linear model is used to estimate
#' generalized propensity scores as the conditional density of treatment given
#' the covariates.
#'
#' ## Binary Treatments
#'
#' For binary treatments, this method estimates the propensity scores using
#' [glm()]. The following estimands are allowed: ATE, ATT, ATC,
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
#' estimand are computed using the standard formulas.
#' Weights can also be computed using marginal mean weighting through
#' stratification for the ATE, ATT, and ATC. See [get_w_from_ps()] for details.
#' Ordinal treatments are treated exactly the same as non-ordinal multi-category
#' treatments except that additional models are available to estimate the
#' generalized propensity score (e.g., ordinal logistic regression).
#'
#' ## Continuous Treatments
#'
#' For continuous treatments, weights are estimated as
#' \eqn{w_i = f_A(a_i) / f_{A|X}(a_i)}, where \eqn{f_A(a_i)} (known as the
#' stabilization factor) is the unconditional density of treatment evaluated the
#' observed treatment value and \eqn{f_{A|X}(a_i)} (known as the generalized
#' propensity score) is the conditional density of treatment given the
#' covariates evaluated at the observed treatment value. The shape of
#' \eqn{f_{A|X}(.)} is controlled by the `density` argument
#' described below (normal distribution by default), and the predicted values
#' used for the mean of the conditional density are estimated using linear
#' regression. \eqn{f_A(.)} is estimated by marginalizing
#' over \eqn{f_{A|X}(.)}. Kernel density estimation can be used
#' instead of assuming a specific density for the denominator by
#' setting `density = "kernel"`. Other arguments to [density()] can be specified
#' to refine the density estimation parameters.
#'
#' ## Multilevel Treatment Models
#'
#' When the model `formula` contains \CRANpkg{lme4}-style random effects terms
#' (e.g., `treat ~ x1 + x2 + (1 | school)`), a multilevel (mixed-effects) model is
#' used to estimate the propensity scores. This can improve balance and overlap
#' when units are clustered (e.g., patients within hospitals or students within
#' schools). The grouping (and any random-slope) variables are taken from `data`.
#' For binary treatments, the model is fit using \pkgfun{lme4}{glmer} with the
#' requested `link` (`"logit"`, `"probit"`, or `"cloglog"`); for continuous
#' treatments, using \pkgfun{lme4}{lmer} or \pkgfun{lme4}{glmer}; and
#' for multi-category treatments, using \pkgfun{mclogit}{mblogit} with its
#' `random` argument (i.e., `multi.method` is set to `"mclogit"`). The propensity
#' scores are the predicted probabilities (or conditional densities) that include
#' the estimated random effects, i.e., cluster-specific predictions. M-estimation
#' is not available for these models; robust (`HC0`) or bootstrap standard errors
#' should be used instead when estimating treatment effects (see the M-estimation
#' section below and [glm_weightit()]).
#'
#' ## Censoring Weights
#'
#' For censoring weights, requested by wrapping the censoring indicator in [`.cens()`][.cens], a single binomial regression model of the probability of being censored is fit, and the weights are `1/P(C = 0 | X)` for the units still under observation and 0 for the censored units. All of the `link` options, `missing` options, and multilevel models described above apply unchanged, and the returned propensity score is the probability of *being censored*.
#'
#' ## Longitudinal Treatments
#'
#' For longitudinal treatments, the weights are the product of the weights
#' estimated at each time point.
#'
#' ## Sampling Weights
#'
#' Sampling weights are supported through `s.weights` in all scenarios except
#' for multi-category treatments with `multi.method = "mnp"` and for binary and
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
#' or `"glm"`. M-estimation is not supported when `subclass` is specified. For continuous treatments, M-estimation is supported when
#' `density` is not `"kernel"`. The conditional treatment variance and
#' unconditional treatment mean and variance are included as parameters to
#' estimate, as these all go into calculation of the weights. For all treatment
#' types, M-estimation is not supported when `missing = "saem"` or when the model
#' `formula` includes random effects terms (i.e., a multilevel model is fit; see
#' *Multilevel Treatment Models* above). See
#' [glm_weightit()] and `vignette("estimating-effects")` for details. For
#' longitudinal treatments, M-estimation is supported whenever the underlying
#' methods are.
#'
#' @section Additional Arguments:
#'
#' For binary treatments, the following additional argument can be specified:
#' \describe{
#'   \item{`link`}{the link used in the generalized linear model for the propensity scores. `link` can be any of those allowed by [binomial()] as well as `"loglog"` and `"clog"`. A `br.` prefix can be added (e.g., `"br.logit"`); this changes the fitting method to the bias-corrected generalized linear models implemented in the \CRANpkg{brglm2} package. `link` can also be either `"flic"` or `"flac"` to fit the corresponding Firth corrected logistic regression models implemented in the \CRANpkg{logistf} package.}
#'   \item{`subclass`}{`integer`; the number of subclasses to use for computing weights using marginal mean weighting through stratification (MMWS). If `NULL`, standard inverse probability weights (and their extensions) will be computed; if a number greater than 1, subclasses will be formed and weights will be computed based on subclass membership. See [get_w_from_ps()] for details and references.}
#' }
#'
#' For multi-category treatments, the following additional arguments can be specified:
#' \describe{
#'   \item{`multi.method`}{the method used to estimate the generalized propensity scores. Allowable options include `"weightit"` (the default) to use multinomial logistic regression implemented in \pkg{WeightIt}, `"glm"` to use a series of binomial models using [glm()], `"mclogit"` to use multinomial logistic regression as implemented in \pkgfun{mclogit}{mblogit}, and `"mnp"` to use Bayesian multinomial probit regression as implemented in \pkgfun{MNP}{MNP}. `"weightit"` and `"mclogit"` should give near-identical results, the main difference being increased robustness and customizability when using `"mclogit"` at the expense of not being able to use M-estimation to compute standard errors after weighting. For ordered treatments, allowable options include `"weightit"` (the default) to use ordinal regression implemented in \pkg{WeightIt} or `"polr"` to use ordinal regression implemented in \pkgfun{MASS}{polr}. Ignored when `missing = "saem"`. Using the defaults allows for the use of M-estimation and requires no additional dependencies, but other packages may provide benefits such as speed and flexibility.
#'
#'   Bias-reduced models are requested by adding a `br.` prefix to `link` rather than by supplying a `multi.method` (see below). Because `MASS::polr()` cannot fit bias-reduced models, supplying a `br.` link with `multi.method = "polr"` uses the \pkg{WeightIt} implementation instead, which also makes M-estimation available.}
#'   \item{`link`}{The link used in the multinomial, binomial, or ordered regression model for the generalized propensity scores depending on the argument supplied to `multi.method`. When `multi.method = "glm"`, `link` can be any of those allowed by [binomial()] as well as `"loglog"` and `"clog"`. When treatment is ordered and `multi.method` is `"weightit"` or `"polr"`, `link` can be any of `"logit"`, `"probit"`, `"loglog"`, `"cloglog"`, or `"cauchit"`. Otherwise, `link` should be `"logit"` or not specified.
#'
#'   For unordered treatments with `multi.method = "weightit"` and for ordered treatments with `multi.method` of `"weightit"` or `"polr"`, a `br.` prefix can be added (e.g., `"br.logit"`) to request mean bias reduction, i.e., to solve the bias-reducing adjusted score equations of Firth (1993) rather than the score equations. See [multinom_weightit()] and [ordinal_weightit()] for details. The estimates are always finite, even when the maximum likelihood estimates are not, which can help when the treatment groups are nearly separated. M-estimation remains available.}
#'   \item{`subclass`}{`integer`; the number of subclasses to use for computing weights using marginal mean weighting through stratification (MMWS). If `NULL`, standard inverse probability weights (and their extensions) will be computed; if a number greater than 1, subclasses will be formed and weights will be computed based on subclass membership. See [get_w_from_ps()] for details and references.}
#' }
#'
#'   For continuous treatments, the following additional arguments may be
#'   supplied:
#'   \describe{
#'     \item{`density`}{A function corresponding the conditional density of the treatment. The standardized residuals of the treatment model will be fed through this function to produce the denominator of the generalized propensity score weights. If blank, [dnorm()] is used as recommended by Robins et al. (2000). This can also be supplied as a string containing the name of the function to be called. If the string contains underscores, the call will be split by the underscores and the latter splits will be supplied as arguments to the second argument and beyond. For example, if `density = "dt_2"` is specified, the density used will be that of a t-distribution with 2 degrees of freedom. Using a t-distribution can be useful when extreme outcome values are observed (Naimi et al., 2014).
#'
#' Can also be `"kernel"` to use kernel density estimation, which calls [density()] to estimate the denominator density for the weights. (This used to be requested by setting `use.kernel = TRUE`, which is now deprecated.)}
#'     \item{`bw`, `adjust`, `kernel`, `n`}{If `density = "kernel"`, the arguments to [density()]. The defaults are the same as those in `density()`.}
#'     \item{`link`}{The link used to fit the linear model for the generalized propensity score. Can be any allowed by [gaussian()].
#'     }
#'   }
#'
#'   Additional arguments to `glm()` can be specified as well when it is used
#'   for fitting. The `method` argument in `glm()` is renamed to `glm.method`.
#'   This can be used to supply alternative fitting functions, such as those
#'   implemented in the \CRANpkg{glm2} package. Other arguments to `weightit()`
#'   are passed to `...` in `glm()`. In the presence of missing data with
#'   `link = "logit"` and `missing = "saem"`, additional arguments are passed to
#'   \pkgfun{misaem}{miss.glm} and \pkgfun{misaem}{predict.miss.glm}, except the
#'   `method` argument in \pkgfun{misaem}{predict.miss.glm} is replaced with
#'   `saem.method`.
#'
#'   For continuous treatments in the presence of missing data with `missing = "saem"`, additional arguments are passed to \pkgfun{misaem}{miss.lm} and \pkgfun{misaem}{predict.miss.lm}.
#'
#'   When the model `formula` includes random effects terms (see *Multilevel Treatment Models* above), additional arguments are passed to the corresponding fitting function: \pkgfun{lme4}{glmer} for binary treatments, \pkgfun{mclogit}{mblogit} for multi-category treatments, and \pkgfun{lme4}{lmer} or \pkgfun{lme4}{glmer} for continuous treatments.
#'
#' @section Additional Outputs:
#' \describe{
#'   \item{`obj`}{When `include.obj = TRUE`, the (generalized) propensity score model fit. For binary treatments, the output of the call to [glm()] or the requested fitting function. For multi-category treatments, the output of the call to the fitting function (or a list thereof if `multi.method = "glm"`). For continuous treatments, the output of the call to `glm()` for the predicted values in the denominator density. When the model `formula` includes random effects terms, the output of the call to \pkgfun{lme4}{glmer}, \pkgfun{lme4}{lmer}, or \pkgfun{mclogit}{mblogit}.
#'   }
#' }
#'
#' @details NULL
#'
#' @seealso [weightit()], [weightitMSM()], [get_w_from_ps()]
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
#' - Bias-reduced regression
#'
#' Firth, D. (1993). Bias reduction of maximum likelihood estimates. *Biometrika*, 80(1), 27–38. \doi{10.1093/biomet/80.1.27}
#'
#' For binary treatments, see also the references for the \CRANpkg{brglm2} package, which does the fitting. For multi-category treatments, the fitting is done by \pkg{WeightIt}; see [multinom_weightit()] and [ordinal_weightit()], which document the adjustments used.
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
#' Li, F., & Li, F. (2019). Propensity score weighting for causal inference with  multiple treatments. *The Annals of Applied Statistics*, 13(4), 2389–2415. \doi{10.1214/19-AOAS1282}
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
#' McCaffrey, D. F., Griffin, B. A., Almirall, D., Slaughter, M. E., Ramchand, R., & Burgette, L. F. (2013). A Tutorial on Propensity Score Estimation for Multiple Treatments Using Generalized Boosted Models. *Statistics in Medicine*, 32(19), 3388–3414. \doi{10.1002/sim.5753}
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
#'
#' summary(W1)
#'
#' bal.tab(W1)
#'
#' #Balancing covariates with respect to race (multi-category)
#' (W2 <- weightit(race ~ age + educ + married +
#'                   nodegree + re74, data = lalonde,
#'                 method = "glm", estimand = "ATE"))
#'
#' summary(W2)
#'
#' bal.tab(W2)
#'
#' #Balancing covariates with respect to re75 (continuous)
#' #with kernel density estimate
#' (W3 <- weightit(re75 ~ age + educ + married +
#'                   nodegree + re74, data = lalonde,
#'                 method = "glm", density = "kernel"))
#'
#' summary(W3)
#'
#' bal.tab(W3)
NULL

weightit2glm <- function(covs, treat, s.weights, subset, estimand, focal,
                         stabilize, missing, .data, verbose, ...) {
  fit.obj <- NULL

  covs <- covs[subset, , drop = FALSE]
  treat <- treat[subset]
  s.weights <- s.weights[subset]

  missing <- .process_missing2(missing, covs)

  re.bars <- ...get(".random")

  if (missing == "saem" && is_not_null(re.bars)) {
    arg::err('random effects are not supported with {.code missing = "saem"}')
  }

  if (missing == "ind") {
    covs <- add_missing_indicators(covs)
  }

  covs <- .make_covs_closer_to_1(covs)

  if (ncol(covs) > 1L) {
    if (missing == "saem") {
      covs0 <- covs
      for (i in colnames(covs)[anyNA_col(covs)]) {
        covs0[is.na(covs0[, i]), i] <- covs0[!is.na(covs0[, i]), i][1L]
      }

      colinear.covs.to.remove <- setdiff(colnames(covs), colnames(make_full_rank(covs0)))
    }
    else {
      colinear.covs.to.remove <- setdiff(colnames(covs), colnames(make_full_rank(covs)))
    }

    covs <- covs[, colnames(covs) %nin% colinear.covs.to.remove, drop = FALSE]
  }

  t.lev <- get_treated_level(treat, estimand, focal)
  treat <- binarize(treat, one = t.lev)

  #Process link
  acceptable.links <- {
    if (missing == "saem") "logit"
    else if (is_not_null(re.bars)) c("logit", "probit", "cloglog", "loglog", "cauchit")
    else c(expand_grid_string(c("", "br."), c("logit", "probit", "cloglog", "loglog", "identity", "log", "clog", "cauchit")),
           "flic", "flac")
  }

  link <- ...get("link")

  if (is_null(link)) {
    link <- acceptable.links[1L]
  }
  else {
    link_match <- pmatch(link, acceptable.links)

    if (anyNA(link_match)) {
      if (missing == "saem") {
        arg::err('only {.val {acceptable.links}} {?is/are} allowed as the link for binary treatments with {.code missing = "saem"}')
      }
      else if (is_not_null(re.bars)) {
        arg::err("only {.val {acceptable.links}} {?is/are} allowed as the link for binary treatments with random effects")
      }
      else {
        arg::err('only {.val {acceptable.links}} {?is/are} allowed as the link for binary treatments')
      }
    }

    link <- acceptable.links[link_match][1L]
  }

  use.br <- startsWith(link, "br.")
  if (use.br) {
    link <- substr(link, 4L, nchar(link))
  }

  use.logistf <- link %in% c("flic", "flac")

  if (missing == "saem") {

    if (!all_the_same(s.weights)) {
      arg::err('sampling weights cannot be used with {.code missing = "saem"}')
    }

    rlang::check_installed("misaem")

    saem.method <- ...get("saem.method", "map")
    control <- ...get("control", list())

    if (is_null(control[["var_cal"]])) control[["var_cal"]] <- TRUE #need TRUE to bypass error in miss.glm.fit()
    if (is_null(control[["ll_obs_cal"]])) control[["ll_obs_cal"]] <- FALSE

    data <- data.frame(treat, covs)

    rlang::try_fetch({
      verbosely({
        fit <- misaem::miss.glm(formula(data), data = data, control = as.list(control))
      }, verbose = verbose)
    },
    warning = function(w) {
      w <- conditionMessage(w)
      if (w != "one argument not used by format '%i '") {
        arg::wrn("(from {.fun misaem::miss.glm}): {w}")
      }
      invokeRestart("muffleWarning")
    },
    error = function(e) {
      e <- conditionMessage(e)
      arg::err("(from {.fun misaem::miss.glm}): {e}")
    })

    p.score <- drop(predict(fit, newdata = covs, method = saem.method))
  }
  else if (use.logistf) {
    rlang::check_installed("logistf")

    fit_fun <- switch(link,
                      flic = logistf::flic,
                      flac = logistf::flac)

    ctrl_fun <- logistf::logistf.control

    control <- do.call(ctrl_fun, c(as.list(...get("control")),
                                   ...mget(setdiff(names(formals(ctrl_fun))[pmatch(...names(), names(formals(ctrl_fun)), 0L)],
                                                   names(...get("control"))))))

    modctrl_fun <- logistf::logistf.mod.control
    modcontrol <- do.call(modctrl_fun, c(as.list(...get("modcontrol")),
                                         ...mget(setdiff(names(formals(modctrl_fun))[pmatch(...names(), names(formals(modctrl_fun)), 0L)],
                                                         names(...get("modctrl_fun"))))))

    rlang::try_fetch({verbosely({
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
      if (w != "non-integer #successes in a binomial glm!") {
        arg::wrn("(from {.fun logistf::{link}}): {w}")
      }
      invokeRestart("muffleWarning")
    },
    error = function(e) {
      arg::err("(from {.fun logistf::{link}}): {conditionMessage(e)}")
    })

    p.score <- fit$predict
  }
  else if (is_not_null(re.bars)) {
    rlang::check_installed("lme4")

    link <- .make_link(link)

    df <- .make_re_data_formula(covs, treat, s.weights, re.bars, .data, subset)

    rlang::try_fetch({verbosely({
      fit <- do.call(lme4::glmer,
                     list(df[["formula"]], data = df[["data"]],
                          family = binomial(link = link),
                          weights = quote(.s.weights)))
    }, verbose = verbose)},
    warning = function(w) {
      w <- conditionMessage(w)
      if (w != "non-integer #successes in a binomial glm!") {
        arg::wrn("(from {.fun lme4::glmer}): {w}")
      }
      invokeRestart("muffleWarning")
    },
    error = function(e) {
      arg::err("(from {.fun lme4::glmer}): {conditionMessage(e)}")
    })

    p.score <- as.numeric(predict(fit, type = "response"))
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

    rlang::try_fetch({verbosely({
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
      if (w != "non-integer #successes in a binomial glm!") {
        arg::wrn("(from {.fun glm}): {w}")
      }

      invokeRestart("muffleWarning")
    },
    error = function(e) {
      arg::err("(from {.fun glm}): {conditionMessage(e)}")
    })

    p.score <- fit$fitted.values

    .psi <- .get_glm_psi(fit)
  }

  if (any(p.score <= 1e-14) || any(p.score >= 1 - 1e-14)) {
    arg::wrn("propensity scores numerically equal to 0 or 1 were estimated, indicating perfect separation. These may yield problems with inference. See {.help [{.code ?method_glm}](WeightIt::method_glm)} for details")
  }

  #Computing weights
  w <- .get_w_from_ps_internal_bin(ps = p.score, treat = treat, estimand,
                                   stabilize = stabilize,
                                   subclass = ...get("subclass"))

  Mparts <- NULL
  if (missing != "saem" && !use.logistf && is_null(...get("subclass")) &&
      !(use.br && identical(fit$type, "correction")) && is_null(re.bars)) {
    Mparts <- list(
      psi_treat = function(Btreat, Xtreat, A, SW) {
        .psi(B = Btreat, X = Xtreat, y = A, weights = SW)
      },
      wfun = function(Btreat, Xtreat, A) {
        XB <- drop(Xtreat %*% Btreat)
        p <- family$linkinv(XB)
        .get_w_from_ps_internal_bin(ps = p, treat = A, estimand,
                                    stabilize = stabilize)
      },
      dw_dBtreat = switch(estimand,
                          ATOS = NULL,
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

  list(w = w, ps = p.score, fit.obj = fit,
       Mparts = Mparts)
}

weightit2glm.cens <- function(covs, treat, s.weights, subset, missing, verbose,
                              estimand = NULL, focal = NULL, stabilize = FALSE, ...) {

  C <- .make_cens_treat(treat)

  out <- .cens_degenerate_out(C[subset])

  if (is_not_null(out)) {
    return(out)
  }

  #The censoring model P(C = 1 | X) is the binary propensity score model with the
  #censored units as the focal group. The GLM score is the same either way, and
  #only the map from the propensity score to the weights differs, so the whole
  #binary machinery (links, `missing = "saem"`, random effects, bias reduction,
  #and the M-estimation parts) is inherited unchanged.
  #
  #`stabilize` must be FALSE: `stabilize_w()` would break the w_ATT + 1 identity
  #that `.att_out_to_cens()` relies on. Stabilized censoring weights are instead
  #formed in `weightit()` from a separate numerator model.
  out <- weightit2glm(covs = covs, treat = C, s.weights = s.weights,
                      subset = subset, estimand = "ATT", focal = 1,
                      stabilize = FALSE, missing = missing,
                      verbose = verbose, ...)

  .att_out_to_cens(out, C[subset])
}

weightit2glm.multi <- function(covs, treat, s.weights, subset, estimand, focal,
                               stabilize, missing, .data, verbose, ...) {
  fit.obj <- NULL

  covs <- covs[subset, , drop = FALSE]
  treat <- factor(treat[subset])
  s.weights <- s.weights[subset]

  missing <- .process_missing2(missing, covs)

  re.bars <- ...get(".random")

  if (missing == "saem" && is_not_null(re.bars)) {
    arg::err('random effects are not supported with {.code missing = "saem"}')
  }

  if (missing == "ind") {
    covs <- add_missing_indicators(covs)
  }

  covs <- .make_covs_closer_to_1(covs)

  if (ncol(covs) > 1L) {
    if (missing == "saem") {
      covs0 <- covs
      for (i in colnames(covs)[anyNA_col(covs)]) {
        covs0[is.na(covs0[, i]), i] <- covs0[!is.na(covs0[, i]), i][1L]
      }

      colinear.covs.to.remove <- setdiff(colnames(covs), colnames(make_full_rank(covs0)))
    }
    else {
      colinear.covs.to.remove <- setdiff(colnames(covs), colnames(make_full_rank(covs)))
    }

    covs <- covs[, colnames(covs) %nin% colinear.covs.to.remove, drop = FALSE]
  }

  for (i in c("use.mlogit", "use.mclogit")) {
    if (is_not_null(...get(i))) {
      arg::wrn('{.arg {i}} is no longer accepted and will be ignored; use {.arg multi.method} instead. See {.help [{.code ?method_glm}](WeightIt::method_glm)} for details')
    }
  }

  ord.treat <- is.ordered(treat) && missing != "saem"

  multi.method <- ...get("multi.method")
  link <- ...get("link")
  link.ignored <- FALSE

  #Random effects (multilevel) PS model: route to mclogit::mblogit()
  if (is_not_null(re.bars)) {
    if (rlang::is_string(multi.method) &&
        tolower(multi.method) %nin% c("mclogit", "mblogit")) {
      arg::wrn('random effects terms in {.arg formula} require {.code multi.method = "mclogit"}; the supplied {.arg multi.method} will be ignored')
    }

    multi.method <- "mclogit"
  }

  #Process multi.method
  if (is_null(multi.method)) {
    if (is_not_null(re.bars)) {
      multi.method <- "mclogit"
    }
    else if (ord.treat) {
      multi.method <- "weightit"
    }
    else {
      multi.method <- {
        if (missing == "saem") "saem"
        else if (identical(link, "bayes.probit")) "mnp"
        else "weightit"
      }
    }
  }
  else if (is_not_null(re.bars)) {
    if (!(rlang::is_string(multi.method) && tolower(multi.method) %in% c("mblogit", "mclogit"))) {
      arg::wrn('random effects terms in {.arg formula} require {.code multi.method = "mclogit"}; the supplied {.arg multi.method} will be ignored')
    }

    multi.method <- "mclogit"
  }
  else if (missing == "saem") {
    if (!identical(multi.method, "saem") &&
        !identical(multi.method, "glm")) {
      arg::wrn('{.arg multi.method} is ignored when {.code missing = "saem"}')
    }

    multi.method <- "saem"
  }
  else {
    arg::arg_string(multi.method)

    if (tolower(multi.method) == "mblogit") {
      multi.method <- "mclogit"
    }

    #Bias reduction is now performed by the in-house fitters, requested with a
    #"br." link, so these no longer name distinct fitting functions
    if (tolower(multi.method) %in% c("brmultinom", "bracl")) {
      arg::wrn(c('{.code multi.method = "{multi.method}"} is deprecated; bias-reduced {if (ord.treat) "ordinal" else "multinomial"} regression is now fit by {.pkg WeightIt} itself, requested with {.code multi.method = "weightit"} and a {.val br.} link.',
                 "i" = if (tolower(multi.method) == "bracl")
                   'This fits a bias-reduced cumulative link model, not the adjacent category logit model {.fun brglm2::bracl} fits, so the estimated propensity scores will differ.'))

      multi.method <- "weightit"

      if (is_null(link)) {
        link <- "br.logit"
      }
    }

    if (ord.treat) {
      allowable.multi.methods <- c("weightit", "polr", "glm", "mclogit", "mnp")
      multi.method <- arg::match_arg(multi.method, allowable.multi.methods)

      ord.treat <- (multi.method %in% c("weightit", "polr"))
    }
    else {
      allowable.multi.methods <- c("weightit", "glm", "mclogit", "mnp")
      multi.method <- arg::match_arg(multi.method, allowable.multi.methods)
    }
  }

  link.ignored <- multi.method %in% c("mclogit", "mnp")

  #Process link. A "br." prefix requests mean bias reduction, available for the
  #in-house multinomial and ordinal fitters.
  acceptable.links <- switch(multi.method,
                             weightit = {
                               if (ord.treat) {
                                 expand_grid_string(c("", "br."),
                                                    c("logit", "probit", "cloglog",
                                                      "loglog", "cauchit"))
                               }
                               else c("logit", "br.logit")
                             },
                             glm = c("logit", "probit", "cloglog", "loglog", "identity",
                                     "log", "clog", "cauchit"),
                             polr = expand_grid_string(c("", "br."),
                                                       c("logit", "probit", "loglog",
                                                         "cloglog", "cauchit")),
                             mnp = c("bayes.probit", "probit"),
                             mclogit = "logit",
                             saem = "logit")

  if (is_null(link)) {
    link <- acceptable.links[1L]
  }
  else if (link.ignored) {
    if (!any_apply(acceptable.links, identical, link)) {
      arg::wrn('{.arg link} is ignored when {.code {if (missing == "saem") "missing" else "multi.method"} = "{multi.method}"}')
    }

    link <- acceptable.links[1L]
  }
  else {
    link_match <- pmatch(link, acceptable.links)

    if (anyNA(link_match)) {
      arg::err('only {.val {acceptable.links}} {?is/are} allowed as the link for {if (ord.treat) "ordinal"} multi-category treatments with {.code multi.method = "{multi.method}"}')
    }

    link <- acceptable.links[link_match][1L]
  }

  use.br <- startsWith(link, "br.")
  if (use.br) {
    link <- substr(link, 4L, nchar(link))

    #`MASS::polr()` cannot do bias reduction, so use the in-house cumulative link
    #fitter, which can (and which additionally supports M-estimation)
    if (multi.method == "polr") {
      multi.method <- "weightit"
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
                                         link = link,
                                         br = use.br)
      }, verbose = verbose)
    }
    else {
      verbosely({
        fit.obj <- .multinom_weightit.fit(x = cbind(`(Intercept)` = 1, covs),
                                          y = treat,
                                          weights = s.weights,
                                          hess = FALSE,
                                          br = use.br)
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

      rlang::try_fetch({verbosely({
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
          arg::wrn("(from {.fun glm}: {w}")
        }

        invokeRestart("muffleWarning")
      },
      error = function(w) {
        arg::err("(from {.fun glm}: {conditionMessage(e)}")
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

    rlang::try_fetch({verbosely({
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
      arg::err("there was a problem fitting the ordinal {link} regressions with {.fun polr}. Try again with an un-ordered treatment.Error message: (from {.fun MASS::polr}):\f{conditionMessage(e)}")
    })

    ps <- fit.obj$fitted.values
  }
  else if (multi.method == "saem") {

    if (!all_the_same(s.weights)) {
      arg::err('sampling weights cannot be used with {.code missing = "saem"}')
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

      rlang::try_fetch({
        verbosely({
          fit.obj[[i]] <- misaem::miss.glm(formula(data_i), data = data_i,
                                           control = as.list(control))
        }, verbose = verbose)
      },
      warning = function(w) {
        w <- conditionMessage(w)
        if (w != "one argument not used by format '%i '") {
          arg::wrn("(from {.fun misaem::miss.glm}): {w}")
        }

        invokeRestart("muffleWarning")
      },
      error = function(e) {
        arg::err("(from {.fun misaem::miss.glm}): {conditionMessage(e)}")
      })

      ps[[i]] <- drop(predict(fit.obj[[i]], newdata = covs, method = saem.method))
    }
  }
  else if (multi.method == "mnp") {
    rlang::check_installed("MNP")

    data <- data.frame(treat, covs)
    formula <- formula(data)

    rlang::try_fetch({verbosely({
      fit.obj <- MNP::mnp(formula, data, verbose = TRUE)
    }, verbose = verbose)},
    error = function(e) {
      arg::err("there was a problem fitting the multinomial regression with {.fun MNP}. Try a different {.arg multi.method}. Error message: (from {.fun MNP::mnp}):\f{conditionMessage(e)}")
    })

    ps <- predict(fit.obj, type = "prob")$p
  }
  else if (multi.method == "mclogit") {
    rlang::check_installed("mclogit")

    #Random effects from formula bars (.random) take precedence over the
    #legacy `random` argument (an mclogit-style formula supplied directly).
    random.arg <- {
      if (is_not_null(re.bars)) .bars_to_mclogit_random(re.bars)
      else ...get("random")
    }

    if (is_not_null(random.arg)) {
      re.vars <- unique(unlist(lapply(
        if (rlang::is_formula(random.arg)) list(random.arg) else as.list(random.arg),
        get_varnames)))

      not.found <- setdiff(re.vars, names(.data))
      if (is_not_null(not.found)) {
        loc <- if (is_not_null(re.bars)) "formula" else "random"

        arg::err("the variable{?s} {.var {not.found}} in the random effects component of {.arg {loc}} {?was/were} not found in the dataset")
      }

      re.data <- as.data.frame(.data)[subset, re.vars, drop = FALSE]

      if (anyNA(re.data)) {
        arg::err("missing values are not allowed in the random effects grouping or slope variables")
      }

      data <- data.frame(re.data, treat = treat, .s.weights = s.weights, covs)
      covnames <- names(data)[-seq_len(ncol(re.data) + 2L)]
      tname <- names(data)[ncol(re.data) + 1L]
      ctrl_fun <- mclogit::mmclogit.control
    }
    else {
      data <- data.frame(treat = treat, .s.weights = s.weights, covs)
      covnames <- names(data)[-(1:2)]
      tname <- names(data)[1L]
      ctrl_fun <- mclogit::mclogit.control
    }

    form <- reformulate(covnames, tname)

    control <- do.call(ctrl_fun,
                       c(as.list(...get("control")),
                         ...mget(setdiff(rlang::fn_fmls_names(ctrl_fun)[pmatch(...names(), rlang::fn_fmls_names(ctrl_fun), 0L)],
                                         names(...get("control"))))))

    rlang::try_fetch({verbosely({
      fit.obj <- do.call(mclogit::mblogit,
                         list(form,
                              data = data,
                              weights = quote(.s.weights),
                              random = random.arg,
                              method = "PQL",
                              estimator = ...get("estimator", eval(formals(mclogit::mclogit)[["estimator"]])),
                              dispersion = ...get("dispersion", eval(formals(mclogit::mclogit)[["dispersion"]])),
                              groups = ...get("groups"),
                              control = control))
    }, verbose = verbose)},
    error = function(e) {
      arg::err("there was a problem fitting the multinomial {link} regression with {.fun mblogit}. Try a different {.arg multi.method}. Error message: (from {.fun mclogit::mblogit}):\f{conditionMessage(e)}")
    })

    #fitted() does not work for random-effects (mmblogit) fits; predict() with
    #conditional = TRUE (the default) returns cluster-specific probabilities.
    ps <- {
      if (is_not_null(random.arg)) predict(fit.obj, type = "response")
      else fitted(fit.obj)
    }

    colnames(ps) <- levels(treat)
  }

  #ps should be matrix of probs for each treat
  #Computing weights
  w <- .get_w_from_ps_internal_multi(ps = ps, treat = treat, estimand, focal = focal,
                                     stabilize = stabilize,
                                     subclass = ...get("subclass"))

  #Get Mparts
  Mparts <- NULL
  if (is_null(...get("subclass"))) {
    if (multi.method == "weightit") {
      Mparts <- list(
        psi_treat = function(Btreat, Xtreat, A, SW) {
          fit.obj$psi(Btreat, Xtreat, A, SW)
        },
        wfun = function(Btreat, Xtreat, A) {
          ps <- fit.obj$get_p(Btreat, Xtreat)
          .get_w_from_ps_internal_multi(ps = ps, treat = A, estimand, focal = focal,
                                        stabilize = stabilize)
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
                                        stabilize = stabilize)
        },
        Xtreat = cbind(`(Intercept)` = 1, covs),
        A = treat,
        btreat = unlist(grab(fit.obj, "coefficients"))
      )
    }
  }

  list(w = w, fit.obj = fit.obj,
       Mparts = Mparts)
}

weightit2glm.cont <- function(covs, treat, s.weights, subset, stabilize, missing, .data, verbose, ...) {

  covs <- covs[subset, , drop = FALSE]
  treat <- treat[subset]
  s.weights <- s.weights[subset]

  missing <- .process_missing2(missing, covs)

  re.bars <- ...get(".random")

  #With nothing to condition on, the conditional density of the treatment is its
  #marginal density and the weights are exactly 1. Fitting the model would instead
  #send them through the numeric marginalization in
  #`.get_w_from_gps_internal_cont()`, which evaluates the density on a grid and so
  #returns values near but not equal to 1. Nothing is estimated, so there are no
  #M-estimation parts either. A formula whose only terms are random effects also has
  #a zero-column `covs`, but is a real model, hence the `re.bars` guard.
  if (ncol(covs) == 0L && is_null(re.bars)) {
    return(list(w = rep_with(1, treat)))
  }

  if (missing == "saem" && is_not_null(re.bars)) {
    arg::err('random effects are not supported with {.code missing = "saem"}')
  }

  if (missing == "ind") {
    covs <- add_missing_indicators(covs)
  }

  covs <- .make_covs_closer_to_1(covs)

  if (ncol(covs) > 1L) {
    if (missing == "saem") {
      covs0 <- covs
      for (i in colnames(covs)[anyNA_col(covs)]) {
        covs0[is.na(covs0[, i]), i] <- covs0[!is.na(covs0[, i]), i][1L]
      }
      colinear.covs.to.remove <- setdiff(colnames(covs), colnames(make_full_rank(covs0)))
    }
    else {
      colinear.covs.to.remove <- setdiff(colnames(covs), colnames(make_full_rank(covs)))
    }

    covs <- covs[, colnames(covs) %nin% colinear.covs.to.remove, drop = FALSE]
  }

  s.weights <- s.weights / mean_fast(s.weights)

  # Process density params
  make_dens_fun <- .get_make_dens_fun(density = ...get("density"),
                                      bw = ...get("bw"),
                                      adjust = ...get("adjust"),
                                      kernel = ...get("kernel"),
                                      n = ...get("n"),
                                      use.kernel = ...get("use.kernel"))

  is_kernel <- isTRUE(attr(make_dens_fun, "is_kernel"))

  #Estimate GPS
  link <- ...get("link", "identity")

  if (is_not_null(re.bars)) {
    if (missing == "saem") {
      arg::err('random effects are not supported with {.code missing = "saem"}')
    }

    rlang::check_installed("lme4")

    acceptable.links <- c("identity", "log", "inverse")
    link_match <- pmatch(link, acceptable.links)

    if (anyNA(link_match)) {
      arg::err("only {.val {acceptable.links}} {?is/are} allowed as the link for continuous treatments with random effects")
    }

    link <- acceptable.links[link_match][1L]

    df <- .make_re_data_formula(covs, treat, s.weights, re.bars, .data, subset)

    if (link == "identity") {
      ctrl_fun <- lme4::lmerControl

      control <- do.call(ctrl_fun,
                         c(as.list(...get("control")),
                           ...mget(setdiff(rlang::fn_fmls_names(ctrl_fun)[pmatch(...names(), rlang::fn_fmls_names(ctrl_fun), 0L)],
                                           names(...get("control"))))))

      rlang::try_fetch({verbosely({
        fit <- do.call(lme4::lmer,
                       list(df[["formula"]], data = df[["data"]],
                            weights = quote(.s.weights),
                            control = control))
      }, verbose = verbose)},
      warning = function(w) {
        arg::wrn("(from {.fun lme4::lmer}): {conditionMessage(w)}")
        invokeRestart("muffleWarning")
      },
      error = function(e) {
        arg::err("(from {.fun lme4::lmer}): {conditionMessage(e)}")
      })

      mu <- as.numeric(predict(fit))
      sd <- sigma(fit)
    }
    else {
      ctrl_fun <- lme4::glmerControl

      control <- do.call(ctrl_fun,
                         c(as.list(...get("control")),
                           ...mget(setdiff(rlang::fn_fmls_names(ctrl_fun)[pmatch(...names(), rlang::fn_fmls_names(ctrl_fun), 0L)],
                                           names(...get("control"))))))

      rlang::try_fetch({verbosely({
        fit <- do.call(lme4::glmer,
                       list(df[["formula"]], data = df[["data"]],
                            family = gaussian(link),
                            weights = quote(.s.weights),
                            control = control))
      }, verbose = verbose)},
      warning = function(w) {
        arg::wrn("(from {.fun lme4::glmer}): {conditionMessage(w)}")
        invokeRestart("muffleWarning")
      },
      error = function(e) {
        arg::err("(from {.fun lme4::glmer}): {conditionMessage(e)}")
      })

      mu <- as.numeric(predict(fit, type = "response"))
      sd <- sigma(fit)
    }
  }
  else if (missing == "saem") {

    if (!all_the_same(s.weights)) {
      arg::err('sampling weights cannot be used with {.code missing = "saem"}')
    }

    rlang::check_installed("misaem")

    acceptable.links <- "identity"

    link_match <- pmatch(link, acceptable.links)

    if (anyNA(link_match)) {
      arg::err('only {.val {acceptable.links}} {?is/are} allowed as the link for continuous treatments with {.code missing = "saem"}')
    }

    link <- acceptable.links[link_match][1L]

    data <- data.frame(treat, covs)
    formula <- formula(data)

    rlang::try_fetch({verbosely({
      fit <- misaem::miss.lm(formula, data = data, control = as.list(...get("control")))
    }, verbose = verbose)},
    warning = function(w) {
      w <- conditionMessage(w)
      if (w != "one argument not used by format '%i '") {
        arg::wrn("(from {.fun misaem::miss.lm}): {w}")
      }

      invokeRestart("muffleWarning")
    },
    error = function(e) {
      arg::err("(from {.fun misaem::miss.lm}): {conditionMessage(e)}")
    })

    saem.method <- ...get("saem.method", "map")

    mu <- drop(predict(fit, newdata = covs, method = saem.method))
    sd <- fit[["s.resid"]]
  }
  else {
    acceptable.links <- c("identity", "log", "inverse")

    link_match <- pmatch(link, acceptable.links)

    if (anyNA(link_match)) {
      arg::err("only {.val {acceptable.links}} {?is/are} allowed as the link for continuous treatments")
    }

    link <- acceptable.links[link_match][1L]

    family <- gaussian(link = link)

    verbosely({
      if (isTRUE(...get("quick"))) {
        fit <- do.call(stats::glm.fit, list(y = treat,
                                            x = cbind(`(Intercept)` = rep_with(1, treat), covs),
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

    mu <- fit$fitted.values

    sd <- NULL #uses formula consistent with M-estimation
    # sd <- sigma(fit) #extracted from model
  }

  #Get weights
  w <- .get_w_from_gps_internal_cont(mu, treat, sd, s.weights,
                                     make_dens_fun)

  Mparts <- NULL
  if (missing != "saem" && !is_kernel && is_null(re.bars)) {
    Mparts <- list(
      psi_treat = function(Btreat, Xtreat, A, SW) {
        s2 <- exp(Btreat[1L])
        lin_pred <- drop(Xtreat %*% Btreat[-1L])
        mu <- family$linkinv(lin_pred)

        SW <- SW / mean_fast(SW)

        cbind(SW * (A - mu)^2 - s2, #conditional variance
              Xtreat * (SW * family$mu.eta(lin_pred) * (A - mu) / family$variance(mu))) #conditional mean
      },
      wfun = function(Btreat, Xtreat, A) {
        sd <- sqrt(exp(Btreat[1L]))
        lin_pred <- drop(Xtreat %*% Btreat[-1L])
        mu <- family$linkinv(lin_pred)

        .get_w_from_gps_internal_cont(mu, A, sd, s.weights,
                                      make_dens_fun)
      },
      Xtreat = cbind(`(Intercept)` = 1, covs),
      A = treat,
      btreat = c("log(s_r^2)" = log(sd^2),
                 fit$coefficients)
    )
  }

  list(w = w, fit.obj = fit,
       Mparts = Mparts)
}

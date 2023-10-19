#' Propensity Score Weighting Using Generalized Linear Models
#' @name method_glm
#' @aliases method_glm method_ps
#' @usage NULL
#'
#' @description
#' This page explains the details of estimating weights from generalized linear model-based propensity scores by setting `method = "glm"` in the call to [weightit()] or [weightitMSM()]. This method can be used with binary, multinomial, and continuous treatments. (This method used to be requested with `method = "ps"`, and this still works.)
#'
#' In general, this method relies on estimating propensity scores with a parametric generalized linear model and then converting those propensity scores into weights using a formula that depends on the desired estimand. For binary and multinomial treatments, a binomial or multinomial regression model is used to estimate the propensity scores as the predicted probability of being in each treatment given the covariates. For ordinal treatments, an ordinal regression model is used to estimate generalized propensity scores. For continuous treatments, a generalized linear model is used to estimate generalized propensity scores as the conditional density of treatment given the covariates.
#'
#' ## Binary Treatments
#'
#' For binary treatments, this method estimates the propensity scores using [glm()]. An additional argument is `link`, which uses the same options as `link` in [family()]. The default link is "logit", but others, including "probit", are allowed. The following estimands are allowed: ATE, ATT, ATC, ATO, ATM, and ATOS. Weights can also be computed using marginal mean weighting through stratification for the ATE, ATT, and ATC. See [get_w_from_ps()] for details.
#'
#' ## Multinomial Treatments
#'
#' For multinomial treatments, the propensity scores are estimated using multinomial regression from one of a few functions depending on the requested link: for logit (`"logit"`) and probit (`"probit"`) links, \pkgfun{mlogit}{mlogit} is used; for the Bayesian probit (`"bayes.probit"`) link, \pkgfun{MNP}{mnp} is used; and for the biased-reduced multinomial logistic regression (`"br.logit"`), \pkgfun{brglm2}{brmultinom} is used. If the treatment variable is an ordered factor, \pkgfun{MASS}{polr} is used to fit ordinal regression unless `link = "br.logit"`, in which case \pkgfun{brglm2}{bracl} is used. Any of the methods allowed in the `method` argument of `polr()` can be supplied to `link`. The following estimands are allowed: ATE, ATT, ATC, ATO, and ATM. The weights for each estimand are computed using the standard formulas or those mentioned above. Weights can also be computed using marginal mean weighting through stratification for the ATE, ATT, and ATC. See [get_w_from_ps()] for details.
#'
#' ## Continuous Treatments
#'
#' For continuous treatments, the generalized propensity score is estimated using linear regression. The conditional density can be specified as normal or another distribution. In addition, kernel density estimation can be used instead of assuming a specific density for the numerator and denominator of the generalized propensity score by setting `use.kernel = TRUE`. Other arguments to [density()] can be specified to refine the density estimation parameters. `plot = TRUE` can be specified to plot the density for the numerator and denominator, which can be helpful in diagnosing extreme weights.
#'
#' ## Longitudinal Treatments
#'
#' For longitudinal treatments, the weights are the product of the weights estimated at each time point.
#'
#' ## Sampling Weights
#'
#' Sampling weights are supported through `s.weights` in all scenarios except for multinomial treatments with `link = "bayes.probit"` and for binary and continuous treatments with `missing = "saem"` (see below). Warning messages may appear otherwise about non-integer successes, and these can be ignored.
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
#' @section Additional Arguments:
#' The following additional arguments can be specified:
#' \describe{
#'   \item{`link`}{The link used in the generalized linear model for the propensity scores. For binary treatments, `link` can be any of those allowed by [binomial()]. A `br.` prefix can be added (e.g., `"br.logit"`); this changes the fitting method to the \pkgfun2{brglm2}{brglmFit}{bias-corrected generalized linear models} implemented in the \CRANpkg{brglm2} package. For multicategory treatments, `link` can be `"logit"`, `"probit"`, `"bayes.probit"`, or `"br.logit"`. For ordered treatments, `link` can be any of those allowed by the `method` argument of \pkgfun{MASS}{polr} or `"br.logit"`. For continuous treatments, `link` can be any of those allowed by [gaussian()].
#'   }
#' }
#' For continuous treatments only, the following arguments may be supplied:
#'   \describe{
#'     \item{`density`}{A function corresponding the conditional density of the treatment. The standardized residuals of the treatment model will be fed through this function to produce the numerator and denominator of the generalized propensity score weights. If blank, [dnorm()] is used as recommended by Robins et al. (2000). This can also be supplied as a string containing the name of the function to be called. If the string contains underscores, the call will be split by the underscores and the latter splits will be supplied as arguments to the second argument and beyond. For example, if `density = "dt_2"` is specified, the density used will be that of a t-distribution with 2 degrees of freedom. Using a t-distribution can be useful when extreme outcome values are observed (Naimi et al., 2014). Ignored if `use.kernel = TRUE` (described below).
#'     }
#'     \item{`use.kernel`}{If `TRUE`, uses kernel density estimation through the [density()] function to estimate the numerator and denominator densities for the weights. If `FALSE`, the argument to the `density` parameter is used instead.
#'     }
#'     \item{`bw`, `adjust`, `kernel`, `n`}{If `use.kernel = TRUE`, the arguments to the [density()] function. The defaults are the same as those in `density` except that `n` is 10 times the number of units in the sample.
#'     }
#'     \item{`plot`}{If `use.kernel = TRUE` with continuous treatments, whether to plot the estimated density.
#'     }
#'   }
#'
#' For binary treatments, additional arguments to `glm()` can be specified as well. The `method` argument in `glm()` is renamed to `glm.method`. This can be used to supply alternative fitting functions, such as those implemented in the \CRANpkg{glm2} package. Other arguments to `weightit()` are passed to `...` in `glm()`. In the presence of missing data with `link = "logit"` and `missing = "saem"`, additional arguments are passed to \pkgfun{misaem}{miss.glm} and \pkgfun{misaem}{predict.miss.glm}, except the `method` argument in \pkgfun{misaem}{predict.miss.glm} is replaced with `saem.method`.
#'
#' For multi-category treatments with `link = "logit"` or `"probit"`, the default is to use multinomial logistic or probit regression using the \CRANpkg{mlogit} package. To request that separate binary logistic or probit regressions are run instead, set `use.mlogit = FALSE`. This can be helpful when `mlogit` is slow or fails to converge. With `link = "logit"`, the option `use.mclogit = TRUE` can be specified to request that \pkgfun{mclogit}{mblogit} from the \CRANpkg{mclogit} package is used instead, which can be faster and is recommended.
#'
#' For continuous treatments in the presence of missing data with `missing = "saem"`, additional arguments are passed to \pkgfun{misaem}{miss.lm} and \pkgfun{misaem}{predict.miss.lm}.
#'
#' @section Additional Outputs:
#' \describe{
#'   \item{`obj`}{When `include.obj = TRUE`, the (generalized) propensity score model fit. For binary treatments, the output of the call to [glm()]. For ordinal treatments, the output of the call to \pkgfun{MASS}{polr}. For multinomial treatments with `link = "logit"` or `"probit"` and `use.mlogit = TRUE`, the output of the call to \pkgfun{mlogit}{mlogit}. For multinomial treatments with `use.mlogit = FALSE`, a list of the `glm()` fits. For multinomial treatments with `link = "br.logit"`, the output of the call to \pkgfun{brglm2}{brmultinom}. For multinomial treatments with `link = "bayes.probit"`, the output of the call to \pkgfun{MNP}{mnp}. For continuous treatments, the output of the call to `glm()` for the predicted values in the denominator density.
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
#' Li, F., Morgan, K. L., & Zaslavsky, A. M. (2018). Balancing covariates via propensity score weighting. Journal of the American Statistical Association, 113(521), 390–400. \doi{10.1080/01621459.2016.1260466}
#'
#' - `estimand = "ATM"`
#'
#' Li, L., & Greene, T. (2013). A Weighting Analogue to Pair Matching in Propensity Score Analysis. The International Journal of Biostatistics, 9(2). \doi{10.1515/ijb-2012-0030}
#'
#' - `estimand = "ATOS"`
#'
#' Crump, R. K., Hotz, V. J., Imbens, G. W., & Mitnik, O. A. (2009). Dealing with limited overlap in estimation of average treatment effects. Biometrika, 96(1), 187–199. \doi{10.1093/biomet/asn055}
#'
#' - Other estimands
#'
#' Austin, P. C. (2011). An Introduction to Propensity Score Methods for Reducing the Effects of Confounding in Observational Studies. Multivariate Behavioral Research, 46(3), 399–424. \doi{10.1080/00273171.2011.568786}
#'
#' - Marginal mean weighting through stratification
#'
#' Hong, G. (2010). Marginal mean weighting through stratification: Adjustment for selection bias in multilevel data. Journal of Educational and Behavioral Statistics, 35(5), 499–531. \doi{10.3102/1076998609359785}
#'
#' - Bias-reduced logistic regression
#'
#' See references for the \pkg{brglm2} \pkgfun2{brglm2}{brglm2}{package}.
#'
#' - SAEM logistic regression for missing data
#'
#' Jiang, W., Josse, J., & Lavielle, M. (2019). Logistic regression with missing covariates — Parameter estimation, model selection and prediction within a joint-modeling framework. Computational Statistics & Data Analysis, 106907. \doi{10.1016/j.csda.2019.106907}
#'
#' ## Multinomial Treatments
#'
#' - `estimand = "ATO"`
#'
#' Li, F., & Li, F. (2019). Propensity score weighting for causal inference with multiple treatments. The Annals of Applied Statistics, 13(4), 2389–2415. \doi{10.1214/19-AOAS1282}
#'
#' - `estimand = "ATM"`
#'
#' Yoshida, K., Hernández-Díaz, S., Solomon, D. H., Jackson, J. W., Gagne, J. J., Glynn, R. J., & Franklin, J. M. (2017). Matching weights to simultaneously compare three treatment groups: Comparison to three-way matching. Epidemiology (Cambridge, Mass.), 28(3), 387–395. \doi{10.1097/EDE.0000000000000627}
#'
#' - Other estimands
#'
#' McCaffrey, D. F., Griffin, B. A., Almirall, D., Slaughter, M. E., Ramchand, R., & Burgette, L. F. (2013). A Tutorial on Propensity Score Estimation for Multiple Treatments Using Generalized Boosted Models. Statistics in Medicine, 32(19), 3388–3414. \doi{10.1002/sim.5753}
#'
#' - Marginal mean weighting through stratification
#'
#' Hong, G. (2012). Marginal mean weighting through stratification: A generalized method for evaluating multivalued and multiple treatments with nonexperimental data. Psychological Methods, 17(1), 44–60. \doi{10.1037/a0024918}
#'
#' ## Continuous treatments
#'
#' Robins, J. M., Hernán, M. Á., & Brumback, B. (2000). Marginal Structural Models and Causal Inference in Epidemiology. Epidemiology, 11(5), 550–560.
#'
#' - Using non-normal conditional densities
#'
#' Naimi, A. I., Moodie, E. E. M., Auger, N., & Kaufman, J. S. (2014). Constructing Inverse Probability Weights for Continuous Exposures: A Comparison of Methods. Epidemiology, 25(2), 292–299. \doi{10.1097/EDE.0000000000000053}
#'
#' - SAEM linear regression for missing data
#'
#' Jiang, W., Josse, J., & Lavielle, M. (2019). Logistic regression with missing covariates — Parameter estimation, model selection and prediction within a joint-modeling framework. Computational Statistics & Data Analysis, 106907. \doi{10.1016/j.csda.2019.106907}
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
#' #Balancing covariates with respect to race (multinomial)
#' (W2 <- weightit(race ~ age + educ + married +
#'                   nodegree + re74, data = lalonde,
#'                 method = "glm", estimand = "ATE",
#'                 use.mlogit = FALSE))
#' summary(W2)
#' bal.tab(W2)
#'
#' #Balancing covariates with respect to re75 (continuous)
#' (W3 <- weightit(re75 ~ age + educ + married +
#'                   nodegree + re74, data = lalonde,
#'                 method = "glm", use.kernel = TRUE))
#' summary(W3)
#' bal.tab(W3)
NULL

weightit2glm <- function(covs, treat, s.weights, subset, estimand, focal,
                         stabilize, subclass, missing, ps, .data, verbose, ...) {
  A <- list(...)

  fit.obj <- NULL

  if (is_null(ps)) {

    covs <- covs[subset, , drop = FALSE]
    treat_sub <- factor(treat[subset])
    s.weights <- s.weights[subset]

    if (missing == "ind") {
      covs <- add_missing_indicators(covs)
    }

    for (i in seq_col(covs)) covs[,i] <- make.closer.to.1(covs[,i])

    if (ncol(covs) > 1) {
      if (missing == "saem") {
        covs0 <- covs
        for (i in colnames(covs)[anyNA_col(covs)]) covs0[is.na(covs0[,i]),i] <- covs0[!is.na(covs0[,i]),i][1]
        colinear.covs.to.remove <- colnames(covs)[colnames(covs) %nin% colnames(make_full_rank(covs0))]
      }
      else colinear.covs.to.remove <- colnames(covs)[colnames(covs) %nin% colnames(make_full_rank(covs))]
      covs <- covs[, colnames(covs) %nin% colinear.covs.to.remove, drop = FALSE]
    }

    #Process link
    acceptable.links <- {
      if (missing == "saem") "logit"
      else expand.grid_string(c("", "br."), c("logit", "probit", "cloglog", "identity", "log", "cauchit"))
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

    if (missing == "saem") {
      rlang::check_installed("misaem")
      if (is_null(A[["saem.method"]])) A[["saem.method"]] <- "map"
      if (is_null(A[["control"]])) A[["control"]] <- list()
      if (is_null(A[["control"]][["var_cal"]])) A[["control"]][["var_cal"]] <- FALSE
      if (is_null(A[["control"]][["ll_obs_cal"]])) A[["control"]][["ll_obs_cal"]] <- FALSE
    }

    t.lev <- get_treated_level(treat_sub)
    c.lev <- setdiff(levels(treat_sub), t.lev)

    ps <- make_df(levels(treat_sub), nrow = length(treat_sub))

    treat <- binarize(treat_sub, one = t.lev)

    if (missing == "saem") {
      data <- data.frame(treat, covs)

      withCallingHandlers({
        verbosely({
          fit <- misaem::miss.glm(formula(data), data = data, control = as.list(A[["control"]]))
        }, verbose = verbose)
      },
      warning = function(w) {
        if (conditionMessage(w) != "one argument not used by format '%i '") .wrn("(from misaem) ", w, tidy = FALSE)
        invokeRestart("muffleWarning")
      })

      ps[[t.lev]] <- p.score <- drop(predict(fit, newdata = covs, method = A[["saem.method"]]))
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
        start <- c(family$linkfun(w.m(treat, s.weights)), rep(0, ncol(covs)))
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
          formula <- formula(data)
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
        if (conditionMessage(w) != "non-integer #successes in a binomial glm!") .wrn("(from `glm()`) ", w, tidy = FALSE)
        invokeRestart("muffleWarning")
      })

      ps[[t.lev]] <- p.score <- fit$fitted.values
      ps[[c.lev]] <- 1 - ps[[t.lev]]
    }

    fit[["call"]] <- NULL
    fit.obj <- fit

  }
  else {
    n <- length(treat)
    p.score <- NULL
    treat <- factor(treat)
    treat_sub <- factor(treat[subset])

    t.lev <- get_treated_level(treat)
    c.lev <- setdiff(levels(treat_sub), t.lev)

    if (is_(ps, c("matrix", "data.frame"))) {
      if (all(dim(ps) == c(n, 2))) {

        if (all(colnames(ps) %in% levels(treat_sub))) {
          ps <- as.data.frame(ps[subset, , drop = FALSE])
        }
        else {
          ps <- as.data.frame(ps[subset, , drop = FALSE])
          names(ps) <- levels(treat_sub)
        }

        p.score <- ps[[t.lev]]
      }
      else if (all(dim(ps) == c(length(treat), 1))) {

        ps <- data.frame(ps[subset,1], 1-ps[subset,1])

        names(ps) <- c(t.lev, c.lev)

        p.score <- ps[[t.lev]]
      }
    }
    else if (is.numeric(ps)) {
      if (length(ps) == n) {
        ps <- data.frame(ps[subset], 1-ps[subset])

        names(ps) <- c(t.lev, c.lev)

        p.score <- ps[[t.lev]]
      }
    }

    if (is_null(p.score)) .err("`ps` must be a numeric vector with a propensity score for each unit")

  }

  #ps should be matrix of probs for each treat
  #Computing weights
  w <- get_w_from_ps(ps = ps, treat = treat_sub, estimand, focal = focal,
                     stabilize = stabilize, subclass = subclass)

  list(w = w, ps = p.score, fit.obj = fit.obj)
}

weightit2glm.multi <- function(covs, treat, s.weights, subset, estimand, focal,
                               stabilize, subclass, missing, ps, .data, verbose, ...) {
  A <- list(...)

  fit.obj <- NULL

  if (is_null(ps)) {

    covs <- covs[subset, , drop = FALSE]
    treat_sub <- factor(treat[subset])
    s.weights <- s.weights[subset]

    ord.treat <- is.ordered(treat_sub)

    if (missing == "ind") {
      covs <- add_missing_indicators(covs)
    }

    for (i in seq_col(covs)) covs[,i] <- make.closer.to.1(covs[,i])

    if (ncol(covs) > 1) {
      if (missing == "saem") {
        covs0 <- covs
        for (i in colnames(covs)[anyNA_col(covs)]) covs0[is.na(covs0[,i]),i] <- covs0[!is.na(covs0[,i]),i][1]
        colinear.covs.to.remove <- colnames(covs)[colnames(covs) %nin% colnames(make_full_rank(covs0))]
      }
      else colinear.covs.to.remove <- colnames(covs)[colnames(covs) %nin% colnames(make_full_rank(covs))]
      covs <- covs[, colnames(covs) %nin% colinear.covs.to.remove, drop = FALSE]
    }

    #Process link
    acceptable.links <- {
      if (ord.treat) c("logit", "probit", "loglog", "cloglog", "cauchit", "br.logit")
      else if (!isFALSE(A$use.mlogit)) c("logit", "probit", "bayes.probit", "br.logit")
      else if (missing == "saem") "logit"
      else expand.grid_string(c("", "br."), c("logit", "probit", "cloglog", "identity", "log", "cauchit"))
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

    if (missing == "saem") {
      rlang::check_installed("misaem")
      if (is_null(A[["saem.method"]])) A[["saem.method"]] <- "map"
      if (is_null(A[["control"]])) A[["control"]] <- list()
      if (is_null(A[["control"]][["var_cal"]])) A[["control"]][["var_cal"]] <- FALSE
      if (is_null(A[["control"]][["ll_obs_cal"]])) A[["control"]][["ll_obs_cal"]] <- FALSE
    }

    if (ord.treat) {
      if (use.br) {
        rlang::check_installed("brglm2")

        ctrl_fun <- brglm2::brglmControl
        control <- do.call(ctrl_fun, c(A[["control"]],
                                       A[setdiff(names(formals(ctrl_fun))[pmatch(names(A), names(formals(ctrl_fun)), 0)],
                                                 names(A[["control"]]))]))

        data <- data.frame(treat_sub, covs)
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
          .err(sprintf("there was a problem with the bias-reduced ordinal logit regression.\n       Try a different link.\n       Error message: %s", conditionMessage(e)), tidy = FALSE)
        })
      }
      else {
        if (A[["link"]] == "logit") A[["link"]] <- "logistic"
        rlang::check_installed("MASS")
        # message(paste("Using ordinal", A$link, "regression."))
        data <- data.frame(treat_sub, covs)
        formula <- formula(data)

        tryCatch({verbosely({
          fit <- do.call(MASS::polr,
                         list(formula,
                              data = data,
                              weights = s.weights,
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
    else {
      if (use.br) {
        rlang::check_installed("brglm2")
        data <- data.frame(treat_sub, covs)
        formula <- formula(data)
        ctrl_fun <- brglm2::brglmControl
        control <- do.call(ctrl_fun, c(A[["control"]],
                                       A[setdiff(names(formals(ctrl_fun))[pmatch(names(A), names(formals(ctrl_fun)), 0)],
                                                 names(A[["control"]]))]))
        tryCatch({verbosely({
          fit <- do.call(brglm2::brmultinom,
                         list(formula, data,
                              weights = s.weights,
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
        data <- data.frame(treat_sub, covs)
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
        ps <- make_df(levels(treat_sub), nrow = length(treat_sub))

        fit.list <- make_list(levels(treat_sub))

        for (i in levels(treat_sub)) {
          t_i <- as.numeric(treat_sub == i)
          data_i <- data.frame(t_i, covs)

          withCallingHandlers({
            verbosely({
              fit.list[[i]] <- misaem::miss.glm(formula(data_i), data = data_i,
                                                control = A[["control"]])
            }, verbose = verbose)
          },
          warning = function(w) {
            if (conditionMessage(w) != "one argument not used by format '%i '") .wrn("(from misaem) ", w, tidy = FALSE)
            invokeRestart("muffleWarning")
          })

          ps[[i]] <- drop(predict(fit, newdata = covs, method = A[["saem.method"]]))
        }
        fit.obj <- fit.list
      }
      else if (isTRUE(A$use.mclogit)) {
        rlang::check_installed("mclogit")

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

        data <- data.frame(treat = treat_sub, .s.weights = s.weights, covs)
        covnames <- names(data)[-c(1,2)]
        tryCatch({verbosely({
          fit <- mlogit::mlogit(as.formula(paste0("treat ~ 1 | ", paste(covnames, collapse = " + "))),
                                data = data,
                                estimate = TRUE,
                                probit = A$link[1] == "probit",
                                weights = .s.weights,
                                varying = NULL,
                                shape = "wide",
                                sep = "",
                                choice = "treat",
                                ...)
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

          t_i <- as.numeric(treat_sub == i)
          data_i <- data.frame(t_i, covs)

          verbosely({
            fit.list[[i]] <- do.call(stats::glm, list(formula(data_i), data = data_i,
                                                      family = quasibinomial(link = A$link),
                                                      weights = s.weights,
                                                      control = control), quote = TRUE)
          }, verbose = verbose)

          ps[[i]] <- fit.list[[i]]$fitted.values

        }
        if (isTRUE(A[["test2"]])) ps <- ps/rowSums(ps)
        fit.obj <- fit.list
      }
    }
    p.score <- NULL
  }
  else {
    n <- length(treat)
    p.score <- NULL
    treat <- factor(treat)
    treat_sub <- factor(treat[subset])
    bin.treat <- is_binary(treat_sub)

    bad.ps <- FALSE
    if (is_(ps, c("matrix", "data.frame"))) {
      if (all(dim(ps) == c(n, nunique(treat)))) {
        ps <- setNames(as.data.frame(ps), levels(treat))[subset, , drop = FALSE]
      }
      else if (all(dim(ps) == c(n, 1))) {
        ps <- setNames(list2DF(lapply(levels(treat), function(x) {
          p_ <- rep(1, length(treat))
          p_[treat == x] <- ps[treat == x, 1]
          p_
        })), levels(treat))[subset, , drop = FALSE]
      }
      else {
        bad.ps <- TRUE
      }
    }
    else if (is.numeric(ps)) {
      if (length(ps) == n) {
        ps <- setNames(list2DF(lapply(levels(treat), function(x) {
          p_ <- rep(1, length(treat))
          p_[treat == x] <- ps[treat == x]
          p_
        })), levels(treat))[subset, , drop = FALSE]
      }
      else {
        bad.ps <- TRUE
      }
    }
    else bad.ps <- TRUE

    if (bad.ps) .err("`ps` must be a numeric vector with a propensity score for each unit or a matrix \n\twith the probability of being in each treatment for each unit")

  }

  #ps should be matrix of probs for each treat
  #Computing weights
  w <- get_w_from_ps(ps = ps, treat = treat_sub, estimand, focal = focal,
                     stabilize = stabilize, subclass = subclass)

  list(w = w, ps = p.score, fit.obj = fit.obj)
}

weightit2glm.cont <- function(covs, treat, s.weights, subset, stabilize, missing, ps, verbose, ...) {
  A <- list(...)

  fit.obj <- NULL

  covs <- covs[subset, , drop = FALSE]
  treat <- treat[subset]
  s.weights <- s.weights[subset]

  if (missing == "ind") {
    covs <- add_missing_indicators(covs)
  }

  for (i in seq_col(covs)) covs[,i] <- make.closer.to.1(covs[,i])

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

  data <- data.frame(treat, covs)
  formula <- formula(data)

  #Process density params
  densfun <- get_dens_fun(use.kernel = isTRUE(A[["use.kernel"]]), bw = A[["bw"]],
                          adjust = A[["adjust"]], kernel = A[["kernel"]],
                          n = A[["n"]], treat = treat, density = A[["density"]])

  #Stabilization - get dens.num
  dens.num <- densfun(treat - mean(treat), s.weights)

  #Estimate GPS
  if (is_null(ps)) {

    if (missing == "saem") {
      rlang::check_installed("misaem")

      withCallingHandlers({verbosely({
        fit <- misaem::miss.lm(formula, data, control = as.list(A[["control"]]))
      }, verbose = verbose)},
      warning = function(w) {
        if (conditionMessage(w) != "one argument not used by format '%i '") {
          .wrn("(from `misaem::miss.lm()`) ", w, tidy = FALSE)
        }
        invokeRestart("muffleWarning")
      })

      if (is_null(A[["saem.method"]])) A[["saem.method"]] <- "map"

      gp.score <- drop(predict(fit, newdata = covs, method = A[["saem.method"]]))
    }
    else {
      if (is_null(A[["link"]])) A[["link"]] <- "identity"
      else {
        if (missing == "saem") acceptable.links <- "identity"
        else acceptable.links <- c("identity", "log", "inverse")

        which.link <- acceptable.links[pmatch(A[["link"]], acceptable.links, nomatch = 0)][1]
        if (is.na(which.link)) {
          .err(sprintf("Only %s allowed as the link for continuous treatments%s",
                       word_list(acceptable.links, quotes = TRUE, is.are = TRUE),
                       if (missing == "saem") ' with missing = "saem"' else ""))
        }

        A[["link"]] <- which.link
      }

      verbosely({
        fit <- do.call("glm", c(list(formula, data = data,
                                     weights = s.weights,
                                     family = gaussian(link = A[["link"]]),
                                     control = as.list(A$control))),
                       quote = TRUE)
      }, verbose = verbose)

      gp.score <- fit$fitted.values
    }

    fit.obj <- fit
  }
  else {
    gp.score <- ps
  }

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
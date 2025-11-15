#' Covariate Balancing Propensity Score Weighting
#' @name method_cbps
#' @usage NULL
#'
#' @description
#' This page explains the details of estimating weights from covariate balancing
#' propensity scores by setting `method = "cbps"` in the call to [weightit()] or
#' [weightitMSM()]. This method can be used with binary, multi-category, and
#' continuous treatments.
#'
#' In general, this method relies on estimating propensity scores using
#' generalized method of moments and then converting those propensity scores
#' into weights using a formula that depends on the desired estimand. This
#' method relies on code written for \pkg{WeightIt} using [optim()].
#'
#' ## Binary Treatments
#'
#' For binary treatments, this method estimates the propensity scores and
#' weights using `optim()` using formulas described by Imai and Ratkovic (2014).
#' The following estimands are allowed: ATE, ATT, ATC, and ATO.
#'
#' ## Multi-Category Treatments
#'
#' For multi-category treatments, this method estimates the generalized
#' propensity scores and weights using `optim()` using formulas described by
#' Imai and Ratkovic (2014). The following estimands are allowed: ATE and ATT.
#'
#' ## Continuous Treatments
#'
#' For continuous treatments, this method estimates the generalized propensity
#' scores and weights using `optim()` using a modification of the formulas
#' described by Fong, Hazlett, and Imai (2018). See Details.
#'
#' ## Longitudinal Treatments
#'
#' For longitudinal treatments, the weights are computed using methods similar
#' to those described by Huffman and van Gameren (2018). This involves
#' specifying moment conditions for the models at each time point as with
#' single-time point treatments but using the product of the time-specific
#' weights as the weights that are used in the balance moment conditions. This
#' yields weights that balance the covariate at each time point. This is not the
#' same implementation as is implemented in `CBPS::CBMSM()`, and results should
#' not be expected to align between the two methods. Any combination of
#' treatment types is supported.
#'
#' For the over-identified version (i.e., setting `over = TRUE`), the empirical
#' variance is used in the objective function, whereas the expected variance
#' averaging over the treatment is used with binary and multi-category point
#' treatments.
#'
#' ## Sampling Weights
#'
#' Sampling weights are supported through `s.weights` in all scenarios.
#'
#' ## Missing Data
#'
#' In the presence of missing data, the following value(s) for `missing` are
#' allowed:
#'     \describe{
#'       \item{`"ind"` (default)}{
#'         First, for each variable with missingness, a new missingness indicator variable is created which takes the value 1 if the original covariate is `NA` and 0 otherwise. The missingness indicators are added to the model formula as main effects. The missing values in the covariates are then replaced with the covariate medians (this value is arbitrary and does not affect estimation). The weight estimation then proceeds with this new formula and set of covariates. The covariates output in the resulting `weightit` object will be the original covariates with the `NA`s.
#'       }
#'     }
#'
#' ## M-estimation
#'
#' M-estimation is supported for the just-identified CBPS (the default, setting
#' `over = FALSE`) for binary and multi-category treatments. Otherwise (i.e.,
#' for continuous or longitudinal treatments or when `over = TRUE`),
#' M-estimation is not supported. See [glm_weightit()] and
#' `vignette("estimating-effects")` for details.
#'
#' @section Additional Arguments:
#'
#' The following additional arguments can be specified:
#' \describe{
#'     \item{`over`}{`logical`; whether to request the over-identified CBPS, which combines the generalized linear model regression score equations (for binary treatments), multinomial logistic regression score equations (for multi-category treatments), or linear regression score equations (for continuous treatments) to the balance moment conditions. Default is `FALSE` to use the just-identified CBPS.
#'     }
#'     \item{`twostep`}{`logical`; when `over = TRUE`, whether to use the two-step approximation to the generalized method of moments variance. Default is `TRUE`. Setting to `FALSE` increases computation time but may improve estimation. Ignored with a warning when `over = FALSE`.
#'     }
#'     \item{`link`}{the link used in the generalized linear model for the propensity scores when treatment is binary. Default is `"logit"` for logistic regression, which is used in the original description of the method by Imai and Ratkovic (2014), but others are allowed, including `"probit"`, `"cauchit"`, `"cloglog"`, `"loglog"`, `"log"`, `"clog"`, and `"identity"`. Note that negative weights are possible with these last three and they should be used with caution. An object of class `"link-glm"` can also be supplied. The argument is passed to [quasibinomial()]. Ignored for multi-category, continuous, and longitudinal treatments.
#'     }
#'     \item{`reltol`}{the relative tolerance for convergence of the optimization. Passed to the `control` argument of `optim()`. Default is `1e-10`.
#'     }
#'     \item{`maxit`}{the maximum number of iterations for convergence of the optimization. Passed to the `control` argument of `optim()`. Default is 1000 for binary and multi-category treatments and 10000 for continuous and longitudinal treatments.
#'     }
#'     \item{`solver`}{the solver to use to estimate the parameters of the just-identified CBPS. Allowable options include `"multiroot"` to use \pkgfun{rootSolve}{multiroot} and `"optim"` to use [stats::optim()]. `"multiroot"` is the default when \pkg{rootSolve} is installed, as it tends to be much faster and more accurate; otherwise, `"optim"` is the default and requires no dependencies. Regardless of `solver`, the output of `optim()` is returned when `include.obj = TRUE` (see below). When `over = TRUE`, the parameter estimates of the just-identified CBPS are used as starting values for the over-identified CBPS.
#'     }
#'     \item{`moments`}{`integer`; the highest power of each covariate to be balanced. For example, if `moments = 3`, each covariate, its square, and its cube will be balanced. Can also be a named vector with a value for each covariate (e.g., `moments = c(x1 = 2, x2 = 4)`). Values greater than 1 for categorical covariates are ignored. Default is 1 to balance covariate means.
#'     }
#'     \item{`int`}{`logical`; whether first-order interactions of the covariates are to be balanced. Default is `FALSE`.
#'     }
#'     \item{`quantile`}{a named list of quantiles (values between 0 and 1) for each continuous covariate, which are used to create additional variables that when balanced ensure balance on the corresponding quantile of the variable. For example, setting `quantile = list(x1 = c(.25, .5. , .75))` ensures the 25th, 50th, and 75th percentiles of `x1` in each treatment group will be balanced in the weighted sample. Can also be a single number (e.g., `.5`) or a vector (e.g., `c(.25, .5, .75)`) to request the same quantile(s) for all continuous covariates. Only allowed with binary and multi-category treatments.
#'     }
#'
#' }
#'
#' @section Additional Outputs:
#' \describe{
#'   \item{`obj`}{When `include.obj = TRUE`, the output of the final call to `optim()` used to produce the model parameters. Note that because of variable transformations, the resulting parameter estimates may not be interpretable.
#'   }
#' }
#'
#' @details
#' CBPS estimates the coefficients of a generalized linear model (for
#' binary treatments), multinomial logistic regression model (for multi-category
#' treatments), or linear regression model (for continuous treatments) that is
#' used to compute (generalized) propensity scores, from which the weights are
#' computed. It involves replacing (or augmenting, in the case of the
#' over-identified version) the standard maximum likelihood score equations with the
#' balance constraints in a generalized method of moments estimation. The idea
#' is to nudge the estimation of the coefficients toward those that produce
#' balance in the weighted sample. The just-identified version (with `over = FALSE`) does away with the maximum likelihood score equations for the coefficients so that only
#' the balance constraints are used, which will therefore
#' produce superior balance on the means (i.e., corresponding to the balance
#' constraints) for binary and multi-category treatments and linear terms for
#' continuous treatments than will the over-identified version.
#'
#' Just-identified CBPS is very similar to entropy balancing and inverse
#' probability tilting. For the ATT, all three methods will yield identical
#' estimates. For other estimands, the results will differ.
#'
#' Note that \pkg{WeightIt} provides different functionality from the \pkg{CBPS}
#' package in terms of the versions of CBPS available; for extensions to CBPS
#' (e.g., optimal CBPS and CBPS for instrumental variables), the \pkg{CBPS}
#' package may be preferred. Note that for longitudinal treatments,
#' `CBPS::CBMSM()` uses different methods and produces different results from
#' `weightitMSM()` called with `method = "cbps"`.
#'
#' @note
#' This method used to rely on functionality in the \pkg{CBPS} package,
#' but no longer does. Slight differences may be found between the two packages
#' in some cases due to numerical imprecision (or, for continuous and
#' longitudinal treatments, due to a difference in the estimator).
#' \pkg{WeightIt} supports arbitrary numbers of groups for the multi-category
#' CBPS and any estimand, whereas \pkg{CBPS} only supports up to four groups and
#' only the ATE. The implementation of the just-identified CBPS for continuous
#' treatments also differs from that of \pkg{CBPS}, and departs slightly from
#' that described by Fong et al. (2018). The treatment mean and treatment
#' variance are treated as random parameters to be estimated and are included in
#' the balance moment conditions. In Fong et al. (2018), the treatment mean and
#' variance are fixed to their empirical counterparts. For continuous treatments
#' with the over-identified CBPS, \pkg{WeightIt} and \pkg{CBPS} use different
#' methods of specifying the GMM variance matrix, which may lead to differing
#' results.
#'
#' Note that the default method differs between the two implementations; by
#' default \pkg{WeightIt} uses the just-identified CBPS, which is faster to fit,
#' yields better balance, and is compatible with M-estimation for estimating the
#' standard error of the treatment effect, whereas \pkg{CBPS} uses the
#' over-identified CBPS by default. However, both the just-identified and
#' over-identified versions are available in both packages.
#'
#' When the \pkg{rootSolve} package is installed, the optimization process will
#' be slightly faster and more accurate because starting values are provided by
#' an initial call to \pkgfun{rootSolve}{multiroot}. However, the package is not
#' required.
#'
#' @seealso [weightit()], [weightitMSM()]
#'
#' [method_ebal] and [method_ipt] for entropy balancing and inverse probability
#' tilting, which work similarly.
#'
#' @references
#' ## Binary treatments
#'
#' Imai, K., & Ratkovic, M. (2014). Covariate balancing propensity score.
#' *Journal of the Royal Statistical Society: Series B (Statistical Methodology)*, 76(1), 243–263.
#'
#' ## Multi-Category treatments
#'
#' Imai, K., & Ratkovic, M. (2014). Covariate balancing propensity score.
#' *Journal of the Royal Statistical Society: Series B (Statistical Methodology)*, 76(1), 243–263.
#'
#' ## Continuous treatments
#'
#' Fong, C., Hazlett, C., & Imai, K. (2018). Covariate balancing propensity
#' score for a continuous treatment: Application to the efficacy of political
#' advertisements. *The Annals of Applied Statistics*, 12(1), 156–177.
#' \doi{10.1214/17-AOAS1101}
#'
#' ## Longitudinal treatments
#'
#' Huffman, C., & van Gameren, E. (2018). Covariate Balancing Inverse
#' Probability Weights for Time-Varying Continuous Interventions. *Journal of Causal Inference*, 6(2). \doi{10.1515/jci-2017-0002}
#'
#' Note: one should not cite Imai & Ratkovic (2015) when using CBPS for
#' longitudinal treatments.
#'
#' Some of the code was inspired by the source code of the \pkg{CBPS} package.
#'
#' @examples
#' data("lalonde", package = "cobalt")
#'
#' #Balancing covariates between treatment groups (binary)
#' (W1a <- weightit(treat ~ age + educ + married +
#'                   nodegree + re74, data = lalonde,
#'                 method = "cbps", estimand = "ATT"))
#'
#' summary(W1a)
#'
#' cobalt::bal.tab(W1a)
#'
#' #Balancing covariates between treatment groups (binary)
#' #using over-identified CBPS with probit link
#' (W1b <- weightit(treat ~ age + educ + married +
#'                   nodegree + re74, data = lalonde,
#'                 method = "cbps", estimand = "ATT",
#'                 over = TRUE, link = "probit"))
#'
#' summary(W1b)
#'
#' cobalt::bal.tab(W1b)
#'
#' #Balancing covariates with respect to race (multi-category)
#' (W2 <- weightit(race ~ age + educ + married +
#'                   nodegree + re74, data = lalonde,
#'                 method = "cbps", estimand = "ATE"))
#'
#' summary(W2)
#'
#' cobalt::bal.tab(W2)
#'
#' #Balancing covariates with respect to re75 (continuous)
#' (W3 <- weightit(re75 ~ age + educ + married +
#'                   nodegree + re74, data = lalonde,
#'                 method = "cbps"))
#'
#' summary(W3)
#'
#' cobalt::bal.tab(W3)
#' \donttest{
#' #Longitudinal treatments
#' data("msmdata")
#' (W4 <- weightitMSM(list(A_1 ~ X1_0 + X2_0,
#'                         A_2 ~ X1_1 + X2_1 +
#'                           A_1 + X1_0 + X2_0),
#'                    data = msmdata,
#'                    method = "cbps"))
#'
#' summary(W4)
#'
#' cobalt::bal.tab(W4)
#' }
NULL

weightit2cbps <- function(covs, treat, s.weights, estimand, focal, subset,
                          stabilize, missing, verbose, ...) {

  covs <- covs[subset, , drop = FALSE]
  treat <- treat[subset]
  s.weights <- s.weights[subset]

  missing <- .process_missing2(missing, covs)

  if (missing == "ind") {
    covs <- add_missing_indicators(covs)
  }

  covs <- covs |>
    .apply_moments_int_quantile(moments = ...get("moments"),
                                int = ...get("int"),
                                quantile = ...get("quantile"),
                                s.weights = s.weights, focal = focal,
                                treat = treat) |>
    .make_covs_full_rank()

  t.lev <- get_treated_level(treat, estimand, focal)
  treat <- binarize(treat, one = t.lev)

  mod_covs <- cbind(`(Intercept)` = 1, scale(svd(covs)$u))
  bal_covs <- mod_covs

  solver <- ...get("solver", NULL)
  if (is_null(solver)) {
    if (rlang::is_installed("rootSolve")) {
      solver <- "multiroot"
    }
    else {
      solver <- "optim"
    }
  }
  else {
    chk::chk_string(solver)
    solver <- match_arg(solver, c("optim", "multiroot"))
  }

  if (solver == "multiroot") {
    rlang::check_installed("rootSolve")
  }

  over <- ...get("over", FALSE)
  chk::chk_flag(over)

  twostep <- ...get("twostep", TRUE)

  if (over) {
    chk::chk_flag(twostep)
  }
  else if (!isTRUE(twostep)) {
    .wrn("`twostep` is ignored when `over = FALSE`")
  }

  reltol <- ...get("reltol", 1e-10)
  chk::chk_number(reltol)

  maxit <- ...get("maxit", 5e3L)
  chk::chk_count(maxit)

  N <- sum(s.weights)

  link <- ...get("link", "logit")

  if (chk::vld_string(link)) {
    chk::chk_subset(link, c("logit", "probit", "cloglog", "loglog", "cauchit",
                            "log", "clog", "identity"))

    link <- .make_link(link)

    if (identical(link, "logit") && estimand == "ATO") {
      over <- FALSE
    }
  }
  else if (inherits(link, "family") && is_not_null(link$linkfun) &&
           is_not_null(link$linkinv) && is_not_null(link$mu.eta) &&
           is_not_null(link$valideta)) {
    link <- list(linkfun = link$linkfun,
                 linkinv = link$linkinv,
                 mu.eta = link$mu.eta,
                 valideta = link$valideta,
                 name = link$link)
    class(link) <- "link-glm"
  }
  else if (!inherits(link, "link-glm")) {
    .err('`link` must be a string or an object of class "link-glm"')
  }

  .fam <- quasibinomial(link)

  # Balance condition for ATT
  psi_bal <- switch(estimand,
                    ATE = function(B, Xm, Xb = Xm, A, SW) {
                      p <- .fam$linkinv(drop(Xm %*% B))
                      SW * (A / p - (1 - A) / (1 - p)) * Xb
                    },
                    ATT = function(B, Xm, Xb = Xm, A, SW) {
                      p <- .fam$linkinv(drop(Xm %*% B))
                      SW * ((A - p) / (1 - p)) * Xb
                    },
                    ATC = function(B, Xm, Xb = Xm, A, SW) {
                      p <- .fam$linkinv(drop(Xm %*% B))
                      SW * ((A - p) / p) * Xb
                    },
                    ATO = function(B, Xm, Xb = Xm, A, SW) {
                      p <- .fam$linkinv(drop(Xm %*% B))
                      SW * (A * (1 - p) - (1 - A) * p) * Xb
                    })

  obj_bal <- function(B, Xm, Xb = Xm, A, SW) {
    gbar <- colMeans(psi_bal(B, Xm, Xb, A, SW))
    sqrt(sum(gbar^2))
  }

  # Initialize coefs using glm
  par_glm <- .get_glm_starting_values(X = mod_covs, Y = treat, w = s.weights,
                                      family = .fam)

  # Slightly improve glm coefs to move closer to optimal
  alpha.func <- function(alpha) obj_bal(par_glm * alpha, mod_covs, bal_covs, treat, s.weights)
  par_alpha <- par_glm * optimize(alpha.func, interval = c(.8, 1.1))$min

  if (solver == "multiroot") {
    out <- suppressWarnings({
      try(verbosely({
        rootSolve::multiroot(f = function(...) colMeans(psi_bal(...)),
                             start = par_alpha,
                             Xm = mod_covs,
                             Xb = bal_covs,
                             A = treat,
                             SW = s.weights,
                             rtol = reltol,
                             atol = reltol,
                             ctol = reltol)
      }, verbose = FALSE), silent = TRUE)
    })

    if (!null_or_error(out) && utils::hasName(out, "root") &&
        utils::hasName(out, "estim.precis") &&
        chk::vld_number(out[["estim.precis"]]) &&
        out[["estim.precis"]] < 1e-5) {
      par_alpha <- out[["root"]]
    }
  }

  # Optimize balance objective
  out <- optim(par = par_alpha,
               fn = obj_bal,
               method = "BFGS",
               control = list(maxit = maxit,
                              reltol = reltol,
                              trace = as.integer(!over && verbose)),
               Xm = mod_covs,
               Xb = bal_covs,
               A = treat,
               SW = s.weights)

  par_out <- out$par

  if (over) {
    #Generalized linear model score
    psi_glm <- function(B, Xm, A, SW) {
      lin_pred <- drop(Xm %*% B)
      p <- .fam$linkinv(lin_pred)
      (SW * (A - p) * .fam$mu.eta(lin_pred) / .fam$variance(p)) * Xm
    }

    # Combine LR and balance
    psi <- function(B, Xm, Xb = Xm, A, SW) {
      cbind(psi_glm(B, Xm, A, SW), psi_bal(B, Xm, Xb, A, SW))
    }

    Sigma <- function(B, Xm, Xb = Xm, A, SW) {
      lp <- drop(Xm %*% B)
      p <- .fam$linkinv(lp)
      g <- .fam$mu.eta(lp) / .fam$variance(p)

      S11 <- crossprod(SW * g * (p * (1 - p)) * Xm, SW * g * Xm)
      S12 <- switch(estimand,
                    ATE = crossprod(SW * g * Xm, SW * Xb),
                    ATT = crossprod(SW * g * Xm, SW * p * Xb),
                    ATC = crossprod(SW * g * Xm, SW * (1 - p) * Xb),
                    ATO = crossprod(SW * g * Xm, SW * p * (1 - p) * Xb))
      S21 <- t(S12)
      S22 <- switch(estimand,
                    ATE = crossprod(SW / (p * (1 - p)) * Xb, SW * Xb),
                    ATT = crossprod(SW * (p / (1 - p)) * Xb, SW * Xb),
                    ATC = crossprod(SW * ((1 - p) / p) * Xb, SW * Xb),
                    ATO = crossprod(SW * (p * (1 - p)) * Xb, SW * Xb))

      rbind(cbind(S11, S12),
            cbind(S21, S22)) / N
    }

    obj <- function(B, Xm, Xb = Xm, A, SW, invS = NULL) {
      if (is_null(invS)) {
        invS <- generalized_inverse(Sigma(B, Xm, Xb, A, SW))
      }

      gbar <- colMeans(psi(B, Xm, Xb, A, SW))

      sqrt(drop(t(gbar) %*% invS %*% gbar))
    }

    invS <- {
      if (twostep) generalized_inverse(Sigma(par_alpha, mod_covs, bal_covs, treat, s.weights))
      else NULL
    }

    start.list <- {
      if (max(abs(par_alpha - par_out)) < 1e-6) list(par_out)
      else list(par_alpha, par_out)
    }

    out <- lapply(start.list, function(par_) {
      optim(par = par_,
            fn = obj,
            method = "BFGS",
            control = list(maxit = maxit,
                           reltol = reltol,
                           trace = as.integer(verbose)),
            Xm = mod_covs,
            Xb = bal_covs,
            A = treat,
            SW = s.weights,
            invS = invS)
    })

    out <- out[[which.min(unlist(grab(out, "value")))]]

    par_out <- out$par
  }

  p.score <- .fam$linkinv(drop(mod_covs %*% par_out))

  if (out$converge != 0) {
    .wrn("the optimization failed to converge; try again with a higher value of `maxit`")
  }

  w <- .get_w_from_ps_internal_bin(p.score, treat, estimand = estimand,
                                   stabilize = stabilize)

  Mparts <- NULL
  if (!over) {
    Mparts <- list(
      psi_treat = function(Btreat, Xtreat, A, SW) {
        psi_bal(Btreat, Xtreat, Xtreat, A, SW)
      },
      wfun = function(Btreat, Xtreat, A) {
        ps <- .fam$linkinv(drop(Xtreat %*% Btreat))
        .get_w_from_ps_internal_bin(ps, A, estimand = estimand,
                                    stabilize = stabilize)
      },
      dw_dBtreat = function(Btreat, Xtreat, A, SW) {
        XB <- drop(Xtreat %*% Btreat)
        ps <- .fam$linkinv(XB)
        .dw_dp_bin(ps, A, estimand = estimand) * .fam$mu.eta(XB) * Xtreat
      },
      hess_treat = function(Btreat, Xtreat, A, SW) {
        XB <- drop(Xtreat %*% Btreat)
        ps <- .fam$linkinv(XB)

        dw <- .dw_dp_bin(ps, A, estimand = estimand) * .fam$mu.eta(XB) * SW
        crossprod(Xtreat, dw * (2 * A - 1) * Xtreat)
      },
      Xtreat = mod_covs,
      A = treat,
      btreat = par_out
    )
  }

  list(w = w, ps = p.score, fit.obj = out,
       Mparts = Mparts)
}

weightit2cbps.multi <- function(covs, treat, s.weights, estimand, focal, subset,
                                stabilize, missing, verbose, ...) {

  covs <- covs[subset, , drop = FALSE]
  treat <- droplevels(treat[subset])
  s.weights <- s.weights[subset]

  missing <- .process_missing2(missing, covs)

  if (missing == "ind") {
    covs <- add_missing_indicators(covs)
  }

  covs <- covs |>
    .apply_moments_int_quantile(moments = ...get("moments"),
                                int = ...get("int"),
                                quantile = ...get("quantile"),
                                s.weights = s.weights, focal = focal,
                                treat = treat) |>
    .make_covs_full_rank()

  mod_covs <- cbind(`(Intercept)` = 1, scale(svd(covs)$u))
  bal_covs <- mod_covs

  solver <- ...get("solver", NULL)
  if (is_null(solver)) {
    if (rlang::is_installed("rootSolve")) {
      solver <- "multiroot"
    }
    else {
      solver <- "optim"
    }
  }
  else {
    chk::chk_string(solver)
    solver <- match_arg(solver, c("optim", "multiroot"))
  }

  if (solver == "multiroot") {
    rlang::check_installed("rootSolve")
  }

  over <- ...get("over", FALSE)
  chk::chk_flag(over)

  twostep <- ...get("twostep", TRUE)

  if (over) {
    chk::chk_flag(twostep)
  }
  else if (!isTRUE(twostep)) {
    .wrn("`twostep` is ignored when `over = FALSE`")
  }

  reltol <- ...get("reltol", 1e-10)
  chk::chk_number(reltol)

  maxit <- ...get("maxit", 1e3L)
  chk::chk_count(maxit)

  N <- sum(s.weights)

  #Function to compute predicted probabilities
  get_pp <- function(B, Xm) {
    qq <- exp(Xm %*% matrix(B, nrow = ncol(Xm)))

    pp <- cbind(1, qq) / (1 + rowSums(qq))

    colnames(pp) <- levels(treat)
    pp
  }

  # Balance measured between all combinations of treatment groups
  combs <- {
    if (over) utils::combn(levels(treat), 2L, simplify = FALSE)
    else lapply(levels(treat)[-1L], function(i) c(levels(treat)[1L], i))
  }

  psi_bal <- switch(estimand,
                    ATE = function(B, Xm, Xb = Xm, A, SW) {
                      pp <- get_pp(B, Xm)

                      do.call("cbind", lapply(combs, function(co) {
                        SW * ((A == co[1L]) / pp[, co[1L]] - (A == co[2L]) / pp[, co[2L]]) * Xb
                      }))
                    },
                    ATO = function(B, Xm, Xb = Xm, A, SW) {
                      pp <- get_pp(B, Xm)

                      do.call("cbind", lapply(combs, function(co) {
                        SW * ((A == co[1L]) / pp[, co[1L]] - (A == co[2L]) / pp[, co[2L]]) * Xb / rowSums(1 / pp)
                      }))
                    },
                    function(B, Xm, Xb = Xm, A, SW) {
                      pp <- get_pp(B, Xm)

                      do.call("cbind", lapply(combs, function(co) {
                        SW * pp[, focal] * ((A == co[1L]) / pp[, co[1L]] - (A == co[2L]) / pp[, co[2L]]) * Xb
                      }))
                    })

  obj_bal <- function(B, Xm, Xb = Xm, A, SW) {
    gbar <- colMeans(psi_bal(B, Xm, Xb, A, SW))
    sqrt(sum(gbar^2))
  }

  # Initialize coefs using multinomial logistic regression
  par_glm <- .multinom_weightit.fit(mod_covs, treat, weights = s.weights,
                                    hess = FALSE)$coefficients

  # Slightly improve glm coefs to move closer to optimal
  alpha.func <- function(alpha) obj_bal(par_glm * alpha, Xm = mod_covs, Xb = bal_covs,
                                        A = treat, SW = s.weights)
  par_alpha <- par_glm * optimize(alpha.func, interval = c(.8, 1.1))$min

  if (solver == "multiroot") {
    out <- suppressWarnings({
      try(verbosely({
        rootSolve::multiroot(f = function(...) colMeans(psi_bal(...)),
                             start = par_alpha,
                             Xm = mod_covs,
                             Xb = bal_covs,
                             A = treat,
                             SW = s.weights,
                             rtol = reltol,
                             atol = reltol,
                             ctol = reltol)
      }, verbose = FALSE), silent = TRUE)
    })

    if (!null_or_error(out) && utils::hasName(out, "root") &&
        utils::hasName(out, "estim.precis") &&
        chk::vld_number(out[["estim.precis"]]) &&
        out[["estim.precis"]] < 1e-5) {
      par_alpha <- out[["root"]]
    }
  }

  # Optimize balance objective
  out <- optim(par = par_alpha,
               fn = obj_bal,
               method = "BFGS",
               control = list(maxit = maxit,
                              reltol = reltol,
                              trace = as.integer(!over && verbose)),
               Xm = mod_covs,
               Xb = bal_covs,
               A = treat,
               SW = s.weights)

  par_out <- out$par

  if (over) {
    #Multinomial logistic regression score
    psi_mlr <- function(B, Xm, A, SW) {
      pp <- get_pp(B, Xm)

      do.call("cbind", lapply(levels(treat), function(i) {
        SW * ((A == i) - pp[, i]) * Xm
      }))
    }

    # Combine LR and balance
    psi <- function(B, Xm, Xb = Xm, A, SW) {
      cbind(psi_mlr(B, Xm, A, SW),
            psi_bal(B, Xm, Xb, A, SW))
    }

    Sigma <- function(B, Xm, Xb = Xm, A, SW) {
      pp <- get_pp(B, Xm)

      if (estimand == "ATO") {
        wden <- rowSums(1 / pp)
      }

      swXmXb <- switch(estimand,
                       ATE = crossprod(SW * Xm, SW * Xb),
                       ATO = crossprod(SW * Xm, SW * Xb / wden),
                       crossprod(SW * Xm, SW * pp[, focal] * Xb))

      S <- list()

      for (i in levels(treat)) {
        for (j in levels(treat)) {
          S[[sprintf("m%s_m%s", i, j)]] <- {
            if      (i == j) crossprod(SW * (pp[, i] * (1 - pp[, i])) *  Xm, SW * Xm)
            else if (i < j) -crossprod(SW * (pp[, i] * pp[, j]) * Xm, SW * Xm)
            else t(S[[sprintf("m%s_m%s", j, i)]])
          }
        }
      }

      for (i in levels(treat)) {
        for (jj in combs) {
          m <- {
            if      (i == jj[1L])  swXmXb
            else if (i == jj[2L]) -swXmXb
            else matrix(0, ncol(Xm), ncol(Xb))
          }
          S[[sprintf("m%s_b%s%s", i, jj[1L], jj[2L])]] <- m
          S[[sprintf("b%s%s_m%s", jj[1L], jj[2L], i)]] <- t(m)
        }
      }

      for (ii in combs) {
        for (jj in combs) {
          m <- switch(estimand,
                      ATE = {
                        if (identical(ii, jj))    crossprod(SW * (1 / pp[, ii[1L]] + 1 / pp[, ii[2L]]) * Xb,
                                                            SW * Xb)
                        else if (ii[1L] == jj[1L])  crossprod(SW * (1 / pp[, ii[1L]]) * Xb, SW * Xb)
                        else if (ii[1L] == jj[2L]) -crossprod(SW * (1 / pp[, ii[1L]]) * Xb, SW * Xb)
                        else if (ii[2L] == jj[1L]) -crossprod(SW * (1 / pp[, ii[2L]]) * Xb, SW * Xb)
                        else if (ii[2L] == jj[2L])  crossprod(SW * (1 / pp[, ii[2L]]) * Xb, SW * Xb)
                        else sq_matrix(0, n = ncol(Xb))
                      },
                      ATO = {
                        if (identical(ii, jj))    crossprod(SW * (1 / pp[, ii[1L]] + 1 / pp[, ii[2L]]) * Xb,
                                                            SW * Xb / wden)
                        else if (ii[1L] == jj[1L])  crossprod(SW * (1 / pp[, ii[1L]]) * Xb, SW * Xb / wden)
                        else if (ii[1L] == jj[2L]) -crossprod(SW * (1 / pp[, ii[1L]]) * Xb, SW * Xb / wden)
                        else if (ii[2L] == jj[1L]) -crossprod(SW * (1 / pp[, ii[2L]]) * Xb, SW * Xb / wden)
                        else if (ii[2L] == jj[2L])  crossprod(SW * (1 / pp[, ii[2L]]) * Xb, SW * Xb / wden)
                        else sq_matrix(0, n = ncol(Xb))
                      },
                      {
                        if (identical(ii, jj)) crossprod(SW * (1 / pp[, ii[1L]] + 1 / pp[, ii[2L]]) * pp[, focal] * Xb,
                                                         SW * pp[, focal] * Xb)
                        else if (ii[1L] == jj[1L])  crossprod(SW * (1 / pp[, ii[1L]]) * pp[, focal] * Xb, SW * pp[, focal] * Xb)
                        else if (ii[1L] == jj[2L]) -crossprod(SW * (1 / pp[, ii[1L]]) * pp[, focal] * Xb, SW * pp[, focal] * Xb)
                        else if (ii[2L] == jj[1L]) -crossprod(SW * (1 / pp[, ii[2L]]) * pp[, focal] * Xb, SW * pp[, focal] * Xb)
                        else if (ii[2L] == jj[2L])  crossprod(SW * (1 / pp[, ii[2L]]) * pp[, focal] * Xb, SW * pp[, focal] * Xb)
                        else sq_matrix(0, n = ncol(Xb))
                      })
          S[[sprintf("b%s%s_b%s%s", ii[1L], ii[2L], jj[1L], jj[2L])]] <- m
          S[[sprintf("b%s%s_b%s%s", jj[1L], jj[2L], ii[1L], ii[2L])]] <- m
        }
      }

      nms <- c(sprintf("m%s", levels(treat)),
               sprintf("b%s", vapply(combs, paste0, character(1L), collapse = "")))

      do.call("rbind", lapply(nms, function(i) {
        do.call("cbind", lapply(nms, function(j) {
          S[[sprintf("%s_%s", i, j)]]
        }))
      })) / N
    }

    obj <- function(B, Xm, Xb = Xm, A, SW, invS = NULL) {
      if (is_null(invS)) {
        invS <- generalized_inverse(Sigma(B, Xm, Xb, A, SW))
      }

      gbar <- colMeans(psi(B, Xm, Xb, A, SW))

      sqrt(drop(t(gbar) %*% invS %*% gbar)) #use sqrt to improve convergence
    }

    invS <- {
      if (twostep) generalized_inverse(Sigma(par_alpha, mod_covs, bal_covs,
                                             treat, s.weights))
      else NULL
    }

    start.list <- {
      if (max(abs(par_alpha - par_out)) < 1e-6) list(par_out)
      else list(par_alpha, par_out)
    }

    out <- lapply(start.list, function(par_) {
      optim(par = par_,
            fn = obj,
            method = "BFGS",
            control = list(maxit = maxit,
                           reltol = reltol,
                           trace = as.integer(verbose)),
            Xm = mod_covs,
            Xb = bal_covs,
            A = treat,
            SW = s.weights,
            invS = invS)
    })

    out <- out[[which.min(unlist(grab(out, "value")))]]

    par_out <- out$par
  }

  if (out$converge != 0) {
    .wrn("the optimization failed to converge; try again with a higher value of `maxit`")
  }

  pp <- get_pp(par_out, mod_covs)

  w <- .get_w_from_ps_internal_multi(pp, treat,
                                     estimand = estimand,
                                     focal = focal,
                                     stabilize = stabilize)

  out$pp <- pp

  Mparts <- NULL
  if (!over) {
    Mparts <- list(
      psi_treat = function(Btreat, Xtreat, A, SW) {
        psi_bal(Btreat, Xtreat, Xtreat, A, SW)
      },
      wfun = function(Btreat, Xtreat, A) {
        ps <- get_pp(Btreat, Xtreat)
        .get_w_from_ps_internal_multi(ps, A, estimand = estimand, focal = focal,
                                      stabilize = stabilize)
      },
      dw_dBtreat = function(Btreat, Xtreat, A, SW) {
        ps <- get_pp(Btreat, Xtreat)
        dw <- .dw_dp_multi(ps, A, estimand = estimand, focal = focal)

        do.call("cbind", lapply(levels(A)[-1L], function(i) {
          Xtreat * Reduce("+", lapply(levels(A), function(j) {
            dw[, j] * ps[, j] * ((i == j) - ps[, i])
          }))
        }))
      },
      hess_treat = function(Btreat, Xtreat, A, SW) {
        ps <- get_pp(Btreat, Xtreat)
        dw <- .dw_dp_multi(ps, A, estimand = estimand, focal = focal)

        dpsi_dw <- do.call("cbind", lapply(levels(A)[-1L], function(i) {
          ((A == levels(A)[1L]) - (A == i)) * SW * Xtreat
        }))

        dw_dB <- do.call("cbind", lapply(levels(A)[-1L], function(i) {
          Xtreat * Reduce("+", lapply(levels(A), function(j) {
            dw[, j] * ps[, j] * ((i == j) - ps[, i])
          }))
        }))

        crossprod(dpsi_dw, dw_dB)
      },
      Xtreat = mod_covs,
      A = treat,
      btreat = par_out
    )
  }

  list(w = w, fit.obj = out,
       Mparts = Mparts)
}

weightit2cbps.cont <- function(covs, treat, s.weights, subset, missing, verbose, ...) {

  covs <- covs[subset, , drop = FALSE]
  treat <- as.numeric(treat[subset])
  s.weights <- s.weights[subset]

  missing <- .process_missing2(missing, covs)

  if (missing == "ind") {
    covs <- add_missing_indicators(covs)
  }

  covs <- covs |>
    .apply_moments_int_quantile(moments = ...get("moments"),
                                int = ...get("int")) |>
    .make_covs_full_rank()

  treat <- scale_w(treat, s.weights)

  for (i in seq_col(covs)) {
    covs[, i] <- scale_w(covs[, i], s.weights)
  }

  mod_covs <- cbind(`(Intercept)` = 1, scale(svd(covs)$u))

  bal_covs <- mod_covs

  solver <- ...get("solver", NULL)
  if (is_null(solver)) {
    if (rlang::is_installed("rootSolve")) {
      solver <- "multiroot"
    }
    else {
      solver <- "optim"
    }
  }
  else {
    chk::chk_string(solver)
    solver <- match_arg(solver, c("optim", "multiroot"))
  }

  if (solver == "multiroot") {
    rlang::check_installed("rootSolve")
  }

  over <- ...get("over", FALSE)
  chk::chk_flag(over)

  twostep <- ...get("twostep", TRUE)

  if (over) {
    chk::chk_flag(twostep)
  }
  else if (!isTRUE(twostep)) {
    .wrn("`twostep` is ignored when `over = FALSE`")
  }

  reltol <- ...get("reltol", 1e-10)
  chk::chk_number(reltol)

  maxit <- ...get("maxit", 1e4L)
  chk::chk_count(maxit)

  s.weights <- s.weights / mean_fast(s.weights)

  squish_tol <- 50

  # Balance condition
  psi_bal <- function(B, Xm, Xb = Xm, A, SW) {
    un_s2 <- exp(B[1L])
    un_p <- B[2L]
    log.dens.num <- squish(dnorm(A, un_p, sqrt(un_s2), log = TRUE),
                           lo = -Inf, hi = squish_tol)

    s2 <- exp(B[3L])
    p <- drop(Xm %*% B[-(1:3)])
    log.dens.denom <- squish(dnorm(A, p, sqrt(s2), log = TRUE),
                             lo = -squish_tol, hi = Inf)

    w <- exp(log.dens.num - log.dens.denom)

    cbind(SW * (A - un_p)^2 - un_s2,
          SW * (A - un_p),
          SW * (A - p)^2 - s2,
          SW * w * A * Xb)
  }

  obj_bal <- function(B, Xm, Xb = Xm, A, SW) {
    gbar <- colMeans(psi_bal(B, Xm, Xb, A, SW))
    sqrt(sum(gbar^2))
  }

  # Initialize coefs using linear regression
  init.fit <- lm.wfit(mod_covs, treat, w = s.weights)
  par_glm <- c(0, 0, log(var(init.fit$residuals)), init.fit$coefficients)
  names(par_glm)[1:3] <- c("log(s^2)", "E[A]", "log(s_r^2)")

  # Slightly improve glm coefs to move closer to optimal
  alpha.func <- function(alpha) obj_bal(par_glm * alpha, mod_covs, bal_covs, treat, s.weights)
  par_alpha <- par_glm * optimize(alpha.func, interval = c(.8, 1.1))$min

  if (solver == "multiroot") {
    out <- suppressWarnings({
      try(verbosely({
        rootSolve::multiroot(f = function(...) colMeans(psi_bal(...)),
                             start = par_alpha,
                             Xm = mod_covs,
                             Xb = bal_covs,
                             A = treat,
                             SW = s.weights,
                             rtol = reltol,
                             atol = reltol,
                             ctol = reltol)
      }, verbose = FALSE), silent = TRUE)
    })

    if (!null_or_error(out) && utils::hasName(out, "root") &&
        utils::hasName(out, "estim.precis") &&
        chk::vld_number(out[["estim.precis"]]) &&
        out[["estim.precis"]] < 1e-5) {
      par_alpha <- out[["root"]]
    }
  }

  # Optimize balance objective
  out <- optim(par = par_alpha,
               fn = obj_bal,
               method = "BFGS",
               control = list(maxit = maxit,
                              reltol = reltol,
                              trace = as.integer(!over && verbose)),
               Xm = mod_covs,
               Xb = bal_covs,
               A = treat,
               SW = s.weights)

  par_out <- out$par

  if (over) {

    #Linear regression score + score for marginal mean + var and conditional var
    psi_lm <- function(B, Xm, A, SW) {
      un_s2 <- exp(B[1L])
      un_p <- B[2L]
      s2 <- exp(B[3L])

      p <- drop(Xm %*% B[-(1:3)])

      cbind(SW * (A - p) * Xm)
    }

    # Combine LR and balance
    psi <- function(B, Xm, Xb = Xm, A, SW) {
      cbind(psi_lm(B, Xm, A, SW), psi_bal(B, Xm, Xb, A, SW))
    }

    Sigma <- function(B, Xm, Xb = Xm, A, SW) {
      crossprod(psi(B, Xm, Xb, A, SW))
    }

    obj <- function(B, Xm, Xb = Xm, A, SW, invS = NULL) {
      psi0 <- psi(B, Xm, Xb, A, SW)

      if (is_null(invS)) {
        invS <- generalized_inverse(crossprod(psi0))
      }

      gbar <- colMeans(psi0)
      sqrt(drop(t(gbar) %*% invS %*% gbar))
    }

    invS <- {
      if (twostep) generalized_inverse(Sigma(par_alpha, mod_covs, bal_covs, treat, s.weights))
      else NULL
    }

    start.list <- {
      if (max(abs(par_alpha - par_out)) < 1e-6) list(par_out)
      else list(par_alpha, par_out)
    }

    out <- lapply(start.list, function(par_) {
      optim(par = par_,
            fn = obj,
            method = "BFGS",
            control = list(maxit = maxit,
                           reltol = reltol,
                           trace = as.integer(verbose)),
            Xm = mod_covs,
            Xb = bal_covs,
            A = treat,
            SW = s.weights,
            invS = invS)
    })

    out <- out[[which.min(unlist(grab(out, "value")))]]

    par_out <- out$par
  }

  if (out$converge != 0) {
    .wrn("the optimization failed to converge; try again with a higher value of `maxit`")
  }

  un_s2 <- exp(par_out[1L])
  un_p <- par_out[2L]
  log.dens.num <- squish(dnorm(treat, un_p, sqrt(un_s2), log = TRUE),
                         lo = -Inf, hi = squish_tol)

  s2 <- exp(par_out[3L])
  p <- drop(mod_covs %*% par_out[-(1:3)])
  log.dens.denom <- squish(dnorm(treat, p, sqrt(s2), log = TRUE),
                           lo = -squish_tol, hi = Inf)

  w <- exp(log.dens.num - log.dens.denom)

  Mparts <- NULL

  list(w = w, fit.obj = out,
       Mparts = Mparts)
}

weightitMSM2cbps <- function(covs.list, treat.list, s.weights, subset, missing, verbose, ...) {

  s.weights <- s.weights[subset]
  treat.types <- character(length(treat.list))

  for (i in seq_along(covs.list)) {
    treat.list[[i]] <- treat.list[[i]][subset]

    if (!has_treat_type(treat.list[[i]])) {
      treat.list[[i]] <- assign_treat_type(treat.list[[i]])
    }
    treat.types[i] <- get_treat_type(treat.list[[i]])

    covs.list[[i]] <- covs.list[[i]][subset, , drop = FALSE]

    if (.process_missing2(missing, covs.list[[i]]) == "ind") {
      covs.list[[i]] <- add_missing_indicators(covs.list[[i]])
    }

    if (treat.types[i] %in% c("binary", "multinomial", "multi-category")) {
      covs.list[[i]] <- covs.list[[i]] |>
        .apply_moments_int_quantile(moments = ...get("moments"),
                                    int = ...get("int"),
                                    quantile = ...get("quantile"),
                                    s.weights = s.weights) |>
        .make_covs_full_rank()
    }
    else {
      covs.list[[i]] <- covs.list[[i]] |>
        .apply_moments_int_quantile(moments = ...get("moments"),
                                    int = ...get("int")) |>
        .make_covs_full_rank()
    }

    treat.list[[i]] <- switch(
      treat.types[i],
      binary = binarize(treat.list[[i]], one = get_treated_level(treat.list[[i]], "ATE")),
      `multi-category` =,
      multinomial = factor(treat.list[[i]]),
      continuous = scale_w(treat.list[[i]], s.weights)
    )

    for (j in seq_col(covs.list[[i]])) {
      covs.list[[i]][, j] <- scale_w(covs.list[[i]][, j], s.weights)
    }

    covs.list[[i]] <- cbind(`(Intercept)` = 1, scale(svd(covs.list[[i]])$u))
  }

  solver <- ...get("solver", NULL)
  if (is_null(solver)) {
    if (rlang::is_installed("rootSolve")) {
      solver <- "multiroot"
    }
    else {
      solver <- "optim"
    }
  }
  else {
    chk::chk_string(solver)
    solver <- match_arg(solver, c("optim", "multiroot"))
  }

  if (solver == "multiroot") {
    rlang::check_installed("rootSolve")
  }

  over <- ...get("over", FALSE)
  chk::chk_flag(over)

  twostep <- ...get("twostep", TRUE)

  if (over) {
    chk::chk_flag(twostep)
  }
  else if (!isTRUE(twostep)) {
    .wrn("`twostep` is ignored when `over = FALSE`")
  }

  reltol <- ...get("reltol", 1e-10)
  chk::chk_number(reltol)

  maxit <- ...get("maxit", 1e4L)
  chk::chk_count(maxit)

  coef_ind <- make_list(length(treat.list))
  for (i in seq_along(treat.list)) {
    coef_ind[[i]] <- length(unlist(coef_ind)) + switch(
      treat.types[i],
      binary = seq_col(covs.list[[i]]),
      `multi-category` =,
      multinomial = seq_len((nlevels(treat.list[[i]]) - 1L) * ncol(covs.list[[i]])),
      continuous = seq_len(3L + ncol(covs.list[[i]]))
    )
  }

  squish_tol <- 50

  get_p <- get_w <- get_psi_bal <- get_par_glm <- make_list(length(treat.list))

  for (i in seq_along(treat.list)) {
    if (treat.types[i] == "binary") {
      get_p[[i]] <- function(B, X, A) {
        plogis(drop(X %*% B))
      }

      get_w[[i]] <- function(p, A, B) {
        w <- numeric(length(A))
        w1 <- which(A == 1)

        w[w1] <- 1 / p[w1]
        w[-w1] <- 1 / (1 - p[-w1])

        w
      }

      get_psi_bal[[i]] <- function(w, B, X, A, SW) {
        SW * w * (A - (1 - A)) * X
      }

      get_par_glm[[i]] <- function(X, A, SW) {
        glm.fit(X, A, family = quasibinomial(),
                weights = SW)$coefficients
      }
    }
    else if (treat.types[i] %in% c("multinomial", "multi-category")) {
      get_p[[i]] <- function(B, X, A) {
        qq <- exp(X %*% matrix(B, nrow = ncol(X)))

        pp <- cbind(1, qq) / (1 + rowSums(qq))

        colnames(pp) <- levels(A)

        pp
      }

      get_w[[i]] <- function(p, A, B) {
        w <- numeric(length(A))
        for (a in levels(A)) {
          wa <- which(A == a)
          w[wa] <- 1 / p[wa, a]
        }
        w
      }

      get_psi_bal[[i]] <- function(w, B, X, A, SW) {
        do.call("cbind", lapply(utils::combn(levels(treat.list[[i]]), 2L, simplify = FALSE), function(co) {
          SW * w * ((A == co[1L]) - (A == co[2L])) * X
        }))
      }

      get_par_glm[[i]] <- function(X, A, SW) {
        .multinom_weightit.fit(X, A, hess = FALSE,
                               weights = SW)$coefficients
      }
    }
    else if (treat.types[i] == "continuous") {
      get_p[[i]] <- function(B, X, A) {
        drop(X %*% B[-(1:3)])
      }

      get_w[[i]] <- function(p, A, B) {
        un_s2 <- exp(B[1L])
        un_p <- B[2L]

        log.dens.num <- squish(dnorm(A, un_p, sqrt(un_s2), log = TRUE),
                               lo = -Inf, hi = squish_tol)

        s2 <- exp(B[3L])
        log.dens.denom <- squish(dnorm(A, p, sqrt(s2), log = TRUE),
                                 lo = -squish_tol, hi = Inf)

        exp(log.dens.num - log.dens.denom)
      }

      get_psi_bal[[i]] <- function(w, B, X, A, SW) {
        un_s2 <- exp(B[1L])
        un_p <- B[2L]
        s2 <- exp(B[3L])
        p <- drop(X %*% B[-(1:3)])

        cbind(
          SW * (A - un_p)^2 - un_s2,
          SW * (A - un_p),
          SW * (A - p)^2 - s2,
          SW * w * A * X)
      }

      get_par_glm[[i]] <- function(X, A, SW) {
        init.fit <- lm.wfit(X, A, w = SW)
        b <- c(0, 0, log(var(init.fit$residuals)), init.fit$coefficients)
        names(b)[1:3] <- c("log(s^2)", "E[A]", "log(s_r^2)")
        b
      }
    }
  }

  # Balance condition
  psi_bal <- function(B, X.list, A.list, SW) {
    w <- Reduce("*", lapply(seq_along(A.list), function(i) {
      Bi <- B[coef_ind[[i]]]
      p <- get_p[[i]](Bi, X.list[[i]], A.list[[i]])
      get_w[[i]](p, A.list[[i]], Bi)
    }), init = 1)

    do.call("cbind", lapply(seq_along(A.list), function(i) {
      get_psi_bal[[i]](w, B[coef_ind[[i]]], X.list[[i]], A.list[[i]], SW)
    }))
  }

  obj_bal <- function(B, X.list, A.list, SW) {
    gbar <- colMeans(psi_bal(B, X.list, A.list, SW))
    sqrt(sum(gbar^2))
  }

  # Initialize coefs using GLMs
  par_glm <- unlist(lapply(seq_along(treat.list), function(i) {
    get_par_glm[[i]](covs.list[[i]], treat.list[[i]], s.weights)
  }))

  # Slightly improve glm coefs to move closer to optimal
  alpha.func <- function(alpha) obj_bal(par_glm * alpha, covs.list, treat.list, s.weights)
  par_alpha <- par_glm * optimize(alpha.func, interval = c(.8, 1.1))$min

  if (solver == "multiroot") {
    out <- suppressWarnings({
      try(verbosely({
        rootSolve::multiroot(f = function(...) colMeans(psi_bal(...)),
                             start = par_alpha,
                             X.list = covs.list,
                             A.list = treat.list,
                             SW = s.weights,
                             rtol = reltol,
                             atol = reltol,
                             ctol = reltol)
      }, verbose = FALSE), silent = TRUE)
    })

    if (!null_or_error(out) && out$estim.precis < 1e-5) {
      par_alpha <- out$root
    }
  }

  # Optimize balance objective
  out <- optim(par = par_alpha,
               fn = obj_bal,
               method = "BFGS",
               control = list(maxit = maxit,
                              reltol = reltol,
                              trace = as.integer(!over && verbose)),
               X.list = covs.list,
               A.list = treat.list,
               SW = s.weights)

  par_out <- out$par

  if (over) {
    get_psi_glm <- make_list(length(treat.list))

    for (i in seq_along(treat.list)) {
      if (treat.types[i] == "binary") {
        get_psi_glm[[i]] <- function(p, X, A, SW) {
          SW * (A - p) * X
        }
      }
      else if (treat.types[i] %in% c("multinomial", "multi-category")) {
        get_psi_glm[[i]] <- function(p, X, A, SW) {
          do.call("cbind", lapply(levels(A), function(i) {
            SW * ((A == i) - p[, i]) * X
          }))
        }
      }
      else if (treat.types[i] == "continuous") {
        get_psi_glm[[i]] <- function(p, X, A, SW) {
          SW * (A - p) * X
        }
      }
    }

    psi <- function(B, X.list, A.list, SW) {
      p.list <- lapply(seq_along(A.list), function(i) {
        get_p[[i]](B[coef_ind[[i]]], X.list[[i]], A.list[[i]])
      })

      w <- Reduce("*", lapply(seq_along(A.list), function(i) {
        get_w[[i]](p.list[[i]], A.list[[i]], B[coef_ind[[i]]])
      }), init = 1)

      cbind(
        do.call("cbind", lapply(seq_along(A.list), function(i) {
          get_psi_bal[[i]](w, B[coef_ind[[i]]], X.list[[i]], A.list[[i]], SW)
        })),
        do.call("cbind", lapply(seq_along(A.list), function(i) {
          get_psi_glm[[i]](p.list[[i]], X.list[[i]], A.list[[i]], SW)
        }))
      )
    }

    Sigma <- function(B, X.list, A.list, SW) {
      crossprod(psi(B, X.list, A.list, SW))
    }

    obj <- function(B, X.list, A.list, SW, invS = NULL) {
      psi0 <- psi(B, X.list, A.list, SW)

      if (is_null(invS)) {
        invS <- generalized_inverse(crossprod(psi0))
      }

      gbar <- colMeans(psi0)
      sqrt(drop(t(gbar) %*% invS %*% gbar))
    }

    invS <- {
      if (twostep) generalized_inverse(Sigma(par_alpha, covs.list, treat.list, s.weights))
      else NULL
    }

    start.list <- {
      if (max(abs(par_alpha - par_out)) < 1e-6) list(par_out)
      else list(par_alpha, par_out)
    }

    out <- lapply(start.list, function(par_) {
      optim(par = par_,
            fn = obj,
            method = "BFGS",
            control = list(maxit = maxit,
                           reltol = reltol,
                           trace = as.integer(verbose)),
            X.list = covs.list,
            A.list = treat.list,
            SW = s.weights,
            invS = invS)
    })

    out <- out[[which.min(unlist(grab(out, "value")))]]

    par_out <- out$par
  }

  if (out$converge != 0) {
    .wrn("the optimization failed to converge; try again with a higher value of `maxit`")
  }

  w <- Reduce("*", lapply(seq_along(treat.list), function(i) {
    Bi <- par_out[coef_ind[[i]]]
    p <- get_p[[i]](Bi, covs.list[[i]], treat.list[[i]])
    get_w[[i]](p, treat.list[[i]], Bi)
  }), init = 1)

  Mparts <- NULL

  list(w = w, fit.obj = out,
       Mparts = Mparts)
}

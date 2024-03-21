#' Covariate Balancing Propensity Score Weighting
#' @name method_cbps
#' @aliases method_cbps
#' @usage NULL
#'
#' @description
#'
#' This page explains the details of estimating weights from covariate balancing propensity scores by setting `method = "cbps"` in the call to [weightit()] or [weightitMSM()]. This method can be used with binary, multi-category, and continuous treatments.
#'
#' In general, this method relies on estimating propensity scores using generalized method of moments and then converting those propensity scores into weights using a formula that depends on the desired estimand. This method relies on code written for \pkg{WeightIt} using [optim()].
#'
#' ## Binary Treatments
#'
#' For binary treatments, this method estimates the propensity scores and weights using `optim()` using formulas described by Imai and Ratkovic (2014). The following estimands are allowed: ATE, ATT, and ATC.
#'
#' ## Multi-Category Treatments
#'
#' For multi-category treatments, this method estimates the generalized propensity scores and weights using `optim()` using formulas described by Imai and Ratkovic (2014). The following estimands are allowed: ATE and ATT.
#'
#' ## Continuous Treatments
#'
#' For continuous treatments, this method estimates the generalized propensity scores and weights using `optim()` using formulas described by Fong, Hazelett, and Imai (2018).
#'
#' ## Longitudinal Treatments
#'
#' For longitudinal treatments, the weights are the product of the weights estimated at each time point. This is not how \pkgfun{CBPS}{CBMSM} in the \pkg{CBPS} package estimates weights for longitudinal treatments.
#'
#' ## Sampling Weights
#'
#' Sampling weights are supported through `s.weights` in all scenarios.
#'
#' ## Missing Data
#'
#' In the presence of missing data, the following value(s) for `missing` are allowed:
#'     \describe{
#'       \item{`"ind"` (default)}{
#'         First, for each variable with missingness, a new missingness indicator variable is created which takes the value 1 if the original covariate is `NA` and 0 otherwise. The missingness indicators are added to the model formula as main effects. The missing values in the covariates are then replaced with the covariate medians (this value is arbitrary and does not affect estimation). The weight estimation then proceeds with this new formula and set of covariates. The covariates output in the resulting `weightit` object will be the original covariates with the `NA`s.
#'       }
#'     }
#'
#' ## M-estimation
#'
#' M-estimation is supported for the just-identified CBPS (the default, setting `over = FALSE`) for all scenarios. See [glm_weightit()] and `vignette("estimating-effects")` for details.
#'
#' @section Additional Arguments:
#'
#' The following additional arguments can be specified:
#'   \describe{
#'     \item{`over`}{`logical`; whether to request the over-identified CBPS, which combines the generalized linear model regression score equations (for binary treatments), multinomial logistic regression score equations (for multi-category treatments), or linear regression score equations (for continuous treatments) to the balance moment conditions. Default is `FALSE` to use the just-identified CBPS.
#'     }
#'     \item{`twostep`}{`logical`; when `over = TRUE`, whether to use the two-step approximation to the generalized method of moments variance. Default is `TRUE`. Ignored when `over = FALSE`.
#'     }
#'     \item{`link`}{`"string"`; the link used in the generalized linear model for the propensity scores when treatment is binary. Default is `"logit"` for logistic regression, which is used in the original description of the method by Imai and Ratkovic (2014), but others are allowed: `"logit"`, `"probit"`, `"cauchit"`, and `"cloglog"` all use the binomial likelihood, `"log"` uses the Poisson likelihood, and `"identity"` uses the Gaussian likelihood (i.e., the linear probability model). Note that negative weights are possible with these last two and they should be used with caution. Ignored for multi-category and continuous treatments.
#'     }
#'     \item{`reltol`}{the relative tolerance for convergence of the optimization. Passed to the `control` argument of `optim()`. Default is `sqrt(.Machine$double.eps)`.
#'     }
#'     \item{`maxit`}{the maximum number of iterations for convergence of the optimization. Passed to the `control` argument of `optim()`. Default is 1000.
#'     }
#' }
#'
#' @section Additional Outputs:
#' \describe{
#'   \item{`obj`}{When `include.obj = TRUE`, the output of the final call to `optim()` used to produce the model parameters. Note that because of variable transformations, the resulting parameter estimates may not be interpretable.
#'   }
#' }
#'
#' @details
#' CBPS estimates the coefficients of a generalized linear model (for binary treatments), multinomial logistic regression model (for multi-category treatments), or linear regression model (for continuous treatments) that is used to compute (generalized) propensity scores, from which the weights are computed. It involves replacing (or augmenting, in the case of the over-identified version) the standard regression score equations with the balance constraints in a generalized method of moments estimation. The idea is to nudge the estimation of the coefficients toward those that produce balance in the weighted sample. The just-identified version (with `exact = FALSE`) does away with the score equations for the coefficients so that only the balance constraints (and the score equation for the variance of the error with a continuous treatment) are used. The just-identified version will therefore produce superior balance on the means (i.e., corresponding to the balance constraints) for binary and multi-category treatments and linear terms for continuous treatments than will the over-identified version.
#'
#' Just-identified CBPS is very similar to entropy balancing and inverse probability tilting. For the ATT, all three methods will yield identical estimates. For other estimands, the results will differ.
#'
#' Note that \pkg{WeightIt} provides different functionality from the \pkg{CBPS} package in terms of the versions of CBPS available; for extensions to CBPS (e.g., optimal CBPS, CBPS for instrumental variables, and jointly estimated CBPS for longitudinal treatments), the \pkg{CBPS} package may be preferred.
#'
#' @note
#' This method used to rely on functionality in the \pkg{CBPS} package, but no longer does. Slight differences may be found between the two packages in some cases due to numerical imprecision. \pkg{WeightIt} supports arbitrary numbers of groups for the multi-category CBPS and any estimand, whereas \pkg{CBPS} only supports up to four groups and only the ATE. For continuous treatments with the over-identified CBPS, \pkg{WeightIt} and \pkg{CBPS} use different methods of specifying the GMM variance matrix, which may lead to differing results. Note that the default method differs between the two implementations; by default \pkg{WeightIt} uses the just-identified CBPS, which is faster to fit, yields better balance, and is compatible with M-estimation for estimating the standard error of the treatment effect, whereas \pkg{CBPS} uses the over-identified CBPS by default. However, both the just-identified and over-identified versions are available in both packages.
#'
#' @seealso
#' [weightit()], [weightitMSM()]
#'
#' [method_ebal] and [method_ipt] for entropy balancing and inverse probability tilting, which work similarly.
#'
#' @references
#' ## Binary treatments
#'
#' Imai, K., & Ratkovic, M. (2014). Covariate balancing propensity score. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 76(1), 243–263.
#'
#' ## Multi-Category Treatments
#'
#' Imai, K., & Ratkovic, M. (2014). Covariate balancing propensity score. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 76(1), 243–263.
#'
#' ## Continuous treatments
#'
#' Fong, C., Hazlett, C., & Imai, K. (2018). Covariate balancing propensity score for a continuous treatment: Application to the efficacy of political advertisements. The Annals of Applied Statistics, 12(1), 156–177. \doi{10.1214/17-AOAS1101}
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
#' summary(W1a)
#' cobalt::bal.tab(W1a)
#'
#' #Balancing covariates between treatment groups (binary)
#' #using over-identified CBPS with probit link
#' (W1b <- weightit(treat ~ age + educ + married +
#'                   nodegree + re74, data = lalonde,
#'                 method = "cbps", estimand = "ATT",
#'                 over = TRUE, link = "probit"))
#' summary(W1b)
#' cobalt::bal.tab(W1b)
#'
#' #Balancing covariates with respect to race (multi-category)
#' (W2 <- weightit(race ~ age + educ + married +
#'                   nodegree + re74, data = lalonde,
#'                 method = "cbps", estimand = "ATE"))
#' summary(W2)
#' cobalt::bal.tab(W2)
#'
#' #Balancing covariates with respect to re75 (continuous)
#' (W3 <- weightit(re75 ~ age + educ + married +
#'                   nodegree + re74, data = lalonde,
#'                 method = "cbps"))
#' summary(W3)
#' cobalt::bal.tab(W3)
NULL

weightit2cbps <- function(covs, treat, s.weights, estimand, focal, subset,
                          stabilize, subclass, missing, moments, int, verbose, ...) {

  A <- list(...)

  covs <- covs[subset, , drop = FALSE]
  treat <- treat[subset]
  s.weights <- s.weights[subset]

  if (missing == "ind") {
    covs <- add_missing_indicators(covs)
  }

  covs <- cbind(covs, int.poly.f(covs, poly = moments, int = int, center = TRUE))

  covs <- cbind(covs, quantile_f(covs, qu = A[["quantile"]], s.weights = s.weights,
                                 focal = focal, treat = treat))

  t.lev <- get_treated_level(treat)
  treat <- binarize(treat, one = t.lev)

  colinear.covs.to.remove <- setdiff(colnames(covs), colnames(make_full_rank(covs)))
  covs <- covs[, colnames(covs) %nin% colinear.covs.to.remove, drop = FALSE]

  mod_covs <- svd(cbind(`(Intercept)` = 1, scale(covs)))$u
  bal_covs <- svd(cbind(`(Intercept)` = 1, scale(covs)))$u

  over <- if (is_null(A$over)) FALSE else A$over
  chk::chk_logical(over)

  twostep <- if (is_null(A$twostep)) TRUE else A$twostep
  if (over) chk::chk_logical(twostep)

  reltol <- if (is_null(A$reltol)) sqrt(.Machine$double.eps) else A$reltol
  chk::chk_number(reltol)

  maxit <- if (is_null(A$maxit)) 1e3 else A$maxit
  chk::chk_count(maxit)

  N <- sum(s.weights)

  if (is_null(A$link)) A$link <- "logit"
  link <- A$link
  chk::chk_string(link)
  link <- match_arg(link, c("logit", "probit", "cauchit", "cloglog", "log", "identity"))

  .fam <- switch(link,
                 "log" = quasipoisson("log"),
                 "identity" = gaussian("identity"),
                 quasibinomial(link))

  #Generalized linear model score
  psi_glm <- function(B, Xm, A, SW) {
    lin_pred <- drop(Xm %*% B)
    p <- .fam$linkinv(lin_pred)
    (SW * (A - p) * .fam$mu.eta(lin_pred) / .fam$variance(p)) * Xm
  }

  if (estimand == "ATE") {
    # Balance condition for ATE
    psi_bal <- function(B, Xm, Xb = Xm, A, SW) {
      p <- .fam$linkinv(drop(Xm %*% B))
      SW * (A/p - (1-A)/(1-p)) * Xb
    }

    obj_bal <- function(B, Xm, Xb = Xm, A, SW) {
      gbar <- colMeans(psi_bal(B, Xm, Xb, A, SW))
      sqrt(sum(gbar^2))
    }

    # Initialize coefs using logistic regression
    par_glm <- glm.fit(mod_covs, treat, family = .fam, weights = s.weights)$coefficients

    # Slightly improve glm coefs to move closer to optimal
    alpha.func <- function(alpha) obj_bal(par_glm * alpha, mod_covs, bal_covs, treat, s.weights)
    par_alpha <- par_glm * optimize(alpha.func, interval = c(.8, 1.1))$min

    # Optimize balance objective
    out <- optim(par = par_alpha,
                 fn = obj_bal,
                 method = "BFGS",
                 control = list(maxit = maxit,
                                reltol = reltol),
                 Xm = mod_covs,
                 Xb = bal_covs,
                 A = treat,
                 SW = s.weights)

    par <- out$par

    if (over) {
      # Combine LR and balance
      psi <- function(B, Xm, Xb = Xm, A, SW) {
        cbind(psi_glm(B, Xm, A, SW), psi_bal(B, Xm, Xb, A, SW))
      }

      Sigma <- function(B, Xm, Xb = Xm, A, SW) {
        lp <- drop(Xm %*% B)
        p <- .fam$linkinv(lp)
        sw.5 <- sqrt(SW)
        g <- .fam$mu.eta(lp) / .fam$variance(p)

        S11 <- crossprod(sw.5 * g * sqrt(p * (1 - p)) * Xm)
        S12 <- crossprod(sw.5 * g * Xm, sw.5 * Xb)
        S21 <- t(S12)
        S22 <- crossprod(sw.5 / sqrt(p * (1 - p)) * Xb)

        rbind(cbind(S11, S12),
              cbind(S21, S22)) / N
      }

      obj <- function(B, Xm, Xb = Xm, A, SW, invS = NULL) {
        if (is.null(invS)) {
          invS <- generalized_inverse(Sigma(B, Xm, Xb, A, SW))
        }

        gbar <- colMeans(psi(B, Xm, Xb, A, SW))
        sqrt(drop(t(gbar) %*% invS %*% gbar))
      }

      invS <- {
        if (twostep) generalized_inverse(Sigma(par_alpha, mod_covs, bal_covs, treat, s.weights))
        else NULL
      }

      out <- lapply(list(par_alpha, par), function(par_) {
        optim(par = par_,
              fn = obj,
              method = "BFGS",
              control = list(maxit = maxit,
                             reltol = reltol),
              Xm = mod_covs,
              Xb = bal_covs,
              A = treat,
              SW = s.weights,
              invS = invS)
      })

      out <- out[[which.min(sapply(out, `[[`, "value"))]]

      par <- out$par
    }

    p.score <- .fam$linkinv(drop(mod_covs %*% par))
  }
  else {
    # Balance condition for ATT
    psi_bal <- switch(estimand,
                      "ATT" = function(B, Xm, Xb = Xm, A, SW) {
                        p <- .fam$linkinv(drop(Xm %*% B))
                        SW * ((A - p) / (1 - p)) * Xb
                      },
                      "ATC" = function(B, Xm, Xb = Xm, A, SW) {
                        p <- .fam$linkinv(drop(Xm %*% B))
                        SW * ((A - p) / p) * Xb
                      })

    obj_bal <- function(B, Xm, Xb = Xm, A, SW) {
      gbar <- colMeans(psi_bal(B, Xm, Xb, A, SW))
      sqrt(sum(gbar^2))
    }

    # Initialize coefs using logistic regression
    par_glm <- glm.fit(mod_covs, treat, family = .fam, weights = s.weights)$coefficients

    # Slightly improve glm coefs to move closer to optimal
    alpha.func <- function(alpha) obj_bal(par_glm * alpha, mod_covs, bal_covs, treat, s.weights)
    par_alpha <- par_glm * optimize(alpha.func, interval = c(.8, 1.1))$min

    # Optimize balance objective
    out <- optim(par = par_alpha,
                 fn = obj_bal,
                 method = "BFGS",
                 control = list(maxit = maxit,
                                reltol = reltol),
                 Xm = mod_covs,
                 Xb = bal_covs,
                 A = treat,
                 SW = s.weights)

    par <- out$par

    if (over) {
      # Combine LR and balance
      psi <- function(B, Xm, Xb = Xm, A, SW) {
        cbind(psi_glm(B, Xm, A, SW), psi_bal(B, Xm, Xb, A, SW))
      }

      Sigma <- function(B, Xm, Xb = Xm, A, SW) {
        lp <- drop(Xm %*% B)
        p <- .fam$linkinv(lp)
        sw.5 <- sqrt(SW)
        g <- .fam$mu.eta(lp) / .fam$variance(p)

        S11 <- crossprod(sw.5 * g * sqrt(p * (1 - p)) * Xm)
        S12 <- switch(estimand,
                      "ATT" = crossprod(sw.5 * g * Xm, sw.5 * p * Xb),
                      "ATC" = crossprod(sw.5 * g * Xm, sw.5 * (1 - p) * Xb))
        S21 <- t(S12)
        S22 <- switch(estimand,
                      "ATT" = crossprod(sw.5 * sqrt(p / (1 - p)) * Xb),
                      "ATC" = crossprod(sw.5 * sqrt((1 - p) / p) * Xb))

        rbind(cbind(S11, S12),
              cbind(S21, S22)) / N
      }

      obj <- function(B, Xm, Xb = Xm, A, SW, invS = NULL) {
        if (is.null(invS)) {
          invS <- generalized_inverse(Sigma(B, Xm, Xb, A, SW))
        }

        gbar <- colMeans(psi(B, Xm, Xb, A, SW))
        sqrt(drop(t(gbar) %*% invS %*% gbar))
      }

      invS <- {
        if (twostep) generalized_inverse(Sigma(par_alpha, mod_covs, bal_covs, treat, s.weights))
        else NULL
      }

      out <- lapply(list(par_alpha, par), function(par_) {
        optim(par = par_,
              fn = obj,
              method = "BFGS",
              control = list(maxit = maxit,
                             reltol = reltol),
              Xm = mod_covs,
              Xb = bal_covs,
              A = treat,
              SW = s.weights,
              invS = invS)
      })

      out <- out[[which.min(sapply(out, `[[`, "value"))]]

      par <- out$par
    }

    p.score <- .fam$linkinv(drop(mod_covs %*% par))
  }

  if (out$converge != 0) {
    .wrn("the optimziation failed to converge; try again with a higher value of `maxit`")
  }

  w <- .get_w_from_ps_internal_bin(p.score, treat, estimand = estimand,
                                   subclass = subclass,
                                   stabilize = stabilize)

  Mparts <- NULL
  if (!over) {
    Mparts <- list(
      psi_treat = function(Btreat, A, Xtreat, SW) {
          psi_bal(Btreat, Xtreat, Xtreat, A, SW)
      },
      wfun = function(Btreat, Xtreat, A) {
        ps <- .fam$linkinv(drop(Xtreat %*% Btreat))
        .get_w_from_ps_internal_bin(ps, A, estimand = estimand,
                                    subclass = subclass, stabilize = stabilize)
      },
      Xtreat = mod_covs,
      A = treat,
      btreat = par
    )
  }

  list(w = w, ps = p.score, fit.obj = out,
       Mparts = Mparts)
}

weightit2cbps.multi <- function(covs, treat, s.weights, estimand, focal, subset,
                                stabilize, subclass, missing, verbose, ...) {

  A <- list(...)

  covs <- covs[subset, , drop = FALSE]
  treat <- factor(treat[subset])
  s.weights <- s.weights[subset]

  if (missing == "ind") {
    covs <- add_missing_indicators(covs)
  }

  colinear.covs.to.remove <- setdiff(colnames(covs), colnames(make_full_rank(covs)))
  covs <- covs[, colnames(covs) %nin% colinear.covs.to.remove, drop = FALSE]

  mod_covs <- svd(cbind(`(Intercept)` = 1, scale(covs)))$u
  bal_covs <- svd(cbind(`(Intercept)` = 1, scale(covs)))$u

  over <- if (is_null(A$over)) FALSE else A$over
  chk::chk_logical(over)

  twostep <- if (is_null(A$twostep)) TRUE else A$twostep
  if (over) chk::chk_logical(twostep)

  reltol <- if (is_null(A$reltol)) sqrt(.Machine$double.eps) else A$reltol
  chk::chk_number(reltol)

  maxit <- if (is_null(A$maxit)) 1e3 else A$maxit
  chk::chk_count(maxit)

  N <- sum(s.weights)

  treat_num <- as.integer(treat)

  K <- nlevels(treat)
  kk <- seq_len(K)

  coef_ind <- setNames(lapply(kk[-K], function(i) {
    (i - 1) * ncol(mod_covs) + seq_len(ncol(mod_covs))
  }), levels(treat)[-K])

  get_pp <- function(B, Xm) {
    qq <- lapply(kk[-K], function(i) {
      exp(drop(Xm %*% B[coef_ind[[i]]]))
    })

    pden <- 1 + rowSums(do.call("cbind", qq))

    out <- do.call("cbind", c(qq, list(1))) / pden
    colnames(out) <- levels(treat)
    out
  }

  #Multinomial logistic regression score
  psi_mlr <- function(B, Xm, A, SW) {
    pp <- get_pp(B, Xm)

    do.call("cbind", lapply(levels(treat), function(i) {
      SW * ((A == i) - pp[,i]) * Xm
    }))
  }

  if (estimand == "ATE") {
    combs <- utils::combn(levels(treat), 2, simplify = FALSE)
    # combs <- lapply(kk, function(i) c(i, i - 1))

    # Balance condition for ATE
    psi_bal <- function(B, Xm, Xb = Xm, A, SW) {
      pp <- get_pp(B, Xm)

      do.call("cbind", lapply(combs, function(co) {
        SW * ((A == co[1]) / pp[,co[1]] - (A == co[2]) / pp[,co[2]]) * Xb
      }))
    }

    obj_bal <- function(B, Xm, Xb = Xm, A, SW) {
      gbar <- colMeans(psi_bal(B, Xm, Xb, A, SW))
      sqrt(sum(gbar^2))
    }

    # Initialize coefs using logistic regression
    par_glm <- unlist(lapply(kk[-K], function(i) {
      glm.fit(mod_covs, treat_num == i, family = binomial(), weights = s.weights)$coefficients
    }))

    # Slightly improve glm coefs to move closer to optimal
    alpha.func <- function(alpha) obj_bal(par_glm * alpha, Xm = mod_covs, Xb = bal_covs,
                                          A = treat, SW = s.weights)
    par_alpha <- par_glm * optimize(alpha.func, interval = c(.8, 1.1))$min

    # Optimize balance objective
    out <- optim(par = par_alpha,
                 fn = obj_bal,
                 method = "BFGS",
                 control = list(maxit = maxit,
                                reltol = reltol),
                 Xm = mod_covs,
                 Xb = bal_covs,
                 A = treat,
                 SW = s.weights)

    par <- out$par

    if (over) {
      # Combine LR and balance
      psi <- function(B, Xm, Xb = Xm, A, SW) {
        cbind(psi_mlr(B, Xm, A, SW), psi_bal(B, Xm, Xb, A, SW))
      }

      Sigma <- function(B, Xm, Xb = Xm, A, SW) {
        pp <- get_pp(B, Xm)
        sw.5 <- sqrt(SW)

        S <- list()

        for (i in levels(treat)) {
          for (j in levels(treat)) {
            S[[sprintf("m%s_m%s", i, j)]] <- {
              if      (i == j) crossprod(sw.5 * sqrt(pp[,i] * (1 - pp[,i])) *  Xm)
              else if (i < j) -crossprod(sw.5 * sqrt(pp[,i] * pp[,j]) * Xm)
              else t(S[[sprintf("m%s_m%s", j, i)]])
            }
          }
        }

        for (i in levels(treat)) {
          for (jj in combs) {
            m <- {
              if      (i == jj[1])  crossprod(sw.5 * Xm, sw.5 * Xb)
              else if (i == jj[2]) -crossprod(sw.5 * Xm, sw.5 * Xb)
              else matrix(0, ncol(Xm), ncol(Xb))
            }
            S[[sprintf("m%s_b%s%s", i, jj[1], jj[2])]] <- m
            S[[sprintf("b%s%s_m%s", jj[1], jj[2], i)]] <- t(m)
          }
        }

        for (ii in combs) {
          for (jj in combs) {
            m <- {
              if (identical(ii, jj))    crossprod(sw.5 * sqrt(1 / pp[,ii[1]] + 1 / pp[,ii[2]]) * Xb)
              else if (ii[1] == jj[1])  crossprod(sw.5 * sqrt(1 / pp[,ii[1]]) * Xb)
              else if (ii[1] == jj[2]) -crossprod(sw.5 * sqrt(1 / pp[,ii[1]]) * Xb)
              else if (ii[2] == jj[1]) -crossprod(sw.5 * sqrt(1 / pp[,ii[2]]) * Xb)
              else if (ii[2] == jj[2])  crossprod(sw.5 * sqrt(1 / pp[,ii[2]]) * Xb)
              else matrix(0, ncol(Xb), ncol(Xb))
            }
            S[[sprintf("b%s%s_b%s%s", ii[1], ii[2], jj[1], jj[2])]] <- m
            S[[sprintf("b%s%s_b%s%s", jj[1], jj[2], ii[1], ii[2])]] <- m
          }
        }

        nms <- c(sprintf("m%s", levels(treat)), sprintf("b%s", vapply(combs, paste0, character(1L), collapse = "")))

        do.call("rbind", lapply(nms, function(i) {
          do.call("cbind", lapply(nms, function(j) {
            S[[sprintf("%s_%s", i, j)]]
          }))
        })) / N
      }

      obj <- function(B, Xm, Xb = Xm, A, SW, invS = NULL) {
        if (is.null(invS)) {
          invS <- generalized_inverse(Sigma(B, Xm, Xb, A, SW))
        }

        gbar <- colMeans(psi(B, Xm, Xb, A, SW))

        sqrt(drop(t(gbar) %*% invS %*% gbar))
      }

      invS <- {
        if (twostep) generalized_inverse(Sigma(par_alpha, mod_covs, bal_covs, treat_num, s.weights))
        else NULL
      }

      out <- lapply(list(par_alpha, par), function(par_) {
        optim(par = par_,
              fn = obj,
              method = "BFGS",
              control = list(maxit = maxit,
                             reltol = reltol),
              Xm = mod_covs,
              Xb = bal_covs,
              A = treat,
              SW = s.weights,
              invS = invS)
      })

      out <- out[[which.min(sapply(out, `[[`, "value"))]]

      par <- out$par
    }
  }
  else {
    # combs <- lapply(setdiff(levels(treat), focal), function(i) c(focal, i))
    combs <- utils::combn(levels(treat), 2, simplify = FALSE)

    # Balance condition for ATT
    psi_bal <- function(B, Xm, Xb = Xm, A, SW) {
      pp <- get_pp(B, Xm)

      do.call("cbind", lapply(combs, function(co) {
        SW * pp[,focal] * ((A == co[1]) / pp[,co[1]] - (A == co[2]) / pp[,co[2]]) * Xb
      }))
    }

    obj_bal <- function(B, Xm, Xb = Xm, A, SW) {
      gbar <- colMeans(psi_bal(B, Xm, Xb, A, SW))
      sqrt(sum(gbar^2)) #use sqrt to improve convergence
    }

    # Initialize coefs using logistic regression
    par_glm <- unlist(lapply(kk[-K], function(i) {
      glm.fit(mod_covs, treat_num == i, family = binomial(), weights = s.weights)$coefficients
    }))

    # Slightly improve glm coefs to move closer to optimal
    alpha.func <- function(alpha) obj_bal(par_glm * alpha, Xm = mod_covs, Xb = bal_covs,
                                          A = treat, SW = s.weights)
    par_alpha <- par_glm * optimize(alpha.func, interval = c(.8, 1.1))$min

    # Optimize balance objective
    out <- optim(par = par_alpha,
                 fn = obj_bal,
                 method = "BFGS",
                 control = list(maxit = maxit,
                                reltol = reltol),
                 Xm = mod_covs,
                 Xb = bal_covs,
                 A = treat,
                 SW = s.weights)

    par <- out$par

    if (over) {
      # Combine LR and balance
      psi <- function(B, Xm, Xb = Xm, A, SW) {
        cbind(psi_mlr(B, Xm, A, SW), psi_bal(B, Xm, Xb, A, SW))
      }

      Sigma <- function(B, Xm, Xb = Xm, A, SW) {
        pp <- get_pp(B, Xm)
        sw.5 <- sqrt(SW)

        S <- list()

        for (i in levels(treat)) {
          for (j in levels(treat)) {
            S[[sprintf("m%s_m%s", i, j)]] <- {
              if      (i == j) crossprod(sw.5 * sqrt(pp[,i] * (1 - pp[,i])) *  Xm)
              else if (i < j) -crossprod(sw.5 * sqrt(pp[,i] * pp[,j]) * Xm)
              else t(S[[sprintf("m%s_m%s", j, i)]])
            }
          }
        }

        for (i in levels(treat)) {
          for (jj in combs) {
            m <- {
              if      (i == jj[1])  crossprod(sw.5 * Xm, sw.5 * pp[,focal] * Xb)
              else if (i == jj[2]) -crossprod(sw.5 * Xm, sw.5 * pp[,focal] * Xb)
              else matrix(0, ncol(Xm), ncol(Xb))
            }
            S[[sprintf("m%s_b%s%s", i, jj[1], jj[2])]] <- m
            S[[sprintf("b%s%s_m%s", jj[1], jj[2], i)]] <- t(m)
          }
        }

        for (ii in combs) {
          for (jj in combs) {
            m <- {
              if (identical(ii, jj))    crossprod(sw.5 * sqrt(1 / pp[,ii[1]] + 1 / pp[,ii[2]]) * pp[,focal] * Xb)
              else if (ii[1] == jj[1])  crossprod(sw.5 * sqrt(1 / pp[,ii[1]]) * pp[,focal] * Xb)
              else if (ii[1] == jj[2]) -crossprod(sw.5 * sqrt(1 / pp[,ii[1]]) * pp[,focal] * Xb)
              else if (ii[2] == jj[1]) -crossprod(sw.5 * sqrt(1 / pp[,ii[2]]) * pp[,focal] * Xb)
              else if (ii[2] == jj[2])  crossprod(sw.5 * sqrt(1 / pp[,ii[2]]) * pp[,focal] * Xb)
              else matrix(0, ncol(Xb), ncol(Xb))
            }
            S[[sprintf("b%s%s_b%s%s", ii[1], ii[2], jj[1], jj[2])]] <- m
            S[[sprintf("b%s%s_b%s%s", jj[1], jj[2], ii[1], ii[2])]] <- m
          }
        }

        nms <- c(sprintf("m%s", levels(treat)), sprintf("b%s", vapply(combs, paste0, character(1L), collapse = "")))

        do.call("rbind", lapply(nms, function(i) {
          do.call("cbind", lapply(nms, function(j) {
            S[[sprintf("%s_%s", i, j)]]
          }))
        })) / N
      }

      obj <- function(B, Xm, Xb = Xm, A, SW, invS = NULL) {
        if (is.null(invS)) {
          invS <- generalized_inverse(Sigma(B, Xm, Xb, A, SW))
        }

        gbar <- colMeans(psi(B, Xm, Xb, A, SW))

        sqrt(drop(t(gbar) %*% invS %*% gbar)) #use sqrt to improve convergence
      }

      invS <- {
        if (twostep) generalized_inverse(Sigma(par_alpha, mod_covs, bal_covs, treat_num, s.weights))
        else NULL
      }

      out <- lapply(list(par_alpha, par), function(par_) {
        optim(par = par_,
              fn = obj,
              method = "BFGS",
              control = list(maxit = maxit,
                             reltol = reltol),
              Xm = mod_covs,
              Xb = bal_covs,
              A = treat,
              SW = s.weights,
              invS = invS)
      })

      out <- out[[which.min(sapply(out, `[[`, "value"))]]

      par <- out$par
    }
  }

  if (out$converge != 0) {
    .wrn("the optimziation failed to converge; try again with a higher value of `maxit`")
  }

  pp <- get_pp(par, mod_covs)

  w <- .get_w_from_ps_internal_multi(pp, treat, estimand = estimand, subclass = subclass,
                                     focal = focal, stabilize = stabilize)

  out$pp <- pp

  Mparts <- NULL
  if (!over) {
    Mparts <- list(
      psi_treat = function(Btreat, A, Xtreat, SW) {
          psi_bal(Btreat, Xtreat, Xtreat, A, SW)
        },
      wfun = function(Btreat, Xtreat, A) {
        ps <- get_pp(Btreat, Xtreat)
        .get_w_from_ps_internal_multi(ps, A, estimand = estimand, focal = focal,
                                      subclass = subclass, stabilize = stabilize)
      },
      Xtreat = mod_covs,
      A = treat,
      btreat = par
    )
  }

  list(w = w, fit.obj = out,
       Mparts = Mparts)
}

weightit2cbps.cont <- function(covs, treat, s.weights, subset, missing, verbose, ...) {

  A <- list(...)

  covs <- covs[subset, , drop = FALSE]
  treat <- treat[subset]
  s.weights <- s.weights[subset]

  if (missing == "ind") {
    covs <- add_missing_indicators(covs)
  }

  colinear.covs.to.remove <- setdiff(colnames(covs), colnames(make_full_rank(covs)))
  covs <- covs[, colnames(covs) %nin% colinear.covs.to.remove, drop = FALSE]

  treat <- (treat - w.m(treat, s.weights)) / sqrt(col.w.v(treat, s.weights))

  mod_covs <- cbind(`(Intercept)` = 1,
                    scale(covs,
                          center = col.w.m(covs, s.weights),
                          scale = sqrt(col.w.v(covs, s.weights))))

  bal_covs <- mod_covs

  over <- if (is_null(A$over)) FALSE else A$over
  chk::chk_logical(over)

  twostep <- if (is_null(A$twostep)) TRUE else A$twostep
  if (over) chk::chk_logical(twostep)

  reltol <- if (is_null(A$reltol)) sqrt(.Machine$double.eps) else A$reltol
  chk::chk_number(reltol)

  maxit <- if (is_null(A$maxit)) 1e3 else A$maxit
  chk::chk_count(maxit)

  dens.num <- dnorm(treat, log = TRUE)

  #Linear regression score; residual variance omitted, included
  #in balance condition
  psi_lm <- function(B, Xm, A, SW) {
    p <- drop(Xm %*% B[-1])

    cbind(SW * (A - p) * Xm)
  }

  # Balance condition
  psi_bal <- function(B, Xm, Xb = Xm, A, SW) {
    s2 <- exp(B[1])
    p <- drop(Xm %*% B[-1])

    dens.denom <- dnorm(A, p, sqrt(s2), log = TRUE)

    w <- exp(dens.num - dens.denom)

    cbind(w * SW * A * Xm,
          SW * (A - p)^2/s2 - 1)
  }

  obj_bal <- function(B, Xm, Xb = Xm, A, SW) {
    gbar <- colMeans(psi_bal(B, Xm, Xb, A, SW))
    sqrt(sum(gbar^2))
  }

  # Initialize coefs using logistic regression
  init.fit <- lm.wfit(mod_covs, treat, w = s.weights)
  par_glm <- c(log(var(init.fit$residuals)), init.fit$coefficients)

  # Slightly improve glm coefs to move closer to optimal
  alpha.func <- function(alpha) obj_bal(par_glm * alpha, mod_covs, bal_covs, treat, s.weights)
  par_alpha <- par_glm * optimize(alpha.func, interval = c(.8, 1.1))$min

  # Optimize balance objective
  out <- optim(par = par_alpha,
               fn = obj_bal,
               method = "BFGS",
               control = list(maxit = maxit,
                              reltol = reltol),
               Xm = mod_covs,
               Xb = bal_covs,
               A = treat,
               SW = s.weights)

  par <- out$par

  if (over) {
    # Combine LR and balance
    psi <- function(B, Xm, Xb = Xm, A, SW) {
      cbind(psi_lm(B, Xm, A, SW), psi_bal(B, Xm, Xb, A, SW))
    }

    Sigma <- function(B, Xm, Xb = Xm, A, SW) {
      crossprod(psi(B, Xm, Xb, A, SW))
    }

    obj <- function(B, Xm, Xb = Xm, A, SW, invS = NULL) {
      if (is.null(invS)) {
        invS <- generalized_inverse(Sigma(B, Xm, Xb, A, SW))
      }

      gbar <- colMeans(psi(B, Xm, Xb, A, SW))
      sqrt(drop(t(gbar) %*% invS %*% gbar))
    }

    invS <- {
      if (twostep) generalized_inverse(Sigma(par_alpha, mod_covs, bal_covs, treat, s.weights))
      else NULL
    }

    out <- lapply(list(par_alpha, par), function(par_) {
      optim(par = par_,
            fn = obj,
            method = "BFGS",
            control = list(maxit = maxit,
                           reltol = reltol),
            Xm = mod_covs,
            Xb = bal_covs,
            A = treat,
            SW = s.weights,
            invS = invS)
    })

    out <- out[[which.min(sapply(out, `[[`, "value"))]]

    par <- out$par
  }

  if (out$converge != 0) {
    .wrn("the optimziation failed to converge; try again with a higher value of `maxit`")
  }

  s2 <- exp(par[1])

  p <- drop(mod_covs %*% par[-1])
  dens.denom <- dnorm(treat, p, sqrt(s2), log = TRUE)

  w <- exp(dens.num - dens.denom)

  Mparts <- NULL
  if (!over) {
    Mparts <- list(
      psi_treat = function(Btreat, A, Xtreat, SW) {
        psi_bal(Btreat, Xtreat, Xtreat, A, SW)
      },
      wfun = function(Btreat, Xtreat, A) {
        s2 <- exp(Btreat[1])

        p <- drop(Xtreat %*% Btreat[-1])
        dens.denom <- dnorm(A, p, sqrt(s2), log = TRUE)
        dens.num <- dnorm(A, log = TRUE)

        exp(dens.num - dens.denom)
      },
      Xtreat = mod_covs,
      A = treat,
      btreat = par
    )
  }

  list(w = w, fit.obj = out,
       Mparts = Mparts)
}

weightitMSM2cbps <- function(covs.list, treat.list, s.weights, subset, missing, verbose, ...) {
  .err("CBMSM doesn't work yet")
}

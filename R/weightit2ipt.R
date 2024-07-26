#' Inverse Probability Tilting
#' @name method_ipt
#' @aliases method_ipt
#' @usage NULL
#'
#' @description
#' This page explains the details of estimating weights using inverse probability tilting by setting `method = "ipt"` in the call to [weightit()] or [weightitMSM()]. This method can be used with binary and multi-category treatments.
#'
#' In general, this method relies on estimating propensity scores using a modification of the usual generalized linear model score equations to enforce balance and then converting those propensity scores into weights using a formula that depends on the desired estimand. This method relies on code written for \pkg{WeightIt} using \pkgfun{rootSolve}{multiroot}.
#'
#' ## Binary Treatments
#'
#' For binary treatments, this method estimates the weights using formulas described by Graham, Pinto, and Egel (2012). The following estimands are allowed: ATE, ATT, and ATC. When the ATE is requested, the optimization is run twice, once for each treatment group.
#'
#' ## Multi-Category Treatments
#'
#' For multi-category treatments, this method estimates the weights using modifications of the formulas described by Graham, Pinto, and Egel (2012). The following estimands are allowed: ATE and ATT. When the ATE is requested, estimation is performed once for each treatment group. When the ATT is requested, estimation is performed once for each non-focal (i.e., control) group.
#'
#' ## Continuous Treatments
#'
#' Inverse probability tilting is not compatible with continuous treatments.
#'
#' ## Longitudinal Treatments
#'
#' For longitudinal treatments, the weights are the product of the weights estimated at each time point. This method is not guaranteed to yield exact balance at each time point. NOTE: the use of inverse probability tilting with longitudinal treatments has not been validated!
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
#' M-estimation is supported for all scenarios. See [glm_weightit()] and `vignette("estimating-effects")` for details.
#'
#' @section Additional Arguments:
#' `moments` and `int` are accepted. See [weightit()] for details.
#'
#' \describe{
#'   \item{`quantile`}{
#'     A named list of quantiles (values between 0 and 1) for each continuous covariate, which are used to create additional variables that when balanced ensure balance on the corresponding quantile of the variable. For example, setting `quantile = list(x1 = c(.25, .5. , .75))` ensures the 25th, 50th, and 75th percentiles of `x1` in each treatment group will be balanced in the weighted sample. Can also be a single number (e.g., `.5`) or an unnamed list of length 1 (e.g., `list(c(.25, .5, .75))`) to request the same quantile(s) for all continuous covariates, or a named vector (e.g., `c(x1 = .5, x2 = .75`) to request one quantile for each covariate.
#'   }
#'   \item{`link`}{`"string"`; the link used to determine the inverse link for computing the (generalized) propensity scores. Default is `"logit"`, which is used in the original description of the method by Graham, Pinto, and Egel (2012), but `"probit"`, `"cauchit"`, and `"cloglog"` are also allowed.
#'   }
#' }
#'
#' The `stabilize` argument is ignored.
#'
#' @section Additional Outputs:
#' \describe{
#'   \item{`obj`}{When `include.obj = TRUE`, the output of the call to [optim()], which contains the coefficient estimates and convergence information. For ATE fits or with multi-category treatments, a list of `rootSolve::multiroot()` outputs, one for each weighted group.
#'   }
#' }
#'
#' @details
#' Inverse probability tilting (IPT) involves specifying estimating equations that fit the parameters of one or more generalized linear models with a modification that ensures exact balance on the covariate means. These estimating equations are solved, and the estimated parameters are used in the (generalized) propensity score, which is used to compute the weights. Conceptually and mathematically, IPT is very similar to entropy balancing and just-identified CBPS. For the ATT and ATC, entropy balancing, just-identified CBPS, and IPT will yield identical results. For the ATE or when `link` is specified as something other than `"logit"`, the three methods differ.
#'
#' Treatment effect estimates for binary treatments are consistent if the true propensity score is a logistic regression or the outcome model is linear in the covariates and their interaction with treatments. For entropy balancing, this is only true for the ATT, and for just-identified CBPS, this is only true if there is no effect modification by covariates. In this way, IPT provides additional theoretical guarantees over the other two methods, though potentially with some cost in precision.
#'
#' @seealso
#' [weightit()], [weightitMSM()]
#'
#' [method_ebal] and [method_cbps] for entropy balancing and CBPS, which work similarly.
#'
#' @references
#' ## `estimand = "ATE"`
#'
#' Graham, B. S., De Xavier Pinto, C. C., & Egel, D. (2012). Inverse Probability Tilting for Moment Condition Models with Missing Data. *The Review of Economic Studies*, 79(3), 1053–1079. \doi{10.1093/restud/rdr047}
#'
#' ## `estimand = "ATT"`
#'
#' Sant'Anna, P. H. C., & Zhao, J. (2020). Doubly robust difference-in-differences estimators. *Journal of Econometrics*, 219(1), 101–122. \doi{10.1016/j.jeconom.2020.06.003}
#'
#' @examples
#' data("lalonde", package = "cobalt")
#'
#' #Balancing covariates between treatment groups (binary)
#' (W1 <- weightit(treat ~ age + educ + married +
#'                   nodegree + re74, data = lalonde,
#'                 method = "ipt", estimand = "ATT"))
#' summary(W1)
#' cobalt::bal.tab(W1)
#'
#' #Balancing covariates with respect to race (multi-category)
#' (W2 <- weightit(race ~ age + educ + married +
#'                   nodegree + re74, data = lalonde,
#'                 method = "ipt", estimand = "ATE"))
#' summary(W2)
#' cobalt::bal.tab(W2)
#'
NULL

weightit2ipt <- function(covs, treat, s.weights, subset, estimand, focal,
                         stabilize, missing, moments, int, verbose, ...) {
  A <- list(...)

  rlang::check_installed("rootSolve")

  covs <- covs[subset, , drop = FALSE]
  treat <- treat[subset]
  s.weights <- s.weights[subset]

  if (missing == "ind") {
    covs <- add_missing_indicators(covs)
  }

  covs <- cbind(covs, .int_poly_f(covs, poly = moments, int = int, center = TRUE))

  covs <- cbind(covs, .quantile_f(covs, qu = A[["quantile"]], s.weights = s.weights,
                                  focal = focal, treat = treat))

  for (i in seq_col(covs)) covs[,i] <- .make_closer_to_1(covs[,i])

  colinear.covs.to.remove <- setdiff(colnames(covs), colnames(make_full_rank(covs)))
  covs <- covs[, colnames(covs) %nin% colinear.covs.to.remove, drop = FALSE]

  C <- cbind(`(Intercept)` = 1, covs)

  t.lev <- get_treated_level(treat)
  treat <- binarize(treat, one = t.lev)

  link <- if_null_then(A$link, "logit")
  chk::chk_string(link)
  chk::chk_subset(link, c("logit", "probit", "cauchit", "cloglog"))

  .fam <- quasibinomial(link)

  groups_to_weight <- switch(estimand,
                             "ATT" = 0,
                             "ATC" = 1,
                             0:1)

  fit.list <- make_list(groups_to_weight)
  par.list <- make_list(groups_to_weight)

  n <- length(treat)
  k <- ncol(C)

  # start <- setNames(rep.int(0, k), colnames(C))
  start <- glm.fit(C, treat, weights = s.weights, family = .fam)$coefficients

  f <- function(B, X, A, SW, .psi) {
    .colMeans(.psi(B, X, A, SW), n, k)
  }

  if (estimand == "ATE") {
    ps <- rep.int(0, n)

    # Control weights
    psi0 <- function(B, X, A, SW) {
      p <- rep.int(0, n)
      p[A == 0] <- .fam$linkinv(drop(X[A == 0,, drop = FALSE] %*% B))
      SW * ((1 - A)/(1 - p) - 1) * X
    }

    verbosely({
      fit.list[["0"]] <- rootSolve::multiroot(f, start = start,
                                              X = C,
                                              A = treat,
                                              SW = s.weights,
                                              .psi = psi0,
                                              rtol = 1e-10, atol = 1e-10, ctol = 1e-10,
                                              verbose = TRUE)
    }, verbose = verbose)

    par.list[["0"]] <- fit.list[["0"]]$root

    ps[treat == 0] <- .fam$linkinv(drop(C[treat == 0,, drop = FALSE] %*% par.list[["0"]]))

    #Treated weights
    psi1 <- function(B, X, A, SW) {
      p <- rep.int(1, n)
      p[A == 1] <- .fam$linkinv(drop(X[A == 1,, drop = FALSE] %*% B))
      SW * (A/p - 1) * X
    }

    verbosely({
      fit.list[["1"]] <- rootSolve::multiroot(f, start = par.list[["0"]],
                                              X = C,
                                              A = treat,
                                              SW = s.weights,
                                              .psi = psi1,
                                              rtol = 1e-10, atol = 1e-10, ctol = 1e-10,
                                              verbose = TRUE)
    }, verbose = verbose)

    par.list[["1"]] <- fit.list[["1"]]$root

    ps[treat == 1] <- .fam$linkinv(drop(C[treat == 1,, drop = FALSE] %*% par.list[["1"]]))
  }
  else {
    psi <- switch(estimand,
                  "ATT" = function(B, X, A, SW) {
                    p <- .fam$linkinv(drop(X %*% B))
                    SW * (A - (1 - A) * p / (1 - p)) * X
                  },
                  "ATC" = function(B, X, A, SW) {
                    p <- .fam$linkinv(drop(X %*% B))
                    SW * (A * (1 - p) / p - (1 - A)) * X
                  })

    verbosely({
      fit.list[[1]] <- rootSolve::multiroot(f, start = start,
                                            A = treat,
                                            X = C,
                                            SW = s.weights,
                                            .psi = psi,
                                            rtol = 1e-10, atol = 1e-10, ctol = 1e-10,
                                            verbose = TRUE)
    }, verbose = verbose)

    par.list[[1]] <- fit.list[[1]]$root

    ps <- .fam$linkinv(drop(C %*% par.list[[1]]))
  }

  w <- .get_w_from_ps_internal_bin(ps, treat, estimand = estimand)

  Mparts <- list(
    psi_treat = switch(
      estimand,
      "ATT" = function(Btreat, A, Xtreat, SW) {
        psi(Btreat, Xtreat, A, SW)
      },
      "ATE" = function(Btreat, A, Xtreat, SW) {
        p0 <- seq_len(length(Btreat) / 2)
        cbind(psi0(Btreat[p0], Xtreat, A, SW),
              psi1(Btreat[-p0], Xtreat, A, SW))
      }),
    wfun = switch(
      estimand,
      "ATE" = function(Btreat, Xtreat, A) {
        p0 <- seq_len(length(Btreat) / 2)

        ps <- numeric(length(A))
        ps[A == 0] <- .fam$linkinv(drop(Xtreat[A == 0,, drop = FALSE] %*% Btreat[p0]))
        ps[A == 1] <- .fam$linkinv(drop(Xtreat[A == 1,, drop = FALSE] %*% Btreat[-p0]))

        .get_w_from_ps_internal_bin(ps, A, estimand = estimand)
      },
      function(Btreat, Xtreat, A) {
        ps <- .fam$linkinv(drop(Xtreat %*% Btreat))
        .get_w_from_ps_internal_bin(ps, A, estimand = estimand)
      }),
    Xtreat = C,
    A = treat,
    btreat = unlist(par.list)
  )

  list(w = w, ps = ps, fit.obj = fit.list,
       Mparts = Mparts)
}

weightit2ipt.multi <- function(covs, treat, s.weights, subset, estimand, focal,
                               stabilize, missing, moments, int, verbose, ...) {
  A <- list(...)

  rlang::check_installed("rootSolve")

  covs <- covs[subset, , drop = FALSE]
  treat <- factor(treat[subset])
  s.weights <- s.weights[subset]

  if (missing == "ind") {
    covs <- add_missing_indicators(covs)
  }

  covs <- cbind(covs, .int_poly_f(covs, poly = moments, int = int, center = TRUE))

  covs <- cbind(covs, .quantile_f(covs, qu = A[["quantile"]], s.weights = s.weights,
                                  focal = focal, treat = treat))

  for (i in seq_col(covs)) covs[,i] <- .make_closer_to_1(covs[,i])

  colinear.covs.to.remove <- setdiff(colnames(covs), colnames(make_full_rank(covs)))
  covs <- covs[, colnames(covs) %nin% colinear.covs.to.remove, drop = FALSE]

  C <- cbind(`(Intercept)` = 1, covs)

  link <- if_null_then(A$link, "logit")
  chk::chk_string(link)
  link <- match_arg(link, c("logit", "probit", "cauchit", "cloglog"))

  .fam <- binomial(link)

  groups_to_weight <- switch(estimand,
                             "ATE" = levels(treat),
                             setdiff(levels(treat), focal))

  w <- rep.int(1, length(treat))

  fit.list <- make_list(groups_to_weight)
  par.list <- make_list(groups_to_weight)

  k <- ncol(C)

  start <- setNames(rep.int(0, k), colnames(C))

  f <- function(B, X, A, SW, .psi) {
    .colMeans(.psi(B, X, A, SW), length(A), k)
  }

  if (estimand == "ATE") {
    psi <- function(B, X, A, SW) {
      p <- .fam$linkinv(drop(X %*% B))
      SW * (A / p - 1) * X
    }

    for (i in groups_to_weight) {
      verbosely({
        fit.list[[i]] <- rootSolve::multiroot(f, start = start,
                                              A = as.numeric(treat == i),
                                              X = C,
                                              SW = s.weights,
                                              .psi = psi,
                                              rtol = 1e-10, atol = 1e-10, ctol = 1e-10,
                                              verbose = TRUE)
      }, verbose = verbose)

      par <- fit.list[[i]]$root

      ps_i <- .fam$linkinv(drop(C[treat == i,, drop = FALSE] %*% par))

      w[treat == i] <- 1 / ps_i

      par.list[[i]] <- par
    }
  }
  else {

    psi <- function(B, X, A, SW) {
      p <- .fam$linkinv(drop(X %*% B))
      SW * (A - (1 - A) * p / (1 - p)) * X
    }

    for (i in groups_to_weight) {

      verbosely({
        fit.list[[i]] <- rootSolve::multiroot(f, start = start,
                                              A = as.numeric(treat[treat %in% c(i, focal)] == focal),
                                              X = C[treat %in% c(i, focal),, drop = FALSE],
                                              SW = s.weights[treat %in% c(i, focal)],
                                              .psi = psi,
                                              rtol = 1e-10, atol = 1e-10, ctol = 1e-10,
                                              verbose = TRUE)
      }, verbose = verbose)
      par <- fit.list[[i]]$root

      ps_i <- .fam$linkinv(drop(C[treat == i,, drop = FALSE] %*% par))

      w[treat == i] <- ps_i / (1 - ps_i)

      par.list[[i]] <- par
    }
  }

  Mparts <- list(
    psi_treat = switch(
      estimand,
      "ATT" = function(Btreat, A, Xtreat, SW) {
        coef_ind <- setNames(lapply(seq_along(groups_to_weight), function(i) {
          (i - 1) * ncol(Xtreat) + seq_col(Xtreat)
        }), groups_to_weight)

        do.call("cbind", lapply(groups_to_weight, function(i) {
          m <- matrix(0, nrow = length(A), ncol = length(Btreat[coef_ind[[i]]]))

          m[A %in% c(i, focal),] <- psi(Btreat[coef_ind[[i]]],
                                        Xtreat[A %in% c(i, focal),, drop = FALSE],
                                        as.numeric(A[A %in% c(i, focal)] == focal),
                                        SW[A %in% c(i, focal)])
          m
        }))
      },
      "ATE" = function(Btreat, A, Xtreat, SW) {
        coef_ind <- setNames(lapply(seq_along(groups_to_weight), function(i) {
          (i - 1) * ncol(Xtreat) + seq_col(Xtreat)
        }), groups_to_weight)

        do.call("cbind", lapply(groups_to_weight, function(i) {
          psi(Btreat[coef_ind[[i]]], Xtreat, as.numeric(A == i), SW)
        }))
      }),
    wfun = switch(
      estimand,
      "ATT" = function(Btreat, Xtreat, A) {
        coef_ind <- setNames(lapply(seq_along(groups_to_weight), function(i) {
          (i - 1) * ncol(Xtreat) + seq_col(Xtreat)
        }), groups_to_weight)

        w <- rep.int(1, length(A))

        for (i in groups_to_weight) {
          ps_i <- .fam$linkinv(drop(Xtreat[A == i,, drop = FALSE] %*% Btreat[coef_ind[[i]]]))
          w[A == i] <- ps_i / (1 - ps_i)
        }
      },
      "ATE" = function(Btreat, Xtreat, A) {
        coef_ind <- setNames(lapply(seq_along(groups_to_weight), function(i) {
          (i - 1) * ncol(Xtreat) + seq_col(Xtreat)
        }), groups_to_weight)

        w <- rep.int(1, length(A))

        for (i in groups_to_weight) {
          ps_i <- .fam$linkinv(drop(Xtreat[A == i,, drop = FALSE] %*% Btreat[coef_ind[[i]]]))
          w[A == i] <- 1 / ps_i
        }
      }),
    Xtreat = C,
    A = treat,
    btreat = unlist(par.list)
  )

  list(w = w, fit.obj = fit.list,
       Mparts = Mparts)
}

weightit2ipt.cont <- function(covs, treat, s.weights, subset, missing, verbose, ...) {
  .err("`method = \"ipt\"` cannot be used with continuous treatments")
}
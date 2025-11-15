#' Inverse Probability Tilting
#' @name method_ipt
#' @usage NULL
#'
#' @description
#' This page explains the details of estimating weights using
#' inverse probability tilting by setting `method = "ipt"` in the call to
#' [weightit()] or [weightitMSM()]. This method can be used with binary and
#' multi-category treatments.
#'
#' In general, this method relies on estimating propensity scores using a
#' modification of the usual generalized linear model score equations to enforce
#' balance and then converting those propensity scores into weights using a
#' formula that depends on the desired estimand. This method relies on code
#' written for \pkg{WeightIt} using \pkgfun{rootSolve}{multiroot}.
#'
#' ## Binary Treatments
#'
#' For binary treatments, this method estimates the weights using formulas
#' described by Graham, Pinto, and Egel (2012). The following estimands are
#' allowed: ATE, ATT, and ATC. When the ATE is requested, the optimization is
#' run twice, once for each treatment group.
#'
#' ## Multi-Category Treatments
#'
#' For multi-category treatments, this method estimates the weights using
#' modifications of the formulas described by Graham, Pinto, and Egel (2012).
#' The following estimands are allowed: ATE and ATT. When the ATE is requested,
#' estimation is performed once for each treatment group. When the ATT is
#' requested, estimation is performed once for each non-focal (i.e., control)
#' group.
#'
#' ## Continuous Treatments
#'
#' Inverse probability tilting is not compatible with continuous treatments.
#'
#' ## Longitudinal Treatments
#'
#' For longitudinal treatments, the weights are the product of the weights
#' estimated at each time point. This method is not guaranteed to yield exact
#' balance at each time point. **NOTE: the use of inverse probability tilting with longitudinal treatments has not been validated!**
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
#' M-estimation is supported for all scenarios. See [glm_weightit()] and
#' `vignette("estimating-effects")` for details.
#'
#' @section Additional Arguments:
#'
#' \describe{
#'   \item{`moments`}{`integer`; the highest power of each covariate to be balanced. For example, if `moments = 3`, each covariate, its square, and its cube will be balanced. Can also be a named vector with a value for each covariate (e.g., `moments = c(x1 = 2, x2 = 4)`). Values greater than 1 for categorical covariates are ignored. Default is 1 to balance covariate means.
#'     }
#'     \item{`int`}{`logical`; whether first-order interactions of the covariates are to be balanced. Default is `FALSE`.
#'     }
#'     \item{`quantile`}{a named list of quantiles (values between 0 and 1) for each continuous covariate, which are used to create additional variables that when balanced ensure balance on the corresponding quantile of the variable. For example, setting `quantile = list(x1 = c(.25, .5. , .75))` ensures the 25th, 50th, and 75th percentiles of `x1` in each treatment group will be balanced in the weighted sample. Can also be a single number (e.g., `.5`) or a vector (e.g., `c(.25, .5, .75)`) to request the same quantile(s) for all continuous covariates.
#'     }
#'   \item{`link`}{the link used to determine the inverse link for computing the (generalized) propensity scores. Default is `"logit"`, which is used in the original description of the method by Graham, Pinto, and Egel (2012), but `"probit"`, `"cauchit"`, `"cloglog"`, `"loglog"`, `"log"`, and `"clog"` are also allowed. Note that negative weights are possible with these last two and they should be used with caution. An object of class `"link-glm"` can also be supplied. The argument is passed to [quasibinomial()].
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
#' Inverse probability tilting (IPT) involves specifying estimating
#' equations that fit the parameters of two or more generalized linear models
#' with a modification that ensures exact balance on the covariate means. These
#' estimating equations are solved, and the estimated parameters are used in the
#' (generalized) propensity score, which is used to compute the weights.
#' Conceptually and mathematically, IPT is very similar to entropy balancing and
#' just-identified CBPS. For the ATT and ATC, entropy balancing, just-identified
#' CBPS, and IPT will yield identical results. For the ATE or when `link` is
#' specified as something other than `"logit"`, the three methods differ.
#'
#' Treatment effect estimates for binary treatments are consistent if the true
#' propensity score is a logistic regression or the outcome model is linear in
#' the covariates and their interaction with treatments. For entropy balancing,
#' this is only true for the ATT, and for just-identified CBPS, this is only
#' true if there is no effect modification by covariates. In this way, IPT
#' provides additional theoretical guarantees over the other two methods, though
#' potentially with some cost in precision.
#'
#' @seealso [weightit()], [weightitMSM()]
#'
#' [method_ebal] and [method_cbps] for entropy balancing and CBPS, which work
#' similarly.
#'
#' @references
#' ## `estimand = "ATE"`
#'
#' Graham, B. S., De Xavier Pinto, C. C., & Egel, D. (2012). Inverse Probability
#' Tilting for Moment Condition Models with Missing Data. *The Review of Economic Studies*, 79(3), 1053–1079. \doi{10.1093/restud/rdr047}
#'
#' ## `estimand = "ATT"`
#'
#' Sant'Anna, P. H. C., & Zhao, J. (2020). Doubly robust difference-in-differences estimators. *Journal of Econometrics*, 219(1), 101–122. \doi{10.1016/j.jeconom.2020.06.003}
#'
#' @examplesIf rlang::is_installed("rootSolve")
#' data("lalonde", package = "cobalt")
#'
#' #Balancing covariates between treatment groups (binary)
#' (W1 <- weightit(treat ~ age + educ + married +
#'                   nodegree + re74, data = lalonde,
#'                 method = "ipt", estimand = "ATT"))
#'
#' summary(W1)
#'
#' cobalt::bal.tab(W1)
#'
#' #Balancing covariates with respect to race (multi-category)
#' (W2 <- weightit(race ~ age + educ + married +
#'                   nodegree + re74, data = lalonde,
#'                 method = "ipt", estimand = "ATE"))
#'
#' summary(W2)
#'
#' cobalt::bal.tab(W2)
NULL

weightit2ipt <- function(covs, treat, s.weights, subset, estimand, focal,
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
    .make_covs_closer_to_1() |>
    .make_covs_full_rank()

  C <- cbind(`(Intercept)` = 1, covs)

  t.lev <- get_treated_level(treat, estimand, focal)
  treat <- binarize(treat, one = t.lev)

  link <- ...get("link", "logit")

  if (chk::vld_string(link)) {
    chk::chk_subset(link, c("logit", "probit", "cloglog", "loglog",
                            "cauchit", "log", "clog"))

    link <- .make_link(link)
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

  groups_to_weight <- switch(estimand,
                             ATT = 0,
                             ATC = 1,
                             0:1)

  fit.list <- par.list <- make_list(groups_to_weight)

  n <- length(treat)
  k <- ncol(C)

  start <- .get_glm_starting_values(X = C, Y = treat, w = s.weights,
                                    family = .fam)

  f <- function(B, X, A, SW, .psi) {
    .colMeans(.psi(B, X, A, SW), n, k)
  }

  if (estimand == "ATE") {
    ps <- rep.int(0, n)

    psi.list <- list(
      "0" = function(B, X, A, SW) {
        p <- rep.int(0, n)
        p[A == 0] <- .fam$linkinv(drop(X[A == 0, , drop = FALSE] %*% B))
        SW * ((1 - A) / (1 - p) - 1) * X
      },
      "1" = function(B, X, A, SW) {
        p <- rep.int(1, n)
        p[A == 1] <- .fam$linkinv(drop(X[A == 1, , drop = FALSE] %*% B))
        SW * (A / p - 1) * X
      })

    for (i in groups_to_weight) {
      ii <- as.character(i)

      verbosely({
        fit.list[[ii]] <- rootSolve::multiroot(f,
                                               start = start,
                                               X = C,
                                               A = treat,
                                               SW = s.weights,
                                               .psi = psi.list[[ii]],
                                               rtol = 1e-10, atol = 1e-10, ctol = 1e-10,
                                               verbose = TRUE)
      }, verbose = verbose)

      par.list[[ii]] <- fit.list[[ii]]$root

      ps[treat == i] <- .fam$linkinv(drop(C[treat == i, , drop = FALSE] %*% par.list[[ii]]))
    }
  }
  else {
    psi <- switch(estimand,
                  ATT = function(B, X, A, SW) {
                    p <- .fam$linkinv(drop(X %*% B))
                    SW * (A - (1 - A) * p / (1 - p)) * X
                  },
                  ATC = function(B, X, A, SW) {
                    p <- .fam$linkinv(drop(X %*% B))
                    SW * (A * (1 - p) / p - (1 - A)) * X
                  })

    verbosely({
      fit.list[[1L]] <- rootSolve::multiroot(f,
                                             start = start,
                                             A = treat,
                                             X = C,
                                             SW = s.weights,
                                             .psi = psi,
                                             rtol = 1e-10, atol = 1e-10, ctol = 1e-10,
                                             verbose = TRUE)
    }, verbose = verbose)

    par.list[[1L]] <- fit.list[[1L]]$root

    ps <- .fam$linkinv(drop(C %*% par.list[[1L]]))
  }

  if (any(unlist(grab(fit.list, "estim.precis")) > 1e-5)) {
    .wrn("the optimization failed to converge; consider using fewer covariates or a different link function")
  }

  w <- .get_w_from_ps_internal_bin(ps, treat, estimand = estimand)

  Mparts <- list(
    psi_treat = switch(
      estimand,
      ATE = function(Btreat, Xtreat, A, SW) {
        p0 <- seq_len(length(Btreat) / 2)
        cbind(psi.list[["0"]](Btreat[p0], Xtreat, A, SW),
              psi.list[["1"]](Btreat[-p0], Xtreat, A, SW))
      },
      function(Btreat, Xtreat, A, SW) {
        psi(Btreat, Xtreat, A, SW)
      }),
    wfun = switch(
      estimand,
      ATE = function(Btreat, Xtreat, A) {
        p0 <- seq_len(length(Btreat) / 2)
        A0 <- A == 0

        ps <- numeric(length(A))
        ps[A0] <- .fam$linkinv(drop(Xtreat[A0, , drop = FALSE] %*% Btreat[p0]))
        ps[!A0] <- .fam$linkinv(drop(Xtreat[!A0, , drop = FALSE] %*% Btreat[-p0]))

        .get_w_from_ps_internal_bin(ps, A, estimand = estimand)
      },
      function(Btreat, Xtreat, A) {
        ps <- .fam$linkinv(drop(Xtreat %*% Btreat))
        .get_w_from_ps_internal_bin(ps, A, estimand = estimand)
      }),
    dw_dBtreat = switch(
      estimand,
      ATE = function(Btreat, Xtreat, A, SW) {
        p0 <- seq_len(length(Btreat) / 2)
        A0 <- A == 0

        XB <- numeric(length(A))

        XB[A0] <- drop(Xtreat[A0, , drop = FALSE] %*% Btreat[p0])
        XB[!A0] <- drop(Xtreat[!A0, , drop = FALSE] %*% Btreat[-p0])

        ps <- .fam$linkinv(XB)

        .dw_dp_bin(ps, A, estimand = estimand) * .fam$mu.eta(XB) * cbind((1 - A) * Xtreat, A * Xtreat)
      },
      function(Btreat, Xtreat, A, SW) {
        XB <- drop(Xtreat %*% Btreat)
        ps <- .fam$linkinv(XB)

        .dw_dp_bin(ps, A, estimand = estimand) * .fam$mu.eta(XB) * Xtreat
      }),
    hess_treat = switch(
      estimand,
      ATE = function(Btreat, Xtreat, A, SW) {
        p0 <- seq_len(length(Btreat) / 2)
        A0 <- A == 0

        XB <- numeric(length(A))

        XB[A0] <- drop(Xtreat[A0, , drop = FALSE] %*% Btreat[p0])
        XB[!A0] <- drop(Xtreat[!A0, , drop = FALSE] %*% Btreat[-p0])

        ps <- .fam$linkinv(XB)

        dw <- .dw_dp_bin(ps, A, estimand = estimand) * .fam$mu.eta(XB) * SW

        .block_diag(list(crossprod(Xtreat[A0, , drop = FALSE],
                                   dw[A0] * Xtreat[A0, , drop = FALSE]),
                         crossprod(Xtreat[!A0, , drop = FALSE],
                                   dw[!A0] * Xtreat[!A0, , drop = FALSE])))
      },
      function(Btreat, Xtreat, A, SW) {
        XB <- drop(Xtreat %*% Btreat)
        ps <- .fam$linkinv(XB)

        dw <- .dw_dp_bin(ps, A, estimand = estimand) * .fam$mu.eta(XB) * SW

        crossprod(Xtreat, dw * (2 * A - 1) * Xtreat)
      }),
    Xtreat = C,
    A = treat,
    btreat = unlist(par.list)
  )

  list(w = w, ps = ps, fit.obj = fit.list,
       Mparts = Mparts)
}

weightit2ipt.multi <- function(covs, treat, s.weights, subset, estimand, focal,
                               stabilize, missing, verbose, ...) {

  covs <- covs[subset, , drop = FALSE]
  treat <- factor(treat[subset])
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
    .make_covs_closer_to_1() |>
    .make_covs_full_rank()

  C <- cbind(`(Intercept)` = 1, covs)

  link <- ...get("link", "logit")

  if (chk::vld_string(link)) {
    chk::chk_subset(link, c("logit", "probit", "cloglog", "loglog", "cauchit", "log", "clog"))

    link <- .make_link(link)
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

  groups_to_weight <- switch(estimand,
                             ATE = levels(treat),
                             setdiff(levels(treat), focal))

  w <- rep_with(1, treat)

  fit.list <- par.list <- make_list(groups_to_weight)

  k <- ncol(C)

  f <- function(B, X, A, SW, .psi) {
    .colMeans(.psi(B, X, A, SW), length(A), k)
  }

  if (estimand == "ATE") {
    psi <- function(B, X, A, SW) {
      p <- .fam$linkinv(drop(X %*% B))
      SW * (A / p - 1) * X
    }

    for (i in groups_to_weight) {
      start <- .get_glm_starting_values(X = C,
                                        Y = as.numeric(treat == i),
                                        w = s.weights,
                                        family = .fam)
      verbosely({
        fit.list[[i]] <- rootSolve::multiroot(f,
                                              start = start,
                                              A = as.numeric(treat == i),
                                              X = C,
                                              SW = s.weights,
                                              .psi = psi,
                                              rtol = 1e-10, atol = 1e-10, ctol = 1e-10,
                                              verbose = TRUE)
      }, verbose = verbose)

      par_out <- fit.list[[i]]$root

      ps_i <- .fam$linkinv(drop(C[treat == i, , drop = FALSE] %*% par_out))

      w[treat == i] <- 1 / ps_i

      par.list[[i]] <- par_out
    }
  }
  else {
    psi <- function(B, X, A, SW) {
      p <- .fam$linkinv(drop(X %*% B))
      SW * (A - (1 - A) * p / (1 - p)) * X
    }

    for (i in groups_to_weight) {
      in_fi <- treat %in% c(i, focal)

      start <- .get_glm_starting_values(X = C[in_fi, , drop = FALSE],
                                        Y = as.numeric(treat[in_fi] == focal),
                                        w = s.weights[in_fi],
                                        family = .fam)

      verbosely({
        fit.list[[i]] <- rootSolve::multiroot(f,
                                              start = start,
                                              A = as.numeric(treat[in_fi] == focal),
                                              X = C[in_fi, , drop = FALSE],
                                              SW = s.weights[in_fi],
                                              .psi = psi,
                                              rtol = 1e-10, atol = 1e-10, ctol = 1e-10,
                                              verbose = TRUE)
      }, verbose = verbose)
      par_out <- fit.list[[i]]$root

      ps_i <- .fam$linkinv(drop(C[treat == i, , drop = FALSE] %*% par_out))

      w[treat == i] <- ps_i / (1 - ps_i)

      par.list[[i]] <- par_out
    }
  }

  if (any(unlist(grab(fit.list, "estim.precis")) > 1e-5)) {
    .wrn("the optimization failed to converge; consider using fewer covariates or a different link function")
  }

  Mparts <- list(
    psi_treat = switch(
      estimand,
      ATE = function(Btreat, Xtreat, A, SW) {
        Bmat <- matrix(Btreat, nrow = ncol(Xtreat),
                       dimnames = list(colnames(Xtreat), groups_to_weight))

        do.call("cbind", lapply(groups_to_weight, function(i) {
          psi(Bmat[, i], Xtreat, as.numeric(A == i), SW)
        }))
      },
      function(Btreat, Xtreat, A, SW) {
        Bmat <- matrix(Btreat, nrow = ncol(Xtreat),
                       dimnames = list(colnames(Xtreat), groups_to_weight))

        do.call("cbind", lapply(groups_to_weight, function(i) {
          m <- matrix(0, nrow = length(A), ncol = length(Bmat[, i]))

          m[A %in% c(i, focal), ] <- psi(Bmat[, i],
                                         Xtreat[A %in% c(i, focal), , drop = FALSE],
                                         as.numeric(A[A %in% c(i, focal)] == focal),
                                         SW[A %in% c(i, focal)])
          m
        }))
      }),
    wfun = switch(
      estimand,
      ATE = function(Btreat, Xtreat, A) {
        Bmat <- matrix(Btreat, nrow = ncol(Xtreat),
                       dimnames = list(colnames(Xtreat), groups_to_weight))

        w <- rep_with(1, A)

        for (i in groups_to_weight) {
          ps_i <- .fam$linkinv(drop(Xtreat[A == i, , drop = FALSE] %*% Bmat[, i]))
          w[A == i] <- 1 / ps_i
        }

        w
      },
      function(Btreat, Xtreat, A) {
        Bmat <- matrix(Btreat, nrow = ncol(Xtreat),
                       dimnames = list(colnames(Xtreat), groups_to_weight))

        w <- rep_with(1, A)

        for (i in groups_to_weight) {
          ps_i <- .fam$linkinv(drop(Xtreat[A == i, , drop = FALSE] %*% Bmat[, i]))
          w[A == i] <- ps_i / (1 - ps_i)
        }

        w
      }),
    dw_dBtreat = switch(
      estimand,
      ATE = function(Btreat, Xtreat, A, SW) {
        Bmat <- matrix(Btreat, nrow = ncol(Xtreat),
                       dimnames = list(colnames(Xtreat),
                                       levels(A)))

        XB <- Xtreat %*% Bmat
        ps <- .fam$linkinv(XB)

        dw <- .dw_dp_multi(ps, A, estimand = estimand, focal = focal) * .fam$mu.eta(XB)

        do.call("cbind", lapply(levels(A), function(i) {
          dw[, i] * (A == i) * Xtreat
        }))
      },
      function(Btreat, Xtreat, A, SW) {
        Bmat <- matrix(Btreat, nrow = ncol(Xtreat),
                       dimnames = list(colnames(Xtreat),
                                       groups_to_weight))

        XB <- Xtreat %*% Bmat

        ps <- .fam$linkinv(XB)

        do.call("cbind", lapply(setdiff(levels(A), focal), function(i) {
          dw2 <- rep_with(0, A)
          in_fi <- A %in% c(i, focal)
          dw2[in_fi] <- .dw_dp_bin(ps[in_fi, i], as.numeric(A[in_fi] == focal), "ATT")

          dw2 * .fam$mu.eta(XB[, i]) * (in_fi) * Xtreat
        }))
      }),
    hess_treat = switch(
      estimand,
      ATE = function(Btreat, Xtreat, A, SW) {
        Bmat <- matrix(Btreat, nrow = ncol(Xtreat),
                       dimnames = list(colnames(Xtreat),
                                       levels(A)))

        XB <- Xtreat %*% Bmat
        ps <- .fam$linkinv(XB)

        dw <- .dw_dp_multi(ps, A, estimand = estimand, focal = focal) * .fam$mu.eta(XB) * SW

        .block_diag(lapply(levels(A), function(i) {
          crossprod(Xtreat, dw[, i] * Xtreat)
        }))
      },
      function(Btreat, Xtreat, A, SW) {
        Bmat <- matrix(Btreat, nrow = ncol(Xtreat),
                       dimnames = list(colnames(Xtreat),
                                       groups_to_weight))

        XB <- Xtreat %*% Bmat

        ps <- .fam$linkinv(XB)

        .block_diag(lapply(setdiff(levels(A), focal), function(i) {
          in_fi <- A %in% c(i, focal)

          A_fi <- as.numeric(A[in_fi] == focal)
          X_fi <- Xtreat[in_fi, , drop = FALSE]

          dw_fi <- .dw_dp_bin(ps[in_fi, i], A_fi, "ATT")

          crossprod(X_fi, dw_fi * .fam$mu.eta(XB[in_fi, i]) * SW[in_fi] * (2 * A_fi - 1) * X_fi)
        }))
      }),
    Xtreat = C,
    A = treat,
    btreat = unlist(par.list)
  )

  list(w = w, fit.obj = fit.list,
       Mparts = Mparts)
}

weightit2ipt.cont <- function(covs, treat, s.weights, subset, missing, verbose, ...) {
  .err('`method = "ipt"` cannot be used with continuous treatments')
}

#Navigated weighting https://doi.org/10.48550/arXiv.2005.10998
weightit2nawt <- function(covs, treat, s.weights, subset, estimand, focal,
                         stabilize, missing, verbose, ...) {

  covs <- covs[subset, , drop = FALSE]
  treat <- treat[subset]
  s.weights <- s.weights[subset]

  missing <- .process_missing2(missing, covs)

  if (missing == "ind") {
    covs <- add_missing_indicators(covs)
  }

  covs <- .apply_moments_int_quantile(covs,
                                      moments = ...get("moments"),
                                      int = ...get("int"),
                                      quantile = ...get("quantile"),
                                      s.weights = s.weights, focal = focal,
                                      treat = treat)

  for (i in seq_col(covs)) {
    covs[, i] <- .make_closer_to_1(covs[, i])
  }

  colinear.covs.to.remove <- setdiff(colnames(covs), colnames(make_full_rank(covs)))
  covs <- covs[, colnames(covs) %nin% colinear.covs.to.remove, drop = FALSE]

  C <- cbind(`(Intercept)` = 1, covs)

  t.lev <- get_treated_level(treat, estimand, focal)
  treat <- binarize(treat, one = t.lev)

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

  alpha <- ...get("alpha", 2)
  chk::chk_number(alpha)
  # chk::chk_gte(alpha, 0)

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
        lin_pred <- drop(X %*% B)
        p <- .fam$linkinv(lin_pred)
        (SW * (A - p) * .fam$mu.eta(lin_pred) / .fam$variance(p)) * X * (1 - p)^alpha
      },
      "1" = function(B, X, A, SW) {
        lin_pred <- drop(X %*% B)
        p <- .fam$linkinv(lin_pred)
        (SW * (A - p) * .fam$mu.eta(lin_pred) / .fam$variance(p)) * X * p^alpha
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
                    lin_pred <- drop(X %*% B)
                    p <- .fam$linkinv(lin_pred)
                    (SW * (A - p) * .fam$mu.eta(lin_pred) / .fam$variance(p)) * X * p^alpha
                  },
                  ATC = function(B, X, A, SW) {
                    lin_pred <- drop(X %*% B)
                    p <- .fam$linkinv(lin_pred)
                    (SW * (A - p) * .fam$mu.eta(lin_pred) / .fam$variance(p)) * X * (1 - p)^alpha
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
    # hess_treat = switch(
    #   estimand,
    #   ATE = function(Btreat, Xtreat, A, SW) {
    #     p0 <- seq_len(length(Btreat) / 2)
    #     A0 <- A == 0
    #
    #     XB <- numeric(length(A))
    #
    #     XB[A0] <- drop(Xtreat[A0, , drop = FALSE] %*% Btreat[p0])
    #     XB[!A0] <- drop(Xtreat[!A0, , drop = FALSE] %*% Btreat[-p0])
    #
    #     ps <- .fam$linkinv(XB)
    #
    #     dw <- .dw_dp_bin(ps, A, estimand = estimand) * .fam$mu.eta(XB) * SW
    #
    #     .block_diag(list(crossprod(Xtreat[A0, , drop = FALSE],
    #                                dw[A0] * Xtreat[A0, , drop = FALSE]),
    #                      crossprod(Xtreat[!A0, , drop = FALSE],
    #                                dw[!A0] * Xtreat[!A0, , drop = FALSE])))
    #   },
    #   function(Btreat, Xtreat, A, SW) {
    #     XB <- drop(Xtreat %*% Btreat)
    #     ps <- .fam$linkinv(XB)
    #
    #     dw <- .dw_dp_bin(ps, A, estimand = estimand) * .fam$mu.eta(XB) * SW
    #
    #     crossprod(Xtreat, dw * (2 * A - 1) * Xtreat)
    #   }),
    Xtreat = C,
    A = treat,
    btreat = unlist(par.list)
  )

  list(w = w, ps = ps, fit.obj = fit.list,
       Mparts = Mparts)
}
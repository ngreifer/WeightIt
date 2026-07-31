#' Fitting (Weighted) Ordinal Regression Models
#'
#' @description
#' `ordinal_weightit()` fits an ordinal regression model with a
#' covariance matrix that accounts for estimation of weights, if supplied. By default, this function uses M-estimation to construct a robust covariance
#' matrix using the estimating equations for the weighting model and the outcome
#' model when available.
#'
#' @inheritParams glm_weightit
#' @param link a string corresponding to the desired link function. Allowable options include `"logit"`, `"probit"`, `"loglog"`, `"cloglog"`, and `"cauchit"`. Default is `"logit"` for ordinal logistic regression.
#' @param br `logical`; whether to use mean bias reduction, i.e., to solve the
#'   bias-reducing adjusted score equations of Kosmidis (2014) rather than the
#'   score equations. This yields estimates with smaller asymptotic bias that are
#'   always finite, even when the maximum likelihood estimates are not (e.g., when
#'   an end category is unobserved). Default is `FALSE`. See Details.
#'
#' @returns
#' An `ordinal_weightit` object.
#'
#' Unless `vcov = "none"`, the `vcov` component contains the covariance matrix
#' adjusted for the estimation of the weights if requested and a compatible
#' `weightit` object was supplied. The `vcov_type` component contains the type
#' of variance matrix requested. If `cluster` is supplied, it will be stored in
#' the `"cluster"` attribute of the output object, even if not used.
#'
#' The `model` component of the output object (also the `model.frame()` output)
#' will include two extra columns when `weightit` is supplied: `(weights)`
#' containing the weights used in the model (the product of the estimated
#' weights and the sampling weights, if any) and `(s.weights)` containing the
#' sampling weights, which will be 1 if `s.weights` is not supplied in the
#' original `weightit()` call.
#'
#' @details
#' `ordinal_weightit()` implements proportional odds ordinal regression using a
#' custom function in \pkg{WeightIt} that optionally computes a coefficient variance matrix that can be adjusted to
#' account for estimation of the weights if a `weightit` or `weightitMSM` object
#' is supplied to the `weightit` argument. Estimation of coefficients should align with that from
#' `MASS::polr()` unless `br = TRUE`.
#'
#' When no argument is supplied to
#' `weightit` or there is no `"Mparts"` attribute in the supplied object, the
#' default variance matrix returned will be the "HC0" sandwich variance matrix,
#' which is robust to misspecification of the outcome family (including
#' heteroscedasticity). Otherwise, the default variance matrix uses M-estimation
#' to additionally adjust for estimation of the weights. When possible, this
#' often yields smaller (and more accurate) standard errors. See the individual
#' methods pages to see whether and when an `"Mparts"` attribute is included in
#' the supplied object. To request that a variance matrix be computed that
#' doesn't account for estimation of the weights even when a compatible
#' `weightit` object is supplied, set `vcov = "HC0"`, which treats the weights
#' as fixed.
#'
#' Bootstrapping can also be used to compute the coefficient variance matrix;
#' when `vcov = "BS"` or `vcov = "FWB"`, which implement the traditional
#' resampling-based and fractional weighted bootstrap, respectively, the entire
#' process of estimating the weights and fitting the outcome model is repeated
#' in bootstrap samples (if a `weightit` object is supplied). This accounts for
#' estimation of the weights and can be used with any weighting method. It is
#' important to set a seed using `set.seed()` to ensure replicability of the
#' results. The fractional weighted bootstrap is more reliable but requires the
#' weighting method to accept sampling weights (which most do, and you'll get an
#' error if it doesn't). Setting `vcov = "FWB"` and supplying `fwb.args = list(wtype = "multinom")`
#' also performs the resampling-based bootstrap but
#' with the additional features \pkg{fwb} provides (e.g., a progress bar and
#' parallelization).
#'
#' ## Bias reduction
#'
#' When `br = TRUE`, the coefficients solve the bias-reducing adjusted score
#' equations for cumulative link models derived by Kosmidis (2014) instead of the
#' score equations, which removes the first-order term in the asymptotic bias of
#' the estimates. Reduction of bias is equivalent to a parameter-dependent additive
#' adjustment of the multinomial counts, and the resulting estimates are always
#' finite and remain invariant to linear transformations of the parameters and to
#' reversal of the order of the response categories. In the special case of two
#' response categories, this reduces to the bias-reduced generalized linear model
#' available with `br = TRUE` in [glm_weightit()]; for a saturated model without
#' covariates and a logit link, it amounts to adding \eqn{1/2} to the counts of the
#' first and last categories.
#'
#' The adjusted score equations are solved by quasi-Fisher scoring started at the
#' maximum likelihood estimates. Two components of `control` govern this: `br.maxit`
#' (the maximum number of iterations, default 100) and `br.tol` (the convergence
#' tolerance for the adjusted score relative to the sum of the weights, default
#' `1e-10`).
#'
#' Weights are treated as multinomial totals, which makes the estimates invariant
#' to whether the data are supplied as individual units or as groups of identical
#' units with weights equal to their counts. As for `br = TRUE` in
#' [glm_weightit()], the reported variance matrix uses the information matrix at
#' the estimates rather than the Jacobian of the adjusted score, i.e., the
#' adjustment is treated as fixed; the two differ by a term that vanishes
#' asymptotically. M-estimation and bootstrapping can be used with `br = TRUE` just
#' as they can without it.
#'
#' Note that \pkgfun{brglm2}{bracl}, which also performs bias reduction for ordinal
#' responses, fits an adjacent category logit model rather than a cumulative link
#' model, so its estimates are not comparable to those from `ordinal_weightit()`.
#' This is true of `parallel = TRUE` as well, which constrains the adjacent
#' category logits rather than the cumulative ones to be parallel and so yields
#' different fitted probabilities.
#'
#' @seealso
#' * [glm_weightit()] for fitting generalized linear models that adjust for estimation of the weights.
#' * [multinom_weightit()] for fitting multinomial regression models that adjust for estimation of the weights.
#' * [coxph_weightit()] for fitting Cox proportional hazards models that adjust for estimation of the weights.
#' * \pkgfun{MASS}{polr} for fitting ordinal regression models that do not account for estimation of the weights.
#'
#' @references
#' Firth, D. (1993). Bias reduction of maximum likelihood estimates.
#' *Biometrika*, 80(1), 27–38. \doi{10.1093/biomet/80.1.27}
#'
#' Kosmidis, I. (2014). Improved estimation in cumulative link models. *Journal of
#' the Royal Statistical Society: Series B (Statistical Methodology)*, 76(1),
#' 169–196. \doi{10.1111/rssb.12025}
#'
#' @examples
#' data("lalonde", package = "cobalt")
#'
#' # Logistic regression ATT weights
#' w.out <- weightit(treat ~ age + educ + married + re74,
#'                   data = lalonde, method = "glm",
#'                   estimand = "ATT")
#'
#' lalonde$re78_3o <- factor(findInterval(lalonde$re78,
#'                                       c(0, 5e3, 1e4)),
#'                          ordered = TRUE)
#'
#'
#' # Ordinal probit regression that adjusts for estimation
#' # of weights
#' fit <- ordinal_weightit(re78_3o ~ treat,
#'                          data = lalonde,
#'                          link = "probit",
#'                          weightit = w.out)
#'
#' summary(fit)
#'
#' # Same model using mean bias reduction
#' fit_br <- ordinal_weightit(re78_3o ~ treat,
#'                            data = lalonde,
#'                            link = "probit",
#'                            weightit = w.out,
#'                            br = TRUE)
#'
#' summary(fit_br)

#' @export
ordinal_weightit <- function(formula, data, link = "logit", weightit = NULL,
                             vcov = NULL, cluster, R = 500L,
                             offset, start = NULL,
                             control = list(...),
                             x = FALSE, y = TRUE,
                             contrasts = NULL, fwb.args = list(),
                             br = FALSE, ...) {

  vcov <- .process_vcov(vcov, weightit, R, fwb.args)

  if (missing(cluster)) {
    cluster <- NULL
  }

  model_call <- match.call()

  if (is_not_null(...get("weights"))) {
    arg::wrn("{.arg weights} is not an allowable argument to {.fun {rlang::call_name(model_call)}} and will be ignored. To fit a weighted model, supply a {.cls weightit} or {.cls weightitMSM} object to the {.arg weightit} argument")

    model_call[["weights"]] <- NULL
  }

  ###
  if (is_not_null(...get("family"))) {
    arg::err(c("{.arg family} cannot be used with {.fun ordinal_weightit}.",
               "i" = "Did you mean to use {.arg link} instead?"))
  }

  arg::arg_flag(br)

  internal_model_call <- .build_internal_model_call(model = "ordinal",
                                                    model_call = model_call,
                                                    weightit = weightit,
                                                    vcov = vcov,
                                                    br = br)

  fit <- .eval_fit(internal_model_call,
                   errors = c("missing values in object" = "missing values are not allowed in the model variables"),
                   from = FALSE)

  fit$family$family <- "ordinal"

  fit$br <- br
  ###

  # Class must be assigned before `.compute_vcov()`, which dispatches on it
  class(fit) <- "ordinal_weightit"

  fit$vcov <- .compute_vcov(fit, weightit, vcov, cluster, model_call, internal_model_call)

  fit <- .process_fit(fit, weightit, vcov, model_call, x, y)

  fit
}

# Expected information matrix of a cumulative link model, F = sum_r Z_r' M_r Z_r,
# where Z_r = [-x_r' | I_q] maps the coefficients (ordered as in WeightIt, with the
# thresholds last) to the q linear predictors and M_r is the information for those
# linear predictors. M_r / w_r is tridiagonal with
#   diagonal      g_s^2 (1 / P_s + 1 / P_{s + 1})
#   off-diagonal  -g_s g_{s + 1} / P_{s + 1}
# (Kosmidis, 2014, Section 4.3). `P` is the n x (q + 1) matrix of category
# probabilities and `g` the n x q matrix of first derivatives of the inverse link
# at the linear predictors.
.ordinal_info <- function(X, P, g, weights) {
  q <- ncol(g)
  p <- ncol(X)

  P0 <- P[, seq_len(q), drop = FALSE]  #P_s,       s = 1, ..., q
  P1 <- P[, -1L, drop = FALSE]         #P_{s + 1}, s = 1, ..., q

  d <- g^2 * (1 / P0 + 1 / P1)

  Faa <- sq_matrix(0, n = q)
  diag(Faa) <- colSums(d * weights)

  M1 <- d #row sums of M_r, i.e., M_r %*% 1

  if (q > 1L) {
    e <- -g[, -q, drop = FALSE] * g[, -1L, drop = FALSE] / P1[, -q, drop = FALSE]

    M1[, -q] <- M1[, -q] + e
    M1[, -1L] <- M1[, -1L] + e

    for (s in seq_len(q - 1L)) {
      Faa[s, s + 1L] <- Faa[s + 1L, s] <- sum(e[, s] * weights)
    }
  }

  M1 <- M1 * weights

  if (p == 0L) {
    return(Faa)
  }

  Fba <- -crossprod(X, M1)

  rbind(cbind(crossprod(X, X * rowSums(M1)), Fba),
        cbind(t(Fba), Faa))
}

# Bias-reducing adjustment to the score contributions of a cumulative link model
# (Kosmidis, 2014, eqns 7-9), returned as the n x (p + q) matrix to be added to the
# maximum likelihood score contributions. Reduction of bias amounts to adjusting the
# multinomial count for category s by c_rs - c_{r,s-1}, where
#     c_rs = m_r g'(eta_rs) v_rss / 2,  c_r0 = c_rk = 0
# and v_rss is the sth diagonal element of Z_r F^{-1} Z_r', the asymptotic
# covariance matrix of the linear predictors for unit r.
.ordinal_br_psi <- function(X, P, g, gp, weights) {
  q <- ncol(g)
  p <- ncol(X)

  Finv <- .solve_info(.ordinal_info(X, P, g, weights))

  a_ind <- p + seq_len(q)

  v <- matrix(diag(Finv)[a_ind], nrow(X), q, byrow = TRUE)

  if (p > 0L) {
    b_ind <- seq_len(p)

    v <- v + rowSums((X %*% Finv[b_ind, b_ind, drop = FALSE]) * X) -
      2 * X %*% Finv[b_ind, a_ind, drop = FALSE]
  }

  cc <- .5 * weights * gp * v

  #Adjustment to the counts, a_rs = c_rs - c_{r,s-1}, for s = 1, ..., k
  adj <- cbind(cc, 0) - cbind(0, cc)

  E <- g * (adj[, seq_len(q), drop = FALSE] / P[, seq_len(q), drop = FALSE] -
              adj[, -1L, drop = FALSE] / P[, -1L, drop = FALSE])

  if (p == 0L) {
    return(E)
  }

  cbind(-X * rowSums(E), E)
}

# Solves the bias-reducing adjusted score equations of a cumulative link model by
# quasi-Fisher scoring, with step halving to keep the thresholds ordered and the
# adjusted score shrinking. `gr()` supplies the adjusted score and `info()` the
# expected information at a set of coefficients; `n_a` is the number of thresholds.
.ordinal_br_solve <- function(start, X, y, weights, offset, gr, info, n_a,
                              maxit = 100L, tol = 1e-10) {
  a_ind <- length(start) - n_a + seq_len(n_a)

  #The score scales with the sum of the weights, so the tolerance should, too
  crit <- tol * max(1, sum(weights))

  theta <- start
  U <- gr(theta, X, y, weights, offset)

  i <- 0L

  while (max(abs(U)) > crit && i < maxit) {
    i <- i + 1L

    step <- try(drop(solve(info(theta), U)), silent = TRUE)

    if (null_or_error(step) || anyNA(step)) break

    U_new <- NULL

    for (h in 0L:25L) {
      theta_new <- theta + step / 2^h

      if (n_a > 1L && any(diff(theta_new[a_ind]) <= 0)) next

      U_try <- try(gr(theta_new, X, y, weights, offset), silent = TRUE)

      if (null_or_error(U_try) || anyNA(U_try)) next

      if (max(abs(U_try)) < max(abs(U))) {
        U_new <- U_try
        break
      }
    }

    if (is_null(U_new)) break

    theta <- theta_new
    U <- U_new
  }

  list(par = theta, converged = max(abs(U)) <= crit, iter = i, score = U)
}

# Ordinal regression
.ordinal_weightit.fit <- function(x, y, weights = NULL, start = NULL, offset = NULL,
                                  link = "logit", hess = TRUE, control = list(),
                                  br = FALSE, ...) {
  arg::arg_atomic(y)
  arg::arg_numeric(x)
  arg::arg_matrix(x)
  arg::arg_flag(br)

  if (rlang::is_string(link)) {
    arg::arg_element(link, c("logit", "probit", "loglog", "cloglog", "cauchit"))

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
    arg::err("{.arg link} must be a string or an object of class {.cls link-glm}")
  }

  fam <- binomial(link)

  if (!is.function(fam$linkinv)) {
    arg::err("the supplied link seems not to create a valid binomial family object")
  }

  .linkfun <- fam$linkfun
  .linkinv <- fam$linkinv
  .mu.eta <- fam$mu.eta
  .d2mu.deta <- {
    if (br) .make_d2mu.deta(link)
    else NULL
  }

  y <- droplevels(as.factor(y))
  n <- length(y)

  if (is_null(weights)) weights <- rep.int(1, n)
  else arg::arg_numeric(weights)

  if (is_null(offset)) offset <- rep.int(0, n)
  else arg::arg_numeric(offset)

  if (!all_the_same(c(length(y), nrow(x), length(weights), length(offset)))) {
    arg::err('{.arg {c("y", "x", "weights", "offset")}} must all have the same number of units')
  }

  x <- x[, colnames(x) != "(Intercept)", drop = FALSE]

  m <- nlevels(y) #num. thresholds
  k0 <- ncol(x) + m - 1L #num. params

  nm <- c(colnames(x), paste(levels(y)[-m], levels(y)[-1L], sep = "|"))

  aliased_X <- !colnames(x) %in% colnames(make_full_rank(x, with.intercept = TRUE))
  aliased_B <- c(aliased_X, rep.int(FALSE, m - 1L))

  x_ <- x[, !aliased_X, drop = FALSE]
  y_ <- as.integer(y)

  no_x <- is_null(x_)

  if (no_x) {
    start <- .linkfun(cumsum(tabulate(y_)[-m] / n))
  }
  else if (is_null(start)) {
    q1 <- floor(median(y_))
    y1 <- as.numeric(y_ > q1)
    coefs <- .get_glm_starting_values(X = cbind(1, x_), Y = y1, w = weights,
                                      family = fam, offset = offset)

    if (m > 2L) {
      spacing <- .linkfun(cumsum(tabulate(y_)[-m] / n))
      start <- c(coefs[-1L], -coefs[1L] + spacing - spacing[q1])
    }
    else {
      start <- c(coefs[-1L], -coefs[1L])
    }
  }
  else {
    arg::arg_numeric(start)
    arg::arg_length(start, k0)

    start <- start[!aliased_B]

    if (any(diff(start[-seq_col(x)]) <= 0)) {
      arg::err("starting values for the thresholds must be in ascending order")
    }
  }

  # Adjust start to use cumsum parameterization
  if (m > 2L) {
    if (no_x) {
      start[-1L] <- log(diff(start))
    }
    else {
      start[-seq_len(ncol(x_) + 1L)] <- log(diff(start[-seq_len(ncol(x_))]))
    }
  }

  names(start) <- nm[!aliased_B]

  ind_mat <- cbind(seq_along(y), y_)

  y_mat <- 1 * vapply(seq_len(m), function(j) y_ == j, logical(length(y_)))
  y_mat0 <- y_mat[, -m, drop = FALSE]
  y_mat1 <- y_mat[, -1L, drop = FALSE]

  #Get predictors on a smaller scale
  if (!no_x) {
    sds <- apply(x_, 2L, sd)
    x_ <- sweep(x_, 2L, sds, "/")
    start <- start * c(sds, rep.int(1, m - 1L))
  }

  # Get predicted probabilities for all units for category y
  get_p <- function(y, xb, a) {
    squish(.linkinv(c(a, Inf)[y] - xb) - .linkinv(c(-Inf, a)[y] - xb),
           lo = 1e-16)
  }

  # Category probabilities and the first two derivatives of the inverse link at all
  # `m - 1` linear predictors, used to compute the bias-reducing adjustment
  # (natural parameterization of `a`)
  get_br_parts <- function(B, X, offset) {
    if (no_x) {
      a <- B
      Xb <- offset
    }
    else {
      a <- B[-seq_col(X)]
      Xb <- offset + drop(X %*% B[seq_col(X)])
    }

    eta <- outer(-Xb, a, "+")
    d <- dim(eta)

    GG <- matrix(.linkinv(eta), d[1L], d[2L])

    list(P = squish(cbind(GG, 1) - cbind(0, GG), lo = 1e-16, hi = Inf),
         g = matrix(.mu.eta(eta), d[1L], d[2L]),
         gp = matrix(.d2mu.deta(eta), d[1L], d[2L]))
  }

  #Ordinal regression LL (cumsum paramaterization)
  ll <- function(B, X, y, weights, offset, .cumsum_param = TRUE) {
    if (no_x) {
      a <- B
      Xb <- offset
    }
    else {
      a <- B[-seq_col(X)]
      b <- B[seq_col(X)]
      Xb <- offset + drop(X %*% b)
    }

    if (.cumsum_param && m > 2L) {
      a <- c(a[1L], a[1L] + cumsum(exp(a[-1L])))
    }

    # Probability of observed outcome
    p <- get_p(y, Xb, a)

    sum(weights * log(p))
  }

  #Controls for solving the adjusted score equations when `br = TRUE`; removed from
  #`control` because `optim()` warns about control components it doesn't recognize
  br_maxit <- control$br.maxit %or% 100L
  br_tol <- control$br.tol %or% 1e-10

  control[c("br.maxit", "br.tol")] <- NULL

  m_control <- list(fnscale = -1, #maximize likelihood; optim() minimizes by default
                    trace = 0,
                    maxit = 1e3,
                    reltol = 1e-12)

  control <- utils::modifyList(m_control, control)

  # Estimate using cumsum parameterization to get estimates
  out0 <- optim(par = start,
                ll,
                X = x_,
                y = y_,
                weights = weights,
                offset = offset,
                .cumsum_param = TRUE,
                method = "BFGS",
                hessian = FALSE,
                control = control)

  theta0 <- out0$par

  # Convert to natural param
  if (m > 2L) {
    if (no_x) {
      theta0[-1L] <- theta0[1L] + cumsum(exp(theta0[-1L]))
    }
    else {
      a1 <- theta0[ncol(x_) + 1L]
      theta0[-seq_len(ncol(x_))] <- c(a1, a1 + cumsum(exp(theta0[-seq_len(ncol(x_) + 1L)])))
    }
  }

  # Psi function and gradient using natural parameterization. `.adjust` controls
  # whether the bias-reducing adjustment is included; it is excluded when computing
  # the Hessian, which should be that of the log-likelihood (i.e., minus the
  # information matrix) even for a bias-reduced fit.
  psi <- function(B, X, y, weights, offset = NULL, .adjust = br) {
    if (is_null(offset)) {
      offset <- rep_with(0, y)
    }

    if (no_x) {
      a <- B
      Xb <- offset
    }
    else {
      a <- B[-seq_col(X)]
      b <- B[seq_col(X)]
      Xb <- offset + drop(X %*% b)
    }

    pp <- get_p(y, Xb, a)

    gj <- .mu.eta(c(a, Inf)[y] - Xb)
    gj1 <- .mu.eta(c(-Inf, a)[y] - Xb)

    .psi_a <- gj * y_mat0 - gj1 * y_mat1

    if (no_x) {
      out <- as.matrix(.psi_a) * (weights / pp)
    }
    else {
      .psi_b <- X * (gj1 - gj)
      out <- cbind(.psi_b, .psi_a) * (weights / pp)
    }

    if (.adjust) {
      parts <- get_br_parts(B, X, offset)

      out <- out + .ordinal_br_psi(X, parts[["P"]], parts[["g"]], parts[["gp"]],
                                   weights)
    }

    colnames(out) <- names(B)
    out
  }

  gr <- function(B, X, y, weights, offset = NULL, .adjust = br) {
    colSums(psi(B, X, y, weights, offset, .adjust))
  }

  gr_ll <- function(B, X, y, weights, offset = NULL) {
    gr(B, X, y, weights, offset, .adjust = FALSE)
  }

  if (br) {
    # Solve the bias-reducing adjusted score equations by quasi-Fisher scoring,
    # starting from the maximum likelihood estimates. Unlike for the multinomial
    # logit model, the adjusted score of a cumulative link model is not the gradient
    # of a penalized likelihood, so there is no objective function to optimize.
    br_out <- .ordinal_br_solve(theta0, x_, y_, weights, offset,
                                gr = gr, info = function(B) {
                                  parts <- get_br_parts(B, x_, offset)

                                  .ordinal_info(x_, parts[["P"]], parts[["g"]],
                                                weights)
                                },
                                n_a = m - 1L,
                                maxit = br_maxit,
                                tol = br_tol)

    if (!br_out[["converged"]]) {
      arg::wrn("the bias-reducing adjusted score equations did not converge; estimates should not be trusted. Try increasing {.code br.maxit} in {.arg control} or simplifying the model")
    }

    theta0 <- br_out[["par"]]
  }

  hessian <- NULL
  if (hess) {
    # Estimate using natural parameterization to get hessian
    hessian <- try(optimHess(par = theta0,
                             function(...) ll(..., .cumsum_param = FALSE),
                             X = x_,
                             y = y_,
                             weights = weights,
                             offset = offset,
                             gr = gr_ll,
                             control = list(fnscale = -1)),
                   silent = TRUE)

    # If optimization fails, use numeric differentiation of gradient to get hessian
    if (null_or_error(hessian)) {
      hessian <- .gradient(gr_ll, theta0,
                           X = x_,
                           y = y_,
                           weights = weights,
                           offset = offset)
    }

    if (!no_x) {
      #Rescale hessian to be on original scale of predictors
      hessian <- hessian * tcrossprod(c(sds, rep.int(1, m - 1L)))
    }

    colnames(hessian) <- rownames(hessian) <- names(theta0)
  }

  theta <- theta0

  grad <- psi(theta, X = x_, y = y_,
              weights = weights, offset = offset)

  # Get predicted probabilities for all units for all categories,
  # natural parameterization of `a`
  get_pp <- function(B, X, offset = NULL) {
    if (is_null(offset)) {
      offset <- rep.int(0, n)
    }

    if (ncol(X) == 0L) {
      a <- B
      Xb <- offset
    }
    else {
      a <- B[-seq_col(X)]
      b <- B[seq_col(X)]
      Xb <- offset + drop(X %*% b)
    }

    GG <- vapply(a, function(a_) .linkinv(a_ - Xb),
                 numeric(length(Xb)))

    pp <- cbind(GG, 1) - cbind(0, GG)

    dimnames(pp) <- list(rownames(X), levels(y))

    pp
  }

  # Fitted values
  pp <- get_pp(theta, x_, offset)

  # Residuals
  res <- setNames(1 - pp[ind_mat], rownames(x))

  # Must be computed here, before theta is put back on the original scale
  # below -- x_ is never de-standardized, so theta and x_ need to be on the
  # same (standardized) scale when this is computed.
  linear.predictors <- offset + drop(x_ %*% theta[seq_col(x_)])

  # Adjust estimates and gradient to be put on original scale
  if (!no_x) {
    theta <- theta / c(sds, rep.int(1, m - 1L))
    grad <- sweep(grad, 2L, c(sds, rep.int(1, m - 1L)), "*")
  }

  coefs <- setNames(rep.int(NA_real_, ncol(x) + m - 1L), nm)
  coefs[!aliased_B] <- theta

  list(coefficients = coefs,
       residuals = res,
       fitted.values = pp,
       family = fam,
       linear.predictors = linear.predictors,
       solve = out0,
       psi = psi,
       f = gr,
       get_p = get_pp,
       df.residual = length(res) - ncol(x_) - (m - 1L),
       x = x,
       y = y,
       weights = weights,
       gradient = grad,
       hessian = hessian,
       br = br,
       varx = cov(x))
}

.ordinal_weightit <- function(formula, data, link = "logit", weights, subset, start = NULL, na.action,
                              hess = TRUE, control = list(), model = TRUE,
                              x = FALSE, y = TRUE, contrasts = NULL, br = FALSE, ...) {
  cal <- match.call()

  arg::arg_supplied(formula)
  arg::arg_string(link)
  arg::arg_flag(hess)
  arg::arg_flag(br)
  arg::arg_flag(model)
  arg::arg_flag(x)
  arg::arg_flag(y)

  if (missing(data)) {
    data <- environment(formula)
  }

  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"),
             names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())

  mt <- .attr(mf, "terms")

  Y <- model.response(mf, "any")
  if (length(dim(Y)) == 1L) {
    nm <- rownames(Y)
    dim(Y) <- NULL
    if (is_not_null(nm)) {
      names(Y) <- nm
    }
  }

  if (is.character(Y) || (is.factor(Y) && !is.ordered(Y))) {
    arg::wrn("the outcome variable is {.type {Y}} without a provided ordering, so the resulting estimates may be nonsensical. Code it as an {.help [ordered factor](base::ordered)} with a verified ordering to prevent invalid results")
  }

  X <- {
    if (is.empty.model(mt)) matrix(NA_real_, nrow = NROW(Y), ncol = 0L)
    else model.matrix(mt, mf, contrasts)
  }

  weights <- as.vector(model.weights(mf))

  if (is_not_null(weights)) {
    arg::arg_numeric(weights)
    arg::arg_gte(weights, 0)
  }

  offset <- as.vector(model.offset(mf))
  if (is_not_null(offset)) {
    arg::arg_numeric(offset)

    if (length(offset) != NROW(Y)) {
      arg::err("number of offsets is {length(offset)}; should equal {NROW(Y)} (number of observations)")
    }
  }

  fit <- eval(call(".ordinal_weightit.fit",
                   x = X, y = Y, link = link, weights = weights,
                   offset = offset, start = start,
                   hess = hess, control = control, br = br))

  if (model) fit$model <- mf
  fit$na.action <- .attr(mf, "na.action")
  if (!x) fit$x <- NULL
  if (!y) fit$y <- NULL

  c(fit,
    list(call = cal, formula = formula, terms = mt,
         data = data, offset = offset,
         contrasts = .attr(X, "contrasts"),
         xlevels = .getXlevels(mt, mf)))
}

.get_hess_ordinal <- function(fit) {
  x <- fit[["x"]] %or% model.matrix(fit)
  y <- fit[["y"]] %or% model.response(model.frame(fit))
  fam <- fit[["family"]]
  weights <- weights(fit)
  offset <- fit$offset
  coefs <- coef(fit)

  .linkinv <- fam$linkinv
  .mu.eta <- fam$mu.eta

  y <- droplevels(as.factor(y))
  n <- length(y)

  if (is_null(weights)) weights <- rep.int(1, n)

  if (is_null(offset)) offset <- rep.int(0, n)

  m <- nlevels(y) #num. thresholds

  aliased_X <- !colnames(x) %in% colnames(make_full_rank(x, with.intercept = TRUE))

  x_ <- x[, !aliased_X, drop = FALSE]
  y_ <- as.integer(y)

  no_x <- is_null(x_)

  y_mat <- 1 * vapply(seq_len(m), function(j) y_ == j, logical(length(y_)))
  y_mat0 <- y_mat[, -m, drop = FALSE]
  y_mat1 <- y_mat[, -1L, drop = FALSE]

  #Get predictors on a smaller scale
  if (!no_x) {
    sds <- apply(x_, 2L, sd)
    x_ <- sweep(x_, 2L, sds, "/")
  }

  # Get predicted probabilities for all units for category y
  get_p <- function(y, xb, a) {
    squish(.linkinv(c(a, Inf)[y] - xb) - .linkinv(c(-Inf, a)[y] - xb),
           lo = 1e-16)
  }

  #Ordinal regression LL (cumsum paramaterization)
  ll <- function(B, X, y, weights, offset, .cumsum_param = TRUE) {
    if (no_x) {
      a <- B
      Xb <- offset
    }
    else {
      a <- B[-seq_col(X)]
      b <- B[seq_col(X)]
      Xb <- offset + drop(X %*% b)
    }

    if (.cumsum_param && m > 2L) {
      a <- c(a[1L], a[1L] + cumsum(exp(a[-1L])))
    }

    # Probability of observed outcome
    p <- get_p(y, Xb, a)

    sum(weights * log(p))
  }

  theta0 <- na.rem(coefs) * c(sds, rep.int(1, m - 1L))

  # Psi function and gradient using natural parameterization
  psi <- function(B, X, y, weights, offset = NULL) {
    if (is_null(offset)) {
      offset <- rep_with(0, y)
    }

    if (no_x) {
      a <- B
      Xb <- offset
    }
    else {
      a <- B[-seq_col(X)]
      b <- B[seq_col(X)]
      Xb <- offset + drop(X %*% b)
    }

    pp <- get_p(y, Xb, a)

    gj <- .mu.eta(c(a, Inf)[y] - Xb)
    gj1 <- .mu.eta(c(-Inf, a)[y] - Xb)

    .psi_a <- gj * y_mat0 - gj1 * y_mat1

    if (no_x) {
      out <- as.matrix(.psi_a) * (weights / pp)
    }
    else {
      .psi_b <- X * (gj1 - gj)
      out <- cbind(.psi_b, .psi_a) * (weights / pp)
    }

    colnames(out) <- names(B)
    out
  }

  gr <- function(B, X, y, weights, offset = NULL) {
    colSums(psi(B, X, y, weights, offset))
  }

  # Estimate using natural parameterization to get hessian
  hessian <- try(optimHess(par = theta0,
                           function(...) ll(..., .cumsum_param = FALSE),
                           X = x_,
                           y = y_,
                           weights = weights,
                           offset = offset,
                           gr = gr,
                           control = list(fnscale = -1)),
                 silent = TRUE)

  # If optimization fails, use numeric differentiation of gradient to get hessian
  if (null_or_error(hessian)) {
    hessian <- .gradient(gr, theta0,
                         X = x_,
                         y = y_,
                         weights = weights,
                         offset = offset)
  }

  if (!no_x) {
    hessian <- hessian * tcrossprod(c(sds, rep.int(1, m - 1L)))
  }

  colnames(hessian) <- rownames(hessian) <- names(theta0)

  hessian
}

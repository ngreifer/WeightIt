#Families for continuous treatments with `method = "glm"` other than the default
#Gaussian one: a Gaussian location-scale model ("gaussian_ls"), Poisson,
#negative binomial ("negbin"), and Gamma. See `weightit2glm.cont()`, which handles
#the Gaussian family itself and calls `.weightit2glm.cont_family()` for the others.

.cont_families <- c("gaussian", "gaussian_ls", "poisson", "negbin", "Gamma")

.cont_family_links <- function(family) {
  switch(family,
         gaussian = ,
         gaussian_ls = c("identity", "log", "inverse"),
         poisson = ,
         negbin = c("log", "identity", "sqrt"),
         Gamma = c("log", "inverse", "identity"))
}

#Processes `family` and `link` for a continuous treatment. `family` is the name of
#one of `.cont_families` (case-insensitive) or a `gaussian()`, `poisson()`, or
#`Gamma()` family object, whose link is used when `link` is not supplied. Bias
#reduction (a `br.` prefix) is not available for continuous treatments.
.process_cont_family <- function(family = NULL, link = NULL) {
  if (is_null(family)) {
    family <- "gaussian"
  }

  if (inherits(family, "family")) {
    if (!family[["family"]] %in% c("gaussian", "poisson", "Gamma")) {
      arg::err(c("{.arg family} must be one of {.or {.val {(.cont_families)}}} for continuous treatments with {.code method = \"glm\"}.",
                 "i" = "A family object can be {.fun gaussian}, {.fun poisson}, or {.fun Gamma}."))
    }

    link <- link %or% family[["link"]]
    family <- family[["family"]]
  }
  else {
    if (!rlang::is_string(family)) {
      arg::err("{.arg family} must be a string or a family object")
    }

    m <- match(tolower(family), tolower(.cont_families))

    if (anyNA(m)) {
      arg::err("{.arg family} must be one of {.or {.val {(.cont_families)}}} for continuous treatments with {.code method = \"glm\"}")
    }

    family <- .cont_families[m]
  }

  links <- .cont_family_links(family)

  if (is_null(link)) {
    link <- links[1L]
  }

  if (!rlang::is_string(link)) {
    arg::err("{.arg link} must be a string")
  }

  link_match <- pmatch(link, links)

  if (anyNA(link_match)) {
    arg::err('with {.code family = "{family}"}, {.arg link} must be one of {.or {.val {links}}}')
  }

  link <- links[link_match]

  list(name = family,
       link = link,
       linkobj = make.link(link))
}

#Checks that the treatment lies in the support of the family
.check_cont_family_support <- function(treat, family) {
  if (family %in% c("poisson", "negbin")) {
    if (any(treat < 0) || any(abs(treat - round(treat)) > sqrt(.Machine$double.eps))) {
      arg::err(c('With {.code family = "{family}"}, the treatment must take only nonnegative integer values.',
                 "i" = 'For other nonnegative treatments, try {.code family = "Gamma"} (for positive values) or {.code family = "gaussian_ls"}.'))
    }
  }
  else if (family == "Gamma") {
    if (any(treat <= 0)) {
      arg::err(c('With {.code family = "Gamma"}, the treatment must take only positive values.',
                 "i" = 'For count treatments, try {.code family = "poisson"} or {.code family = "negbin"}; otherwise, try {.code family = "gaussian_ls"}.'))
    }
  }
}

#Variance function (up to the dispersion) for the families with one. `disp` is
#theta for the negative binomial and is unused otherwise.
.cont_family_variance <- function(family) {
  switch(family,
         poisson = function(mu, disp) mu,
         negbin = function(mu, disp) mu + mu^2 / disp,
         Gamma = function(mu, disp) mu^2)
}

#Conditional density of the treatment. `disp` is theta for the negative binomial,
#the shape for Gamma, and the SD for the Gaussian families.
.cont_family_density <- function(family) {
  switch(family,
         gaussian = ,
         gaussian_ls = function(a, mu, disp) dnorm(a, mu, disp),
         poisson = function(a, mu, disp) dpois(a, mu),
         negbin = function(a, mu, disp) dnbinom(a, size = disp, mu = mu),
         Gamma = function(a, mu, disp) dgamma(a, shape = disp, rate = disp / mu))
}

#ML score for the log of the Gamma shape `k` given the means, per unit
.gamma_shape_psi <- function(k, A, mu, SW) {
  SW * (log(k) + 1 - digamma(k) + log(A / mu) - A / mu)
}

#ML score for the log of the negative binomial theta given the means, per unit
.negbin_theta_psi <- function(theta, A, mu, SW) {
  SW * (digamma(A + theta) - digamma(theta) + log(theta) + 1 -
          log(theta + mu) - (A + theta) / (theta + mu))
}

#ML estimate of the Gamma shape given the means. Its score is strictly decreasing
#in the shape (its derivative is 1/k - trigamma(k) < 0), so the root is unique.
.gamma_shape_ml <- function(A, mu, SW) {
  f <- function(logk) sum(.gamma_shape_psi(exp(logk), A, mu, SW))

  exp(uniroot(f, c(-5, 5), extendInt = "downX", tol = 1e-12)$root)
}

#ML estimate of the negative binomial theta given the means
.negbin_theta_ml <- function(A, mu, SW, limit = 50L) {
  rlang::check_installed("MASS")

  theta <- rlang::try_fetch({
    MASS::theta.ml(A, mu, weights = SW, limit = limit)
  },
  warning = function(w) {
    arg::wrn("(from {.fun MASS::theta.ml}): {conditionMessage(w)}")
    invokeRestart("muffleWarning")
  })

  if (!is.finite(theta) || theta <= 0) {
    arg::err(c("The negative binomial dispersion parameter could not be estimated.",
               "i" = 'The treatment may not be overdispersed relative to a Poisson distribution; try {.code family = "poisson"}.'))
  }

  as.numeric(theta)
}

#Score for the coefficients of the mean model of a Poisson, negative binomial, or
#Gamma model, per unit. The dispersion (for the negative binomial, theta) enters
#only through the variance function; for Gamma, the shape cancels.
.cont_family_psi_coef <- function(B, X, A, SW, fam, disp = NULL) {
  eta <- drop(X %*% B)
  mu <- fam$linkobj$linkinv(eta)
  d <- fam$linkobj$mu.eta(eta)
  V <- .cont_family_variance(fam$name)(mu, disp)

  X * (SW * d * (A - mu) / V)
}

#Fit of a Gaussian location-scale model (a double GLM): the mean is linear on the
#`link` scale, and the log of the variance is linear in the same covariates. The
#coefficients solve the normal score equations, by alternating a weighted
#`glm.fit()` for the mean (with weights `SW / s2`) with Fisher scoring for the
#variance coefficients (whose score is `X' SW ((A - mu)^2 / s2 - 1) / 2`).
.fit_gaussian_ls <- function(X, A, SW, fam, control = list(), maxit = 100L,
                             tol = 1e-10) {
  family <- gaussian(link = fam$link)

  fit_mean <- function(w, start = NULL) {
    fit <- suppressWarnings({
      glm.fit(x = X, y = A, weights = w, family = family, start = start,
              control = do.call(glm.control, as.list(control)))
    })

    if (!isTRUE(fit$converged)) {
      arg::err('the mean model of the {.code family = "gaussian_ls"} treatment model failed to converge; try another {.arg link}')
    }

    fit
  }

  fit <- fit_mean(SW)
  B <- fit$coefficients

  r2 <- (A - fit$fitted.values)^2
  G <- c(log(w.m(r2, SW)), rep.int(0, ncol(X) - 1L))

  #Normal log-likelihood of the variance coefficients given the squared residuals
  ll_var <- function(G) {
    xg <- drop(X %*% G)
    -sum(SW * (xg + r2 * exp(-xg))) / 2
  }

  converged <- FALSE

  for (iter in seq_len(maxit)) {
    #Variance coefficients given the means: Fisher scoring to convergence, halving
    #the step until the log-likelihood does not decrease (a full step can overshoot
    #when the squared residuals are heavy-tailed)
    for (j in seq_len(100L)) {
      s2 <- exp(drop(X %*% G))
      step <- drop(solve(crossprod(sqrt(SW) * X), colSums(X * (SW * (r2 / s2 - 1)))))

      ll0 <- ll_var(G)
      a <- 1

      repeat {
        G_new <- G + a * step
        ll1 <- ll_var(G_new)

        if (is.finite(ll1) && ll1 >= ll0 - 1e-10 * abs(ll0)) {
          break
        }

        a <- a / 2

        if (a < 1e-10) {
          G_new <- G
          break
        }
      }

      G <- G_new

      if (max(abs(step)) < tol || a < 1e-10) break
    }

    s2 <- exp(drop(X %*% G))

    fit <- fit_mean(SW / s2, start = B)

    change <- max(abs(fit$coefficients - B))

    B <- fit$coefficients
    r2 <- (A - fit$fitted.values)^2

    if (change < tol) {
      converged <- TRUE
      break
    }
  }

  if (!converged) {
    arg::wrn('the {.code family = "gaussian_ls"} treatment model did not converge in {maxit} iterations')
  }

  s2 <- exp(drop(X %*% G))

  names(G) <- names(B) <- colnames(X)

  list(coefficients = B,
       var.coefficients = G,
       fitted.values = fit$fitted.values,
       sd = sqrt(s2),
       mean.fit = fit,
       iter = iter,
       converged = converged)
}

#Derivative of the weights of a Gaussian location-scale model with a normal
#density with respect to (G, B), the coefficients of the log variance and of the
#mean, as an n x (length(G) + length(B)) matrix. With
#  w_i = N_i / D_i,  N_i = sum_j c_j phi(z_ij) / s_j,  D_i = phi(z_ii) / s_i,
#  z_ij = (A_i - mu_j) / s_j,  s_j = exp(x_j'G / 2),  c_j = SW_j / sum(SW),
#dw/dtheta = w * (dN/dtheta / N - dlog(D)/dtheta), where
#  dlog(D_i)/dB = z_ii mu'_i x_i / s_i,       dlog(D_i)/dG = (z_ii^2 - 1) x_i / 2,
#  dN_i/dB = sum_j c_j phi(z_ij) z_ij mu'_j x_j / s_j^2,
#  dN_i/dG = sum_j c_j phi(z_ij) (z_ij^2 - 1) x_j / (2 s_j).
#The sums over units are computed at each distinct treatment value in chunks, as
#in `.get_w_from_gps_internal_cont_ls()`; with a numeric derivative, the O(n^2)
#numerator would instead be recomputed for every parameter and step.
.dw_dB_gaussian_ls <- function(G, B, X, A, s.weights, linkobj, chunk = 2e6) {
  eta <- drop(X %*% B)
  mu <- linkobj$linkinv(eta)
  dmu <- linkobj$mu.eta(eta)
  s <- exp(drop(X %*% G) / 2)
  cw <- s.weights / sum(s.weights)

  ua <- unique(A)
  rows <- max(1L, floor(chunk / length(mu)))

  N <- numeric(length(ua))
  dN_dB <- dN_dG <- matrix(0, nrow = length(ua), ncol = ncol(X))

  for (start in seq(1L, length(ua), by = rows)) {
    i <- start:min(length(ua), start + rows - 1L)
    s_rep <- rep(s, each = length(i))

    Z <- outer(ua[i], mu, "-") / s_rep
    P <- matrix(dnorm(Z) / s_rep, nrow = length(i))

    N[i] <- drop(P %*% cw)
    dN_dB[i, ] <- (P * Z) %*% (X * (cw * dmu / s))
    dN_dG[i, ] <- (P * (Z^2 - 1)) %*% (X * (cw / 2))
  }

  m <- match(A, ua)
  z <- (A - mu) / s
  w <- N[m] / (dnorm(z) / s)

  w * cbind(dN_dG[m, , drop = FALSE] / N[m] - X * ((z^2 - 1) / 2),
            dN_dB[m, , drop = FALSE] / N[m] - X * (z * dmu / s))
}

#Weights `f(A) / f(A | X)` for the Poisson, negative binomial, and Gamma families,
#where `f(A)` is the sampling-weighted average over units of the conditional
#density. For the discrete families, it is computed exactly at each distinct
#treatment value. For Gamma, it is computed on a grid evenly spaced on the log
#scale of the treatment and interpolated there, with the log density (a spline
#through a sharply peaked density can dip below 0). A grid evenly spaced on the
#treatment's own scale puts few points where a right-skewed treatment has most of
#its values and fails near 0 when the shape is below 1; on 2000 units, the log
#grid of 200 points gives weights within 5e-4 (relative) of the exact average for
#shapes from 0.3 to 200.
.get_w_from_cont_family <- function(mu, A, disp, family, s.weights, n_grid = 200L) {
  dfun <- .cont_family_density(family)

  dens_denom <- dfun(A, mu, disp)

  if (family %in% c("poisson", "negbin")) {
    ua <- sort(unique(A))

    dens_num <- vapply(ua, function(a) w.m(dfun(a, mu, disp), s.weights),
                       numeric(1L))[match(A, ua)]
  }
  else {
    log_grid <- seq(log(min(A)), log(max(A)), length.out = n_grid)

    dens_grid <- vapply(exp(log_grid), function(a) w.m(dfun(a, mu, disp), s.weights),
                        numeric(1L))

    dens_num <- spline(x = log_grid, y = log(dens_grid), xout = log(A))[["y"]] |>
      exp()
  }

  if (any(!is.finite(dens_denom)) || any(dens_denom <= 0) ||
      any(!is.finite(dens_num)) || any(dens_num <= 0)) {
    arg::err("the estimated density of the treatment is not positive and finite for all units. Try another {.arg family} or {.arg link}")
  }

  dens_num / dens_denom
}

#Fits the treatment model for a continuous treatment with a family other than the
#Gaussian one and returns the weights, the fit, and the M-estimation components.
#`covs` has already been subset and processed by `weightit2glm.cont()`, and
#`s.weights` normalized to a mean of 1.
.weightit2glm.cont_family <- function(covs, treat, s.weights, fam, make_dens_fun,
                                      re.bars, .data, subset, verbose, ...) {

  .check_cont_family_support(treat, fam$name)

  is_kernel <- isTRUE(attr(make_dens_fun, "is_kernel"))

  if (is_not_null(re.bars)) {
    #A location-scale model would need a model for the variance with random
    #effects, which lme4 does not fit
    if (fam$name == "gaussian_ls") {
      arg::err('random effects are not available with {.code family = "gaussian_ls"}')
    }

    if (fam$name == "negbin" && fam$link != "log") {
      arg::err('with random effects and {.code family = "negbin"}, {.arg link} must be {.val log}')
    }

    rlang::check_installed("lme4")

    df <- .make_re_data_formula(covs, treat, s.weights, re.bars, .data, subset)

    ctrl_fun <- lme4::glmerControl

    control <- do.call(ctrl_fun,
                       c(as.list(...get("control")),
                         ...mget(setdiff(rlang::fn_fmls_names(ctrl_fun)[pmatch(...names(), rlang::fn_fmls_names(ctrl_fun), 0L)],
                                         names(...get("control"))))))

    fit_fun_name <- if (fam$name == "negbin") "glmer.nb" else "glmer"
    fit_fun <- if (fam$name == "negbin") lme4::glmer.nb else lme4::glmer

    fit_args <- list(df[["formula"]], data = df[["data"]],
                     weights = quote(.s.weights),
                     control = control)

    if (fam$name != "negbin") {
      fit_args[["family"]] <- {
        if (fam$name == "poisson") poisson(link = fam$link)
        else Gamma(link = fam$link)
      }
    }

    rlang::try_fetch({verbosely({
      fit <- do.call(fit_fun, fit_args)
    }, verbose = verbose)},
    warning = function(w) {
      arg::wrn("(from {.fun lme4::{fit_fun_name}}): {conditionMessage(w)}")
      invokeRestart("muffleWarning")
    },
    error = function(e) {
      arg::err("(from {.fun lme4::{fit_fun_name}}): {conditionMessage(e)}")
    })

    mu <- as.numeric(predict(fit, type = "response"))

    #For Gamma, the shape is estimated by maximum likelihood given the fitted
    #conditional means, as without random effects; `1 / sigma(fit)^2` is a
    #moment-type estimate that underestimates small shapes
    disp <- switch(fam$name,
                   negbin = lme4::getME(fit, "glmer.nb.theta"),
                   Gamma = .gamma_shape_ml(treat, mu, s.weights),
                   NULL)

    w <- .get_w_from_cont_family(mu, treat, disp, fam$name, s.weights)

    return(list(w = w, fit.obj = fit, Mparts = NULL))
  }

  Xtreat <- cbind(`(Intercept)` = 1, covs)

  if (fam$name == "gaussian_ls") {
    verbosely({
      fit <- .fit_gaussian_ls(Xtreat, treat, s.weights, fam,
                              control = ...get("control", list()))
    }, verbose = verbose)

    w <- .get_w_from_gps_internal_cont(fit$fitted.values, treat, fit$sd, s.weights,
                                       make_dens_fun)

    Mparts <- NULL

    if (!is_kernel) {
      linkobj <- fam$linkobj
      k <- ncol(Xtreat)

      Mparts <- list(
        psi_treat = function(Btreat, Xtreat, A, SW) {
          G <- Btreat[seq_len(k)]
          B <- Btreat[-seq_len(k)]

          s2 <- exp(drop(Xtreat %*% G))
          eta <- drop(Xtreat %*% B)
          mu <- linkobj$linkinv(eta)

          SW <- SW / mean_fast(SW)

          cbind(Xtreat * (SW * ((A - mu)^2 / s2 - 1) / 2), #log conditional variance
                Xtreat * (SW * linkobj$mu.eta(eta) * (A - mu) / s2)) #conditional mean
        },
        wfun = function(Btreat, Xtreat, A) {
          G <- Btreat[seq_len(k)]
          B <- Btreat[-seq_len(k)]

          sd <- sqrt(exp(drop(Xtreat %*% G)))
          mu <- linkobj$linkinv(drop(Xtreat %*% B))

          .get_w_from_gps_internal_cont(mu, A, sd, s.weights,
                                        make_dens_fun)
        },
        dw_dBtreat = {
          #Closed form for the normal density with the exact numerator; otherwise
          #the numeric derivative is used
          if (isTRUE(attr(make_dens_fun, "is_dnorm")) && length(unique(treat)) <= 2e4) {
            function(Btreat, Xtreat, A, SW) {
              .dw_dB_gaussian_ls(G = Btreat[seq_len(k)], B = Btreat[-seq_len(k)],
                                 X = Xtreat, A = A, s.weights = s.weights,
                                 linkobj = linkobj)
            }
          }
          else NULL
        },
        Xtreat = Xtreat,
        A = treat,
        btreat = c(setNames(fit$var.coefficients,
                            paste0("log(s_r^2):", names(fit$var.coefficients))),
                   fit$coefficients)
      )
    }

    return(list(w = w, fit.obj = fit, Mparts = Mparts))
  }

  #Poisson, negative binomial, and Gamma
  if (fam$name == "negbin") {
    rlang::check_installed("MASS")

    data <- data.frame(treat, covs)
    formula <- if (ncol(covs) > 0L) formula(data) else treat ~ 1

    control <- do.call(glm.control, as.list(...get("control", list())))

    #`glm.nb()` takes `link` unevaluated, so it is inserted into the call as a name
    nb_call <- as.call(list(quote(MASS::glm.nb), formula,
                            data = quote(data),
                            weights = quote(s.weights),
                            control = quote(control),
                            link = as.name(fam$link)))

    rlang::try_fetch({verbosely({
      fit <- eval(nb_call)
    }, verbose = verbose)},
    warning = function(w) {
      arg::wrn("(from {.fun MASS::glm.nb}): {conditionMessage(w)}")
      invokeRestart("muffleWarning")
    },
    error = function(e) {
      arg::err("(from {.fun MASS::glm.nb}): {conditionMessage(e)}")
    })

    mu <- fit$fitted.values
    coefs <- fit$coefficients

    #`glm.nb()` alternates the coefficients and theta only to a tolerance, so theta
    #is re-solved given the final means, which is what the estimating equations
    #below require
    disp <- .negbin_theta_ml(treat, mu, s.weights)
  }
  else {
    family <- {
      if (fam$name == "poisson") poisson(link = fam$link)
      else Gamma(link = fam$link)
    }

    #With a link other than the log, `glm()`'s default starting values can give
    #means outside the support; starting from the intercept-only model avoids that
    start <- {
      if (fam$link == "log") NULL
      else c(fam$linkobj$linkfun(w.m(treat, s.weights)), rep.int(0, ncol(covs)))
    }

    rlang::try_fetch({verbosely({
      if (isTRUE(...get("quick"))) {
        fit <- do.call(stats::glm.fit, list(y = treat,
                                            x = Xtreat,
                                            weights = s.weights,
                                            family = family,
                                            start = start,
                                            control = as.list(...get("control"))),
                       quote = TRUE)
      }
      else {
        data <- data.frame(treat, covs)
        formula <- if (ncol(covs) > 0L) formula(data) else treat ~ 1

        fit <- do.call(stats::glm, list(formula, data = data,
                                        weights = s.weights,
                                        family = family,
                                        start = start,
                                        control = as.list(...get("control"))),
                       quote = TRUE)
      }
    }, verbose = verbose)},
    warning = function(w) {
      arg::wrn("(from {.fun glm}): {conditionMessage(w)}")
      invokeRestart("muffleWarning")
    },
    error = function(e) {
      arg::err("(from {.fun glm}): {conditionMessage(e)}")
    })

    mu <- fit$fitted.values
    coefs <- fit$coefficients

    disp <- if (fam$name == "Gamma") .gamma_shape_ml(treat, mu, s.weights) else NULL
  }

  w <- .get_w_from_cont_family(mu, treat, disp, fam$name, s.weights)

  coefs <- setNames(coefs, colnames(Xtreat))

  has_disp <- fam$name %in% c("negbin", "Gamma")
  disp_psi <- switch(fam$name,
                     negbin = .negbin_theta_psi,
                     Gamma = .gamma_shape_psi,
                     NULL)
  linkinv <- fam$linkobj$linkinv

  Mparts <- list(
    psi_treat = function(Btreat, Xtreat, A, SW) {
      SW <- SW / mean_fast(SW)

      if (has_disp) {
        disp <- exp(Btreat[1L])
        B <- Btreat[-1L]
        mu <- linkinv(drop(Xtreat %*% B))

        cbind(disp_psi(disp, A, mu, SW),
              .cont_family_psi_coef(B, Xtreat, A, SW, fam, disp))
      }
      else {
        .cont_family_psi_coef(Btreat, Xtreat, A, SW, fam)
      }
    },
    wfun = function(Btreat, Xtreat, A) {
      if (has_disp) {
        disp <- exp(Btreat[1L])
        B <- Btreat[-1L]
      }
      else {
        disp <- NULL
        B <- Btreat
      }

      mu <- linkinv(drop(Xtreat %*% B))

      .get_w_from_cont_family(mu, A, disp, fam$name, s.weights)
    },
    Xtreat = Xtreat,
    A = treat,
    btreat = {
      if (has_disp) c(setNames(log(disp),
                               if (fam$name == "negbin") "log(theta)" else "log(shape)"),
                      coefs)
      else coefs
    }
  )

  list(w = w, fit.obj = fit, Mparts = Mparts)
}

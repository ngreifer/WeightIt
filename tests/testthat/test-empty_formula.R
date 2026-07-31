# An empty model formula (`A ~ 1`, `.cens(C) ~ 1`) leaves nothing for any method to
# model or balance, and every method's target is met by the same weights. Most
# methods cannot run at all on a zero-column `covs`, so `weightit.fit()` computes
# those weights with an intercept-only GLM instead, dropping the method-specific
# arguments. That is an implementation shortcut, not a change of method: the
# requested `method` is what gets reported. See `.route_empty_covs()`.

skip_on_cran()

eps <- if (capabilities("long.double")) 1e-5 else 1e-3

test_data <- readRDS(test_path("fixtures", "test_data.rds"))

make_data <- function(seed = 123L) {
  set.seed(seed)
  d <- test_data
  d$C <- rbinom(nrow(d), 1L, prob = plogis(-0.9 + 0.8 * d$X1 - 0.5 * d$X3))
  d$G <- factor(d$X5, labels = c("g0", "g1"))
  d
}

# The installed methods reachable with a given treatment type, named as
# `.weightit_methods` names them
methods_for <- function(treat.type) {
  Filter(function(m) {
    treat.type %in% .weightit_methods[[m]]$treat_type &&
      all(vapply(.weightit_methods[[m]]$packages_needed, rlang::is_installed,
                 logical(1L)))
  }, names(.weightit_methods))
}

# Arguments some methods require, so that the only thing under test is the empty
# formula
extra_args <- list(gbm = list(criterion = "smd.mean"),
                   super = list(SL.library = c("SL.glm", "SL.mean")))

fit_empty <- function(f, m, data, ...) {
  do.call("weightit", c(list(f, data = data, method = m),
                        extra_args[[m]], list(...)))
}

# ---- Every method routes, and they all agree -------------------------------

test_that("every method gives the glm weights for an empty formula", {
  d <- make_data()

  for (tt in c("binary", "multinomial", "continuous", "censoring")) {
    f <- switch(tt,
                binary = A ~ 1,
                multinomial = Am ~ 1,
                continuous = Ac ~ 1,
                censoring = .cens(C) ~ 1)

    W_glm <- weightit(f, data = d, method = "glm")

    for (m in setdiff(methods_for(tt), "glm")) {
      W <- fit_empty(f, m, d)

      # The shortcut is invisible: `method` is reported as supplied
      expect_identical(as.character(W$method), m,
                       label = sprintf("reported method for %s/%s", tt, m))
      expect_equal(unname(W$weights), unname(W_glm$weights), tolerance = eps)
    }
  }
})

test_that("the shortcut is silent", {
  d <- make_data()

  # Computing the weights with a GLM rather than the requested method is an
  # implementation detail: the answer is the same either way
  expect_no_message(weightit(A ~ 1, data = d, method = "ebal"))
  expect_no_message(weightit(.cens(C) ~ 1, data = d, method = "energy"))
  expect_no_message(weightit(A ~ 1, data = d, method = "glm"))

  # A non-empty formula is untouched
  expect_no_message(W <- weightit(A ~ X1, data = d, method = "ebal"))
  expect_identical(as.character(W$method), "ebal")
})

test_that("method-specific arguments are dropped rather than passed on", {
  d <- make_data()

  # Arguments that tune a method that is no longer used, or that describe
  # covariates that do not exist, cannot change the result -- and must not error
  W <- weightit(A ~ 1, data = d, method = "glm")

  expect_equal(unname(fit_empty(A ~ 1, "optweight", d, tols = 0.001)$weights),
               unname(W$weights), tolerance = eps)
  expect_equal(unname(fit_empty(A ~ 1, "energy", d, moments = 2L,
                                int = TRUE)$weights),
               unname(W$weights), tolerance = eps)
  expect_equal(unname(fit_empty(A ~ 1, "cbps", d, over = TRUE)$weights),
               unname(W$weights), tolerance = eps)

  # Including glm's own arguments: an empty formula always gives the same weights
  expect_equal(unname(weightit(A ~ 1, data = d, method = "glm",
                               link = "probit")$weights),
               unname(W$weights), tolerance = eps)
})

test_that("estimand, focal, by, and s.weights survive the routing", {
  d <- make_data()

  W_att <- fit_empty(A ~ 1, "ebal", d, estimand = "ATT")
  expect_identical(W_att$estimand, "ATT")
  expect_equal(unname(W_att$weights),
               unname(weightit(A ~ 1, data = d, method = "glm",
                               estimand = "ATT")$weights),
               tolerance = eps)

  W_by <- fit_empty(A ~ 1, "energy", d, by = ~G)
  expect_equal(unname(W_by$weights),
               unname(weightit(A ~ 1, data = d, method = "glm", by = ~G)$weights),
               tolerance = eps)

  W_sw <- fit_empty(A ~ 1, "cbps", d, s.weights = "SW")
  expect_equal(unname(W_sw$weights),
               unname(weightit(A ~ 1, data = d, method = "glm",
                               s.weights = "SW")$weights),
               tolerance = eps)
})

test_that("M-estimation is available after routing", {
  skip_if_not_installed("rootSolve")

  d <- make_data()

  # The weights come from a glm, so the glm M-estimation parts apply even for
  # methods that supply none of their own
  for (m in intersect(c("energy", "cfd"), methods_for("censoring"))) {
    W <- fit_empty(.cens(C) ~ 1, m, d)

    expect_false(is_null(attr(W, "Mparts", exact = TRUE)))
    expect_M_parts_okay(W, tolerance = eps)
  }
})

test_that("an empty continuous formula gives weights of exactly 1", {
  d <- make_data()

  # With nothing to condition on, the conditional density of the treatment is its
  # marginal density, so the weights are 1 rather than merely close to it: the
  # numeric marginalization that would otherwise compute them is approximate
  for (m in methods_for("continuous")) {
    W <- fit_empty(Ac ~ 1, m, d)

    expect_identical(unname(W$weights), rep.int(1, nrow(d)),
                     label = sprintf("weights for method = \"%s\"", m))

    # Nothing was estimated, so there is nothing for M-estimation to account for
    expect_null(attr(W, "Mparts", exact = TRUE))
    expect_null(attr(W, "Mparts.list", exact = TRUE))
  }

  # Constant weights of 1 are the right answer here, not a sign of failure
  expect_no_warning(weightit(Ac ~ 1, data = d, method = "glm"))

  # And an outcome model still gets its standard errors, treating the (unestimated)
  # weights as fixed
  W <- weightit(Ac ~ 1, data = d, method = "glm")
  fit <- glm_weightit(Y_C ~ Ac, data = d, weightit = W)
  expect_true(all(is.finite(sqrt(diag(vcov(fit))))))
})

# ---- What is deliberately not routed ---------------------------------------

test_that("method = NULL is not routed", {
  d <- make_data()

  # `NULL` estimates nothing to begin with
  expect_no_message(W <- weightit(A ~ 1, data = d, method = NULL))
  expect_null(W$method)
  expect_true(all(W$weights == 1))

  expect_no_message(Wc <- weightit(.cens(C) ~ 1, data = d, method = NULL))
  expect_equal(unname(Wc$weights), 1 - d$C, tolerance = eps)
})

test_that("a user-defined method is not routed", {
  d <- make_data()

  # The user's function may handle a covariate-free problem deliberately
  my.fun <- function(treat, covs, ...) {
    2 + treat
  }

  expect_no_message(W <- weightit(A ~ 1, data = d, method = my.fun))
  expect_equal(unname(W$weights), 2 + d$A, tolerance = eps)
})

test_that("a supplied ps is not routed", {
  d <- make_data()

  ps <- plogis(-0.9 + 0.8 * d$X1)

  expect_no_message(W <- weightit(A ~ 1, data = d, ps = ps))
  expect_equal(unname(W$weights), unname(get_w_from_ps(ps, d$A)),
               tolerance = eps)
})

test_that("a formula with only random effects terms is not routed", {
  skip_if_not_installed("lme4")

  d <- make_data()

  # `A ~ (1 | G)` also has a zero-column fixed-effects `covs`, but it is a real
  # model: the propensity scores must be cluster-specific, not constant
  expect_no_message(W <- weightit(A ~ (1 | G), data = d, method = "glm"))

  expect_gt(nunique(W$ps), 1L)

  expect_not_equal(unname(W$weights),
                   unname(weightit(A ~ 1, data = d, method = "glm")$weights))
})

# ---- weightit.fit() and weightitMSM() -------------------------------------

test_that("weightit.fit() accepts a zero-column covs matrix", {
  d <- make_data()

  W <- weightit.fit(matrix(0, nrow = nrow(d), ncol = 0L), treat = d$A,
                    method = "ebal")

  expect_identical(as.character(W$method), "ebal")
  expect_equal(unname(W$weights),
               unname(weightit(A ~ 1, data = d, method = "glm")$weights),
               tolerance = eps)

  # The requested method's package is not needed, since it is never called
  expect_no_error(weightit.fit(matrix(0, nrow = nrow(d), ncol = 0L),
                               treat = d$A, method = "optweight"))
})

make_msm_data <- function(seed = 1234L) {
  set.seed(seed)
  d <- msmdata
  d$C_2 <- rbinom(nrow(d), 1L, prob = plogis(-2 + 0.2 * d$X1_1 - 0.3 * d$A_1))
  is.na(d[d$C_2 == 1L, c("X1_2", "X2_2", "A_3")]) <- TRUE
  d
}

msm_empty_formulas <- list(A_1 ~ X1_0 + X2_0,
                           A_2 ~ X1_1 + X2_1 + A_1,
                           .cens(C_2) ~ 1,
                           A_3 ~ X1_2 + X2_2 + A_2)

test_that("weightitMSM() routes each empty formula independently", {
  skip_if_not_installed("rootSolve")

  d <- make_msm_data()
  cens <- d$C_2 == 1L

  # `cbps` cannot run on a zero-column `covs` at all, so without the shortcut this
  # errors. Only the empty time point takes it; the others really are fit with
  # `method`, which the fit objects show directly.
  # (`is.MSM.method = FALSE` messages that cbps could fit one model instead)
  W <- suppressMessages(weightitMSM(msm_empty_formulas, data = d,
                                    method = "cbps", is.MSM.method = FALSE,
                                    include.obj = TRUE))

  expect_identical(as.character(W$method), "cbps")
  expect_named(W$obj, c("A_1", "A_2", "C_2", "A_3"))
  expect_s3_class(W$obj$C_2, "glm")
  expect_false(inherits(W$obj$A_1, "glm"))

  w1 <- weightit(A_1 ~ X1_0 + X2_0, data = d, method = "cbps")$weights
  w2 <- weightit(A_2 ~ X1_1 + X2_1 + A_1, data = d, method = "cbps")$weights
  cw <- ifelse(cens, 0, 1 / mean(!cens))
  w3 <- rep.int(1, nrow(d))
  w3[!cens] <- weightit(A_3 ~ X1_2 + X2_2 + A_2, data = d[!cens, ],
                        method = "cbps")$weights

  expect_equal(unname(W$weights), unname(w1 * w2 * cw * w3), tolerance = eps)
  expect_true(all(W$weights[cens] == 0))
})

test_that("simultaneous CBPS handles a time point with no covariates", {
  skip_if_not_installed("rootSolve")

  d <- make_msm_data()

  # `is.MSM.method = TRUE` estimates one set of weights for all time points at
  # once, so a single empty formula cannot change the method. Such a time point
  # instead contributes only the intercept condition, which `svd()` cannot be
  # asked for.
  W <- weightitMSM(msm_empty_formulas, data = d, method = "cbps",
                   is.MSM.method = TRUE, include.obj = TRUE)

  expect_true(all(W$weights[d$C_2 == 1L] == 0))
  expect_true(all(W$weights[d$C_2 == 0L] > 0))
  expect_true(all(is.finite(W$weights)))

  # All the moment conditions, treatment and censoring, are still solved jointly
  expect_lt(W$obj$value, 1e-8)

  # And each treatment condition holds within its own risk set
  for (tp in c("A_1", "A_2", "A_3")) {
    cvs <- switch(tp, A_1 = c("X1_0", "X2_0"), A_2 = c("X1_1", "X2_1"),
                  A_3 = c("X1_2", "X2_2"))
    A <- d[[tp]]
    X <- as.matrix(d[cvs])

    ok <- W$at.risk[, tp] & !is.na(A)
    t1 <- ok & A == 1
    t0 <- ok & A == 0

    m1 <- colSums(X[t1, , drop = FALSE] * W$weights[t1]) / sum(W$weights[t1])
    m0 <- colSums(X[t0, , drop = FALSE] * W$weights[t0]) / sum(W$weights[t0])

    expect_lt(max(abs(m1 - m0)), 1e-6)
  }

  # An empty treatment formula also runs
  expect_no_error(weightitMSM(list(A_1 ~ 1, A_2 ~ X1_1 + A_1), data = msmdata,
                              method = "cbps", is.MSM.method = TRUE))
})

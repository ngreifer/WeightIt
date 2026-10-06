test_that("Binary treatment", {
  skip_on_cran()
  skip_if_not_installed("rootSolve")
  skip_if_not_installed("cobalt")
  skip_if_not_installed("patrick")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  set.seed(123)
  base_weights <- runif(nrow(test_data))

  expect_no_condition({
    W0 <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                   data = test_data, method = "ebal", estimand = "ATE",
                   include.obj = TRUE, solver = "optim")
  })

  expect_M_parts_okay(W0, tolerance = eps)

  # Weights from the configurations already run, so each new one can be checked
  # against all of them
  seen <- new.env()

  patrick::with_parameters_test_that(
    "Ebal: sw = {sw}, bw = {bw}, estimand = {estimand}",
    {
      W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                    data = test_data, method = "ebal", estimand = estimand,
                    s.weights = if (sw) "SW" else NULL,
                    base.weights = if (bw) base_weights else NULL,
                    include.obj = TRUE, solver = "multiroot")

      expect_M_parts_okay(W, tolerance = eps)
      expect_equal(cobalt::col_w_smd(W$covs, W$treat, W$weights,
                                     s.weights = W$s.weights),
                   0 * cobalt::col_w_smd(W$covs, W$treat,
                                         s.weights = W$s.weights),
                   expected.label = "all 0s",
                   tolerance = eps)

      expect_true(is_null(W$ps))
      expect_false(is_null(W$obj))

      if (estimand %in% c("ATT", "ATC")) {
        expect_ATT_weights_okay(W, tolerance = eps)
      }

      for (i in 0:1) {
        e <- {
          if (estimand == "ATT" && i == 1) expect_equal
          else if (estimand == "ATC" && i == 0) expect_equal
          else expect_not_equal
        }

        e(unname(W$weights[W$treat == i]),
          rep(1, sum(W$treat == i)),
          label = sprintf("%s weights", i),
          expected.label = "all 1s",
          tolerance = eps)
      }

      for (other in ls(seen)) {
        expect_not_equal(unname(W$weights), seen[[other]],
                         expected.label = sprintf("weights for %s", other),
                         tolerance = eps)
      }

      seen[[sprintf("sw = %s, bw = %s, estimand = %s", sw, bw, estimand)]] <- unname(W$weights)
    },
    .cases = expand.grid(estimand = c("ATE", "ATT", "ATC"),
                         bw = c(FALSE, TRUE),
                         sw = c(FALSE, TRUE),
                         stringsAsFactors = FALSE)
  )

  # Estimands
  expect_error({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "ebal", estimand = "ATO")
  }, "not an allowable estimand", ignore.case = TRUE)

  #Non-full rank
  W <- expect_no_condition({
    weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 +
               I(1 - X5) + I(X9 * 2),
             data = test_data, method = "ebal", estimand = "ATE",
             include.obj = TRUE, solver = "optim")
  })

  expect_M_parts_okay(W, tolerance = eps)
  expect_equal(W$weights, W0$weights, tolerance = eps)

  # All categorical covariates (issue #86)
  W <- expect_no_condition({
    weightit(A ~ cut(X1, 3) + cut(X2, 3) + cut(X3, 3),
             data = test_data, method = "ebal", estimand = "ATE",
             include.obj = TRUE, solver = "optim", reltol = 1e-12)
  })

  expect_M_parts_okay(W, tolerance = eps)

  # tols > 0
  W <- expect_no_condition({
    weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
              data = test_data, method = "ebal", estimand = "ATE",
              include.obj = TRUE, tols = .05)
  })

  expect_not_equal(W$weights, W0$weights)

  expect_true(all(abs(cobalt::bal.tab(W)$Balance$Diff.Adj) <= .05 + eps)) #None worse than tols
  expect_true(any(abs(abs(cobalt::bal.tab(W)$Balance$Diff.Adj) - .05) <= eps)) #Some exactly tols
  expect_true(any(abs(cobalt::bal.tab(W)$Balance$Diff.Adj) > eps)) #Some worse than 0

  W <- expect_no_condition({
    weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
             data = test_data, method = "ebal", estimand = "ATE",
             include.obj = TRUE, tols = .05, s.weights = "SW")
  })

  expect_true(all(abs(cobalt::bal.tab(W)$Balance$Diff.Adj) <= .05 + eps)) #None worse than tols
  expect_true(any(abs(abs(cobalt::bal.tab(W)$Balance$Diff.Adj) - .05) <= eps)) #Some exactly tols
  expect_true(any(abs(cobalt::bal.tab(W)$Balance$Diff.Adj) > eps)) #Some worse than 0

  #Should be equivalent to CBPS and IPT with logit link for ATT
  patrick::with_parameters_test_that(
    "Ebal matches CBPS/IPT for ATT: sw = {sw}",
    {
      W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                    data = test_data, method = "ebal", estimand = "ATT",
                    s.weights = if (sw) "SW" else NULL,
                    include.obj = TRUE, solver = "optim")

      Wcbps <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                        data = test_data, method = "cbps", estimand = "ATT",
                        s.weights = if (sw) "SW" else NULL,
                        link = "logit", solver = "multiroot",
                        include.obj = TRUE)

      expect_equal(ESS(W$weights[W$treat == 0] * W$s.weights[W$treat == 0]),
                   ESS(Wcbps$weights[Wcbps$treat == 0] * Wcbps$s.weights[Wcbps$treat == 0]),
                   expected.label = "ESS for CBPS",
                   tolerance = .01)

      Wipt <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                       data = test_data, method = "ipt", estimand = "ATT",
                       s.weights = if (sw) "SW" else NULL,
                       link = "logit",
                       include.obj = TRUE)

      expect_equal(ESS(W$weights[W$treat == 0] * W$s.weights[W$treat == 0]),
                   ESS(Wipt$weights[Wipt$treat == 0] * Wipt$s.weights[Wipt$treat == 0]),
                   expected.label = "ESS for IPT",
                   tolerance = .01)
    },
    sw = c(FALSE, TRUE)
  )
})

test_that("Multi-category treatment", {
  skip_on_cran()
  skip_if_not_installed("rootSolve")
  skip_if_not_installed("cobalt")
  skip_if_not_installed("patrick")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  set.seed(123)
  base_weights <- runif(nrow(test_data))

  expect_no_condition({
    W0 <- weightit(Am ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                   data = test_data, method = "ebal", estimand = "ATE",
                   include.obj = TRUE, solver = "optim")
  })

  # Weights from the configurations already run, so each new one can be checked
  # against all of them
  seen <- new.env()

  patrick::with_parameters_test_that(
    "Ebal: sw = {sw}, bw = {bw}, estimand = {estimand}",
    {
      W <- weightit(Am ~ X1 + X2 + X3 + X4 + X5,
                    data = test_data, method = "ebal", estimand = estimand,
                    focal = if (estimand == "ATE") NULL else "T",
                    s.weights = if (sw) "SW" else NULL,
                    base.weights = if (bw) base_weights else NULL,
                    include.obj = TRUE, solver = "multiroot")

      expect_M_parts_okay(W, tolerance = eps)
      for (tt in combn(levels(W$treat), 2, simplify = FALSE)) {
        in_tt <- W$treat %in% tt
        expect_equal(cobalt::col_w_smd(W$covs[in_tt,], W$treat[in_tt], W$weights[in_tt],
                                       s.weights = W$s.weights[in_tt]),
                     0 * cobalt::col_w_smd(W$covs[in_tt,], W$treat[in_tt],
                                           s.weights = W$s.weights[in_tt]),
                     label = sprintf("SMDs for %s", paste(tt, collapse = " vs. ")),
                     expected.label = "all 0s",
                     tolerance = eps)
      }

      expect_true(is_null(W$ps))
      expect_false(is_null(W$obj))

      if (estimand %in% c("ATT", "ATC")) {
        expect_ATT_weights_okay(W, tolerance = eps)
      }

      for (i in levels(W$treat)) {
        e <- {
          if (estimand == "ATT" && i == W$focal) expect_equal
          else expect_not_equal
        }

        e(unname(W$weights[W$treat == i]),
          rep(1, sum(W$treat == i)),
          label = sprintf("%s weights", i),
          expected.label = "all 1s",
          tolerance = eps)
      }

      for (other in ls(seen)) {
        expect_not_equal(unname(W$weights), seen[[other]],
                         expected.label = sprintf("weights for %s", other),
                         tolerance = eps)
      }

      seen[[sprintf("sw = %s, bw = %s, estimand = %s", sw, bw, estimand)]] <- unname(W$weights)
    },
    .cases = expand.grid(estimand = c("ATE", "ATT"),
                         bw = c(FALSE, TRUE),
                         sw = c(FALSE, TRUE),
                         stringsAsFactors = FALSE)
  )

  # tols > 0
  W <- expect_no_condition({
    weightit(Am ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
             data = test_data, method = "ebal", estimand = "ATE",
             include.obj = TRUE, tols = .05)
  })

  expect_not_equal(W$weights, W0$weights)

  expect_true(all(abs(cobalt::bal.tab(W)$Balance$Max.Diff.Adj) <= .05 + eps)) #None worse than tols
  expect_true(any(abs(abs(cobalt::bal.tab(W)$Balance$Max.Diff.Adj) - .05) <= eps)) #Some exactly tols
  expect_true(any(abs(cobalt::bal.tab(W)$Balance$Max.Diff.Adj) > eps)) #Some worse than 0

  W <- expect_no_condition({
    weightit(Am ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
             data = test_data, method = "ebal", estimand = "ATE",
             include.obj = TRUE, tols = .05, s.weights = "SW")
  })

  expect_true(all(abs(cobalt::bal.tab(W)$Balance$Max.Diff.Adj) <= .05 + eps)) #None worse than tols
  expect_true(any(abs(abs(cobalt::bal.tab(W)$Balance$Max.Diff.Adj) - .05) <= eps)) #Some exactly tols
  expect_true(any(abs(cobalt::bal.tab(W)$Balance$Max.Diff.Adj) > eps)) #Some worse than 0
})

test_that("Continuous treatment", {
  skip_on_cran()
  skip_if_not_installed("rootSolve")
  skip_if_not_installed("cobalt")
  skip_if_not_installed("patrick")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  set.seed(123)
  base_weights <- runif(nrow(test_data))

  expect_no_condition({
    W0 <- weightit(Ac ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                   data = test_data, method = "ebal",
                   include.obj = TRUE, solver = "optim")
  })

  expect_M_parts_okay(W0, tolerance = eps)

  # Weights from the configurations already run, so each new one can be checked
  # against all of them
  seen <- new.env()

  patrick::with_parameters_test_that(
    "Ebal: sw = {sw}, bw = {bw}, d.moments = {d.moments}",
    {
      W <- weightit(Ac ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                    data = test_data, method = "ebal",
                    d.moments = d.moments,
                    s.weights = if (sw) "SW" else NULL,
                    base.weights = if (bw) base_weights else NULL,
                    include.obj = TRUE, solver = "multiroot")

      expect_M_parts_okay(W, tolerance = eps)
      expect_equal(cobalt::col_w_cov(W$covs, W$treat, W$weights, std = TRUE,
                                     s.weights = W$s.weights),
                   0 * cobalt::col_w_cov(W$covs, W$treat, std = TRUE,
                                         s.weights = W$s.weights),
                   expected.label = "all 0s",
                   tolerance = eps)

      expect_equal(cobalt::col_w_mean(cbind(poly(W$treat, d.moments), W$covs), W$weights,
                                      s.weights = W$s.weights),
                   cobalt::col_w_mean(cbind(poly(W$treat, d.moments), W$covs),
                                      s.weights = W$s.weights),
                   expected.label = "unweighted means",
                   tolerance = eps)

      expect_true(is_null(W$ps))
      expect_false(is_null(W$obj))

      for (other in ls(seen)) {
        expect_not_equal(unname(W$weights), seen[[other]],
                         expected.label = sprintf("weights for %s", other),
                         tolerance = eps)
      }

      seen[[sprintf("sw = %s, bw = %s, d.moments = %s", sw, bw, d.moments)]] <- unname(W$weights)
    },
    .cases = expand.grid(d.moments = c(1, 3),
                         bw = c(FALSE, TRUE),
                         sw = c(FALSE, TRUE))
  )

  #Non-full rank
  W <- expect_no_condition({
    weightit(Ac ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 +
               I(1 - X5) + I(X9 * 2),
             data = test_data, method = "ebal",
             include.obj = TRUE, solver = "optim")
  })

  expect_M_parts_okay(W, tolerance = eps)
  expect_equal(W$weights, W0$weights, tolerance = eps)

  # tols > 0
  W <- expect_no_condition({
    weightit(Ac ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
             data = test_data, method = "ebal", estimand = "ATE",
             include.obj = TRUE, tols = .05)
  })

  expect_not_equal(W$weights, W0$weights)

  expect_true(all(abs(cobalt::bal.tab(W)$Balance$Corr.Adj) <= .05 + eps)) #None worse than tols
  expect_true(any(abs(abs(cobalt::bal.tab(W)$Balance$Corr.Adj) - .05) <= eps)) #Some exactly tols
  expect_true(any(abs(cobalt::bal.tab(W)$Balance$Corr.Adj) > eps)) #Some worse than 0

  W <- expect_no_condition({
    weightit(Ac ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
             data = test_data, method = "ebal", estimand = "ATE",
             include.obj = TRUE, tols = .05, s.weights = "SW")
  })

  expect_true(all(abs(cobalt::bal.tab(W)$Balance$Corr.Adj) <= .05 + eps)) #None worse than tols
  expect_true(any(abs(abs(cobalt::bal.tab(W)$Balance$Corr.Adj) - .05) <= eps)) #Some exactly tols
  expect_true(any(abs(cobalt::bal.tab(W)$Balance$Corr.Adj) > eps)) #Some worse than 0
})

test_that("Continuous treatment: M-estimation parts assemble into a variance", {
  skip_on_cran()
  skip_if_not_installed("rootSolve")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  # The continuous parts supply `dw_dBtreat` but no `hess_treat`, so assembling them
  # takes the mixed analytic/numeric branch in `.get_glm_weightit_vcov()`. The parts
  # themselves were checked above; nothing verified that the assembly produces a usable
  # variance, which is where a wrong standard error would appear without any error.
  W <- weightit(Ac ~ X1 + X2 + X3, data = test_data, method = "ebal")

  expect_M_parts_okay(W, tolerance = eps)
  expect_false(is_null(attr(W, "Mparts")$dw_dBtreat))
  expect_null(attr(W, "Mparts")$hess_treat)

  fit <- lm_weightit(Y_C ~ Ac, data = test_data, weightit = W)

  se_asympt <- sqrt(diag(vcov(fit)))
  se_hc0 <- sqrt(diag(vcov(fit, vcov = "HC0")))

  expect_true(all(is.finite(se_asympt)))
  expect_true(all(se_asympt > 0))

  # Accounting for estimation of the weights changes the answer, so the parts are
  # actually contributing rather than being dropped
  expect_not_equal(unname(se_asympt), unname(se_hc0))
})

test_that("Continuous treatment: d.moments is a floor on the target moments", {
  skip_on_cran()

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  covs <- c("X1", "X2", "X3", "X4")

  # How many leading moments of `v` the weights hold at their unweighted values
  n_held <- function(W, v) {
    x <- test_data[[v]]
    held <- vapply(1:3, function(k) {
      abs(weighted.mean(x^k, W$weights) - mean(x^k)) / sd(x^k) < 1e-6
    }, logical(1L))

    as.integer(sum(cumprod(held)))
  }

  fit <- function(...) {
    weightit(Ac ~ X1 + X2 + X3 + X4, data = test_data, method = "ebal", ...)
  }

  # `moments` alone: each covariate's target moments are its own entry, and the
  # covariates not named get the default of 1
  W <- fit(moments = c(X1 = 2, X2 = 3))
  expect_identical(vapply(covs, n_held, 0L, W = W),
                   c(X1 = 2L, X2 = 3L, X3 = 1L, X4 = 1L))

  # `d.moments` at or below every entry of `moments` changes nothing for the
  # covariates -- it is a floor, not a replacement
  W2 <- fit(moments = c(X1 = 2, X2 = 3), d.moments = 2)
  expect_identical(vapply(covs, n_held, 0L, W = W2),
                   c(X1 = 2L, X2 = 3L, X3 = 1L, X4 = 1L))

  # Raising it above an entry raises that covariate, and reaches the covariates
  # `moments` does not name. This used to collapse every covariate to its mean,
  # i.e. asking for more moments produced fewer.
  W3 <- fit(moments = c(X1 = 2, X2 = 3), d.moments = 3)
  expect_identical(vapply(covs, n_held, 0L, W = W3),
                   c(X1 = 3L, X2 = 3L, X3 = 3L, X4 = 3L))

  # `d.moments` governs the treatment on its own; `moments` does not raise it
  n_held_treat <- function(W) {
    a <- test_data$Ac
    as.integer(sum(cumprod(vapply(1:3, function(k) {
      abs(weighted.mean(a^k, W$weights) - mean(a^k)) / sd(a^k) < 1e-6
    }, logical(1L)))))
  }

  expect_identical(n_held_treat(W), 1L)
  expect_identical(n_held_treat(W2), 2L)
  expect_identical(n_held_treat(W3), 3L)

  # The treatment-covariate correlations follow `moments`, not `d.moments`: X1's
  # cube is held to its target above but its correlation with the treatment is
  # not driven to zero
  wcor <- function(W, v, k) {
    abs(cov.wt(cbind(test_data$Ac, test_data[[v]]^k), wt = W$weights, cor = TRUE)$cor[1L, 2L])
  }

  expect_lt(wcor(W3, "X1", 2L), 1e-6)
  expect_gt(wcor(W3, "X1", 3L), 1e-3)
})

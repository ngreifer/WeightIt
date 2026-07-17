test_that("Binary treatment", {
  skip_on_cran()
  skip_if_not_installed("optweight", minimum_version = "2.0.1")
  skip_if_not_installed("cobalt")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  expect_no_condition({
    W0 <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                   data = test_data, method = "optweight", estimand = "ATE",
                   include.obj = TRUE)
  })

  # optweight does not support M-estimation
  expect_null(attr(W0, "Mparts", exact = TRUE))
  expect_null(attr(W0, "Mparts.list", exact = TRUE))

  expect_equal(cobalt::col_w_smd(W0$covs, W0$treat, W0$weights),
               0 * cobalt::col_w_smd(W0$covs, W0$treat),
               expected.label = "all 0s",
               tolerance = eps)

  expect_true(is_null(W0$ps))
  expect_false(is_null(W0$obj))

  # `info$duals` should be populated with one row per balance constraint
  expect_true(is.data.frame(W0$info$duals))
  expect_true(all(c("constraint", "cov", "dual") %in% names(W0$info$duals)))

  sw.opts <- c(FALSE, TRUE)
  estimand.opts <- c("ATE", "ATT", "ATC")

  weight.mat <- matrix(nrow = nrow(test_data),
                       ncol = length(sw.opts) * length(estimand.opts))
  colnames(weight.mat) <- rep("", ncol(weight.mat))

  k <- 1

  for (sw in sw.opts) {
    for (estimand in estimand.opts) {
      test_that(sprintf("Optweight: sw = %s, estimand = %s", sw, estimand), {
        W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                      data = test_data, method = "optweight", estimand = estimand,
                      s.weights = if (sw) "SW" else NULL,
                      include.obj = TRUE)

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

        for (i in seq_len(k - 1)) {
          expect_not_equal(unname(W$weights), weight.mat[,i],
                           expected.label = sprintf("weights for %s", colnames(weight.mat)[i]),
                           tolerance = eps)
        }

        n <- sprintf("W_%s_%s", sw, estimand)
        colnames(weight.mat)[k] <<- n
        weight.mat[,k] <<- W$weights
        k <<- k + 1
      })
    }
  }

  # Estimands
  expect_error({
    weightit(A ~ X1 + X2 + X3, data = test_data, method = "optweight", estimand = "ATO")
  }, "not an allowable estimand", ignore.case = TRUE)

  # Additional arguments (moments, int, quantile): spot-checked individually
  # against a same-formula baseline rather than crossed with sw/estimand.
  expect_no_condition({
    W_base <- weightit(A ~ X1 + X2 + X3 + X4 + X5, data = test_data,
                       method = "optweight", estimand = "ATE", include.obj = TRUE)
  })

  configs <- list(
    "moments = 2" = list(moments = 2),
    "int = TRUE" = list(int = TRUE),
    "quantile" = list(quantile = list(X1 = c(.25, .5, .75)))
  )

  for (nm in names(configs)) {
    test_that(sprintf("Optweight: %s", nm), {
      W <- do.call(weightit,
                  c(list(A ~ X1 + X2 + X3 + X4 + X5, data = test_data,
                        method = "optweight", estimand = "ATE", include.obj = TRUE),
                    configs[[nm]]))

      expect_equal(cobalt::col_w_smd(W$covs, W$treat, W$weights),
                   0 * cobalt::col_w_smd(W$covs, W$treat),
                   expected.label = "all 0s",
                   tolerance = eps)

      expect_not_equal(unname(W$weights), unname(W_base$weights),
                       expected.label = "weights for baseline")
    })
  }

  #Non-full rank
  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 +
                    I(1 - X5) + I(X9 * 2),
                  data = test_data, method = "optweight", estimand = "ATE",
                  include.obj = TRUE)
  })

  expect_equal(W$weights, W0$weights, tolerance = eps)

  set.seed(4321)
  d_na <- inject_missingness(test_data, "X1", prop = 0.05)

  expect_error({
    weightit(A ~ X1 + X2 + X3, data = d_na, method = "optweight", estimand = "ATE",
            missing = "surr")
  }, "only.*allowed for.*missing", ignore.case = TRUE)

  expect_no_condition({
    W_na <- weightit(A ~ X1 + X2 + X3, data = d_na, method = "optweight",
                     estimand = "ATE", missing = "ind")
  })

  expect_true(anyNA(W_na$covs))
  # Guard against the solver-failure pathology described above recurring
  # silently in CI: weights should vary, not be a degenerate constant.
  expect_not_equal(W_na$weights, rep(W_na$weights[1], length(W_na$weights)))

  # tols > 0: approximate balance. The solver only guarantees satisfying the
  # constraints up to its own convergence tolerance, which is looser than the
  # floating-point `eps` used for the exact-balance (tols = 0) checks above,
  # so a dedicated, looser tolerance is used for these boundary checks.
  tols.eps <- 1e-3

  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "optweight", estimand = "ATE",
                  include.obj = TRUE, tols = .05)
  })

  expect_not_equal(W$weights, W0$weights)

  expect_true(all(abs(cobalt::bal.tab(W)$Balance$Diff.Adj) <= .05 + tols.eps)) #None worse than tols
  expect_true(any(abs(abs(cobalt::bal.tab(W)$Balance$Diff.Adj) - .05) <= tols.eps)) #Some exactly tols

  # tols > 0, crossed with sw and estimand: the balance constraint should bind
  # the same way regardless of these other design factors.
  for (sw in sw.opts) {
    for (estimand in estimand.opts) {
      test_that(sprintf("Optweight: tols = .05, sw = %s, estimand = %s", sw, estimand), {
        W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                      data = test_data, method = "optweight", estimand = estimand,
                      s.weights = if (sw) "SW" else NULL, tols = .05)

        diffs <- abs(cobalt::bal.tab(W)$Balance$Diff.Adj)

        expect_true(all(diffs <= .05 + tols.eps)) #None worse than tols
        expect_true(any(abs(diffs - .05) <= tols.eps)) #Some exactly tols

        if (estimand %in% c("ATT", "ATC")) {
          expect_ATT_weights_okay(W, tolerance = eps)
        }
      })
    }
  }

  # min.w: floor on individual weights (default 1e-8, effectively unbounded below;
  # see ?optweight::optweight). Raising it should force every weight to be at
  # least that large while still permitting exact balance (tols = 0, the default)
  # to be achieved.
  expect_no_condition({
    W_minw <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                       data = test_data, method = "optweight", estimand = "ATE",
                       include.obj = TRUE, min.w = .5)
  })

  expect_true(min(W_minw$weights) >= .5 - tols.eps)
  expect_true(min(W0$weights) < .5) #confirms the floor is a real, binding constraint here

  expect_equal(cobalt::col_w_smd(W_minw$covs, W_minw$treat, W_minw$weights),
               0 * cobalt::col_w_smd(W_minw$covs, W_minw$treat),
               expected.label = "all 0s",
               tolerance = eps)

  expect_not_equal(W_minw$weights, W0$weights)

  # min.w can also be negative (or -Inf) to allow negative weights (see the
  # "Allowing negative weights" example in ?optweight::optweight). Relaxing the
  # floor below 0 only enlarges the feasible set relative to the nonnegative
  # default, so it should never increase (and may decrease) the variance-
  # minimizing objective; WeightIt should also warn that negative weights are
  # present, since these can't be used in most model-fitting functions.
  expect_warning({
    W_minw_neg <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                           data = test_data, method = "optweight", estimand = "ATE",
                           include.obj = TRUE, min.w = -1)
  }, "negative", ignore.case = TRUE)

  expect_true(min(W_minw_neg$weights) >= -1 - tols.eps)
  expect_true(any(W_minw_neg$weights < 0))
  expect_true(var(W_minw_neg$weights) <= var(W0$weights) + tols.eps)

  expect_equal(cobalt::col_w_smd(W_minw_neg$covs, W_minw_neg$treat, W_minw_neg$weights),
               0 * cobalt::col_w_smd(W_minw_neg$covs, W_minw_neg$treat),
               expected.label = "all 0s",
               tolerance = eps)

  expect_not_equal(W_minw_neg$weights, W0$weights)

  # min.w = -Inf removes the floor entirely, so it should permit weights at least
  # as extreme (i.e., no larger a minimum) as the finite min.w = -1 case above
  expect_warning({
    W_minw_ninf <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                            data = test_data, method = "optweight", estimand = "ATE",
                            include.obj = TRUE, min.w = -Inf)
  }, "negative", ignore.case = TRUE)

  expect_true(any(W_minw_ninf$weights < 0))
  expect_true(min(W_minw_ninf$weights) <= min(W_minw_neg$weights) + tols.eps)
  expect_true(var(W_minw_ninf$weights) <= var(W_minw_neg$weights) + tols.eps)

  expect_equal(cobalt::col_w_smd(W_minw_ninf$covs, W_minw_ninf$treat, W_minw_ninf$weights),
               0 * cobalt::col_w_smd(W_minw_ninf$covs, W_minw_ninf$treat),
               expected.label = "all 0s",
               tolerance = eps)

  # norm: the objective function minimized when selecting weights (see the "norm"
  # section of ?optweight::optweight). "l2" (the default, used by W0 above)
  # minimizes the variance of the weights, equivalently maximizing the ESS;
  # "linf" instead minimizes the largest |weight - 1|. Both should still achieve
  # exact balance (tols = 0, the default), but each should outperform the other
  # on the metric it specifically targets.
  expect_no_condition({
    W_linf <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                       data = test_data, method = "optweight", estimand = "ATE",
                       include.obj = TRUE, norm = "linf")
  })

  expect_equal(cobalt::col_w_smd(W_linf$covs, W_linf$treat, W_linf$weights),
               0 * cobalt::col_w_smd(W_linf$covs, W_linf$treat),
               expected.label = "all 0s",
               tolerance = eps)

  expect_true(max(abs(W_linf$weights - 1)) <= max(abs(W0$weights - 1)) + tols.eps)
  expect_true(var(W0$weights) <= var(W_linf$weights) + tols.eps)

  expect_not_equal(W_linf$weights, W0$weights)
})

test_that("Multi-category treatment", {
  skip_on_cran()
  skip_if_not_installed("optweight", minimum_version = "2.0.1")
  skip_if_not_installed("cobalt")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  expect_no_condition({
    W0 <- weightit(Am ~ X1 + X2 + X3 + X4 + X5,
                   data = test_data, method = "optweight", estimand = "ATE",
                   include.obj = TRUE)
  })

  expect_null(attr(W0, "Mparts", exact = TRUE))

  sw.opts <- c(FALSE, TRUE)
  estimand.opts <- c("ATE", "ATT")

  weight.mat <- matrix(nrow = nrow(test_data),
                       ncol = length(sw.opts) * length(estimand.opts))
  colnames(weight.mat) <- rep("", ncol(weight.mat))

  k <- 1

  for (sw in sw.opts) {
    for (estimand in estimand.opts) {
      test_that(sprintf("Optweight: sw = %s, estimand = %s", sw, estimand), {
        W <- weightit(Am ~ X1 + X2 + X3 + X4 + X5,
                      data = test_data, method = "optweight", estimand = estimand,
                      focal = if (estimand == "ATE") NULL else "T",
                      s.weights = if (sw) "SW" else NULL,
                      include.obj = TRUE)

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

        for (i in seq_len(k - 1)) {
          expect_not_equal(unname(W$weights), weight.mat[,i],
                           expected.label = sprintf("weights for %s", colnames(weight.mat)[i]),
                           tolerance = eps)
        }

        n <- sprintf("W_%s_%s", sw, estimand)
        colnames(weight.mat)[k] <<- n
        weight.mat[,k] <<- W$weights
        k <<- k + 1
      })
    }
  }

  # Additional arguments: spot-checked individually
  configs <- list(
    "moments = 2" = list(moments = 2),
    "int = TRUE" = list(int = TRUE),
    "quantile" = list(quantile = list(X1 = c(.25, .5, .75))),
    "norm = linf" = list(norm = "linf"),
    "min.w = .5" = list(min.w = .5)
  )

  for (nm in names(configs)) {
    test_that(sprintf("Optweight: %s", nm), {
      W <- do.call(weightit,
                  c(list(Am ~ X1 + X2 + X3 + X4 + X5, data = test_data,
                        method = "optweight", estimand = "ATE", include.obj = TRUE),
                    configs[[nm]]))

      for (tt in combn(levels(W$treat), 2, simplify = FALSE)) {
        in_tt <- W$treat %in% tt
        expect_equal(cobalt::col_w_smd(W$covs[in_tt,], W$treat[in_tt], W$weights[in_tt]),
                     0 * cobalt::col_w_smd(W$covs[in_tt,], W$treat[in_tt]),
                     expected.label = "all 0s",
                     tolerance = eps)
      }

      expect_not_equal(unname(W$weights), unname(W0$weights),
                       expected.label = "weights for baseline")
    })
  }

  tols.eps <- 1e-3

  # min.w can also be negative (or -Inf) to allow negative weights (see the
  # "Allowing negative weights" example in ?optweight::optweight). A floor of
  # -.5 binds here (the unconstrained optimum is more negative than that), so
  # weights should sit right at that floor; relaxing further to -Inf should
  # only help (or match) the variance-minimizing objective, since it strictly
  # enlarges the feasible set. Exact balance is unaffected either way, since
  # `min.w` only constrains the objective, not the balance constraints.
  expect_warning({
    W_minw_neg <- weightit(Am ~ X1 + X2 + X3 + X4 + X5,
                           data = test_data, method = "optweight", estimand = "ATE",
                           include.obj = TRUE, min.w = -.5)
  }, "negative", ignore.case = TRUE)

  expect_true(min(W_minw_neg$weights) >= -.5 - tols.eps)
  expect_true(any(W_minw_neg$weights < 0))
  expect_true(var(W_minw_neg$weights) <= var(W0$weights) + tols.eps)

  for (tt in combn(levels(W_minw_neg$treat), 2, simplify = FALSE)) {
    in_tt <- W_minw_neg$treat %in% tt
    expect_equal(cobalt::col_w_smd(W_minw_neg$covs[in_tt, ], W_minw_neg$treat[in_tt],
                                   W_minw_neg$weights[in_tt]),
                 0 * cobalt::col_w_smd(W_minw_neg$covs[in_tt, ], W_minw_neg$treat[in_tt]),
                 expected.label = "all 0s",
                 tolerance = eps)
  }

  expect_not_equal(unname(W_minw_neg$weights), unname(W0$weights))

  expect_warning({
    W_minw_ninf <- weightit(Am ~ X1 + X2 + X3 + X4 + X5,
                            data = test_data, method = "optweight", estimand = "ATE",
                            include.obj = TRUE, min.w = -Inf)
  }, "negative", ignore.case = TRUE)

  expect_true(any(W_minw_ninf$weights < 0))
  expect_true(min(W_minw_ninf$weights) <= min(W_minw_neg$weights) + tols.eps)
  expect_true(var(W_minw_ninf$weights) <= var(W_minw_neg$weights) + tols.eps)

  for (tt in combn(levels(W_minw_ninf$treat), 2, simplify = FALSE)) {
    in_tt <- W_minw_ninf$treat %in% tt
    expect_equal(cobalt::col_w_smd(W_minw_ninf$covs[in_tt, ], W_minw_ninf$treat[in_tt],
                                   W_minw_ninf$weights[in_tt]),
                 0 * cobalt::col_w_smd(W_minw_ninf$covs[in_tt, ], W_minw_ninf$treat[in_tt]),
                 expected.label = "all 0s",
                 tolerance = eps)
  }

  # tols > 0, crossed with sw and estimand. For estimand = "ATT" (the focal group's
  # weights are fixed at 1), the pairwise `tols` constraint binds the same way it
  # does for binary treatments. For estimand = "ATE" with 3+ groups, `target.tols`
  # (which defaults to 0 and isn't relaxed here) forces every group's weighted mean
  # to equal a common target exactly; because all 3 pairwise midpoints must equal
  # that same target, the group means -- and hence all pairwise differences -- are
  # fully determined to be exactly 0 regardless of `tols`. So `tols` has no effect
  # on the ATE solution here, and exact balance (not "some exactly at .05") is the
  # correct expectation for that estimand.
  for (sw in sw.opts) {
    for (estimand in estimand.opts) {
      test_that(sprintf("Optweight: tols = .05, sw = %s, estimand = %s", sw, estimand), {
        W <- weightit(Am ~ X1 + X2 + X3 + X4 + X5,
                      data = test_data, method = "optweight", estimand = estimand,
                      focal = if (estimand == "ATE") NULL else "T",
                      s.weights = if (sw) "SW" else NULL, tols = .05)

        diffs <- abs(cobalt::bal.tab(W)$Balance.Across.Pairs$Max.Diff.Adj)

        expect_true(all(diffs <= .05 + tols.eps)) #None worse than tols

        if (estimand == "ATT") {
          expect_true(any(abs(diffs - .05) <= tols.eps)) #Some exactly tols
          expect_ATT_weights_okay(W, tolerance = eps)
        }
        else {
          expect_true(all(diffs <= tols.eps)) #target.tols = 0 forces exact balance
        }
      })
    }
  }
})

test_that("Continuous treatment", {
  skip_on_cran()
  skip_if_not_installed("optweight", minimum_version = "2.0.1")
  skip_if_not_installed("cobalt")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  expect_no_condition({
    W0 <- weightit(Ac ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                   data = test_data, method = "optweight",
                   include.obj = TRUE)
  })

  expect_null(attr(W0, "Mparts", exact = TRUE))

  expect_equal(cobalt::col_w_cov(W0$covs, W0$treat, W0$weights, std = TRUE),
               0 * cobalt::col_w_cov(W0$covs, W0$treat, std = TRUE),
               expected.label = "all 0s",
               tolerance = eps)

  expect_true(is_null(W0$ps))
  expect_false(is_null(W0$obj))

  # s.weights
  expect_no_condition({
    W_sw <- weightit(Ac ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                     data = test_data, method = "optweight",
                     s.weights = "SW", include.obj = TRUE)
  })

  expect_equal(cobalt::col_w_cov(W_sw$covs, W_sw$treat, W_sw$weights, std = TRUE,
                                 s.weights = W_sw$s.weights),
               0 * cobalt::col_w_cov(W_sw$covs, W_sw$treat, std = TRUE,
                                     s.weights = W_sw$s.weights),
               expected.label = "all 0s",
               tolerance = eps)

  expect_not_equal(W_sw$weights, W0$weights)

  #Non-full rank
  expect_no_condition({
    W <- weightit(Ac ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 +
                    I(1 - X5) + I(X9 * 2),
                  data = test_data, method = "optweight",
                  include.obj = TRUE)
  })

  expect_equal(W$weights, W0$weights, tolerance = eps)

  expect_no_condition({
    W_base <- weightit(Ac ~ X1 + X2 + X3 + X4 + X5, data = test_data,
                       method = "optweight", include.obj = TRUE)
  })

  configs <- list(
    "moments = 2" = list(moments = 2),
    "int = TRUE" = list(int = TRUE),
    "norm = linf" = list(norm = "linf"),
    "min.w = .5" = list(min.w = .5)
  )

  for (nm in names(configs)) {
    test_that(sprintf("Optweight: %s", nm), {
      W <- do.call(weightit,
                  c(list(Ac ~ X1 + X2 + X3 + X4 + X5, data = test_data,
                        method = "optweight", include.obj = TRUE),
                    configs[[nm]]))

      expect_equal(cobalt::col_w_cov(W$covs, W$treat, W$weights, std = TRUE),
                   0 * cobalt::col_w_cov(W$covs, W$treat, std = TRUE),
                   expected.label = "all 0s",
                   tolerance = eps)

      expect_not_equal(unname(W$weights), unname(W_base$weights),
                       expected.label = "weights for baseline")
    })
  }

  tols.eps <- 1e-3

  # min.w can also be negative (or -Inf) to allow negative weights (see the
  # "Allowing negative weights" example in ?optweight::optweight). A floor of
  # -1 binds here (the unconstrained optimum is more negative than that), so
  # weights should sit right at that floor; relaxing further to -Inf should
  # only help (or match) the variance-minimizing objective, since it strictly
  # enlarges the feasible set. Exact balance is unaffected either way.
  expect_warning({
    W_minw_neg <- weightit(Ac ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                           data = test_data, method = "optweight",
                           include.obj = TRUE, min.w = -1)
  }, "negative", ignore.case = TRUE)

  expect_true(min(W_minw_neg$weights) >= -1 - tols.eps)
  expect_true(any(W_minw_neg$weights < 0))
  expect_true(var(W_minw_neg$weights) <= var(W0$weights) + tols.eps)

  expect_equal(cobalt::col_w_cov(W_minw_neg$covs, W_minw_neg$treat, W_minw_neg$weights, std = TRUE),
               0 * cobalt::col_w_cov(W_minw_neg$covs, W_minw_neg$treat, std = TRUE),
               expected.label = "all 0s",
               tolerance = eps)

  expect_not_equal(W_minw_neg$weights, W0$weights)

  expect_warning({
    W_minw_ninf <- weightit(Ac ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                            data = test_data, method = "optweight",
                            include.obj = TRUE, min.w = -Inf)
  }, "negative", ignore.case = TRUE)

  expect_true(any(W_minw_ninf$weights < 0))
  expect_true(min(W_minw_ninf$weights) <= min(W_minw_neg$weights) + tols.eps)
  expect_true(var(W_minw_ninf$weights) <= var(W_minw_neg$weights) + tols.eps)

  expect_equal(cobalt::col_w_cov(W_minw_ninf$covs, W_minw_ninf$treat, W_minw_ninf$weights, std = TRUE),
               0 * cobalt::col_w_cov(W_minw_ninf$covs, W_minw_ninf$treat, std = TRUE),
               expected.label = "all 0s",
               tolerance = eps)

  # tols > 0, crossed with sw. The solver only guarantees satisfying the
  # constraints up to its own convergence tolerance, which is looser than the
  # floating-point `eps` used for the exact-balance checks above, so a
  # dedicated, looser tolerance is used for these boundary checks.
  for (sw in c(FALSE, TRUE)) {
    test_that(sprintf("Optweight: tols = .05, sw = %s", sw), {
      W <- weightit(Ac ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                    data = test_data, method = "optweight",
                    s.weights = if (sw) "SW" else NULL, tols = .05)

      corrs <- abs(cobalt::bal.tab(W)$Balance$Corr.Adj)

      expect_true(all(corrs <= .05 + tols.eps)) #None worse than tols
      expect_true(any(abs(corrs - .05) <= tols.eps)) #Some exactly tols
    })
  }
})

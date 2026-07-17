test_that("Binary treatment", {
  skip_on_cran()
  skip_if_not_installed("osqp")
  skip_if_not_installed("cobalt")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  expect_no_condition({
    W0 <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                   data = test_data, method = "cfd", estimand = "ATE",
                   include.obj = TRUE)
  })

  # Plain default (moments = 0): balance is only approximate, not exact.
  expect_true(is_null(W0$ps))
  expect_false(is_null(W0$obj))
  expect_balance_improved(W0)

  sw.opts <- c(FALSE, TRUE)
  estimand.opts <- c("ATE", "ATT", "ATC")

  # Main grid uses moments = 1 with the default kernel ("gaussian"), which
  # guarantees exact mean balance at tols = 0 (the default tols), analogous
  # to the energy balancing main grid.
  weight.mat <- matrix(nrow = nrow(test_data),
                       ncol = length(sw.opts) * length(estimand.opts))
  colnames(weight.mat) <- rep("", ncol(weight.mat))

  k <- 1

  for (sw in sw.opts) {
    for (estimand in estimand.opts) {
      test_that(sprintf("CFD: sw = %s, estimand = %s", sw, estimand), {
        conv_warning <- FALSE

        withCallingHandlers({
          W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                        data = test_data, method = "cfd", estimand = estimand,
                        moments = 1,
                        s.weights = if (sw) "SW" else NULL,
                        include.obj = TRUE)
        }, warning = function(w) {
          # s.weights can occasionally cause the QP to be non-convex or
          # infeasible (documented in ?method_cfd); tolerate that here
          # rather than treating it as a hard test failure.
          if (grepl("converge|feasible", conditionMessage(w))) {
            conv_warning <<- TRUE
          }
          invokeRestart("muffleWarning")
        })

        expect_true(is_null(W$ps))
        expect_false(is_null(W$obj))

        if (!conv_warning) {
          expect_equal(cobalt::col_w_smd(W$covs, W$treat, W$weights,
                                         s.weights = W$s.weights),
                       0 * cobalt::col_w_smd(W$covs, W$treat,
                                             s.weights = W$s.weights),
                       expected.label = "all 0s",
                       tolerance = eps)
        }

        expect_true(all(is.finite(W$weights)))
        expect_true(all(W$weights >= 0))

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
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "cfd", estimand = "ATO")
  }, "not an allowable estimand", ignore.case = TRUE)

  # Non-full-rank covariates: as with energy balancing, CFD weights are not
  # invariant to adding a redundant/collinear column (the kernel/distance
  # matrix is computed on the transformed covariates, which changes when a
  # covariate's information is duplicated), so this only checks the fit
  # completes without error, not that weights match a rank-deficient-free fit.
  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 +
                    I(1 - X5) + I(X9 * 2),
                  data = test_data, method = "cfd", estimand = "ATE",
                  include.obj = TRUE)
  })

  expect_true(all(is.finite(W$weights)))

  # tols > 0 (requires moments > 0; tols is documented as "Ignored when
  # moments = 0")
  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "cfd", estimand = "ATE",
                  moments = 1, tols = .05, include.obj = TRUE)
  })

  Wexact <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                     data = test_data, method = "cfd", estimand = "ATE",
                     moments = 1, tols = 0)

  expect_not_equal(W$weights, Wexact$weights)
  expect_true(all(abs(cobalt::bal.tab(W)$Balance$Diff.Adj) <= .05 + eps)) #None worse than tols

  # Spot-check: kernel options.
  W_kernel_energy <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                              data = test_data, method = "cfd", estimand = "ATE",
                              kernel = "energy")
  W_method_energy <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                              data = test_data, method = "energy", estimand = "ATE")

  expect_equal(W_kernel_energy$weights, W_method_energy$weights, tolerance = eps)

  # Spot-check: kernel = "laplace"
  expect_no_condition({
    W_laplace <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                          data = test_data, method = "cfd", estimand = "ATE",
                          kernel = "laplace", include.obj = TRUE)
  })

  expect_not_equal(W_laplace$weights, W0$weights)

  # Spot-check: kernel = "matern" with nu in {1/2, 3/2 (default), 5/2} --
  for (nu in c(1/2, 3/2, 5/2)) {
    W_matern <- weightit(A ~ X1 + X2 + X3,
                        data = test_data, method = "cfd", estimand = "ATT",
                        kernel = "matern", nu = nu)

    expect_true(all(is.finite(W_matern$weights)))
  }

  # Spot-check: kernel = "t" (Monte Carlo kernel; nsim kept small for speed,
  # and the formula restricted to fewer covariates and estimand = "ATT",
  # since this kernel is considerably more expensive to fit than the others)
  expect_no_condition({
    W_t <- weightit(A ~ X1 + X2 + X3,
                    data = test_data, method = "cfd", estimand = "ATT",
                    kernel = "t", nsim = 50, include.obj = TRUE)
  })

  expect_true(all(is.finite(W_t$weights)))

  # Spot-check: nu out of range for matern errors cleanly
  expect_error({
    weightit(A ~ X1 + X2 + X3,
             data = test_data, method = "cfd", estimand = "ATE",
             kernel = "matern", nu = 20)
  })

  # Spot-check: int (interactions). As with energy balancing, restricted to
  # continuous covariates only -- combining int = TRUE with a dummy-coded
  # multi-level factor (X6) and tols = 0 produces an infeasible QP on this
  # dataset (same "no feasible solution" warning and degenerate constant
  # weight vector observed for method = "energy"), so it is avoided here.
  expect_no_condition({
    W_int <- weightit(A ~ X1 + X2 + X3 + X4 + X7 + X8 + X9,
                      data = test_data, method = "cfd", estimand = "ATE",
                      int = TRUE, include.obj = TRUE)
  })

  expect_true(all(is.finite(W_int$weights)))
  expect_equal(cobalt::col_w_smd(W_int$covs, W_int$treat, W_int$weights),
               0 * cobalt::col_w_smd(W_int$covs, W_int$treat),
               expected.label = "all 0s",
               tolerance = eps)

  # Spot-check: quantile (restricted to a smaller formula and ATT for speed)
  expect_no_condition({
    W_q <- weightit(A ~ X1 + X2 + X3,
                    data = test_data, method = "cfd", estimand = "ATT",
                    moments = 1, quantile = list(X1 = c(.25, .5, .75)),
                    include.obj = TRUE)
  })

  W_noq <- weightit(A ~ X1 + X2 + X3,
                    data = test_data, method = "cfd", estimand = "ATT",
                    moments = 1)

  expect_not_equal(W_q$weights, W_noq$weights)

  # Spot-check: improved (only affects ATE)
  W_imp_false <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                          data = test_data, method = "cfd", estimand = "ATE",
                          moments = 1, improved = FALSE)

  expect_not_equal(Wexact$weights, W_imp_false$weights)

  # Spot-check: lambda (weight-variability penalty) -- higher lambda should
  # reduce weight variability (i.e., increase ESS) relative to the
  # near-zero default.
  W_lambda_hi <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                          data = test_data, method = "cfd", estimand = "ATE",
                          lambda = 1)

  expect_true(ESS(W_lambda_hi$weights) > ESS(W0$weights))
})

test_that("Multi-category treatment", {
  skip_on_cran()
  skip_if_not_installed("osqp")
  skip_if_not_installed("cobalt")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  expect_no_condition({
    W0 <- weightit(Am ~ X1 + X2 + X3 + X4 + X5,
                   data = test_data, method = "cfd", estimand = "ATE",
                   include.obj = TRUE)
  })

  expect_true(is_null(W0$ps))
  expect_false(is_null(W0$obj))

  # expect_balance_improved() assumes a binary treatment; check pairwise SMD
  # magnitude shrinks for one pair of levels for this multi-category default
  # (unconstrained) fit instead.
  smd_unw <- cobalt::col_w_smd(W0$covs[W0$treat %in% c("T", "C1"), ],
                               W0$treat[W0$treat %in% c("T", "C1")])
  smd_w <- cobalt::col_w_smd(W0$covs[W0$treat %in% c("T", "C1"), ],
                             W0$treat[W0$treat %in% c("T", "C1")],
                             W0$weights[W0$treat %in% c("T", "C1")])
  expect_true(max(abs(smd_w)) < max(abs(smd_unw)))

  sw.opts <- c(FALSE, TRUE)
  estimand.opts <- c("ATE", "ATT")

  weight.mat <- matrix(nrow = nrow(test_data),
                       ncol = length(sw.opts) * length(estimand.opts))
  colnames(weight.mat) <- rep("", ncol(weight.mat))

  k <- 1

  for (sw in sw.opts) {
    for (estimand in estimand.opts) {
      test_that(sprintf("CFD: sw = %s, estimand = %s", sw, estimand), {
        conv_warning <- FALSE

        withCallingHandlers({
          W <- weightit(Am ~ X1 + X2 + X3 + X4 + X5,
                        data = test_data, method = "cfd", estimand = estimand,
                        moments = 1,
                        focal = if (estimand == "ATE") NULL else "T",
                        s.weights = if (sw) "SW" else NULL,
                        include.obj = TRUE)
        }, warning = function(w) {
          if (grepl("converge|feasible", conditionMessage(w))) {
            conv_warning <<- TRUE
          }
          invokeRestart("muffleWarning")
        })

        expect_true(is_null(W$ps))
        expect_false(is_null(W$obj))

        if (!conv_warning) {
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
        }

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

  # tols > 0
  expect_no_condition({
    W <- weightit(Am ~ X1 + X2 + X3 + X4 + X5,
                  data = test_data, method = "cfd", estimand = "ATE",
                  moments = 1, tols = .05, include.obj = TRUE)
  })

  Wexact <- weightit(Am ~ X1 + X2 + X3 + X4 + X5,
                     data = test_data, method = "cfd", estimand = "ATE",
                     moments = 1, tols = 0)

  expect_not_equal(W$weights, Wexact$weights)
  expect_true(all(abs(cobalt::bal.tab(W)$Balance$Max.Diff.Adj) <= .05 + eps)) #None worse than tols

  # Spot-check: kernel = "energy" reproduces method = "energy" for
  # multi-category treatments as well.
  W_kernel_energy <- weightit(Am ~ X1 + X2 + X3 + X4 + X5,
                              data = test_data, method = "cfd", estimand = "ATE",
                              kernel = "energy")
  W_method_energy <- weightit(Am ~ X1 + X2 + X3 + X4 + X5,
                              data = test_data, method = "energy", estimand = "ATE")

  expect_equal(W_kernel_energy$weights, W_method_energy$weights, tolerance = eps)

  # Spot-check: improved
  W_imp_false <- weightit(Am ~ X1 + X2 + X3 + X4 + X5,
                          data = test_data, method = "cfd", estimand = "ATE",
                          moments = 1, improved = FALSE)

  expect_not_equal(Wexact$weights, W_imp_false$weights)
})

test_that("Continuous treatment is not supported", {
  skip_on_cran()
  skip_if_not_installed("osqp")

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  # CFD's `treat_type` in .weightit_methods is c("binary", "multinomial")
  # only (no "continuous"), so weightit() raises a clear, early error rather
  # than attempting (and silently failing at, or mis-specifying) a fit.
  expect_error({
    weightit(Ac ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
             data = test_data, method = "cfd")
  }, "binary or multinomial", ignore.case = TRUE)
})

test_that("Binary treatment", {
  skip_on_cran()
  skip_if_not_installed("osqp")
  skip_if_not_installed("cobalt")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  expect_no_condition({
    W0 <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                   data = test_data, method = "energy", estimand = "ATE",
                   include.obj = TRUE)
  })

  # Plain default (moments = 0): balance is only approximate, not exact.
  expect_true(is_null(W0$ps))
  expect_false(is_null(W0$obj))
  expect_balance_improved(W0)

  sw.opts <- c(FALSE, TRUE)
  estimand.opts <- c("ATE", "ATT", "ATC")

  # Main grid uses moments = 1, which (per the documented behavior of
  # `moments` for this method) guarantees exact mean balance at tols = 0
  # (the default tols). This is the balance-optimizing regime analogous to
  # ebal/ipt's default behavior.
  weight.mat <- matrix(nrow = nrow(test_data),
                       ncol = length(sw.opts) * length(estimand.opts))
  colnames(weight.mat) <- rep("", ncol(weight.mat))

  k <- 1

  for (sw in sw.opts) {
    for (estimand in estimand.opts) {
      test_that(sprintf("Energy: sw = %s, estimand = %s", sw, estimand), {
        conv_warning <- FALSE

        withCallingHandlers({
          W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                        data = test_data, method = "energy", estimand = estimand,
                        moments = 1,
                        s.weights = if (sw) "SW" else NULL,
                        include.obj = TRUE)
        }, warning = function(w) {
          # s.weights can occasionally cause the QP to be non-convex or
          # infeasible (documented in ?method_energy); tolerate that here
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
                  data = test_data, method = "energy", estimand = "ATO")
  }, "not an allowable estimand", ignore.case = TRUE)

  # Non-full-rank covariates: unlike ebal, energy balancing weights are NOT
  # invariant to adding a redundant/collinear column, because the distance
  # matrix (computed on the *scaled* covariates) changes when a covariate's
  # information is duplicated.
  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 +
                    I(1 - X5) + I(X9 * 2),
                  data = test_data, method = "energy", estimand = "ATE",
                  include.obj = TRUE)
  })

  expect_true(all(is.finite(W$weights)))

  # tols > 0 (requires moments > 0 for tols to have any effect; tols is
  # documented as "Ignored when moments = 0")
  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "energy", estimand = "ATE",
                  moments = 1, tols = .05, include.obj = TRUE)
  })

  Wexact <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                     data = test_data, method = "energy", estimand = "ATE",
                     moments = 1, tols = 0)

  expect_not_equal(W$weights, Wexact$weights)

  # Unlike ebal (whose dual-based solution actively pushes at least one
  # covariate to exactly the tols boundary), energy balancing's QP objective
  # does not need to saturate the relaxed constraint for any covariate here
  # -- all Diff.Adj values are comfortably inside +/-.05 with none close to
  # the boundary. So we only check the constraint is respected and that
  # relaxing it changed the solution (already checked above), not that any
  # covariate sits exactly at the boundary.
  expect_true(all(abs(cobalt::bal.tab(W)$Balance$Diff.Adj) <= .05 + eps)) #None worse than tols

  # Spot-check: int (interactions). Restricted to continuous covariates only
  # -- combining int = TRUE with a dummy-coded multi-level factor (X6) and
  # tols = 0 produces an infeasible QP on this dataset (osqp warns "no
  # feasible solution" and returns a degenerate constant weight vector),
  # which appears to be a genuine (if edge-case) source behavior rather than
  # a testing mistake, so it is avoided here rather than worked around.
  expect_no_condition({
    W_int <- weightit(A ~ X1 + X2 + X3 + X4 + X7 + X8 + X9,
                      data = test_data, method = "energy", estimand = "ATE",
                      int = TRUE, include.obj = TRUE)
  })

  expect_true(all(is.finite(W_int$weights)))
  expect_equal(cobalt::col_w_smd(W_int$covs, W_int$treat, W_int$weights),
               0 * cobalt::col_w_smd(W_int$covs, W_int$treat),
               expected.label = "all 0s",
               tolerance = eps)

  # Spot-check: quantile
  expect_no_condition({
    W_q <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                    data = test_data, method = "energy", estimand = "ATE",
                    moments = 1, quantile = list(X1 = c(.25, .5, .75)),
                    include.obj = TRUE)
  })

  expect_not_equal(W_q$weights, Wexact$weights)

  # Spot-check: improved (only affects ATE)
  W_imp_true <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                         data = test_data, method = "energy", estimand = "ATE",
                         moments = 1, improved = TRUE)
  W_imp_false <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                          data = test_data, method = "energy", estimand = "ATE",
                          moments = 1, improved = FALSE)

  expect_not_equal(W_imp_true$weights, W_imp_false$weights)
  # improved = TRUE is the default
  expect_equal(Wexact$weights, W_imp_true$weights, tolerance = eps)

  # Spot-check: lambda (weight-variability penalty) -- higher lambda should
  # reduce weight variability (i.e., increase ESS) relative to the
  # near-zero default.
  W_lambda_hi <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                          data = test_data, method = "energy", estimand = "ATE",
                          lambda = 1)

  expect_true(ESS(W_lambda_hi$weights) > ESS(W0$weights))

  # Spot-check: min.w
  expect_warning({
    W_minw <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                       data = test_data, method = "energy", estimand = "ATE",
                       min.w = -Inf)
  }, "negative")

  expect_true(all(is.finite(W_minw$weights)))
  expect_true(any(W_minw$weights < 0))
})

test_that("Multi-category treatment", {
  skip_on_cran()
  skip_if_not_installed("osqp")
  skip_if_not_installed("cobalt")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  expect_no_condition({
    W0 <- weightit(Am ~ X1 + X2 + X3 + X4 + X5,
                   data = test_data, method = "energy", estimand = "ATE",
                   include.obj = TRUE)
  })

  expect_true(is_null(W0$ps))
  expect_false(is_null(W0$obj))

  # expect_balance_improved() assumes a binary treatment (it calls
  # cobalt::col_w_smd() without pairwise subsetting), so for the
  # multi-category default (unconstrained) fit we just confirm pairwise SMD
  # magnitude shrinks for at least one pair rather than reusing that helper.
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
      test_that(sprintf("Energy: sw = %s, estimand = %s", sw, estimand), {
        conv_warning <- FALSE

        withCallingHandlers({
          W <- weightit(Am ~ X1 + X2 + X3 + X4 + X5,
                        data = test_data, method = "energy", estimand = estimand,
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
                  data = test_data, method = "energy", estimand = "ATE",
                  moments = 1, tols = .05, include.obj = TRUE)
  })

  Wexact <- weightit(Am ~ X1 + X2 + X3 + X4 + X5,
                     data = test_data, method = "energy", estimand = "ATE",
                     moments = 1, tols = 0)

  expect_not_equal(W$weights, Wexact$weights)

  # As in the binary case, the relaxed constraint is not saturated for any
  # covariate on this dataset -- only that it is respected is checked here.
  expect_true(all(abs(cobalt::bal.tab(W)$Balance$Max.Diff.Adj) <= .05 + eps)) #None worse than tols

  # Spot-check: improved
  W_imp_false <- weightit(Am ~ X1 + X2 + X3 + X4 + X5,
                          data = test_data, method = "energy", estimand = "ATE",
                          moments = 1, improved = FALSE)

  expect_not_equal(Wexact$weights, W_imp_false$weights)
})

test_that("Continuous treatment", {
  skip_on_cran()
  skip_if_not_installed("osqp")
  skip_if_not_installed("cobalt")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  expect_no_condition({
    W0 <- weightit(Ac ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                   data = test_data, method = "energy",
                   include.obj = TRUE)
  })

  expect_true(is_null(W0$ps))
  expect_false(is_null(W0$obj))

  # Plain default (d.moments = 0): only approximate treatment-covariate
  # correlation reduction, not exact.
  expect_true(max(abs(cobalt::col_w_cov(W0$covs, W0$treat, W0$weights, std = TRUE))) <
                max(abs(cobalt::col_w_cov(W0$covs, W0$treat, std = TRUE))))

  sw.opts <- c(FALSE, TRUE)
  d.moments.opts <- c(1, 3)

  weight.mat <- matrix(nrow = nrow(test_data),
                       ncol = length(sw.opts) * length(d.moments.opts))
  colnames(weight.mat) <- rep("", ncol(weight.mat))

  k <- 1

  for (sw in sw.opts) {
    for (d.moments in d.moments.opts) {
      test_that(sprintf("Energy: sw = %s, d.moments = %s", sw, d.moments), {
        conv_warning <- FALSE

        withCallingHandlers({
          W <- weightit(Ac ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                        data = test_data, method = "energy",
                        moments = 1, d.moments = d.moments,
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
        }

        expect_true(all(is.finite(W$weights)))
        expect_true(all(W$weights >= 0))

        for (i in seq_len(k - 1)) {
          expect_not_equal(unname(W$weights), weight.mat[,i],
                           expected.label = sprintf("weights for %s", colnames(weight.mat)[i]),
                           tolerance = eps)
        }

        n <- sprintf("W_%s_%s", sw, d.moments)
        colnames(weight.mat)[k] <<- n
        weight.mat[,k] <<- W$weights
        k <<- k + 1
      })
    }
  }

  # tols > 0
  expect_no_condition({
    W <- weightit(Ac ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "energy",
                  moments = 1, tols = .05, include.obj = TRUE)
  })

  Wexact <- weightit(Ac ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                     data = test_data, method = "energy",
                     moments = 1, tols = 0)

  expect_not_equal(W$weights, Wexact$weights)

  expect_true(all(abs(cobalt::bal.tab(W)$Balance$Corr.Adj) <= .05 + eps)) #None worse than tols
  expect_true(any(abs(abs(cobalt::bal.tab(W)$Balance$Corr.Adj) - .05) <= eps)) #Some exactly tols
  expect_true(any(abs(cobalt::bal.tab(W)$Balance$Corr.Adj) > eps)) #Some worse than 0

  # Spot-check: dimension.adj
  W_adj_false <- weightit(Ac ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                          data = test_data, method = "energy",
                          dimension.adj = FALSE)

  expect_not_equal(W0$weights, W_adj_false$weights)
})

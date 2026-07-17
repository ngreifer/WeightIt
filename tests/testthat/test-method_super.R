test_that("Binary treatment", {
  skip_on_cran()
  skip_if_not_installed("SuperLearner")
  skip_if_not_installed("cobalt")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  # Small/fast library throughout
  SL.lib <- c("SL.mean", "SL.glm", "SL.step.interaction")

  # `expect_no_error()` (not `expect_no_condition()`) because the very first
  # SuperLearner fit in a session emits a `packageStartupMessage` when it
  # lazily loads `nnls` (used by the default `SL.method = "method.NNLS"`).
  set.seed(123)
  suppressPackageStartupMessages({
    W0 <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                   data = test_data, method = "super", estimand = "ATE",
                   SL.library = SL.lib, include.obj = TRUE)
  })

  expect_true(is.numeric(W0$ps))
  expect_true(all(W0$ps > 0 & W0$ps < 1))
  expect_false(is_null(W0$obj))

  # SuperLearner-based weighting does not support M-estimation
  expect_null(attr(W0, "Mparts", exact = TRUE))

  sw.opts <- c(FALSE, TRUE)
  estimand.opts <- c("ATE", "ATT", "ATC", "ATO", "ATM", "ATOS")

  weight.mat <- matrix(nrow = nrow(test_data),
                       ncol = length(sw.opts) * length(estimand.opts))
  colnames(weight.mat) <- rep("", ncol(weight.mat))

  k <- 1

  for (sw in sw.opts) {
    for (estimand in estimand.opts) {
      test_that(sprintf("Super: sw = %s, estimand = %s", sw, estimand), {
        set.seed(123)
        W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                      data = test_data, method = "super", estimand = estimand,
                      s.weights = if (sw) "SW" else NULL,
                      SL.library = SL.lib,
                      include.obj = TRUE)

        expect_true(is.numeric(W$ps))
        expect_true(all(W$ps > 0 & W$ps < 1))
        expect_false(is_null(W$obj))
        expect_true(all(is.finite(W$weights) & W$weights >= 0))

        # SuperLearner is a machine-learning PS method; it approximates but
        # does not solve exactly for balance, so we check improvement over
        # the unweighted sample rather than exact-zero SMDs.
        expect_balance_improved(W)

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

  # `SL.library` is required (no default)
  expect_error({
    weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
             data = test_data, method = "super", estimand = "ATE")
  }, "SL.library", fixed = TRUE)

  # Estimands
  expect_error({
    weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
             data = test_data, method = "super", estimand = "XYZ",
             SL.library = SL.lib)
  }, "not an allowable estimand", ignore.case = TRUE)

  # Non-full rank
  set.seed(123)
  expect_no_error({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 +
                    I(1 - X5) + I(X9 * 2),
                  data = test_data, method = "super", estimand = "ATE",
                  SL.library = SL.lib)
  })

  expect_true(is.numeric(W$ps))
  expect_true(all(W$ps > 0 & W$ps < 1))

  # `discrete = TRUE` selects the single best-performing learner rather than
  # an ensemble; spot-check rather than cross into the main grid.
  set.seed(123)
  expect_no_error({
    W.discrete <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                           data = test_data, method = "super", estimand = "ATE",
                           SL.library = SL.lib, discrete = TRUE)
  })

  expect_true(is.numeric(W.discrete$ps))
  expect_true(all(W.discrete$ps > 0 & W.discrete$ps < 1))
  expect_balance_improved(W.discrete)
  expect_not_equal(W.discrete$weights, W0$weights)

  # `SL.method = "method.balance"` uses Balance SuperLearner (Pirracchio &
  # Carone, 2018) with a `criterion` for choosing predictions
  set.seed(123)
    expect_no_error({
      W.bal <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                        data = test_data, method = "super", estimand = "ATT",
                        SL.library = SL.lib, SL.method = "method.balance",
                        criterion = "smd.mean")
    })

  expect_true(is.numeric(W.bal$ps))
  expect_true(all(W.bal$ps > 0 & W.bal$ps < 1))
  expect_ATT_weights_okay(W.bal, tolerance = eps)

  # `subclass` computes MMWS weights instead of standard IPW; restricted to
  # ATE/ATT/ATC. Spot-check via the sharp drop in the number of unique
  # weights rather than relying on a `"subclass"` attribute (which is not
  # retained on the final `weightit` object's weights).
  set.seed(123)
  expect_no_error({
    W.sub <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                      data = test_data, method = "super", estimand = "ATE",
                      SL.library = SL.lib, subclass = 10)
  })

  expect_true(length(unique(W.sub$weights)) < nrow(test_data) / 10)
  expect_not_equal(W.sub$weights, W0$weights)

  # `subclass` is incompatible with estimands other than ATE/ATT/ATC
  expect_error({
    weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
             data = test_data, method = "super", estimand = "ATO",
             SL.library = SL.lib, subclass = 10)
  }, "ATE, ATT, and ATC", fixed = TRUE)

  # Missing data (`missing = "ind"` is the only allowed value)
  data_na <- test_data
  set.seed(4321)
  is.na(data_na$X1) <- sample(nrow(data_na), round(.05 * nrow(data_na)))

  set.seed(123)
  expect_no_condition({
    W.na <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                     data = data_na, method = "super", estimand = "ATE",
                     missing = "ind",
                     SL.library = SL.lib)
  })

  expect_true(anyNA(W.na$covs$X1))
  expect_true(is.numeric(W.na$ps))
  expect_true(all(W.na$ps[!is.na(W.na$ps)] > 0 & W.na$ps[!is.na(W.na$ps)] < 1))

  expect_error({
    weightit(A ~ X1 + X2 + X3, data = data_na, method = "super",
             missing = "surr", SL.library = SL.lib)
  }, "only.*allowed for.*missing", ignore.case = TRUE)
})

test_that("Multi-category treatment", {
  skip_on_cran()
  skip_if_not_installed("SuperLearner")
  skip_if_not_installed("cobalt")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  SL.lib <- c("SL.mean", "SL.glm", "SL.step.interaction")

  set.seed(123)
  suppressPackageStartupMessages({
    W0 <- weightit(Am ~ X1 + X2 + X3 + X4 + X5,
                   data = test_data, method = "super", estimand = "ATE",
                   SL.library = SL.lib, include.obj = TRUE)
  })

  expect_true(is_null(W0$ps)) #ps not returned for multi-category super
  expect_false(is_null(W0$obj))

  sw.opts <- c(FALSE, TRUE)
  estimand.opts <- c("ATE", "ATT", "ATO", "ATM")

  weight.mat <- matrix(nrow = nrow(test_data),
                       ncol = length(sw.opts) * length(estimand.opts))
  colnames(weight.mat) <- rep("", ncol(weight.mat))

  k <- 1

  for (sw in sw.opts) {
    for (estimand in estimand.opts) {
      test_that(sprintf("Super: sw = %s, estimand = %s", sw, estimand), {
        set.seed(123)
        # `suppressWarnings()`: see note above on "non-integer #successes".
        suppressWarnings({
          W <- weightit(Am ~ X1 + X2 + X3 + X4 + X5,
                        data = test_data, method = "super", estimand = estimand,
                        focal = if (estimand == "ATT") "T" else NULL,
                        s.weights = if (sw) "SW" else NULL,
                        SL.library = SL.lib,
                        include.obj = TRUE)
        })

        expect_true(is_null(W$ps))
        expect_false(is_null(W$obj))
        expect_true(all(is.finite(W$weights) & W$weights >= 0))

        # Only check pairs involving "T": by construction (see
        # fixtures/make_test_data.R), "C1" and "C2" are randomly assigned
        # among untreated units and share no true relationship with the
        # covariates, so a C1-vs-C2 "balance improved" check would be
        # comparing pure noise and can fail by chance.
        for (tt in combn(levels(W$treat), 2, simplify = FALSE)) {
          if ("T" %nin% tt) next

          in_tt <- W$treat %in% tt
          W_sub <- list(covs = W$covs[in_tt, , drop = FALSE],
                        treat = factor(W$treat[in_tt]),
                        weights = W$weights[in_tt],
                        s.weights = W$s.weights[in_tt])
          expect_balance_improved(W_sub,
                                  label = sprintf("SMDs for %s", paste(tt, collapse = " vs. ")))
        }

        if (estimand == "ATT") {
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

  # Documented behavior noted above: ATC == ATT for multi-category treatments
  # with the same focal, for any method that routes through
  # `.get_w_from_ps_internal_multi()` (super, bart, glm, cbps).
  set.seed(123)
  W.att <- weightit(Am ~ X1 + X2 + X3 + X4 + X5,
                    data = test_data, method = "super", estimand = "ATT",
                    focal = "T", SL.library = SL.lib)
  set.seed(123)
  W.atc <- weightit(Am ~ X1 + X2 + X3 + X4 + X5,
                    data = test_data, method = "super", estimand = "ATC",
                    focal = "T", SL.library = SL.lib)

  expect_equal(W.att$weights, W.atc$weights)

  # `ATOS` is not an allowable estimand for multi-category treatments
  expect_error({
    weightit(Am ~ X1 + X2 + X3 + X4 + X5,
             data = test_data, method = "super", estimand = "ATOS",
             SL.library = SL.lib)
  }, "not an allowable estimand", ignore.case = TRUE)

  # `SL.method = "method.balance"` is not supported for multi-category
  expect_error({
    weightit(Am ~ X1 + X2 + X3 + X4 + X5,
             data = test_data, method = "super", estimand = "ATE",
             SL.library = SL.lib, SL.method = "method.balance")
  }, "method.balance", fixed = TRUE)
})

test_that("Continuous treatment", {
  skip_on_cran()
  skip_if_not_installed("SuperLearner")
  skip_if_not_installed("cobalt")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  SL.lib <- c("SL.mean", "SL.glm", "SL.step.interaction")

  set.seed(123)
  suppressPackageStartupMessages({
    W0 <- weightit(Ac ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                   data = test_data, method = "super",
                   SL.library = SL.lib, include.obj = TRUE)
  })

  expect_false(is_null(W0$obj))
  expect_true(all(is.finite(W0$weights) & W0$weights >= 0))

  sw.opts <- c(FALSE, TRUE)
  density.opts <- c("dnorm", "kernel")

  weight.mat <- matrix(nrow = nrow(test_data),
                       ncol = length(sw.opts) * length(density.opts))
  colnames(weight.mat) <- rep("", ncol(weight.mat))

  k <- 1

  for (sw in sw.opts) {
    for (density in density.opts) {
      test_that(sprintf("Super: sw = %s, density = %s", sw, density), {
        set.seed(123)
        # `suppressWarnings()`: kernel density estimation with `s.weights`
        # produces a benign "Selecting bandwidth *not* using 'weights'"
        # warning from `stats::density()`.
        suppressWarnings({
          W <- weightit(Ac ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                        data = test_data, method = "super",
                        density = density,
                        s.weights = if (sw) "SW" else NULL,
                        SL.library = SL.lib,
                        include.obj = TRUE)
        })

        expect_false(is_null(W$obj))
        expect_true(all(is.finite(W$weights) & W$weights >= 0))

        # ML-based GPS is approximate; check improvement in weighted
        # treatment-covariate correlation rather than exact-zero.
        weighted <- abs(cobalt::col_w_cov(W$covs, W$treat, W$weights, std = TRUE,
                                          s.weights = W$s.weights))
        unweighted <- abs(cobalt::col_w_cov(W$covs, W$treat, std = TRUE,
                                            s.weights = W$s.weights))
        expect_true(max(weighted) < max(unweighted))

        for (i in seq_len(k - 1)) {
          expect_not_equal(unname(W$weights), weight.mat[,i],
                           expected.label = sprintf("weights for %s", colnames(weight.mat)[i]),
                           tolerance = eps)
        }

        n <- sprintf("W_%s_%s", sw, density)
        colnames(weight.mat)[k] <<- n
        weight.mat[,k] <<- W$weights
        k <<- k + 1
      })
    }
  }

  # Non-full rank
  set.seed(123)
  expect_no_error({
    W <- weightit(Ac ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 +
                    I(1 - X5) + I(X9 * 2),
                  data = test_data, method = "super",
                  SL.library = SL.lib)
  })

  expect_true(all(is.finite(W$weights) & W$weights >= 0))
})

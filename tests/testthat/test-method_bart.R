test_that("Binary treatment", {
  skip_on_cran()
  skip_if_not_installed("dbarts")
  skip_if_not_installed("cobalt")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  # Small tree/posterior-draw settings throughout to keep fits fast (dbarts's
  # `bart2()` defaults are n.trees = 75, n.samples = 500, n.burn = 500,
  # n.chains = 4). `n.chains = 1, n.threads = 1` plus `set.seed()` before
  # every fit ensures reproducibility (see *Reproducibility* section of
  # `?method_bart`: multi-threaded fits are not guaranteed reproducible
  # across machines even with `seed` set, so we single-thread throughout the
  # main grid).
  n.trees <- 20
  n.samples <- 50
  n.burn <- 50

  set.seed(123)
  expect_no_condition({
    W0 <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                   data = test_data, method = "bart", estimand = "ATE",
                   n.trees = n.trees, n.samples = n.samples, n.burn = n.burn,
                   n.chains = 1, n.threads = 1,
                   include.obj = TRUE)
  })

  expect_true(is.numeric(W0$ps))
  expect_true(all(W0$ps > 0 & W0$ps < 1))
  expect_false(is_null(W0$obj))

  # BART-based weighting does not support M-estimation, unlike
  # cbps/ebal/ipt/glm -- confirm no Mparts are attached.
  expect_null(attr(W0, "Mparts", exact = TRUE))

  # `set.seed()` + `n.threads = 1` reproduces exactly
  set.seed(123)
  W0b <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "bart", estimand = "ATE",
                  n.trees = n.trees, n.samples = n.samples, n.burn = n.burn,
                  n.chains = 1, n.threads = 1)

  expect_equal(W0$weights, W0b$weights)

  # `s.weights` is documented as unsupported for BART (`s.weights_ok = FALSE`
  # in `.weightit_methods`). Confirmed behavior: `weightit()` errors *before*
  # ever reaching `weightit2bart()`, via `.check_method_s.weights()`, with a
  # clear message -- it does not silently ignore or mishandle `s.weights`.
  expect_error({
    weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
             data = test_data, method = "bart", estimand = "ATE",
             s.weights = "SW",
             n.trees = n.trees, n.samples = n.samples, n.burn = n.burn,
             n.chains = 1, n.threads = 1)
  }, "sampling weights cannot be used", ignore.case = TRUE)

  # A constant "s.weights" (i.e., no actual sampling weighting) is allowed,
  # since `.check_method_s.weights()` only errors when the weights actually
  # vary.
  set.seed(123)
  expect_no_condition({
    weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
             data = test_data, method = "bart", estimand = "ATE",
             s.weights = rep(1, nrow(test_data)),
             n.trees = n.trees, n.samples = n.samples, n.burn = n.burn,
             n.chains = 1, n.threads = 1)
  })

  estimand.opts <- c("ATE", "ATT", "ATC", "ATO", "ATM", "ATOS")

  weight.mat <- matrix(nrow = nrow(test_data), ncol = length(estimand.opts))
  colnames(weight.mat) <- rep("", ncol(weight.mat))

  k <- 1

  for (estimand in estimand.opts) {
    test_that(sprintf("BART: estimand = %s", estimand), {
      set.seed(123)
      W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                    data = test_data, method = "bart", estimand = estimand,
                    n.trees = n.trees, n.samples = n.samples, n.burn = n.burn,
                    n.chains = 1, n.threads = 1,
                    include.obj = TRUE)

      expect_true(is.numeric(W$ps))
      expect_true(all(W$ps > 0 & W$ps < 1))
      expect_false(is_null(W$obj))
      expect_true(all(is.finite(W$weights) & W$weights >= 0))

      # BART is a machine-learning PS method; it approximates but does not
      # solve exactly for balance, so we check improvement over the
      # unweighted sample rather than exact-zero SMDs.
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

      n <- sprintf("W_%s", estimand)
      colnames(weight.mat)[k] <<- n
      weight.mat[,k] <<- W$weights
      k <<- k + 1
    })
  }

  # Estimands
  expect_error({
    weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
             data = test_data, method = "bart", estimand = "XYZ",
             n.trees = n.trees, n.samples = n.samples, n.burn = n.burn)
  }, "not an allowable estimand", ignore.case = TRUE)

  # Non-full rank
  set.seed(123)
  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 +
                    I(1 - X5) + I(X9 * 2),
                  data = test_data, method = "bart", estimand = "ATE",
                  n.trees = n.trees, n.samples = n.samples, n.burn = n.burn,
                  n.chains = 1, n.threads = 1)
  })

  expect_true(is.numeric(W$ps))
  expect_true(all(W$ps > 0 & W$ps < 1))

  # `subclass` computes MMWS weights instead of standard IPW; restricted to
  # ATE/ATT/ATC. Spot-check via the sharp drop in the number of unique
  # weights.
  set.seed(123)
  expect_no_condition({
    W.sub <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                      data = test_data, method = "bart", estimand = "ATE",
                      n.trees = n.trees, n.samples = n.samples, n.burn = n.burn,
                      n.chains = 1, n.threads = 1, subclass = 10)
  })

  expect_true(length(unique(W.sub$weights)) < nrow(test_data) / 10)
  expect_not_equal(W.sub$weights, W0$weights)

  # `n.threads` > 1: spot-check that multi-threaded fitting still produces
  # valid, structurally sound weights. We don't assert exact
  # cross-run reproducibility here since the docs explicitly note that is
  # only guaranteed with `n.threads = 1`.
  set.seed(123)
  expect_no_condition({
    W.mt <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                     data = test_data, method = "bart", estimand = "ATE",
                     n.trees = n.trees, n.samples = n.samples, n.burn = n.burn,
                     n.chains = 2, n.threads = 2, seed = 1)
  })

  expect_true(is.numeric(W.mt$ps))
  expect_true(all(W.mt$ps > 0 & W.mt$ps < 1))
  expect_true(all(is.finite(W.mt$weights) & W.mt$weights >= 0))

  # `seed` argument (passed to `dbarts::bart2()`) reproduces exactly,
  # independent of `set.seed()`
  W.seed1 <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                      data = test_data, method = "bart", estimand = "ATE",
                      n.trees = n.trees, n.samples = n.samples, n.burn = n.burn,
                      n.chains = 1, n.threads = 1, seed = 99)
  W.seed2 <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                      data = test_data, method = "bart", estimand = "ATE",
                      n.trees = n.trees, n.samples = n.samples, n.burn = n.burn,
                      n.chains = 1, n.threads = 1, seed = 99)

  expect_equal(W.seed1$weights, W.seed2$weights)

  # Missing data (`missing = "ind"` is the only allowed value)
  data_na <- test_data
  set.seed(4321)
  is.na(data_na$X1) <- sample(nrow(data_na), round(.05 * nrow(data_na)))

  set.seed(123)
  expect_no_condition({
    W.na <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                     data = data_na, method = "bart", estimand = "ATE",
                     missing = "ind",
                     n.trees = n.trees, n.samples = n.samples, n.burn = n.burn,
                     n.chains = 1, n.threads = 1)
  })

  expect_true(anyNA(W.na$covs$X1))
  expect_true(is.numeric(W.na$ps))
  expect_true(all(W.na$ps > 0 & W.na$ps < 1))

  expect_error({
    weightit(A ~ X1 + X2 + X3, data = data_na, method = "bart",
             missing = "surr",
             n.trees = n.trees, n.samples = n.samples, n.burn = n.burn)
  }, "only.*allowed for.*missing", ignore.case = TRUE)
})

test_that("Multi-category treatment", {
  skip_on_cran()
  skip_if_not_installed("dbarts")
  skip_if_not_installed("cobalt")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  n.trees <- 20
  n.samples <- 50
  n.burn <- 50

  set.seed(123)
  expect_no_condition({
    W0 <- weightit(Am ~ X1 + X2 + X3 + X4 + X5,
                   data = test_data, method = "bart", estimand = "ATE",
                   n.trees = n.trees, n.samples = n.samples, n.burn = n.burn,
                   n.chains = 1, n.threads = 1,
                   include.obj = TRUE)
  })

  expect_true(is_null(W0$ps)) #ps not directly returned for multi-category bart
  expect_false(is_null(W0$obj))

  estimand.opts <- c("ATE", "ATT", "ATO", "ATM")

  weight.mat <- matrix(nrow = nrow(test_data), ncol = length(estimand.opts))
  colnames(weight.mat) <- rep("", ncol(weight.mat))

  k <- 1

  for (estimand in estimand.opts) {
    test_that(sprintf("BART: estimand = %s", estimand), {
      set.seed(123)
      W <- weightit(Am ~ X1 + X2 + X3 + X4 + X5,
                    data = test_data, method = "bart", estimand = estimand,
                    focal = if (estimand == "ATT") "T" else NULL,
                    n.trees = n.trees, n.samples = n.samples, n.burn = n.burn,
                    n.chains = 1, n.threads = 1,
                    include.obj = TRUE)

      expect_true(is_null(W$ps))
      expect_false(is_null(W$obj))
      expect_true(all(is.finite(W$weights) & W$weights >= 0))

      # Only check pairs involving "T": by construction (see
      # fixtures/make_test_data.R), "C1" and "C2" are randomly assigned
      # among untreated units and share no true relationship with the
      # covariates, so a C1-vs-C2 "balance improved" check would compare
      # pure noise and can fail by chance.
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

      n <- sprintf("W_%s", estimand)
      colnames(weight.mat)[k] <<- n
      weight.mat[,k] <<- W$weights
      k <<- k + 1
    })
  }

  # Documented behavior noted above: ATC == ATT for multi-category treatments
  # with the same focal.
  set.seed(123)
  W.att <- weightit(Am ~ X1 + X2 + X3 + X4 + X5,
                    data = test_data, method = "bart", estimand = "ATT",
                    focal = "T",
                    n.trees = n.trees, n.samples = n.samples, n.burn = n.burn,
                    n.chains = 1, n.threads = 1)
  set.seed(123)
  W.atc <- weightit(Am ~ X1 + X2 + X3 + X4 + X5,
                    data = test_data, method = "bart", estimand = "ATC",
                    focal = "T",
                    n.trees = n.trees, n.samples = n.samples, n.burn = n.burn,
                    n.chains = 1, n.threads = 1)

  expect_equal(W.att$weights, W.atc$weights)

  # `ATOS` is not an allowable estimand for multi-category treatments
  expect_error({
    weightit(Am ~ X1 + X2 + X3 + X4 + X5,
             data = test_data, method = "bart", estimand = "ATOS",
             n.trees = n.trees, n.samples = n.samples, n.burn = n.burn)
  }, "not an allowable estimand", ignore.case = TRUE)

  # `s.weights` errors for multi-category treatments too
  expect_error({
    weightit(Am ~ X1 + X2 + X3 + X4 + X5,
             data = test_data, method = "bart", estimand = "ATE",
             s.weights = "SW",
             n.trees = n.trees, n.samples = n.samples, n.burn = n.burn)
  }, "sampling weights cannot be used", ignore.case = TRUE)
})

test_that("Continuous treatment", {
  skip_on_cran()
  skip_if_not_installed("dbarts")
  skip_if_not_installed("cobalt")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  n.trees <- 20
  n.samples <- 50
  n.burn <- 50

  set.seed(123)
  expect_no_condition({
    W0 <- weightit(Ac ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                   data = test_data, method = "bart",
                   n.trees = n.trees, n.samples = n.samples, n.burn = n.burn,
                   n.chains = 1, n.threads = 1,
                   include.obj = TRUE)
  })

  expect_false(is_null(W0$obj))
  expect_true(all(is.finite(W0$weights) & W0$weights >= 0))

  density.opts <- c("dnorm", "kernel")

  weight.mat <- matrix(nrow = nrow(test_data), ncol = length(density.opts))
  colnames(weight.mat) <- rep("", ncol(weight.mat))

  k <- 1

  for (density in density.opts) {
    test_that(sprintf("BART: density = %s", density), {
      set.seed(123)
      W <- weightit(Ac ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                    data = test_data, method = "bart",
                    density = density,
                    n.trees = n.trees, n.samples = n.samples, n.burn = n.burn,
                    n.chains = 1, n.threads = 1,
                    include.obj = TRUE)

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

      n <- sprintf("W_%s", density)
      colnames(weight.mat)[k] <<- n
      weight.mat[,k] <<- W$weights
      k <<- k + 1
    })
  }

  # Non-full rank
  set.seed(123)
  expect_no_condition({
    W <- weightit(Ac ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 +
                    I(1 - X5) + I(X9 * 2),
                  data = test_data, method = "bart",
                  n.trees = n.trees, n.samples = n.samples, n.burn = n.burn,
                  n.chains = 1, n.threads = 1)
  })

  expect_true(all(is.finite(W$weights) & W$weights >= 0))

  # `s.weights` errors for continuous treatments too
  expect_error({
    weightit(Ac ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
             data = test_data, method = "bart", s.weights = "SW",
             n.trees = n.trees, n.samples = n.samples, n.burn = n.burn)
  }, "sampling weights cannot be used", ignore.case = TRUE)
})

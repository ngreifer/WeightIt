test_that("Binary treatment", {
  skip_if_not_installed("rootSolve")
  skip_if_not_installed("cobalt")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-1

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  set.seed(123)
  base_weights <- runif(nrow(test_data))

  expect_no_condition({
    W0 <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                   data = test_data, method = "ebal", estimand = "ATE",
                   include.obj = TRUE, solver = "optim")
  })

  expect_M_parts_okay(W0, tolerance = eps)

  sw.opts <- c(FALSE, TRUE)
  bw.opts <- c(FALSE, TRUE)
  estimand.opts <- c("ATE", "ATT", "ATC")

  weight.mat <- matrix(nrow = nrow(test_data),
                       ncol = length(sw.opts) *
                         length(bw.opts) * length(estimand.opts))
  colnames(weight.mat) <- rep("", ncol(weight.mat))

  k <- 1

  for (sw in sw.opts) {
    for (bw in bw.opts) {
      for (estimand in estimand.opts) {
        expect_no_condition({
          W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                        data = test_data, method = "ebal", estimand = estimand,
                        s.weights = if (sw) "SW" else NULL,
                        base.weights = if (bw) base_weights else NULL,
                        include.obj = TRUE, solver = "multiroot")
        })

        n <- sprintf("W_%s_%s_%s", sw, bw, estimand)
        assign(n, W)

        expect_M_parts_okay(W, tolerance = eps)
        expect_equal(cobalt::col_w_smd(W$covs, W$treat, W$weights,
                                       s.weights = W$s.weights),
                     0 * cobalt::col_w_smd(W$covs, W$treat,
                                           s.weights = W$s.weights),
                     label = sprintf("SMDs for %s", n),
                     expected.label = "all 0s",
                     tolerance = eps)

        expect_true(is_null((!!{{ str2lang(n) }})$ps))
        expect_false(is_null((!!{{ str2lang(n) }})$obj))

        for (i in 0:1) {
          e <- {
            if (estimand == "ATT" && i == 1) expect_equal
            else if (estimand == "ATC" && i == 0) expect_equal
            else expect_not_equal
          }

          e(W$weights[W$treat == i],
            rep(1, sum(W$treat == i)),
            label = sprintf("%s weights for %s", i, n),
            expected.label = "all 1s",
            tolerance = eps)
        }

        for (i in seq_len(k - 1)) {
          expect_not_equal(W$weights, weight.mat[,i],
                           label = sprintf("Weights for %s", n),
                           expected.label = sprintf("weights for %s", colnames(weight.mat)[i]),
                           tolerance = eps)
        }

        colnames(weight.mat)[k] <- n
        weight.mat[,k] <- W$weights
        k <- k + 1

        rm(list = n)
      }
    }
  }

  # Estimands
  expect_error({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "ebal", estimand = "ATO")
  }, "not an allowable estimand")

  #Non-full rank
  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 +
                    I(1 - X5) + I(X9 * 2),
                  data = test_data, method = "ebal", estimand = "ATE",
                  include.obj = TRUE, solver = "optim")
  })

  expect_M_parts_okay(W, tolerance = eps)
  expect_equal(W$weights, W0$weights, tolerance = eps)

  # All categorical covariates (issue #86)
  expect_no_condition({
    W <- weightit(A ~ cut(X1, 3) + cut(X2, 3) + cut(X3, 3),
                  data = test_data, method = "ebal", estimand = "ATE",
                  include.obj = TRUE, solver = "optim", reltol = 1e-12)
  })

  expect_M_parts_okay(W, tolerance = eps)

  #Should be equivalent to CBPS and IPT with logit link for ATT
  for (sw in sw.opts) {
    expect_no_condition({
      W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                    data = test_data, method = "ebal", estimand = "ATT",
                    s.weights = if (sw) "SW" else NULL,
                    include.obj = TRUE, solver = "optim")
    })

    n <- sprintf("W_%s", sw)

    expect_no_condition({
      Wcbps <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                        data = test_data, method = "cbps", estimand = "ATT",
                        s.weights = if (sw) "SW" else NULL,
                        link = "logit", solver = "multiroot",
                        include.obj = TRUE)
    })

    ncbps <- sprintf("Wcbps_%s", sw)

    expect_equal(ESS(W$weights[W$treat == 0] * W$s.weights[W$treat == 0]),
                 ESS(Wcbps$weights[Wcbps$treat == 0] * Wcbps$s.weights[Wcbps$treat == 0]),
                 label = sprintf("ESS for %s", n),
                 expected.label = sprintf("ESS for %s", ncbps),
                 tolerance = .01)

    expect_no_condition({
      Wipt <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                       data = test_data, method = "ipt", estimand = "ATT",
                       s.weights = if (sw) "SW" else NULL,
                       link = "logit",
                       include.obj = TRUE)
    })

    nipt <- sprintf("Wipt_%s", sw)

    expect_equal(ESS(W$weights[W$treat == 0] * W$s.weights[W$treat == 0]),
                 ESS(Wipt$weights[Wipt$treat == 0] * Wipt$s.weights[Wipt$treat == 0]),
                 label = sprintf("ESS for %s", n),
                 expected.label = sprintf("ESS for %s", nipt),
                 tolerance = .01)
  }
})

test_that("Multi-category treatment", {
  skip_if_not_installed("rootSolve")
  skip_if_not_installed("cobalt")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-1

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  set.seed(123)
  base_weights <- runif(nrow(test_data))

  expect_no_condition({
    W0 <- weightit(Am ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                   data = test_data, method = "ebal", estimand = "ATE",
                   include.obj = TRUE, solver = "optim")
  })

  sw.opts <- c(FALSE, TRUE)
  bw.opts <- c(FALSE, TRUE)
  estimand.opts <- c("ATE", "ATT")

  weight.mat <- matrix(nrow = nrow(test_data),
                       ncol = length(sw.opts) *
                         length(bw.opts) * length(estimand.opts))
  colnames(weight.mat) <- rep("", ncol(weight.mat))

  k <- 1

  for (sw in sw.opts) {
    for (bw in bw.opts) {
      for (estimand in estimand.opts) {
        expect_no_condition({
          W <- weightit(Am ~ X1 + X2 + X3 + X4 + X5,
                        data = test_data, method = "ebal", estimand = estimand,
                        focal = if (estimand == "ATE") NULL else "T",
                        s.weights = if (sw) "SW" else NULL,
                        base.weights = if (bw) base_weights else NULL,
                        include.obj = TRUE, solver = "multiroot")
        })

        n <- sprintf("W_%s_%s_%s", sw, bw, estimand)
        assign(n, W)

        expect_M_parts_okay(W, tolerance = eps)
        for (tt in combn(levels(W$treat), 2, simplify = FALSE)) {
          in_tt <- W$treat %in% tt
          expect_equal(cobalt::col_w_smd(W$covs[in_tt,], W$treat[in_tt], W$weights[in_tt],
                                         s.weights = W$s.weights[in_tt]),
                       0 * cobalt::col_w_smd(W$covs[in_tt,], W$treat[in_tt],
                                             s.weights = W$s.weights[in_tt]),
                       label = sprintf("SMDs for %s", n),
                       expected.label = "all 0s",
                       tolerance = eps)
        }

        expect_true(is_null((!!{{ str2lang(n) }})$ps))
        expect_false(is_null((!!{{ str2lang(n) }})$obj))

        for (i in levels(W$treat)) {
          e <- {
            if (estimand == "ATT" && i == W$focal) expect_equal
            else expect_not_equal
          }

          e(W$weights[W$treat == i],
            rep(1, sum(W$treat == i)),
            label = sprintf("%s weights for %s", i, n),
            expected.label = "all 1s",
            tolerance = eps)
        }

        for (i in seq_len(k - 1)) {
          expect_not_equal(W$weights, weight.mat[,i],
                           label = sprintf("Weights for %s", n),
                           expected.label = sprintf("weights for %s", colnames(weight.mat)[i]),
                           tolerance = eps)
        }

        colnames(weight.mat)[k] <- n
        weight.mat[,k] <- W$weights
        k <- k + 1

        rm(list = n)
      }
    }
  }
})

test_that("Continuous treatment", {
  skip_if_not_installed("rootSolve")
  skip_if_not_installed("cobalt")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-1

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  set.seed(123)
  base_weights <- runif(nrow(test_data))

  expect_no_condition({
    W0 <- weightit(Ac ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                   data = test_data, method = "ebal",
                   include.obj = TRUE, solver = "optim")
  })

  expect_M_parts_okay(W0, tolerance = eps)

  sw.opts <- c(FALSE, TRUE)
  bw.opts <- c(FALSE, TRUE)
  d.moments.opts <- c(1, 3)

  weight.mat <- matrix(nrow = nrow(test_data),
                       ncol = length(sw.opts) *
                         length(bw.opts) * length(d.moments.opts))
  colnames(weight.mat) <- rep("", ncol(weight.mat))

  k <- 1

  for (sw in sw.opts) {
    for (bw in bw.opts) {
      for (d.moments in d.moments.opts) {
        expect_no_condition({
          W <- weightit(Ac ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                        data = test_data, method = "ebal",
                        d.moments = d.moments,
                        s.weights = if (sw) "SW" else NULL,
                        base.weights = if (bw) base_weights else NULL,
                        include.obj = TRUE, solver = "multiroot")
        })

        n <- sprintf("W_%s_%s_%s", sw, bw, d.moments)
        assign(n, W)

        expect_M_parts_okay(W, tolerance = eps)
        expect_equal(cobalt::col_w_cov(W$covs, W$treat, W$weights, std = TRUE,
                                       s.weights = W$s.weights),
                     0 * cobalt::col_w_cov(W$covs, W$treat, std = TRUE,
                                           s.weights = W$s.weights),
                     label = sprintf("A-X correlations for %s", n),
                     expected.label = "all 0s",
                     tolerance = eps)

        expect_equal(cobalt::col_w_mean(cbind(poly(W$treat, d.moments), W$covs), W$weights,
                                        s.weights = W$s.weights),
                     cobalt::col_w_mean(cbind(poly(W$treat, d.moments), W$covs),
                                            s.weights = W$s.weights),
                     label = sprintf("A+X weighted means for %s", n),
                     expected.label = "unweighted means",
                     tolerance = eps)

        expect_true(is_null((!!{{ str2lang(n) }})$ps))
        expect_false(is_null((!!{{ str2lang(n) }})$obj))

        for (i in seq_len(k - 1)) {
          expect_not_equal(W$weights, weight.mat[,i],
                           label = sprintf("Weights for %s", n),
                           expected.label = sprintf("weights for %s", colnames(weight.mat)[i]),
                           tolerance = eps)
        }

        colnames(weight.mat)[k] <- n
        weight.mat[,k] <- W$weights
        k <- k + 1

        rm(list = n)
      }
    }
  }

  #Non-full rank
  expect_no_condition({
    W <- weightit(Ac ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 +
                    I(1 - X5) + I(X9 * 2),
                  data = test_data, method = "ebal",
                  include.obj = TRUE, solver = "optim")
  })

  expect_M_parts_okay(W, tolerance = eps)
  expect_equal(W$weights, W0$weights, tolerance = eps)
})

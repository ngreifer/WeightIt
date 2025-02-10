test_that("Binary treatment", {
  skip_if_not_installed("rootSolve")
  skip_if_not_installed("cobalt")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-1

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  expect_no_condition({
    W0 <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                   data = test_data, method = "ipt", estimand = "ATE",
                   include.obj = TRUE)
  })

  expect_M_parts_okay(W0, tolerance = eps)

  sw.opts <- c(FALSE, TRUE)
  estimand.opts <- c("ATE", "ATT", "ATC")
  link.opts <- c("logit", "probit", "loglog", "cauchit")

  weight.mat <- matrix(nrow = nrow(test_data),
                       ncol = length(sw.opts) *
                         length(estimand.opts) * length(link.opts))
  colnames(weight.mat) <- rep("", ncol(weight.mat))

  k <- 1

  for (sw in sw.opts) {
    for (estimand in estimand.opts) {
      for (link in link.opts) {
        expect_no_condition({
          W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                        data = test_data, method = "ipt", estimand = estimand,
                        link = link,
                        s.weights = if (sw) "SW" else NULL,
                        include.obj = TRUE)
        })

        n <- sprintf("W_%s_%s_%s", sw, estimand, link)
        assign(n, W)

        expect_M_parts_okay(W, tolerance = eps)
        expect_equal(cobalt::col_w_smd(W$covs, W$treat, W$weights,
                                       s.weights = W$s.weights),
                     0 * cobalt::col_w_smd(W$covs, W$treat,
                                           s.weights = W$s.weights),
                     label = sprintf("SMDs for %s", n),
                     expected.label = "all 0s",
                     tolerance = eps)

        expect_true(is.numeric((!!{{ str2lang(n) }})$ps))
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
                  data = test_data, method = "ipt", estimand = "ATO")
  }, "not an allowable estimand")

  #Non-full rank
  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 +
                    I(1 - X5) + I(X9 * 2),
                  data = test_data, method = "ipt", estimand = "ATE",
                  include.obj = TRUE)
  })

  expect_M_parts_okay(W, tolerance = eps)
  expect_equal(W$weights, W0$weights, tolerance = eps)

  #Should be equivalent to CBPS for ATT
  for (sw in sw.opts) {
    for (link in link.opts) {
      expect_no_condition({
        W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                      data = test_data, method = "ipt", estimand = "ATT",
                      s.weights = if (sw) "SW" else NULL,
                      link = link,
                      include.obj = TRUE)
      })

      n <- sprintf("W_%s_%s", sw, link)

      expect_no_condition({
        Wcbps <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                          data = test_data, method = "cbps", estimand = "ATT",
                          s.weights = if (sw) "SW" else NULL,
                          link = link, solver = "multiroot",
                          include.obj = TRUE)
      })

      ncbps <- sprintf("Wcbps_%s_%s", sw, link)

      expect_equal(ESS(W$weights[W$treat == 0] * W$s.weights[W$treat == 0]),
                   ESS(Wcbps$weights[Wcbps$treat == 0] * Wcbps$s.weights[Wcbps$treat == 0]),
                   label = sprintf("ESS for %s", n),
                   expected.label = sprintf("ESS for %s", ncbps),
                   tolerance = .01)
    }
  }
})

test_that("Multi-category treatment", {
  skip_if_not_installed("rootSolve")
  skip_if_not_installed("cobalt")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-1

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  expect_no_condition({
    W0 <- weightit(Am ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                   data = test_data, method = "ipt", estimand = "ATE",
                   include.obj = TRUE)
  })

  sw.opts <- c(FALSE, TRUE)
  estimand.opts <- c("ATE", "ATT")
  link.opts <- c("logit", "probit", "loglog", "cauchit")

  weight.mat <- matrix(nrow = nrow(test_data),
                       ncol = length(sw.opts) *
                         length(estimand.opts) * length(link.opts))

  colnames(weight.mat) <- rep("", ncol(weight.mat))

  k <- 1

  for (sw in sw.opts) {
    for (estimand in estimand.opts) {
      for (link in link.opts) {
        expect_no_condition({
          W <- weightit(Am ~ X1 + X2 + X3 + X4 + X5,
                        data = test_data, method = "ipt", estimand = estimand,
                        focal = if (estimand == "ATE") NULL else "T",
                        s.weights = if (sw) "SW" else NULL,
                        link = link,
                        include.obj = TRUE)
        })

        n <- sprintf("W_%s_%s_%s", sw, estimand, link)
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

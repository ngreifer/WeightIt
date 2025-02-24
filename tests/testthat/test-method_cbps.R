test_that("Binary treatment", {
  skip_on_cran() #many tests, take too long
  skip_if_not_installed("rootSolve")
  skip_if_not_installed("cobalt")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-1

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  expect_no_condition({
    W0 <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                   data = test_data, method = "cbps", estimand = "ATE",
                   include.obj = TRUE, solver = "optim")
  })

  sw.opts <- c(FALSE, TRUE)
  over.opts <- c("exact", "twostep", "cont")
  estimand.opts <- c("ATE", "ATT", "ATC", "ATO")
  link.opts <- c("logit", "probit", "loglog", "cauchit")

  weight.mat <- matrix(nrow = nrow(test_data),
                       ncol = length(sw.opts) * length(over.opts) *
                         length(estimand.opts) * length(link.opts))
  colnames(weight.mat) <- rep("", ncol(weight.mat))
  k <- 1

  for (sw in sw.opts) {
    for (over in over.opts) {
      for (estimand in estimand.opts) {
        for (link in link.opts) {
          expect_no_condition({
            suppressWarnings({
            W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                          data = test_data, method = "cbps", estimand = estimand,
                          over = over != "exact", twostep = if (over == "exact") NULL else over == "twostep",
                          link = link,
                          s.weights = if (sw) "SW" else NULL,
                          include.obj = TRUE, solver = "multiroot")
            })
          })

          n <- sprintf("W_%s_%s_%s_%s", sw, over, estimand, link)
          assign(n, W)

          if (over == "exact") {
            expect_M_parts_okay(W, tolerance = eps)
            expect_equal(cobalt::col_w_smd(W$covs, W$treat, W$weights,
                                           s.weights = W$s.weights),
                         0 * cobalt::col_w_smd(W$covs, W$treat,
                                               s.weights = W$s.weights),
                         label = sprintf("SMDs for %s", n),
                         expected.label = "all 0s",
                         tolerance = eps)
          }

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

          if (link != "logit" || estimand != "ATO") {
            for (i in seq_len(k)) {
              expect_not_equal(W$weights, weight.mat[,i],
                               label = sprintf("Weights for %s", n),
                               expected.label = sprintf("weights for model %s", colnames(weight.mat)[i]),
                               tolerance = min(1e-8, eps))
            }
          }

          weight.mat[,k] <- W$weights
          colnames(weight.mat)[k] <- n
          k <- k + 1
        }
      }
    }
  }

  # Estimands
  expect_error({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "cbps", estimand = "ATM")
  }, "not an allowable estimand")

  #Non-full rank
  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 +
                    I(1 - X5) + I(X9 * 2),
                  data = test_data, method = "cbps", estimand = "ATE",
                  include.obj = TRUE)
  })

  expect_equal(W$weights, W0$weights, tolerance = eps)
})

test_that("Multi-category treatment", {
  skip_on_cran() #many tests, take too long
  skip_if_not_installed("rootSolve")
  skip_if_not_installed("cobalt")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-1

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  expect_no_condition({
    W0 <- weightit(Am ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                   data = test_data, method = "cbps", estimand = "ATE",
                   include.obj = TRUE, solver = "optim")
  })

  sw.opts <- c(FALSE, TRUE)
  over.opts <- c("exact", "twostep", "cont")
  estimand.opts <- c("ATE", "ATT", "ATO")

  weight.mat <- matrix(nrow = nrow(test_data),
                       ncol = length(sw.opts) * length(over.opts) *
                         length(estimand.opts))
  colnames(weight.mat) <- rep("", ncol(weight.mat))
  k <- 1

  for (sw in sw.opts) {
    for (over in over.opts) {
      for (estimand in estimand.opts) {
        expect_no_error({
          suppressWarnings({
          W <- weightit(Am ~ X1 + X2 + X3 + X4 + X5,
                        data = test_data, method = "cbps", estimand = estimand,
                        over = over != "exact", twostep = if (over == "exact") NULL else over == "twostep",
                        focal = if (estimand == "ATT") "T" else NULL,
                        s.weights = if (sw) "SW" else NULL,
                        include.obj = TRUE, solver = "multiroot")
          })
        })

        n <- sprintf("W_%s_%s_%s", sw, over, estimand)
        assign(n, W)

        if (over == "exact") {
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
                           tolerance = min(1e-8, eps))
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
  skip_on_cran() #many tests, take too long
  skip_if_not_installed("rootSolve")
  skip_if_not_installed("cobalt")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-1

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  expect_no_condition({
    W0 <- weightit(Ac ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                   data = test_data, method = "cbps",
                   include.obj = TRUE, solver = "optim")
  })

  sw.opts <- c(FALSE, TRUE)
  over.opts <- c("exact", "twostep", "cont")

  weight.mat <- matrix(nrow = nrow(test_data),
                       ncol = length(sw.opts) *
                         length(over.opts))
  colnames(weight.mat) <- rep("", ncol(weight.mat))

  k <- 1

  for (sw in sw.opts) {
    for (over in over.opts) {
      expect_no_error({
        suppressWarnings({
        W <- weightit(Ac ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                      data = test_data, method = "cbps",
                      over = over != "exact", twostep = if (over == "exact") NULL else over == "twostep",
                      s.weights = if (sw) "SW" else NULL,
                      include.obj = TRUE, solver = "multiroot")
        })
      })

      n <- sprintf("W_%s_%s", sw, over)
      assign(n, W)

      if (over == "exact") {
        expect_equal(cobalt::col_w_cov(W$covs, W$treat, W$weights, std = TRUE,
                                       s.weights = W$s.weights),
                     0 * cobalt::col_w_cov(W$covs, W$treat, std = TRUE,
                                           s.weights = W$s.weights),
                     label = sprintf("A-X correlations for %s", n),
                     expected.label = "all 0s",
                     tolerance = eps)
      }

      expect_true(is_null((!!{{ str2lang(n) }})$ps))
      expect_false(is_null((!!{{ str2lang(n) }})$obj))

      for (i in seq_len(k - 1)) {
        expect_not_equal(W$weights, weight.mat[,i],
                         label = sprintf("Weights for %s", n),
                         expected.label = sprintf("weights for %s", colnames(weight.mat)[i]),
                         tolerance = min(1e-8, eps))
      }

      colnames(weight.mat)[k] <- n
      weight.mat[,k] <- W$weights
      k <- k + 1

      rm(list = n)
    }
  }

  #Non-full rank
  expect_no_condition({
    W <- weightit(Ac ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 +
                    I(1 - X5) + I(X9 * 2),
                  data = test_data, method = "cbps",
                  include.obj = TRUE, solver = "optim")
  })

  expect_equal(W$weights, W0$weights, tolerance = eps)
})

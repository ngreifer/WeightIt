test_that("Binary treatment", {
  skip_on_cran()
  skip_if_not_installed("CBPS")
  skip_if_not_installed("cobalt")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  expect_no_condition({
    W0 <- weightit(A ~ X1 + X2 + X3,
                   data = test_data, method = "npcbps", estimand = "ATE",
                   include.obj = TRUE)
  })

  # npCBPS does not support M-estimation (no Mparts/Mparts.list attribute)
  expect_null(attr(W0, "Mparts", exact = TRUE))
  expect_null(attr(W0, "Mparts.list", exact = TRUE))

  expect_true(is_null(W0$ps))
  expect_false(is_null(W0$obj))

  expect_balance_improved(W0)

  # Additional arguments (moments, int, quantile) are spot-checked
  # individually rather than crossed, since only ATE is allowed and there is
  # no s.weights/estimand grid to combine them with.
  configs <- list(
    "moments = 2" = list(moments = 2),
    "int = TRUE" = list(int = TRUE),
    "quantile" = list(quantile = list(X1 = c(.25, .5, .75)))
  )

  weight.mat <- matrix(nrow = nrow(test_data), ncol = length(configs))
  colnames(weight.mat) <- rep("", ncol(weight.mat))

  k <- 1

  for (nm in names(configs)) {
    test_that(sprintf("npCBPS: %s", nm), {
      W <- do.call(weightit,
                  c(list(A ~ X1 + X2 + X3,
                        data = test_data, method = "npcbps", estimand = "ATE",
                        include.obj = TRUE),
                    configs[[nm]]))

      expect_true(is_null(W$ps))
      expect_false(is_null(W$obj))
      expect_balance_improved(W)

      expect_not_equal(unname(W$weights), unname(W0$weights),
                       expected.label = "weights for baseline")

      for (i in seq_len(k - 1)) {
        expect_not_equal(unname(W$weights), weight.mat[,i],
                         expected.label = sprintf("weights for %s", colnames(weight.mat)[i]),
                         tolerance = eps)
      }

      colnames(weight.mat)[k] <<- nm
      weight.mat[,k] <<- W$weights
      k <<- k + 1
    })
  }

  # Estimands: only ATE is allowed for npcbps
  for (estimand in c("ATT", "ATC", "ATO")) {
    expect_error({
      weightit(A ~ X1 + X2 + X3, data = test_data, method = "npcbps",
              estimand = estimand)
    }, "not an allowable estimand", ignore.case = TRUE)
  }

  #Non-full rank
  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + I(1 - X2) + I(X3 * 2),
                  data = test_data, method = "npcbps", estimand = "ATE",
                  include.obj = TRUE)
  })

  expect_equal(W$weights, W0$weights, tolerance = eps)

  d_na <- inject_missingness(test_data, "X1", prop = 0.05)

  expect_error({
    weightit(A ~ X1 + X2 + X3, data = d_na, method = "npcbps", estimand = "ATE",
            missing = "surr")
  }, "only.*allowed for.*missing", ignore.case = TRUE)

  expect_no_condition({
    W_na <- weightit(A ~ X1 + X2 + X3, data = d_na, method = "npcbps",
                     estimand = "ATE", missing = "ind")
  })

  expect_true(anyNA(W_na$covs))

  # s.weights
  expect_error({
    weightit(A ~ X1 + X2 + X3, data = test_data, method = "npcbps",
            estimand = "ATE", s.weights = "SW")
  }, "sampling weights cannot be used", ignore.case = TRUE)
})

test_that("Multi-category treatment", {
  skip_on_cran()
  skip_if_not_installed("CBPS")
  skip_if_not_installed("cobalt")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  expect_no_condition({
    W0 <- weightit(Am ~ X1 + X2 + X3,
                   data = test_data, method = "npcbps", estimand = "ATE",
                   include.obj = TRUE)
  })

  expect_null(attr(W0, "Mparts", exact = TRUE))

  expect_true(is_null(W0$ps))
  expect_false(is_null(W0$obj))

  for (tt in combn(levels(W0$treat), 2, simplify = FALSE)) {
    in_tt <- W0$treat %in% tt
    weighted <- abs(cobalt::col_w_smd(W0$covs[in_tt,], W0$treat[in_tt], W0$weights[in_tt]))
    unweighted <- abs(cobalt::col_w_smd(W0$covs[in_tt,], W0$treat[in_tt]))
    expect_true(max(weighted) < max(unweighted))
  }

  # Estimands: only ATE is allowed
  for (estimand in c("ATT", "ATC")) {
    expect_error({
      weightit(Am ~ X1 + X2 + X3, data = test_data, method = "npcbps",
              estimand = estimand)
    }, "not an allowable estimand", ignore.case = TRUE)
  }

  # Additional arguments: spot-checked individually
  configs <- list(
    "moments = 2" = list(moments = 2),
    "int = TRUE" = list(int = TRUE),
    "quantile" = list(quantile = list(X1 = c(.25, .5, .75)))
  )

  weight.mat <- matrix(nrow = nrow(test_data), ncol = length(configs))
  colnames(weight.mat) <- rep("", ncol(weight.mat))

  k <- 1

  for (nm in names(configs)) {
    test_that(sprintf("npCBPS: %s", nm), {
      W <- do.call(weightit,
                  c(list(Am ~ X1 + X2 + X3,
                        data = test_data, method = "npcbps", estimand = "ATE",
                        include.obj = TRUE),
                    configs[[nm]]))

      expect_true(is_null(W$ps))
      expect_false(is_null(W$obj))

      expect_not_equal(unname(W$weights), unname(W0$weights),
                       expected.label = "weights for baseline")

      for (i in seq_len(k - 1)) {
        expect_not_equal(unname(W$weights), weight.mat[,i],
                         expected.label = sprintf("weights for %s", colnames(weight.mat)[i]),
                         tolerance = eps)
      }

      colnames(weight.mat)[k] <<- nm
      weight.mat[,k] <<- W$weights
      k <<- k + 1
    })
  }
})

test_that("Continuous treatment", {
  skip_on_cran()
  skip_if_not_installed("CBPS")
  skip_if_not_installed("cobalt")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  expect_no_condition({
    W0 <- weightit(Ac ~ X1 + X2 + X3,
                   data = test_data, method = "npcbps",
                   include.obj = TRUE)
  })

  expect_null(attr(W0, "Mparts", exact = TRUE))

  expect_true(is_null(W0$ps))
  expect_false(is_null(W0$obj))

  weighted <- abs(cobalt::col_w_cov(W0$covs, W0$treat, W0$weights, std = TRUE))
  unweighted <- abs(cobalt::col_w_cov(W0$covs, W0$treat, std = TRUE))
  expect_true(max(weighted) < max(unweighted))

  #Non-full rank
  expect_no_condition({
    W <- weightit(Ac ~ X1 + X2 + X3 + I(1 - X2) + I(X3 * 2),
                  data = test_data, method = "npcbps",
                  include.obj = TRUE)
  })

  expect_equal(W$weights, W0$weights, tolerance = eps)

  # moments/int spot-checks
  configs <- list(
    "moments = 2" = list(moments = 2),
    "int = TRUE" = list(int = TRUE)
  )

  weight.mat <- matrix(nrow = nrow(test_data), ncol = length(configs))
  colnames(weight.mat) <- rep("", ncol(weight.mat))

  k <- 1

  for (nm in names(configs)) {
    test_that(sprintf("npCBPS: %s", nm), {
      W <- do.call(weightit,
                  c(list(Ac ~ X1 + X2 + X3,
                        data = test_data, method = "npcbps",
                        include.obj = TRUE),
                    configs[[nm]]))

      expect_true(is_null(W$ps))
      expect_false(is_null(W$obj))

      expect_not_equal(unname(W$weights), unname(W0$weights),
                       expected.label = "weights for baseline")

      for (i in seq_len(k - 1)) {
        expect_not_equal(unname(W$weights), weight.mat[,i],
                         expected.label = sprintf("weights for %s", colnames(weight.mat)[i]),
                         tolerance = eps)
      }

      colnames(weight.mat)[k] <<- nm
      weight.mat[,k] <<- W$weights
      k <<- k + 1
    })
  }
})

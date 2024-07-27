test_that("Binary treatment", {
  skip_if_not_installed("rootSolve")
  skip_if_not_installed("cobalt")
  skip_if_not_installed("brglm2")
  skip_if_not_installed("logistf")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-1

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  expect_no_condition({
    W0 <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                   data = test_data, method = "glm", estimand = "ATE",
                   include.obj = TRUE)
  })

  expect_M_parts_okay(W0)

  expect_true(is.numeric(W0$ps))

  # quick

  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "glm", estimand = "ATE",
                  quick = TRUE, include.obj = TRUE)
  })

  expect_M_parts_okay(W)

  expect_equal(W$weights, W0$weights)

  expect_false(is_null(W$obj))
  expect_false(is_null(W0$obj))
  expect_failure(expect_equal(W$obj, W0$obj))

  # Estimands

  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "glm", estimand = "ATT")
  })

  expect_M_parts_okay(W)

  expect_equal(W$weights[W$treat == 1], rep(1, sum(W$treat == 1)),
               tolerance = eps)

  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "glm", estimand = "ATC")
  })

  expect_M_parts_okay(W)

  expect_equal(W$weights[W$treat == 0], rep(1, sum(W$treat == 0)),
               tolerance = eps)

  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "glm", estimand = "ATO")
  })

  expect_M_parts_okay(W)

  expect_equal(unname(cobalt::col_w_smd(W$covs, W$treat, W$weights)),
               rep(0, 12),
               tolerance = eps)

  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "glm", estimand = "ATM")
  })

  expect_M_parts_okay(W)

  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "glm", estimand = "ATOS")
  })

  expect_M_parts_okay(W)

  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "glm", estimand = "ATE",
                  link = "probit")
  })

  expect_M_parts_okay(W)

  # brglm2
  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "glm", estimand = "ATE",
                  link = "br.logit", epsilon = 1e-10)
  })

  expect_M_parts_okay(W)

  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "glm", estimand = "ATE",
                  link = "br.probit", type = "AS_median", epsilon = 1e-10)
  })

  expect_M_parts_okay(W)

  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "glm", estimand = "ATE",
                  link = "br.logit", type = "correction", epsilon = 1e-10)
  })

  expect_null(attr(W, "Mparts"))

  # logistf
  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "glm", estimand = "ATE",
                  link = "flic")
  })

  expect_null(attr(W, "Mparts"))

  # s.weights

  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "glm", estimand = "ATE",
                  s.weights = "SW")
  })

  expect_equal(test_data$SW, W$s.weights)

  expect_M_parts_okay(W)

  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X5 + X6,
                  data = test_data, method = "glm", estimand = "ATE",
                  s.weights = "SW", link = "log")
  })

  expect_M_parts_okay(W)

  # No warning for non-integer #successes

  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "glm", estimand = "ATE",
                  link = "br.logit", s.weights = "SW", epsilon = 1e-10)
  })
  expect_M_parts_okay(W)

  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "glm", estimand = "ATE",
                  link = "flac", s.weights = "SW")
  })

  expect_equal(test_data$SW, W$s.weights)

  #Stabilization

  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "glm", estimand = "ATE",
                  include.obj = TRUE, stabilize = TRUE)
  })

  expect_M_parts_okay(W)

  expect_null(attr(W, "Mparts", exact = TRUE))
  expect_false(is_null(attr(W, "Mparts.list", exact = TRUE)))

  expect_equal(cobalt::col_w_smd(W$covs, W$treat, W$weights),
               cobalt::col_w_smd(W0$covs, W0$treat, W0$weights),
               tolerance = eps)

  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "glm", estimand = "ATE",
                  include.obj = TRUE, stabilize = ~X1)
  })

  expect_M_parts_okay(W)

  expect_failure(expect_equal(cobalt::col_w_smd(W$covs, W$treat, W$weights),
                              cobalt::col_w_smd(W0$covs, W0$treat, W0$weights)))

  #Non-full rank
  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 +
                    I(1 - X5) + I(X9 * 2),
                  data = test_data, method = "glm", estimand = "ATE",
                  include.obj = TRUE)
  })

  expect_equal(W$weights, W0$weights)

  # Separation
  set.seed(123)
  test_data$Xx <- rbinom(nrow(test_data), 1, .01 + .99 * test_data$A)

  expect_warning({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 +
                    Xx,
                  data = test_data, method = "glm", estimand = "ATE",
                  include.obj = TRUE)
  }, "Propensity scores numerically equal to 0 or 1 were estimated")

  # expect_failure(expect_M_parts_okay(W))

  test_data$Xx <- NULL
})

test_that("Ordinal treatment", {
  skip_if_not_installed("rootSolve")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-1

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))
  test_data$Ao <- ordered(findInterval(test_data$Ac, quantile(test_data$Ac, seq(0, 1, length.out = 5)),
                                       all.inside = TRUE))

  expect_no_condition({
    W0 <- weightit(Ao ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "glm", estimand = "ATE",
                  include.obj = TRUE)
  })

  expect_M_parts_okay(W0)

  # expect_no_condition({
  #   W <- weightit(Ao ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
  #                  data = test_data, method = "glm", estimand = "ATE",
  #                 link = "br.logit", parallel = TRUE,
  #                 include.obj = TRUE)
  # })
  # expect_failure(expect_M_parts_okay(W))

  expect_no_condition({
    W <- weightit(Ao ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "glm", estimand = "ATE",
                  multi.method = "weightit",
                  include.obj = TRUE)
  })

  expect_M_parts_okay(W)
})

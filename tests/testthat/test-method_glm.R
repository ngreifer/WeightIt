test_that("Binary treatment", {
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

  expect_equal(W$weights[W$treat == 1], rep(1, sum(W$treat == 1)))

  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "glm", estimand = "ATC")
  })

  expect_M_parts_okay(W)

  expect_equal(W$weights[W$treat == 0], rep(1, sum(W$treat == 0)))

  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "glm", estimand = "ATO")
  })

  expect_M_parts_okay(W)

  expect_equal(unname(cobalt::col_w_smd(W$covs, W$treat, W$weights)),
               rep(0, 12))

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
                  link = "br.logit")
  })

  expect_null(attr(W, "Mparts"))

  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "glm", estimand = "ATE",
                  link = "flic")
  })

  expect_null(attr(W, "Mparts"))

  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "glm", estimand = "ATE",
                  link = "probit")
  })

  expect_M_parts_okay(W)

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
                  link = "br.logit", s.weights = "SW")
  })

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
               cobalt::col_w_smd(W0$covs, W0$treat, W0$weights))

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

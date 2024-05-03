test_that("Binary treatment", {
  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                   data = test_data, method = "glm", estimand = "ATE",
                   include.obj = TRUE)
  })

  expect_M_parts_okay(W)

  expect_no_condition({
    fit0 <- glm_weightit(Y_C ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                        data = test_data, weightit = W)
  })

  #M-estimation for glm
  expect_no_condition({
    fit <- glm_weightit(Y_C ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                        data = test_data, weightit = W, vcov = "asympt")
  })

  expect_equal(coef(fit0), coef(fit))
  expect_equal(vcov(fit0), vcov(fit))

  expect_no_condition({
    fit <- glm_weightit(Y_C ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                        data = test_data, weightit = W, vcov = "HC0")
  })

  expect_equal(coef(fit0), coef(fit))
  expect_failure(expect_equal(vcov(fit0), vcov(fit)))

  set.seed(123)
  expect_no_condition({
    fit <- glm_weightit(Y_C ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                        data = test_data, weightit = W, vcov = "FWB", R = 50)
  })

  expect_equal(coef(fit0), coef(fit))
  expect_failure(expect_equal(vcov(fit0), vcov(fit)))

  set.seed(123)
  expect_no_condition({
    fit_ <- glm_weightit(Y_C ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                        data = test_data, weightit = W, vcov = "FWB", R = 50,
                        fwb.args = list(wtype = "mammen"))
  })

  expect_equal(coef(fit), coef(fit_))
  expect_failure(expect_equal(vcov(fit), vcov(fit_)))

  expect_no_condition({
    fit <- glm_weightit(Y_C ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                        data = test_data, weightit = W, vcov = "BS", R = 50)
  })

  expect_equal(coef(fit0), coef(fit))
  expect_failure(expect_equal(vcov(fit0), vcov(fit)))


})
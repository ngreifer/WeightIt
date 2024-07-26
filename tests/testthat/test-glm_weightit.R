test_that("No weights", {
  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  expect_no_condition({
    fit0 <- glm_weightit(Y_B ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                         data = test_data, family = binomial)
  })

  #M-estimation for glm
  expect_no_condition({
    fit <- glm_weightit(Y_B ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                        data = test_data, family = binomial, vcov = "HC0")
  })

  expect_equal(coef(fit0), coef(fit))
  expect_equal(vcov(fit0), vcov(fit))

  fit_g <- glm(Y_B ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
               data = test_data, family = binomial)
  expect_equal(coef(fit0), coef(fit_g))
  expect_equal(vcov(fit0), sandwich::sandwich(fit_g),
               tolerance = 1e-6)

  #Offset
  set.seed(123)
  off <- runif(nrow(test_data))

  expect_no_condition({
    fit <- glm_weightit(Y_B ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9) + offset(off),
                        data = test_data, family = binomial)
  })

  expect_failure(expect_equal(coef(fit0), coef(fit)))

  fit_g <- glm(Y_B ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9) + offset(off),
               data = test_data, family = binomial)
  expect_equal(coef(fit), coef(fit_g))
  expect_equal(vcov(fit), sandwich::sandwich(fit_g),
               tolerance = 1e-6)

  #Cluster-robust SEs
  clus <- sample(1:50, nrow(test_data), replace = TRUE)

  expect_no_condition({
    fit <- glm_weightit(Y_B ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                        data = test_data, family = binomial, cluster = clus)
  })

  expect_equal(coef(fit0), coef(fit))

  fit_g <- glm(Y_B ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
               data = test_data, family = binomial)
  expect_equal(coef(fit), coef(fit_g))
  expect_equal(vcov(fit), sandwich::vcovCL(fit_g, cluster = clus),
               tolerance = 1e-5)

  #BR
  expect_no_condition({
    fit <- glm_weightit(Y_B ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                        data = test_data, family = binomial("probit"), br = TRUE)
  })

  fit_g <- glm(Y_B ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
               data = test_data, family = binomial("probit"),
               method = brglm2::brglmFit)

  expect_equal(coef(fit), coef(fit_g))

})

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
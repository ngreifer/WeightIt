test_that("No weights", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("sandwich")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-1

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))
  test_data$Y_O <- with(test_data, factor(findInterval(Y_C, quantile(Y_C, seq(0, 1, length = 5)),
                                                       all.inside = TRUE), ordered = TRUE))
  set.seed(123)
  test_data$off <- runif(nrow(test_data))
  test_data$clus <- sample(1:50, nrow(test_data), replace = TRUE)

  expect_no_condition({
    fit0 <- ordinal_weightit(Y_O ~ A * (X1 + X2 + X3 + X4 + X5),
                            data = test_data)
  })

  #M-estimation for polr
  expect_no_condition({
    fit <- ordinal_weightit(Y_O ~ A * (X1 + X2 + X3 + X4 + X5),
                           data = test_data, vcov = "HC0")
  })

  expect_equal(coef(fit0), coef(fit), tolerance = eps)
  expect_equal(vcov(fit0), vcov(fit), tolerance = eps)

  fit_g <- MASS::polr(Y_O ~ A * (X1 + X2 + X3 + X4 + X5),
                      data = test_data, Hess = TRUE,
                      control = list(reltol = 1e-12))

  .coef <- function(x) {
    c(x$coefficients, x$zeta)
  }

  expect_equal(coef(fit), .coef(fit_g),
               tolerance = eps)
  expect_equal(vcov(fit), sandwich::sandwich(fit_g),
               tolerance = eps)

  expect_no_condition({
    fit <- ordinal_weightit(Y_O ~ A * (X1 + X2 + X3 + X4 + X5), cluster = ~clus,
                         data = test_data)
  })

  expect_equal(coef(fit0), coef(fit), tolerance = eps)

  #Cluster-robust SEs
  expect_equal(vcov(fit),
               sandwich::vcovCL(fit_g, cluster = ~clus),
               tolerance = eps)

  #Offset
  expect_no_condition({
    fit <- ordinal_weightit(Y_O ~ A * (X1 + X2 + X3 + X4 + X5) + offset(off),
                           data = test_data)
  })

  expect_failure(expect_equal(coef(fit0), coef(fit)))

  fit_g <- MASS::polr(Y_O ~ A * (X1 + X2 + X3 + X4 + X5) + offset(off),
                      data = test_data, Hess = TRUE,
                      control = list(reltol = 1e-12))

  expect_equal(coef(fit), .coef(fit_g),
               tolerance = eps)
  # expect_equal(vcov(fit), sandwich::sandwich(fit_g),
  #              tolerance = eps)

  #Probit
  expect_no_condition({
    fit <- ordinal_weightit(Y_O ~ A  * (X1 + X2 + X3 + X4 + X5),
                         data = test_data, vcov = "HC0",
                         link = "probit")
  })

  suppressWarnings({
    fit_g <- MASS::polr(Y_O ~ A * (X1 + X2 + X3 + X4 + X5),
                        data = test_data, Hess = TRUE,
                        control = list(reltol = 1e-12),
                        method = "probit")
  })

  expect_equal(coef(fit), .coef(fit_g),
               tolerance = eps)
  expect_equal(vcov(fit), sandwich::sandwich(fit_g),
               tolerance = eps)
})

test_that("Binary treatment", {
  skip_if_not_installed("MASS")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-1

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))
  test_data$Y_O <- with(test_data, factor(findInterval(Y_C, quantile(Y_C, seq(0, 1, length = 5)),
                                                       all.inside = TRUE), ordered = TRUE))
  set.seed(123)
  test_data$off <- runif(nrow(test_data))
  test_data$clus <- sample(1:50, nrow(test_data), replace = TRUE)

  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "glm", estimand = "ATE",
                  include.obj = TRUE)
  })

  expect_no_condition({
    fit0 <- ordinal_weightit(Y_O ~ A * (X1 + X2 + X3 + X4 + X5),
                            data = test_data, weightit = W)
  })

  #M-estimation for polr
  expect_no_condition({
    fit <- ordinal_weightit(Y_O ~ A * (X1 + X2 + X3 + X4 + X5),
                           data = test_data,  weightit = W, vcov = "asympt")
  })

  expect_equal(coef(fit0), coef(fit), tolerance = eps)
  expect_equal(vcov(fit0), vcov(fit), tolerance = eps)

  expect_no_condition({
    fit <- ordinal_weightit(Y_O ~ A  * (X1 + X2 + X3 + X4 + X5),
                           data = test_data, weightit = W, vcov = "HC0")
  })

  suppressWarnings({
    fit_g <- MASS::polr(Y_O ~ A * (X1 + X2 + X3 + X4 + X5),
                        data = test_data, Hess = TRUE,
                        weights = W$weights,
                        control = list(reltol = 1e-12))
  })

  .coef <- function(x) {
    c(x$coefficients, x$zeta)
  }

  expect_equal(coef(fit), .coef(fit_g),
               tolerance = eps)
  # expect_equal(vcov(fit), sandwich::sandwich(fit_g),
  #              tolerance = eps)

  expect_no_condition({
    fit <- ordinal_weightit(Y_O ~ A * (X1 + X2 + X3 + X4 + X5), cluster = ~clus,
                         data = test_data, weightit = W)
  })

  expect_equal(coef(fit0), coef(fit), tolerance = eps)

  #Cluster-robust SEs
  # expect_equal(vcov(fit),
  #              sandwich::vcovCL(fit_g, cluster = ~clus),
  #              tolerance = eps)

  #Offset
  expect_no_condition({
    fit <- ordinal_weightit(Y_O ~ A * (X1 + X2 + X3 + X4 + X5) + offset(off),
                           data = test_data, weightit = W)
  })

  expect_failure(expect_equal(coef(fit0), coef(fit)))

  suppressWarnings({
    fit_g <- MASS::polr(Y_O ~ A * (X1 + X2 + X3 + X4 + X5) + offset(off),
                        data = test_data, Hess = TRUE,
                        weights = W$weights,
                        control = list(reltol = 1e-12))
  })

  expect_equal(coef(fit), .coef(fit_g),
               tolerance = eps)

  #Probit
  expect_no_condition({
    fit <- ordinal_weightit(Y_O ~ A  * (X1 + X2 + X3 + X4 + X5),
                         data = test_data, vcov = "HC0",
                         link = "probit", weightit = W)
  })

  suppressWarnings({
    fit_g <- MASS::polr(Y_O ~ A * (X1 + X2 + X3 + X4 + X5),
                        data = test_data, Hess = TRUE,
                        control = list(reltol = 1e-12),
                        method = "probit", weights = W$weights)
  })

  expect_equal(coef(fit), .coef(fit_g),
               tolerance = eps)
  # expect_equal(vcov(fit), sandwich::sandwich(fit_g),
  #              tolerance = eps)

})

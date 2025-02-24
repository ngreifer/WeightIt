test_that("No weights", {
  skip_if_not_installed("sandwich")
  skip_if_not_installed("mlogit")
  skip_if_not_installed("dfidx")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-1

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))
  test_data$Y_M <- with(test_data, factor(findInterval(Y_C, quantile(Y_C, seq(0, 1, length = 5)),
                                                       all.inside = TRUE)))
  set.seed(123)
  test_data$off <- runif(nrow(test_data))
  test_data$clus <- sample(1:50, nrow(test_data), replace = TRUE)

  expect_no_condition({
    fit0 <- multinom_weightit(Y_M ~ A * (X1 + X2 + X3 + X4 + X5),
                         data = test_data)
  })

  #M-estimation for mlogit
  expect_no_condition({
    fit <- multinom_weightit(Y_M ~ A * (X1 + X2 + X3 + X4 + X5),
                        data = test_data, vcov = "HC0")
  })

  expect_equal(coef(fit0), coef(fit), tolerance = eps)
  expect_equal(vcov(fit0), vcov(fit), tolerance = eps)

  fit_g <- mlogit::mlogit(Y_M ~ 0 | A * (X1 + X2 + X3 + X4 + X5),
                          data = dfidx::dfidx(test_data, choice = "Y_M", shape = "wide"))

  ind <- unlist(split(seq_along(coef(fit0)), rep(seq_len(nlevels(test_data$Y_M) - 1),
                                                 length(coef(fit0))/(nlevels(test_data$Y_M) - 1))))

  expect_equal(unname(coef(fit0)), unname(coef(fit_g)[ind]),
               tolerance = eps)
  expect_equal(unname(vcov(fit0)), unname(sandwich::sandwich(fit_g)[ind, ind]),
               tolerance = eps)

  expect_no_condition({
    fit <- multinom_weightit(Y_M ~ A * (X1 + X2 + X3 + X4 + X5), cluster = ~clus,
                           data = test_data)
  })

  expect_equal(coef(fit0), coef(fit), tolerance = eps)

  #Cluster-robust SEs
  expect_equal(unname(vcov(fit)),
               unname(sandwich::vcovCL(fit_g, cluster = test_data$clus)[ind, ind]),
               tolerance = eps)

  #Offset
  expect_no_condition({
    fit <- multinom_weightit(Y_M ~ A * (X1 + X2 + X3 + X4 + X5) + offset(off),
                        data = test_data)
  })

  expect_not_equal(coef(fit0), coef(fit), tolerance = eps)

  #Test using sandwich functions
  expect_no_condition({
    fit0 <- glm_weightit(Y_B ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                         data = test_data, family = binomial)
  })

  expect_equal(vcov(fit0), sandwich::sandwich(fit0),
               tolerance = eps)
})

test_that("Binary treatment", {
  skip_if_not_installed("mlogit")
  skip_if_not_installed("dfidx")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-1

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))
  test_data$Y_M <- with(test_data, factor(findInterval(Y_C, quantile(Y_C, seq(0, 1, length = 5)),
                                                       all.inside = TRUE)))
  set.seed(123)
  test_data$off <- runif(nrow(test_data))
  test_data$clus <- sample(1:50, nrow(test_data), replace = TRUE)

  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "glm", estimand = "ATE",
                  include.obj = TRUE)
  })

  expect_no_condition({
    fit0 <- multinom_weightit(Y_M ~ A * (X1 + X2 + X3 + X4 + X5),
                            data = test_data, weightit = W)
  })

  #M-estimation for mlogit
  expect_no_condition({
    fit <- multinom_weightit(Y_M ~ A * (X1 + X2 + X3 + X4 + X5),
                           data = test_data,  weightit = W, vcov = "asympt")
  })

  expect_equal(coef(fit0), coef(fit), tolerance = eps)
  expect_equal(vcov(fit0), vcov(fit), tolerance = eps)

  expect_no_condition({
    fit <- multinom_weightit(Y_M ~ A  * (X1 + X2 + X3 + X4 + X5),
                           data = test_data, weightit = W, vcov = "HC0")
  })

  mlogit_data <- dfidx::dfidx(transform(test_data, .weights = W$weights),
                              choice = "Y_M", shape = "wide")

  fit_g <- mlogit::mlogit(Y_M ~ 0 | A  * (X1 + X2 + X3 + X4 + X5),
                          data = mlogit_data,
                          weights = .weights, tol = 1e-12, ftol = 1e-12)

  ind <- unlist(split(seq_along(coef(fit)), rep(seq_len(nlevels(test_data$Y_M) - 1),
                                                 length(coef(fit))/(nlevels(test_data$Y_M) - 1))))

  expect_equal(unname(coef(fit)), unname(coef(fit_g)[ind]),
               tolerance = eps)
  # expect_equal(unname(vcov(fit)), unname(sandwich::sandwich(fit_g)[ind, ind]),
  #              tolerance = eps)

  expect_no_condition({
    fit <- multinom_weightit(Y_M ~ A * (X1 + X2 + X3 + X4 + X5), cluster = ~clus,
                           data = test_data, weightit = W, vcov = "HC0")
  })

  expect_equal(coef(fit0), coef(fit), tolerance = eps)

  # #Cluster-robust SEs
  # expect_equal(unname(vcov(fit)),
  #              unname(sandwich::vcovCL(fit_g, cluster = test_data$clus, cadjust = T)[ind, ind]),
  #              tolerance = eps)

  #Offset
  expect_no_condition({
    fit <- multinom_weightit(Y_M ~ A * (X1 + X2 + X3 + X4 + X5) + offset(off),
                           data = test_data)
  })

  expect_not_equal(coef(fit0), coef(fit), tolerance = eps)

  #Test using sandwich functions
  expect_no_condition({
    fit0 <- multinom_weightit(Y_M ~ A * (X1 + X2 + X3 + X4 + X5),
                              data = test_data, weightit = W)
  })

  expect_equal(vcov(fit0),
               sandwich::sandwich(fit0),
               tolerance = eps)

  expect_equal(vcov(fit0, type = "HC0"),
               sandwich::sandwich(fit0, asympt = FALSE),
               tolerance = eps)
})

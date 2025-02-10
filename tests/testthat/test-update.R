test_that("update.glm_weightit() works", {
  skip_if(!capabilities("long.double"))
  eps <- if (capabilities("long.double")) 1e-5 else 1e-1

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  expect_no_condition({
    fit0 <- glm_weightit(Y_B ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                         data = test_data, family = binomial)
  })

  #Updating formula
  expect_no_condition({
    fit1 <- glm_weightit(Y_B ~ A,
                         data = test_data, family = binomial)
  })

  expect_equal(update(fit0, formula = . ~ A),
               fit1, tolerance = eps,
               ignore_attr = "class",
               ignore_formula_env = TRUE)

  #Updating dataset
  sub <- test_data$X3 > 0
  test_data_s <- test_data[sub,]

  expect_no_condition({
    fit1 <- glm_weightit(Y_B ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                         data = test_data_s, family = binomial)
  })

  expect_equal(update(fit0, data = test_data_s),
               fit1, tolerance = eps,
               ignore_attr = "class",
               ignore_formula_env = TRUE)

  expect_equal(update(fit0, subset = sub)$coefficients,
               fit1$coefficients, tolerance = eps,
               ignore_attr = "class",
               ignore_formula_env = TRUE)

  #Updating vcov only
  expect_no_condition({
    fit1 <- glm_weightit(Y_B ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                         data = test_data, family = binomial, vcov = "const")
  })

  expect_equal(update(fit0, vcov = "const"),
               fit1, tolerance = eps,
               ignore_attr = "class",
               ignore_formula_env = TRUE)

  expect_not_equal(vcov(fit0), vcov(fit1))

  #Model should not be refit when only vcov is changed
  clus <- sample(1:50, nrow(test_data), replace = TRUE)
  i <- FALSE
  expect_no_condition({
    fit2 <- glm_weightit(Y_B ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                         data = test_data, family = binomial,
                         control = if (i) list(stop("bad error")) else list())
  })

  i <- TRUE
  expect_no_condition({
    update(fit2, vcov = "const")
  })

  expect_error({
    update(fit2, family = binomial("probit"))
  }, "bad error")

  expect_error({
    update(fit2)
  }, "bad error")

  expect_no_condition({
    fit2c <- update(fit2, cluster = ~clus)
  })

  expect_equal(vcov(fit2c),
               vcov(update(fit0, cluster = ~clus)),
               tolerance = eps)

  expect_not_equal(vcov(fit2c),
                   vcov(fit0))

  #Updating s.weights without weightit object should add one
  expect_no_condition({
    suppressWarnings({
      fitg <- glm(Y_B ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                  data = test_data, family = binomial,
                  weights = SW)
    })
  })

  expect_no_condition({
    fit0s <- update(fit0, s.weights = "SW")
  })

  expect_null(fit0$weightit)
  expect_equal(fit0s$weightit,
               list(s.weights = test_data$SW,
                    weights = rep(1, nrow(test_data)),
                    method = NULL),
               ignore_attr = c("names", "class"))
  expect_equal(fit0s$weightit$s.weights,
               fit0s$prior.weights,
               ignore_attr = "names")

  expect_equal(coef(fit0s),
               coef(fitg),
               tolerance = eps)

})

test_that("update.weightit() works", {
  skip_if(!capabilities("long.double"))
  eps <- if (capabilities("long.double")) 1e-5 else 1e-1

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  expect_no_condition({
    W0 <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                   data = test_data)
  })

  #Updating formula
  expect_no_condition({
    W1 <- weightit(A ~ X1 + X2 + X3 + X4,
                   data = test_data)
  })

  expect_equal(update(W0, formula = . ~ X1 + X2 + X3 + X4),
               W1, tolerance = eps,
               ignore_attr = "class",
               ignore_formula_env = TRUE)

  #Updating dataset
  sub <- test_data$X3 > 0
  test_data_s <- test_data[sub,]

  expect_no_condition({
    W1 <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                   data = test_data_s)
  })

  expect_equal(update(W0, data = test_data_s),
               W1, tolerance = eps,
               ignore_attr = "class",
               ignore_formula_env = TRUE)

  #Updating method and estimand
  expect_no_condition({
    W1 <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                   data = test_data, method = "ebal", estimand = "ATT")
  })

  expect_equal(update(W0, method = "ebal", estimand = "ATT"),
               W1, tolerance = eps,
               ignore_attr = "class",
               ignore_formula_env = TRUE)

  #Updating s.weights
  expect_no_condition({
    W1 <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                   data = test_data, s.weights = "SW")
  })

  expect_equal(update(W0, s.weights = "SW"),
               W1, tolerance = eps,
               ignore_attr = "class",
               ignore_formula_env = TRUE)

})

test_that("update.glm_weightit() works with weightit", {
  skip_if(!capabilities("long.double"))
  eps <- if (capabilities("long.double")) 1e-5 else 1e-1

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  expect_no_condition({
    W0 <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                   data = test_data)
  })

  expect_no_condition({
    fit0 <- glm_weightit(Y_B ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                         data = test_data, family = binomial, weightit = W0)
  })

  #Updating weightit
  expect_no_condition({
    fit1 <- glm_weightit(Y_B ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                         data = test_data, family = binomial)
  })

  expect_equal(update(fit0, weightit = NULL),
               fit1, tolerance = eps,
               ignore_attr = "class",
               ignore_formula_env = TRUE)

  #Updating dataset
  sub <- test_data$X3 > 0
  test_data_s <- test_data[sub,]

  expect_no_condition({
    W1 <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                   data = test_data_s)
  })

  expect_no_condition({
    fit1 <- glm_weightit(Y_B ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                         data = test_data_s, family = binomial, weightit = W1)
  })

  expect_equal(update(fit0, data = test_data_s)$coefficients,
               fit1$coefficients, tolerance = eps,
               ignore_attr = "class",
               ignore_formula_env = TRUE)

  #Updating vcov only
  expect_no_condition({
    fit1 <- glm_weightit(Y_B ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                         data = test_data, family = binomial, weightit = W0,
                         vcov = "HC0")
  })

  expect_equal(update(fit0, vcov = "HC0"),
               fit1, tolerance = eps,
               ignore_attr = "class",
               ignore_formula_env = TRUE)

  expect_not_equal(vcov(fit0), vcov(fit1))

  #Model should not be refit when only vcov is changed
  clus <- sample(1:50, nrow(test_data), replace = TRUE)
  i <- FALSE
  expect_no_condition({
    fit2 <- glm_weightit(Y_B ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                         data = test_data, family = binomial, weightit = W0,
                         control = if (i) list(stop("bad error")) else list())
  })

  i <- TRUE
  expect_no_condition({
    update(fit2, vcov = "HC0")
  })

  expect_error({
    update(fit2, family = binomial("probit"))
  }, "bad error")

  expect_error({
    update(fit2)
  }, "bad error")

  expect_no_condition({
    fit2c <- update(fit2, cluster = ~clus)
  })

  expect_equal(vcov(fit2c),
               vcov(update(fit0, cluster = ~clus)),
               tolerance = eps)

  expect_not_equal(vcov(fit2c),
                   vcov(fit0))

  #Updating s.weights updates weightit
  expect_no_condition({
    W1 <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                   data = test_data, s.weights = "SW")
  })

  expect_no_condition({
    fit1 <- glm_weightit(Y_B ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                         data = test_data, family = binomial, weightit = W1)
  })

  expect_no_condition({
    fit1s <- glm_weightit(Y_B ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                          data = test_data, family = binomial,
                          weightit = update(object = W0, s.weights = "SW"))
  })

  expect_no_condition({
    fit0s <- update(fit0, s.weights = "SW")
  })

  expect_equal(fit0s$coefficients,
               fit1$coefficients,
               ignore_attr = c("names", "class"),
               tolerance = eps)

  expect_equal(fit0s$coefficients,
               fit1s$coefficients,
               ignore_attr = c("names", "class"),
               tolerance = eps)

  expect_equal(fit0s,
               fit1s,
               tolerance = eps)

  #Mimic bootstrapping
  boot_dat <- test_data[sample(nrow(test_data), replace = TRUE),]

  expect_no_condition({
    W1 <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                   data = boot_dat)
  })

  expect_no_condition({
    fit1 <- glm_weightit(Y_B ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                         data = boot_dat, family = binomial, weightit = W1)
  })

  expect_equal(update(fit0, data = boot_dat, vcov = "none")$coefficients,
               fit1$coefficients, tolerance = eps,
               ignore_attr = "class",
               ignore_formula_env = TRUE)
})

test_that("update.multinom_weightit() works with weightit", {
  skip_if(!capabilities("long.double"))
  eps <- if (capabilities("long.double")) 1e-5 else 1e-1

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))
  test_data$Y_M <- with(test_data, factor(findInterval(Y_C, quantile(Y_C, seq(0, 1, length = 5)),
                                                       all.inside = TRUE)))

  expect_no_condition({
    W0 <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                   data = test_data)
  })

  expect_no_condition({
    fit0 <- multinom_weightit(Y_M ~ A * (X1 + X2 + X3 + X4 + X5),
                              data = test_data, weightit = W0)
  })

  #Updating weightit
  expect_no_condition({
    fit1 <- multinom_weightit(Y_M ~ A * (X1 + X2 + X3 + X4 + X5),
                              data = test_data)
  })

  expect_equal(update(fit0, weightit = NULL),
               fit1, tolerance = eps,
               ignore_attr = "class",
               ignore_formula_env = TRUE)

  #Updating dataset
  sub <- test_data$X3 > 0
  test_data_s <- test_data[sub,]

  expect_no_condition({
    W1 <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                   data = test_data_s)
  })

  expect_no_condition({
    fit1 <- multinom_weightit(Y_M ~ A * (X1 + X2 + X3 + X4 + X5),
                              data = test_data_s, weightit = W1)
  })

  expect_equal(update(fit0, data = test_data_s)$coefficients,
               fit1$coefficients, tolerance = eps,
               ignore_attr = "class",
               ignore_formula_env = TRUE)

  #Updating vcov only
  expect_no_condition({
    fit1 <- multinom_weightit(Y_M ~ A * (X1 + X2 + X3 + X4 + X5),
                              data = test_data, weightit = W0,
                              vcov = "HC0")
  })

  expect_equal(update(fit0, vcov = "HC0"),
               fit1, tolerance = eps,
               ignore_attr = "class",
               ignore_formula_env = TRUE)

  expect_not_equal(vcov(fit0), vcov(fit1))

  #Model should not be refit when only vcov is changed
  clus <- sample(1:50, nrow(test_data), replace = TRUE)
  i <- FALSE
  expect_no_condition({
    fit2 <- multinom_weightit(Y_M ~ A * (X1 + X2 + X3 + X4 + X5),
                              data = test_data, weightit = W0,
                              control = if (i) list(stop("bad error")) else list())
  })

  i <- TRUE
  expect_no_condition({
    update(fit2, vcov = "HC0")
  })

  expect_error({
    update(fit2)
  }, "bad error")

  expect_no_condition({
    fit2c <- update(fit2, cluster = ~clus)
  })

  expect_equal(vcov(fit2c),
               vcov(update(fit0, cluster = ~clus)),
               tolerance = eps)

  expect_not_equal(vcov(fit2c),
                   vcov(fit0))

  #Updating s.weights updates weightit
  expect_no_condition({
    W1 <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                   data = test_data, s.weights = "SW")
  })

  expect_no_condition({
    fit1 <- multinom_weightit(Y_M ~ A * (X1 + X2 + X3 + X4 + X5),
                              data = test_data, weightit = W1)
  })

  expect_no_condition({
    fit1s <- multinom_weightit(Y_M ~ A * (X1 + X2 + X3 + X4 + X5),
                               data = test_data,
                               weightit = update(object = W0, s.weights = "SW"))
  })

  expect_no_condition({
    fit0s <- update(fit0, s.weights = "SW")
  })

  expect_equal(fit0s$coefficients,
               fit1$coefficients,
               ignore_attr = c("names", "class"),
               tolerance = eps)

  expect_equal(fit0s$coefficients,
               fit1s$coefficients,
               ignore_attr = c("names", "class"),
               tolerance = eps)

  expect_equal(fit0s,
               fit1s,
               tolerance = eps)

  #Mimic bootstrapping
  boot_dat <- test_data[sample(nrow(test_data), replace = TRUE),]

  expect_no_condition({
    W1 <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                   data = boot_dat)
  })

  expect_no_condition({
    fit1 <- multinom_weightit(Y_M ~ A * (X1 + X2 + X3 + X4 + X5),
                              data = boot_dat, weightit = W1, vcov = "none")
  })

  expect_equal(update(fit0, data = boot_dat, vcov = "none")$coefficients,
               fit1$coefficients, tolerance = eps,
               ignore_attr = "class",
               ignore_formula_env = TRUE)
})

test_that("update.ordinal_weightit() works with weightit", {
  skip_if(!capabilities("long.double"))
  eps <- if (capabilities("long.double")) 1e-5 else 1e-1

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))
  test_data$Y_O <- with(test_data, factor(findInterval(Y_C, quantile(Y_C, seq(0, 1, length = 5)),
                                                       all.inside = TRUE)))

  expect_no_condition({
    W0 <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                   data = test_data)
  })

  expect_no_condition({
    fit0 <- ordinal_weightit(Y_O ~ A * (X1 + X2 + X3 + X4 + X5),
                             data = test_data, weightit = W0)
  })

  #Updating weightit
  expect_no_condition({
    fit1 <- ordinal_weightit(Y_O ~ A * (X1 + X2 + X3 + X4 + X5),
                             data = test_data)
  })

  expect_equal(update(fit0, weightit = NULL),
               fit1, tolerance = eps,
               ignore_attr = "class",
               ignore_formula_env = TRUE)

  #Updating dataset
  sub <- test_data$X3 > 0
  test_data_s <- test_data[sub,]

  expect_no_condition({
    W1 <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                   data = test_data_s)
  })

  expect_no_condition({
    fit1 <- ordinal_weightit(Y_O ~ A * (X1 + X2 + X3 + X4 + X5),
                             data = test_data_s, weightit = W1)
  })

  expect_equal(update(fit0, data = test_data_s)$coefficients,
               fit1$coefficients, tolerance = eps,
               ignore_attr = "class",
               ignore_formula_env = TRUE)

  #Updating vcov only
  expect_no_condition({
    fit1 <- ordinal_weightit(Y_O ~ A * (X1 + X2 + X3 + X4 + X5),
                             data = test_data, weightit = W0,
                             vcov = "HC0")
  })

  expect_equal(update(fit0, vcov = "HC0"),
               fit1, tolerance = eps,
               ignore_attr = "class",
               ignore_formula_env = TRUE)

  expect_not_equal(vcov(fit0), vcov(fit1))

  #Model should not be refit when only vcov is changed
  clus <- sample(1:50, nrow(test_data), replace = TRUE)
  i <- FALSE
  expect_no_condition({
    fit2 <- ordinal_weightit(Y_O ~ A * (X1 + X2 + X3 + X4 + X5),
                             data = test_data, weightit = W0,
                             control = if (i) list(stop("bad error")) else list())
  })

  i <- TRUE
  expect_no_condition({
    update(fit2, vcov = "HC0")
  })

  expect_error({
    update(fit2)
  }, "bad error")

  expect_no_condition({
    fit2c <- update(fit2, cluster = ~clus)
  })

  expect_equal(vcov(fit2c),
               vcov(update(fit0, cluster = ~clus)),
               tolerance = eps)

  expect_not_equal(vcov(fit2c),
                   vcov(fit0))

  #Updating s.weights updates weightit
  expect_no_condition({
    W1 <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                   data = test_data, s.weights = "SW")
  })

  expect_no_condition({
    fit1 <- ordinal_weightit(Y_O ~ A * (X1 + X2 + X3 + X4 + X5),
                             data = test_data, weightit = W1)
  })

  expect_no_condition({
    fit1s <- ordinal_weightit(Y_O ~ A * (X1 + X2 + X3 + X4 + X5),
                              data = test_data,
                              weightit = update(object = W0, s.weights = "SW"))
  })

  expect_no_condition({
    fit0s <- update(fit0, s.weights = "SW")
  })

  expect_equal(fit0s$coefficients,
               fit1$coefficients,
               ignore_attr = c("names", "class"),
               tolerance = eps)

  expect_equal(fit0s$coefficients,
               fit1s$coefficients,
               ignore_attr = c("names", "class"),
               tolerance = eps)

  expect_equal(fit0s,
               fit1s,
               tolerance = eps)

  #Mimic bootstrapping
  boot_dat <- test_data[sample(nrow(test_data), replace = TRUE),]

  expect_no_condition({
    W1 <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                   data = boot_dat)
  })

  expect_no_condition({
    fit1 <- ordinal_weightit(Y_O ~ A * (X1 + X2 + X3 + X4 + X5),
                             data = boot_dat, weightit = W1, vcov = "none")
  })

  expect_equal(update(fit0, data = boot_dat, vcov = "none")$coefficients,
               fit1$coefficients, tolerance = eps,
               ignore_attr = "class",
               ignore_formula_env = TRUE)
})

test_that("update.coxph_weightit() works with weightit", {
  skip_if(!capabilities("long.double"))
  skip_if_not_installed("survival")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-1

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  expect_no_condition({
    W0 <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                   data = test_data)
  })

  expect_no_condition({
    fit0 <- coxph_weightit(survival::Surv(Y_S) ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                         data = test_data, weightit = W0)
  })

  #Updating weightit
  expect_no_condition({
    fit1 <- coxph_weightit(survival::Surv(Y_S) ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                         data = test_data)
  })

  expect_equal(update(fit0, weightit = NULL),
               fit1, tolerance = eps,
               ignore_attr = "class",
               ignore_formula_env = TRUE)

  #Updating dataset
  sub <- test_data$X3 > 0
  test_data_s <- test_data[sub,]

  expect_no_condition({
    W1 <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                   data = test_data_s)
  })

  expect_no_condition({
    fit1 <- coxph_weightit(survival::Surv(Y_S) ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                         data = test_data_s, weightit = W1)
  })

  expect_equal(update(fit0, data = test_data_s)$coefficients,
               fit1$coefficients, tolerance = eps,
               ignore_attr = "class",
               ignore_formula_env = TRUE)

  #Updating vcov only
  clus <- sample(1:50, nrow(test_data), replace = TRUE)

  expect_no_condition({
    fit1 <- coxph_weightit(survival::Surv(Y_S) ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                         data = test_data, weightit = W0,
                         cluster = ~clus)
  })

  expect_equal(update(fit0, cluster = ~clus),
               fit1, tolerance = eps,
               ignore_attr = "class",
               ignore_formula_env = TRUE)

  expect_not_equal(vcov(fit0), vcov(fit1))

  #Model should not be refit when only vcov is changed
  i <- FALSE
  expect_no_condition({
    fit2 <- coxph_weightit(survival::Surv(Y_S) ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                         data = test_data, weightit = W0,
                         control = if (i) list(stop("bad error")) else list())
  })

  i <- TRUE
  expect_error({
    update(fit2, ties = "breslow")
  }, "bad error")

  expect_error({
    update(fit2)
  }, "bad error")

  expect_no_condition({
    fit2c <- update(fit2, cluster = ~clus)
  })

  expect_equal(vcov(fit2c),
               vcov(update(fit0, cluster = ~clus)),
               tolerance = eps)

  expect_not_equal(vcov(fit2c),
                   vcov(fit0))

  #Updating s.weights updates weightit
  expect_no_condition({
    W1 <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                   data = test_data, s.weights = "SW")
  })

  expect_no_condition({
    fit1 <- coxph_weightit(survival::Surv(Y_S) ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                         data = test_data, weightit = W1)
  })

  expect_no_condition({
    fit1s <- coxph_weightit(survival::Surv(Y_S) ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                          data = test_data,
                          weightit = update(object = W0, s.weights = "SW"))
  })

  expect_no_condition({
    fit0s <- update(fit0, s.weights = "SW")
  })

  expect_equal(fit0s$coefficients,
               fit1$coefficients,
               ignore_attr = c("names", "class"),
               tolerance = eps)

  expect_equal(fit0s$coefficients,
               fit1s$coefficients,
               ignore_attr = c("names", "class"),
               tolerance = eps)

  expect_equal(fit0s,
               fit1s,
               tolerance = eps)

  #Mimic bootstrapping
  boot_dat <- test_data[sample(nrow(test_data), replace = TRUE),]

  expect_no_condition({
    W1 <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                   data = boot_dat)
  })

  expect_no_condition({
    fit1 <- coxph_weightit(survival::Surv(Y_S) ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                         data = boot_dat, weightit = W1)
  })

  expect_equal(update(fit0, data = boot_dat, vcov = "none")$coefficients,
               fit1$coefficients, tolerance = eps,
               ignore_attr = "class",
               ignore_formula_env = TRUE)
})
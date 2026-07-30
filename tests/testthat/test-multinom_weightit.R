test_that("No weights", {
  skip_on_cran()
  skip_if_not_installed("sandwich")
  skip_if_not_installed("mlogit")
  skip_if_not_installed("dfidx")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))
  test_data$Y_M <- factor(test_data$Y_O, ordered = FALSE)

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
  skip_on_cran()
  skip_if_not_installed("mlogit")
  skip_if_not_installed("dfidx")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))
  test_data$Y_M <- factor(test_data$Y_O, ordered = FALSE)

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

test_that("Bootstrap vcov (BS, FWB)", {
  skip_on_cran()
  skip_if_not_installed("fwb")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))
  test_data$Y_M <- factor(test_data$Y_O, ordered = FALSE)

  #Small/reduced formula for speed
  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "glm", estimand = "ATE",
                  include.obj = TRUE)
  })

  expect_no_condition({
    fit0 <- multinom_weightit(Y_M ~ A + X1,
                              data = test_data, weightit = W)
  })

  set.seed(123)
  expect_no_condition({
    fit_bs <- multinom_weightit(Y_M ~ A + X1,
                                data = test_data, weightit = W,
                                vcov = "BS", R = 30)
  })

  expect_equal(coef(fit0), coef(fit_bs), tolerance = eps)
  expect_not_equal(vcov(fit0), vcov(fit_bs), tolerance = eps)

  set.seed(123)
  expect_no_condition({
    fit_fwb <- multinom_weightit(Y_M ~ A + X1,
                                 data = test_data, weightit = W,
                                 vcov = "FWB", R = 30)
  })

  expect_equal(coef(fit0), coef(fit_fwb), tolerance = eps)
  expect_not_equal(vcov(fit0), vcov(fit_fwb), tolerance = eps)
})

test_that("vcov = 'const' returns the model-based variance", {
  skip_on_cran()
  skip_if_not_installed("mlogit")
  skip_if_not_installed("dfidx")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))
  test_data$Y_M <- factor(test_data$Y_O, ordered = FALSE)

  #The model-based variance is the inverse of the Hessian, which is what
  #mlogit() reports
  expect_no_condition({
    fit <- multinom_weightit(Y_M ~ A * (X1 + X2), data = test_data,
                             vcov = "const")
  })

  expect_identical(fit$vcov_type, "const")
  expect_identical(colnames(vcov(fit)), names(coef(fit)))

  fit_g <- mlogit::mlogit(Y_M ~ 0 | A * (X1 + X2),
                          data = dfidx::dfidx(test_data, choice = "Y_M",
                                              shape = "wide"))

  ind <- unlist(split(seq_along(coef(fit)),
                      rep(seq_len(nlevels(test_data$Y_M) - 1),
                          length(coef(fit)) / (nlevels(test_data$Y_M) - 1))))

  expect_equal(unname(vcov(fit)), unname(vcov(fit_g)[ind, ind]),
               tolerance = eps)

  #Requesting it after the fact from an object fit with a different vcov gives
  #the same matrix, and differs from that object's own (HC0) variance
  expect_no_condition({
    fit0 <- multinom_weightit(Y_M ~ A * (X1 + X2), data = test_data)
  })

  expect_no_condition({
    V_const <- vcov(fit0, vcov = "const")
  })

  expect_equal(V_const, vcov(fit), ignore_attr = TRUE, tolerance = eps)
  expect_not_equal(vcov(fit0), V_const)

  #Point estimates are unaffected by the choice of vcov
  expect_equal(coef(fit), coef(fit0), tolerance = eps)

  #With a 2-level outcome, a multinomial logit is a reparametrization of a
  #binomial logit, so the model-based variance matches glm_weightit()'s
  test_data$Y_B2 <- factor(test_data$Y_B, labels = c("no", "yes"))

  expect_no_condition({
    fit_b <- multinom_weightit(Y_B2 ~ A + X1 + X2 + X3, data = test_data,
                               vcov = "const")
  })

  fit_glm <- glm_weightit(Y_B ~ A + X1 + X2 + X3, data = test_data,
                          family = binomial, vcov = "const")

  expect_equal(unname(vcov(fit_b)), unname(vcov(fit_glm)), tolerance = eps)

  #Weighted
  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "glm", estimand = "ATE",
                  include.obj = TRUE)
  })

  #`vcov = "const"` ignores estimation of the weights, so it warns
  expect_warning({
    fit_b_w <- multinom_weightit(Y_B2 ~ A + X1 + X2 + X3, data = test_data,
                                 weightit = W, vcov = "const")
  }, "const.*should not be used", ignore.case = TRUE)

  expect_identical(fit_b_w$vcov_type, "const")

  expect_warning({
    fit_glm_w <- glm_weightit(Y_B ~ A + X1 + X2 + X3, data = test_data,
                              weightit = W, family = binomial, vcov = "const")
  }, "const.*should not be used", ignore.case = TRUE)

  expect_equal(unname(vcov(fit_b_w)), unname(vcov(fit_glm_w)), tolerance = eps)
})

test_that("Binary (2-level factor) outcome reduces to binomial logit", {
  skip_on_cran()

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))
  test_data$Y_B2 <- factor(test_data$Y_B, labels = c("no", "yes"))

  #No weightit object -- coefficients and HC0 vcov should match glm_weightit()
  #with a binomial family exactly (up to optimizer tolerance), since a
  #2-category multinomial logit is a reparametrization of binomial logit.
  expect_no_condition({
    fit_m <- multinom_weightit(Y_B2 ~ A + X1 + X2 + X3, data = test_data, vcov = "HC0")
  })

  fit_glm <- glm_weightit(Y_B ~ A + X1 + X2 + X3, data = test_data, family = binomial)

  expect_equal(unname(coef(fit_m)), unname(coef(fit_glm)), tolerance = eps)
  expect_equal(unname(vcov(fit_m)), unname(vcov(fit_glm)), tolerance = eps)

  #With a weightit object
  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "glm", estimand = "ATE",
                  include.obj = TRUE)
  })

  expect_no_condition({
    fit_m_w <- multinom_weightit(Y_B2 ~ A + X1 + X2 + X3, data = test_data,
                                 weightit = W, vcov = "HC0")
  })

  fit_glm_w <- glm_weightit(Y_B ~ A + X1 + X2 + X3, data = test_data,
                            weightit = W, family = binomial, vcov = "HC0")

  expect_equal(unname(coef(fit_m_w)), unname(coef(fit_glm_w)), tolerance = eps)
  expect_equal(unname(vcov(fit_m_w)), unname(vcov(fit_glm_w)), tolerance = eps)
})

test_that("collinear covariates give NA for the aliased coefficients", {
  skip_on_cran()

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))
  test_data$Y_M <- factor(test_data$Y_O, ordered = FALSE)
  test_data$X1b <- test_data$X1

  expect_no_condition({
    fit <- multinom_weightit(Y_M ~ A + X1 + X1b + X2, data = test_data)
  })

  #One aliased coefficient per outcome equation
  nq <- nlevels(test_data$Y_M) - 1L
  expect_identical(sum(is.na(coef(fit))), nq)
  expect_identical(names(coef(fit))[is.na(coef(fit))],
                   paste0(levels(test_data$Y_M)[-1L], "~X1b"))

  V_est <- vcov(fit, complete = FALSE)
  expect_identical(dim(V_est), rep(4L * nq, 2L))
  expect_false(anyNA(V_est))

  #`complete = TRUE` expands to the full set with NAs
  expect_identical(dim(vcov(fit)), rep(5L * nq, 2L))

  #summary() reports the estimable coefficients with standard errors aligned to
  #them (they were previously all NA)
  expect_no_condition(s <- summary(fit))
  expect_identical(nrow(s$coefficients), 4L * nq)
  expect_false(anyNA(s$coefficients))
  expect_equal(unname(s$coefficients[, "Std. Error"]),
               unname(sqrt(diag(V_est))), tolerance = eps)

  #With a 2-level outcome this reduces to a binomial logit, so the estimable
  #coefficients and standard errors match glm_weightit()'s
  test_data$Y_B2 <- factor(test_data$Y_B, labels = c("no", "yes"))

  expect_no_condition({
    fit_b <- multinom_weightit(Y_B2 ~ A + X1 + X1b + X2, data = test_data)
  })

  fit_glm <- glm_weightit(Y_B ~ A + X1 + X1b + X2, data = test_data,
                          family = binomial)

  expect_equal(unname(coef(fit_b)), unname(coef(fit_glm)), tolerance = eps)
  expect_equal(unname(vcov(fit_b, complete = FALSE)),
               unname(vcov(fit_glm, complete = FALSE)), tolerance = eps)
})

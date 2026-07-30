test_that("No weights", {
  skip_on_cran()
  skip_if_not_installed("sandwich")
  skip_if_not_installed("brglm2")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

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

  expect_equal(coef(fit0), coef(fit), tolerance = eps)
  expect_equal(vcov(fit0), vcov(fit), tolerance = eps)

  fit_g <- glm(Y_B ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
               data = test_data, family = binomial)
  expect_equal(coef(fit0), coef(fit_g), tolerance = eps)
  expect_equal(vcov(fit0), sandwich::sandwich(fit_g),
               tolerance = eps)

  #Offset
  set.seed(123)
  off <- runif(nrow(test_data))

  expect_no_condition({
    fit <- glm_weightit(Y_B ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9) + offset(off),
                        data = test_data, family = binomial)
  })

  expect_not_equal(coef(fit0), coef(fit), tolerance = eps)

  fit_g <- glm(Y_B ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9) + offset(off),
               data = test_data, family = binomial)
  expect_equal(coef(fit), coef(fit_g), tolerance = eps)
  expect_equal(vcov(fit), sandwich::sandwich(fit_g),
               tolerance = eps)

  #Cluster-robust SEs
  clus <- sample(1:50, nrow(test_data), replace = TRUE)

  expect_no_condition({
    fit <- glm_weightit(Y_B ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                        data = test_data, family = binomial, cluster = clus)
  })

  expect_equal(coef(fit0), coef(fit), tolerance = eps)

  fit_g <- glm(Y_B ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
               data = test_data, family = binomial)
  expect_equal(coef(fit), coef(fit_g), tolerance = eps)
  expect_equal(vcov(fit), sandwich::vcovCL(fit_g, cluster = clus, type = "HC0"),
               tolerance = eps)

  #BR
  expect_no_condition({
    fit <- glm_weightit(Y_B ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                        data = test_data, family = binomial("probit"), br = TRUE)
  })

  fit_g <- glm(Y_B ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
               data = test_data, family = binomial("probit"),
               method = brglm2::brglmFit)

  expect_equal(coef(fit), coef(fit_g), tolerance = eps)
  expect_equal(vcov(fit), sandwich::sandwich(fit_g),
               tolerance = eps * 1e3)

  expect_equal(vcov(fit, cluster = clus),
               sandwich::vcovCL(fit_g, cluster = clus, type = "HC0"),
               tolerance = eps * 1e3)

  #Test error for missingness
  expect_error({
    fit0 <- glm_weightit(Y_B ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                         data = transform(test_data, Y_B = c(NA, Y_B[-1L])),
                         family = binomial)
  }, "missing values", ignore.case = TRUE)

  expect_error({
    fit0 <- glm_weightit(Y_B ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                         data = transform(test_data, X1 = c(NA, X1[-1L])),
                         family = binomial)
  }, "missing values", ignore.case = TRUE)

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
  skip_if_not_installed("fwb")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                   data = test_data, method = "glm", estimand = "ATE",
                   include.obj = TRUE)
  })

  expect_no_condition({
    fit0 <- glm_weightit(Y_C ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                        data = test_data, weightit = W)
  })

  #M-estimation for glm
  expect_no_condition({
    fit <- glm_weightit(Y_C ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                        data = test_data, weightit = W,
                        vcov = "asympt")
  })

  expect_equal(coef(fit0), coef(fit), tolerance = eps)
  expect_equal(vcov(fit0), vcov(fit), tolerance = eps)

  expect_no_condition({
    fit <- glm_weightit(Y_C ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                        data = test_data, weightit = W,
                        vcov = "HC0")
  })

  expect_equal(coef(fit0), coef(fit), tolerance = eps)
  expect_not_equal(vcov(fit0), vcov(fit), tolerance = eps)

  set.seed(123)
  expect_no_condition({
    fit <- glm_weightit(Y_C ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                        data = test_data, weightit = W,
                        vcov = "FWB", R = 50)
  })

  expect_equal(coef(fit0), coef(fit), tolerance = eps)
  expect_not_equal(vcov(fit0), vcov(fit), tolerance = eps)

  set.seed(123)
  expect_no_condition({
    fit_ <- glm_weightit(Y_C ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                        data = test_data, weightit = W,
                        vcov = "FWB", R = 50,
                        fwb.args = list(wtype = "mammen"))
  })

  expect_equal(coef(fit), coef(fit_), tolerance = eps)
  expect_not_equal(vcov(fit), vcov(fit_), tolerance = eps)

  expect_no_condition({
    fit <- glm_weightit(Y_C ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                        data = test_data, weightit = W,
                        vcov = "BS", R = 50)
  })

  expect_equal(coef(fit0), coef(fit), tolerance = eps)
  expect_not_equal(vcov(fit0), vcov(fit), tolerance = eps)

  #Test error for missingness
  expect_error({
    fit0 <- glm_weightit(Y_B ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                         data = transform(test_data, Y_B = c(NA, Y_B[-1L])),
                         family = binomial)
  }, "missing values", ignore.case = TRUE)

  expect_error({
    fit0 <- glm_weightit(Y_B ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                         data = transform(test_data, X1 = c(NA, X1[-1L])),
                         family = binomial)
  }, "missing values", ignore.case = TRUE)

  #Test using sandwich functions
  expect_no_condition({
    fit0 <- glm_weightit(Y_B ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                         data = test_data, weightit = W, family = binomial)
  })

  expect_equal(vcov(fit0),
               sandwich::sandwich(fit0),
               tolerance = eps)

  expect_equal(vcov(fit0, type = "HC0"),
               sandwich::sandwich(fit0, asympt = FALSE),
               tolerance = eps)
})

test_that("Gaussian and Poisson families", {
  skip_on_cran()
  skip_if_not_installed("sandwich")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  #Count outcome constructed locally for this test only (not part of the fixture)
  set.seed(123)
  test_data$Y_P <- rpois(nrow(test_data), lambda = exp(.3 + .05 * test_data$X1))

  ##No weightit object -- Gaussian
  expect_no_condition({
    fit0 <- glm_weightit(Y_C ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                         data = test_data, family = gaussian())
  })

  fit_g <- glm(Y_C ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
              data = test_data, family = gaussian())

  expect_equal(coef(fit0), coef(fit_g), tolerance = eps)
  expect_equal(vcov(fit0), sandwich::sandwich(fit_g), tolerance = eps)

  ##No weightit object -- Poisson
  expect_no_condition({
    fit0_p <- glm_weightit(Y_P ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                           data = test_data, family = poisson())
  })

  fit_g_p <- glm(Y_P ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                data = test_data, family = poisson())

  expect_equal(coef(fit0_p), coef(fit_g_p), tolerance = eps)
  expect_equal(vcov(fit0_p), sandwich::sandwich(fit_g_p), tolerance = eps)

  ##With a weightit object
  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "glm", estimand = "ATE",
                  include.obj = TRUE)
  })

  #Gaussian, M-estimation vs. HC0 (should differ but agree on coefficients)
  expect_no_condition({
    fit <- glm_weightit(Y_C ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                        data = test_data, weightit = W, family = gaussian())
  })

  expect_no_condition({
    fit_hc0 <- glm_weightit(Y_C ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                            data = test_data, weightit = W, family = gaussian(),
                            vcov = "HC0")
  })

  expect_equal(coef(fit), coef(fit_hc0), tolerance = eps)
  expect_not_equal(vcov(fit), vcov(fit_hc0), tolerance = eps)
  expect_equal(vcov(fit_hc0), sandwich::sandwich(fit_hc0, asympt = FALSE), tolerance = eps)

  #Poisson, M-estimation vs. HC0
  expect_no_condition({
    fit_p <- glm_weightit(Y_P ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                          data = test_data, weightit = W, family = poisson())
  })

  expect_no_condition({
    fit_p_hc0 <- glm_weightit(Y_P ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                              data = test_data, weightit = W, family = poisson(),
                              vcov = "HC0")
  })

  expect_equal(coef(fit_p), coef(fit_p_hc0), tolerance = eps)
  expect_not_equal(vcov(fit_p), vcov(fit_p_hc0), tolerance = eps)
  expect_equal(vcov(fit_p_hc0), sandwich::sandwich(fit_p_hc0, asympt = FALSE), tolerance = eps)
})

test_that("family = 'multinomial' is rejected", {
  skip_on_cran()

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))
  test_data$Y_M <- factor(test_data$Y_O, ordered = FALSE)

  #Previously this warned and routed the input to multinom_weightit()
  expect_error(glm_weightit(Y_M ~ A + X1 + X2, data = test_data,
                            family = "multinomial"),
               "not allowed")

  expect_error(glm_weightit(Y_M ~ A + X1 + X2, data = test_data,
                            family = "multinomial"),
               "multinom_weightit")
})

test_that("collinear covariates give NA for the aliased coefficients", {
  skip_on_cran()
  skip_if_not_installed("sandwich")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  #An exact duplicate of an existing covariate
  test_data$X1b <- test_data$X1

  expect_no_condition({
    fit <- glm_weightit(Y_C ~ A + X1 + X1b + X2, data = test_data)
  })

  fit_g <- glm(Y_C ~ A + X1 + X1b + X2, data = test_data)

  #stats::glm() also returns NA for the redundant coefficient
  expect_equal(coef(fit), coef(fit_g), tolerance = eps)
  expect_identical(is.na(coef(fit)),
                   c("(Intercept)" = FALSE, A = FALSE, X1 = FALSE,
                     X1b = TRUE, X2 = FALSE))

  #`complete = TRUE` (the default) expands the variance to the full set of
  #coefficients with NAs; `complete = FALSE` covers the estimable ones only
  V <- vcov(fit)
  expect_identical(dim(V), c(5L, 5L))
  expect_true(all(is.na(V["X1b", ])))

  V_est <- vcov(fit, complete = FALSE)
  expect_identical(dim(V_est), c(4L, 4L))
  expect_false(anyNA(V_est))
  expect_identical(colnames(V_est), c("(Intercept)", "A", "X1", "X2"))

  #The default HC0 variance matches the sandwich variance of the glm() fit
  expect_equal(unname(V_est), unname(sandwich::sandwich(fit_g)), tolerance = eps)

  #print() and summary() report only the estimable coefficients, with standard
  #errors aligned to them
  expect_no_condition(invisible(capture.output(print(fit))))
  expect_no_condition(s <- summary(fit))
  expect_identical(rownames(s$coefficients), c("(Intercept)", "A", "X1", "X2"))
  expect_false(anyNA(s$coefficients))
  expect_equal(unname(s$coefficients[, "Std. Error"]),
               unname(sqrt(diag(sandwich::sandwich(fit_g)))), tolerance = eps)

  #sandwich's estfun()/bread() are also restricted to the estimable coefficients
  expect_identical(ncol(sandwich::estfun(fit)), 4L)
  expect_identical(dim(sandwich::bread(fit)), c(4L, 4L))
  expect_equal(unname(sandwich::sandwich(fit)), unname(V_est), tolerance = eps)

  #vcov = "const" is the model-based variance over the estimable coefficients
  expect_no_condition({
    fit_const <- glm_weightit(Y_C ~ A + X1 + X1b + X2, data = test_data,
                              vcov = "const")
  })

  expect_equal(unname(vcov(fit_const, complete = FALSE)),
               unname(vcov(fit_g, complete = FALSE)), tolerance = eps)

  #lm_weightit() behaves the same way
  expect_no_condition({
    fit_lm <- lm_weightit(Y_C ~ A + X1 + X1b + X2, data = test_data)
  })

  expect_equal(coef(fit_lm), coef(fit), tolerance = eps)
  expect_equal(vcov(fit_lm, complete = FALSE), V_est, tolerance = eps)

  #A redundant factor: two aliased coefficients at once
  test_data$G <- factor(rep(c("a", "b", "c"), length.out = nrow(test_data)))
  test_data$Gb <- test_data$G

  expect_no_condition({
    fit_f <- glm_weightit(Y_B ~ A + G + Gb, data = test_data, family = binomial)
  })

  fit_f_g <- glm(Y_B ~ A + G + Gb, data = test_data, family = binomial)

  expect_equal(coef(fit_f), coef(fit_f_g), tolerance = eps)
  expect_identical(sum(is.na(coef(fit_f))), 2L)
  expect_equal(unname(vcov(fit_f, complete = FALSE)),
               unname(sandwich::sandwich(fit_f_g)), tolerance = eps)
})

test_that("dropping a collinear covariate doesn't change the estimable results", {
  skip_on_cran()

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))
  test_data$X1b <- test_data$X1

  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3, data = test_data, method = "glm",
                  estimand = "ATE")
  })

  #A perfectly collinear column carries no information, so fitting with it must
  #give exactly what fitting without it does. This covers the M-estimation
  #(`asympt`) path, where the aliased columns must be dropped from the model
  #matrix consistently with the NA coefficients.
  for (v in c("const", "HC0")) {
    fit_full <- glm_weightit(Y_B ~ A + X1 + X1b + X2, data = test_data,
                             family = binomial, vcov = v)
    fit_red <- glm_weightit(Y_B ~ A + X1 + X2, data = test_data,
                            family = binomial, vcov = v)

    expect_equal(coef(fit_full)[!is.na(coef(fit_full))], coef(fit_red),
                 tolerance = eps)
    expect_equal(vcov(fit_full, complete = FALSE), vcov(fit_red),
                 ignore_attr = TRUE, tolerance = eps)
  }

  for (v in c("asympt", "HC0")) {
    fit_full <- glm_weightit(Y_B ~ A + X1 + X1b + X2, data = test_data,
                             family = binomial, weightit = W, vcov = v)
    fit_red <- glm_weightit(Y_B ~ A + X1 + X2, data = test_data,
                            family = binomial, weightit = W, vcov = v)

    expect_identical(fit_full$vcov_type, v)
    expect_equal(coef(fit_full)[!is.na(coef(fit_full))], coef(fit_red),
                 tolerance = eps)
    expect_equal(vcov(fit_full, complete = FALSE), vcov(fit_red),
                 ignore_attr = TRUE, tolerance = eps)
  }

  #M-estimation still yields smaller SEs than treating the weights as fixed
  fit_asympt <- glm_weightit(Y_B ~ A + X1 + X1b + X2, data = test_data,
                             family = binomial, weightit = W, vcov = "asympt")
  fit_hc0 <- glm_weightit(Y_B ~ A + X1 + X1b + X2, data = test_data,
                          family = binomial, weightit = W, vcov = "HC0")

  expect_true(all(diag(vcov(fit_asympt, complete = FALSE)) <=
                    diag(vcov(fit_hc0, complete = FALSE))))
})

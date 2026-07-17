# Y_S in the fixture is an uncensored (fully observed) event time. Build a
# right-censored version once per test_that() call (tests should remain
# self-sufficient) by censoring at the 80th percentile.
.censor_Y_S <- function(data) {
  cutoff <- quantile(data$Y_S, .8)
  data$event <- as.numeric(data$Y_S < cutoff)
  data$time <- pmin(data$Y_S, cutoff)
  data
}

test_that("Unweighted fit matches survival::coxph()", {
  skip_on_cran()
  skip_if_not_installed("survival")
  skip_if_not_installed("sandwich")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))
  test_data <- .censor_Y_S(test_data)

  expect_no_condition({
    fit0 <- coxph_weightit(survival::Surv(time, event) ~ A + X1 + X2 + X3,
                           data = test_data)
  })

  fit_ref <- survival::coxph(survival::Surv(time, event) ~ A + X1 + X2 + X3,
                             data = test_data, robust = TRUE)

  expect_equal(unname(coef(fit0)), unname(coef(fit_ref)), tolerance = eps)

  # No `weightit` object supplied -> default vcov is the HC0 sandwich variance
  # (per the "Details" section of ?coxph_weightit)
  expect_identical(fit0$vcov_type, "HC0")

  expect_equal(unname(vcov(fit0)), unname(vcov(fit_ref)),
               tolerance = eps)

  expect_equal(unname(vcov(fit0)), unname(sandwich::sandwich(fit_ref)),
               tolerance = eps)

  # Also works with an uncensored Surv(time) (single-argument) response, as
  # used elsewhere in the package's tests (e.g., test-update.R)
  expect_no_condition({
    fit1 <- coxph_weightit(survival::Surv(Y_S) ~ A + X1 + X2 + X3,
                           data = test_data)
  })

  fit1_ref <- survival::coxph(survival::Surv(Y_S) ~ A + X1 + X2 + X3,
                              data = test_data)

  expect_equal(unname(coef(fit1)), unname(coef(fit1_ref)), tolerance = eps)
})

test_that("Weighted fit (weightit object) matches survival::coxph() with weights", {
  skip_on_cran()
  skip_if_not_installed("survival")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))
  test_data <- .censor_Y_S(test_data)

  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "glm", estimand = "ATE")
  })

  expect_no_condition({
    fit <- coxph_weightit(survival::Surv(time, event) ~ A + X1 + X2 + X3,
                          data = test_data, weightit = W, vcov = "HC0")
  })

  fit_ref <- survival::coxph(survival::Surv(time, event) ~ A + X1 + X2 + X3,
                             data = test_data, weights = W$weights, robust = TRUE)

  expect_equal(unname(coef(fit)), unname(coef(fit_ref)), tolerance = eps)

  # vcov = "HC0" treats the weights as fixed, which matches coxph()'s own
  # robust ("dfbeta") sandwich variance when given the same case weights.
  expect_equal(unname(vcov(fit)), unname(vcov(fit_ref)), tolerance = eps)
  expect_equal(unname(vcov(fit)), unname(sandwich::sandwich(fit_ref)), tolerance = eps)
})

test_that("weightit object with M-estimation parts defaults to vcov = 'asympt'", {
  skip_on_cran()
  skip_if_not_installed("survival")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))
  test_data <- .censor_Y_S(test_data)

  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "glm", estimand = "ATE")
  })

  expect_true(is_not_null(attr(W, "Mparts")))

  expect_no_condition({
    fit_default <- coxph_weightit(survival::Surv(time, event) ~ A + X1 + X2 + X3,
                                  data = test_data, weightit = W)
  })

  expect_identical(fit_default$vcov_type, "asympt")

  fit_asympt <- coxph_weightit(survival::Surv(time, event) ~ A + X1 + X2 + X3,
                               data = test_data, weightit = W, vcov = "asympt")

  expect_equal(vcov(fit_default), vcov(fit_asympt), tolerance = eps)

  fit_hc0 <- coxph_weightit(survival::Surv(time, event) ~ A + X1 + X2 + X3,
                            data = test_data, weightit = W, vcov = "HC0")

  # Same point estimates regardless of vcov type
  expect_equal(coef(fit_asympt), coef(fit_hc0), tolerance = eps)

  # M-estimation (accounting for estimation of the weights) yields different
  # (here, uniformly smaller) SEs than treating the weights as fixed, as
  # described in the "Details" section of ?coxph_weightit
  expect_not_equal(vcov(fit_asympt), vcov(fit_hc0))
  expect_true(all(diag(vcov(fit_asympt)) <= diag(vcov(fit_hc0))))
})

test_that("vcov = 'none' produces no variance matrix", {
  skip_on_cran()
  skip_if_not_installed("survival")

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))
  test_data <- .censor_Y_S(test_data)

  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5, data = test_data,
                  method = "glm", estimand = "ATE")
  })

  expect_no_condition({
    fit_none <- coxph_weightit(survival::Surv(time, event) ~ A + X1 + X2,
                               data = test_data, weightit = W, vcov = "none")
  })

  expect_identical(fit_none$vcov_type, "none")

  expect_warning(v <- vcov(fit_none), "vcov.*none", ignore.case = TRUE)
  expect_null(v)
})

test_that("cluster argument applies a small-sample-corrected cluster-robust variance", {
  skip_on_cran()
  skip_if_not_installed("survival")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))
  test_data <- .censor_Y_S(test_data)
  test_data$clus <- rep(1:20, length.out = nrow(test_data))

  g <- length(unique(test_data$clus))
  adj <- g / (g - 1)

  # Unweighted
  expect_no_condition({
    fit0 <- coxph_weightit(survival::Surv(time, event) ~ A + X1 + X2 + X3,
                           data = test_data, cluster = ~clus)
  })

  fit0_ref <- survival::coxph(survival::Surv(time, event) ~ A + X1 + X2 + X3,
                              data = test_data, cluster = clus)

  expect_equal(unname(vcov(fit0)), unname(fit0_ref$var) * adj, tolerance = eps)

  # Weighted, vcov = "HC0" (weights treated as fixed)
  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "glm", estimand = "ATE")
  })

  expect_no_condition({
    fit <- coxph_weightit(survival::Surv(time, event) ~ A + X1 + X2 + X3,
                          data = test_data, weightit = W, vcov = "HC0",
                          cluster = ~clus)
  })

  fit_ref <- survival::coxph(survival::Surv(time, event) ~ A + X1 + X2 + X3,
                             data = test_data, weights = W$weights, cluster = clus)

  expect_equal(unname(vcov(fit)), unname(fit_ref$var) * adj, tolerance = eps)

  # Cluster attribute is stored even when not used in the variance calculation
  # (a warning fires because cluster is ignored for vcov = "none")
  expect_warning({
    fit_stored <- coxph_weightit(survival::Surv(time, event) ~ A + X1 + X2 + X3,
                                 data = test_data, weightit = W, vcov = "none",
                                 cluster = ~clus)
  }, "cluster.*not used", ignore.case = TRUE)

  expect_equal(attr(fit_stored, "cluster"), ~clus, ignore_attr = TRUE)
})

test_that("vcov = 'BS' and vcov = 'FWB' produce well-formed, reproducible variance matrices", {
  skip_on_cran()
  skip_if_not_installed("survival")
  skip_if_not_installed("fwb")

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))
  test_data <- .censor_Y_S(test_data)

  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5, data = test_data,
                  method = "glm", estimand = "ATE")
  })

  set.seed(1234)
  expect_no_condition({
    fit_bs <- coxph_weightit(survival::Surv(time, event) ~ A + X1 + X2,
                             data = test_data, weightit = W, vcov = "BS", R = 30)
  })

  expect_equal(fit_bs$vcov_type, "BS", ignore_attr = TRUE)
  V_bs <- vcov(fit_bs)
  expect_true(isSymmetric(unname(V_bs)))
  expect_true(all(diag(V_bs) > 0))

  set.seed(1234)
  expect_no_condition({
    fit_fwb <- coxph_weightit(survival::Surv(time, event) ~ A + X1 + X2,
                              data = test_data, weightit = W, vcov = "FWB", R = 30)
  })

  expect_equal(fit_fwb$vcov_type, "FWB", ignore_attr = TRUE)
  V_fwb <- vcov(fit_fwb)
  expect_true(isSymmetric(unname(V_fwb)))
  expect_true(all(diag(V_fwb) > 0))

  # Reproducible given the same seed
  set.seed(1234)
  fit_bs2 <- coxph_weightit(survival::Surv(time, event) ~ A + X1 + X2,
                            data = test_data, weightit = W, vcov = "BS", R = 30)
  expect_equal(vcov(fit_bs), vcov(fit_bs2))

  set.seed(1234)
  fit_fwb2 <- coxph_weightit(survival::Surv(time, event) ~ A + X1 + X2,
                             data = test_data, weightit = W, vcov = "FWB", R = 30)
  expect_equal(vcov(fit_fwb), vcov(fit_fwb2))

  # Point estimates are unaffected by the choice of vcov
  expect_equal(coef(fit_bs), coef(fit_fwb), tolerance = 1e-3)
})

test_that("Weighting method 'ebal' also works as input to coxph_weightit()", {
  skip_on_cran()
  skip_if_not_installed("survival")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))
  test_data <- .censor_Y_S(test_data)

  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "ebal", estimand = "ATE")
  })

  expect_true(is_not_null(attr(W, "Mparts")))

  expect_no_condition({
    fit <- coxph_weightit(survival::Surv(time, event) ~ A + X1 + X2 + X3,
                          data = test_data, weightit = W)
  })

  expect_identical(fit$vcov_type, "asympt")

  fit_ref <- survival::coxph(survival::Surv(time, event) ~ A + X1 + X2 + X3,
                             data = test_data, weights = W$weights)

  expect_equal(unname(coef(fit)), unname(coef(fit_ref)), tolerance = eps)
})

test_that("strata() cannot be used in the formula", {
  skip_on_cran()
  skip_if_not_installed("survival")

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))
  test_data <- .censor_Y_S(test_data)

  expect_error(
    coxph_weightit(survival::Surv(time, event) ~ A + survival::strata(X6),
                   data = test_data),
    "strata.*cannot be used", ignore.case = TRUE
  )
})

test_that("cluster() cannot be used as a formula term (must use cluster= argument)", {
  skip_on_cran()
  skip_if_not_installed("survival")

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))
  test_data <- .censor_Y_S(test_data)
  test_data$clus <- rep(1:20, length.out = nrow(test_data))

  expect_error(
    coxph_weightit(survival::Surv(time, event) ~ A + survival::cluster(clus),
                   data = test_data),
    "cluster.*cannot be used in the model formula", ignore.case = TRUE
  )
})

test_that("frailty() cannot be used in the formula", {
  skip_on_cran()
  skip_if_not_installed("survival")

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))
  test_data <- .censor_Y_S(test_data)
  test_data$clus <- rep(1:20, length.out = nrow(test_data))

  expect_error(
    coxph_weightit(survival::Surv(time, event) ~ A + survival::frailty(clus),
                   data = test_data),
    "frailty.*cannot be used", ignore.case = TRUE
  )
})

test_that("pspline() cannot be used in the formula", {
  skip_on_cran()
  skip_if_not_installed("survival")

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))
  test_data <- .censor_Y_S(test_data)

  expect_error(
    coxph_weightit(survival::Surv(time, event) ~ A + survival::pspline(X1),
                   data = test_data),
    "pspline.*cannot be used", ignore.case = TRUE
  )
})

test_that("tt() cannot be used in the formula", {
  skip_on_cran()
  skip_if_not_installed("survival")

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))
  test_data <- .censor_Y_S(test_data)

  expect_error(
    coxph_weightit(survival::Surv(time, event) ~ survival::tt(A),
                   data = test_data),
    "tt.*cannot be used", ignore.case = TRUE
  )
})

test_that("ridge() cannot be used in the formula", {
  skip_on_cran()
  skip_if_not_installed("survival")

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))
  test_data <- .censor_Y_S(test_data)

  expect_error(
    coxph_weightit(survival::Surv(time, event) ~ A + survival::ridge(X1, X2),
                   data = test_data),
    "ridge.*cannot be used", ignore.case = TRUE
  )
})

test_that("non-right-censored Surv() types are rejected", {
  skip_on_cran()
  skip_if_not_installed("survival")

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))
  test_data <- .censor_Y_S(test_data)

  # Counting-process (start, stop] format
  expect_error(
    coxph_weightit(survival::Surv(time, time + 1, event, type = "counting") ~ A + X1,
                   data = test_data),
    "only supports right-censoring", ignore.case = TRUE
  )

  # Interval censoring
  expect_error(
    coxph_weightit(survival::Surv(time, time + 1, type = "interval2") ~ A + X1,
                   data = test_data),
    "only supports right-censoring", ignore.case = TRUE
  )
})

test_that("a non-Surv response (including Surv2) is rejected", {
  skip_on_cran()
  skip_if_not_installed("survival")

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))
  test_data <- .censor_Y_S(test_data)

  expect_error(
    coxph_weightit(A ~ X1 + X2, data = test_data),
    "must be a survival.*object", ignore.case = TRUE
  )

  # Surv2 (time-varying / (id, time) long format) is also rejected -- it
  # inherits from "Surv" but is explicitly excluded
  expect_error(
    coxph_weightit(survival::Surv2(time, event) ~ A + X1, data = test_data),
    "must be a survival.*object", ignore.case = TRUE
  )
})

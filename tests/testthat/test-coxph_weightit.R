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

test_that("vcov = 'const' returns the model-based variance", {
  skip_on_cran()
  skip_if_not_installed("survival")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))
  test_data <- .censor_Y_S(test_data)

  # Unweighted: the model-based ("naive") variance is exactly the variance
  # coxph() returns when it isn't asked for a robust one
  expect_no_condition({
    fit0 <- coxph_weightit(survival::Surv(time, event) ~ A + X1 + X2 + X3,
                           data = test_data, vcov = "const")
  })

  expect_identical(fit0$vcov_type, "const")

  fit0_ref <- survival::coxph(survival::Surv(time, event) ~ A + X1 + X2 + X3,
                              data = test_data)

  expect_equal(unname(vcov(fit0)), unname(vcov(fit0_ref)), tolerance = eps)

  # Requesting it after the fact from an object fit with a different vcov gives
  # the same matrix (and differs from that object's own HC0 variance)
  fit_hc0 <- coxph_weightit(survival::Surv(time, event) ~ A + X1 + X2 + X3,
                            data = test_data)

  expect_no_condition({
    V_const <- vcov(fit_hc0, vcov = "const")
  })

  expect_equal(V_const, vcov(fit0), ignore_attr = TRUE, tolerance = eps)
  expect_not_equal(vcov(fit_hc0), V_const)

  # Weighted: matches the model-based variance from a coxph() fit given the
  # same case weights. coxph() switches to its own robust variance when given
  # non-integer weights, so `robust = FALSE` is needed to get the naive one.
  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5, data = test_data,
                  method = "glm", estimand = "ATE")
  })

  # `vcov = "const"` ignores estimation of the weights, so it warns
  expect_warning({
    fit_w <- coxph_weightit(survival::Surv(time, event) ~ A + X1 + X2 + X3,
                            data = test_data, weightit = W, vcov = "const")
  }, "const.*should not be used", ignore.case = TRUE)

  expect_identical(fit_w$vcov_type, "const")

  fit_w_ref <- survival::coxph(survival::Surv(time, event) ~ A + X1 + X2 + X3,
                               data = test_data, weights = W$weights,
                               robust = FALSE)

  expect_equal(unname(vcov(fit_w)), unname(vcov(fit_w_ref)), tolerance = eps)

  # Point estimates are unaffected by the choice of vcov
  expect_equal(coef(fit0), coef(fit_hc0), tolerance = eps)
})

test_that("an intercept-only (null) model is supported", {
  skip_on_cran()
  skip_if_not_installed("survival")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))
  test_data <- .censor_Y_S(test_data)

  # A null model has no coefficients, so score residuals (and therefore the
  # gradient and variance matrix) are undefined for it
  expect_no_condition({
    fit_null <- coxph_weightit(survival::Surv(time, event) ~ 1, data = test_data)
  })

  expect_s3_class(fit_null, "coxph_weightit")
  expect_length(coef(fit_null), 0L)
  expect_identical(dim(vcov(fit_null)), c(0L, 0L))
  expect_identical(fit_null$vcov_type, "HC0")

  fit_null_ref <- survival::coxph(survival::Surv(time, event) ~ 1,
                                  data = test_data)

  expect_equal(fit_null$loglik, fit_null_ref$loglik, tolerance = eps)

  # print() and summary() work on a coefficient-free model
  expect_no_condition(invisible(capture.output(print(fit_null))))
  expect_no_condition(invisible(capture.output(print(summary(fit_null)))))

  # Every variance type is a no-op for a model with no coefficients
  for (v in c("const", "HC0")) {
    expect_no_condition({
      fit_v <- coxph_weightit(survival::Surv(time, event) ~ 1,
                              data = test_data, vcov = v)
    })

    expect_identical(dim(vcov(fit_v)), c(0L, 0L))
    expect_identical(fit_v$vcov_type, v)
  }

  # Also works when weights are supplied, including with M-estimation
  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5, data = test_data,
                  method = "glm", estimand = "ATE")
  })

  expect_no_condition({
    fit_null_w <- coxph_weightit(survival::Surv(time, event) ~ 1,
                                 data = test_data, weightit = W)
  })

  expect_identical(fit_null_w$vcov_type, "asympt")
  expect_length(coef(fit_null_w), 0L)
  expect_identical(dim(vcov(fit_null_w)), c(0L, 0L))

  fit_null_w_ref <- survival::coxph(survival::Surv(time, event) ~ 1,
                                    data = test_data, weights = W$weights)

  expect_equal(fit_null_w$loglik, fit_null_w_ref$loglik, tolerance = eps)
})

test_that("a rank-deficient fit gives estimable coefficients and NA for the rest", {
  skip_on_cran()
  skip_if_not_installed("survival")
  skip_if_not_installed("sandwich")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))
  test_data <- .censor_Y_S(test_data)

  # X1dup is collinear with X1, so one of the two is not estimable
  test_data$X1dup <- test_data$X1

  expect_no_condition({
    fit <- coxph_weightit(survival::Surv(time, event) ~ A + X1 + X1dup + X2,
                          data = test_data)
  })

  fit_ref <- survival::coxph(survival::Surv(time, event) ~ A + X1 + X1dup + X2,
                             data = test_data, robust = TRUE)

  # `survival::coxph()` also returns NA for the redundant coefficient
  expect_equal(coef(fit), coef(fit_ref), tolerance = eps)
  expect_identical(is.na(coef(fit)),
                   c(A = FALSE, X1 = FALSE, X1dup = TRUE, X2 = FALSE))

  # `complete = TRUE` (the default) expands the variance to the full set of
  # coefficients, giving the aliased term an NA row and column rather than a
  # bogus variance; `complete = FALSE` covers the estimable coefficients only
  V <- vcov(fit)
  expect_identical(dim(V), c(4L, 4L))
  expect_true(all(is.na(V["X1dup", ])))
  expect_true(all(is.na(V[, "X1dup"])))

  V_est <- vcov(fit, complete = FALSE)
  expect_identical(dim(V_est), c(3L, 3L))
  expect_false(anyNA(V_est))
  expect_identical(colnames(V_est), c("A", "X1", "X2"))

  # The estimable standard errors match coxph()'s robust ones
  expect_equal(unname(sqrt(diag(V_est))),
               unname(sqrt(diag(vcov(fit_ref)))[!is.na(coef(fit_ref))]),
               tolerance = eps)

  # A perfectly collinear column carries no information, so dropping it changes
  # nothing about the estimable results
  expect_no_condition({
    fit_red <- coxph_weightit(survival::Surv(time, event) ~ A + X1 + X2,
                              data = test_data)
  })

  expect_equal(coef(fit)[!is.na(coef(fit))], coef(fit_red), tolerance = eps)
  expect_equal(V_est, vcov(fit_red), ignore_attr = TRUE, tolerance = eps)

  # No `gradient` is stored, so the variance comes from `psi`, evaluated against
  # the model matrix already reduced to its estimable columns
  expect_null(fit$gradient)
  expect_null(fit_red$gradient)

  # print() and summary() report only the estimable coefficients
  expect_no_condition(invisible(capture.output(print(fit))))
  expect_no_condition(s <- summary(fit))
  expect_identical(rownames(s$coefficients), c("A", "X1", "X2"))
  expect_false(anyNA(s$coefficients))

  # sandwich's estfun()/bread() are also restricted to the estimable coefficients
  expect_identical(ncol(sandwich::estfun(fit)), 3L)
  expect_identical(dim(sandwich::bread(fit)), c(3L, 3L))
  expect_equal(unname(sandwich::sandwich(fit)), unname(V_est), tolerance = eps)

  # vcov = "const" is the model-based variance over the estimable coefficients
  expect_no_condition({
    fit_const <- coxph_weightit(survival::Surv(time, event) ~ A + X1 + X1dup + X2,
                                data = test_data, vcov = "const")
  })

  fit_const_ref <- survival::coxph(survival::Surv(time, event) ~ A + X1 + X1dup + X2,
                                   data = test_data)

  expect_equal(unname(vcov(fit_const, complete = FALSE)),
               unname(vcov(fit_const_ref)[!is.na(coef(fit_const_ref)),
                                          !is.na(coef(fit_const_ref))]),
               tolerance = eps)

  # A redundant factor: two aliased coefficients at once
  test_data$G <- factor(rep(c("a", "b", "c"), length.out = nrow(test_data)))
  test_data$Gdup <- test_data$G

  expect_no_condition({
    fit_f <- coxph_weightit(survival::Surv(time, event) ~ A + G + Gdup,
                            data = test_data)
  })

  fit_f_ref <- survival::coxph(survival::Surv(time, event) ~ A + G + Gdup,
                               data = test_data, robust = TRUE)

  expect_equal(coef(fit_f), coef(fit_f_ref), tolerance = eps)
  expect_identical(sum(is.na(coef(fit_f))), 2L)

  expect_equal(unname(sqrt(diag(vcov(fit_f, complete = FALSE)))),
               unname(sqrt(diag(vcov(fit_f_ref)))[!is.na(coef(fit_f_ref))]),
               tolerance = eps)

  # Weights don't change any of this; M-estimation still accounts for them
  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5, data = test_data,
                  method = "glm", estimand = "ATE")
  })

  expect_no_condition({
    fit_w <- coxph_weightit(survival::Surv(time, event) ~ A + X1 + X1dup + X2,
                            data = test_data, weightit = W)
  })

  expect_identical(fit_w$vcov_type, "asympt")
  expect_identical(is.na(coef(fit_w)),
                   c(A = FALSE, X1 = FALSE, X1dup = TRUE, X2 = FALSE))
  expect_false(anyNA(vcov(fit_w, complete = FALSE)))

  expect_no_condition({
    fit_w_hc0 <- coxph_weightit(survival::Surv(time, event) ~ A + X1 + X1dup + X2,
                                data = test_data, weightit = W, vcov = "HC0")
  })

  fit_w_ref <- survival::coxph(survival::Surv(time, event) ~ A + X1 + X1dup + X2,
                               data = test_data, weights = W$weights,
                               robust = TRUE)

  expect_equal(coef(fit_w), coef(fit_w_ref), tolerance = eps)
  expect_equal(unname(vcov(fit_w_hc0, complete = FALSE)),
               unname(vcov(fit_w_ref)[!is.na(coef(fit_w_ref)),
                                      !is.na(coef(fit_w_ref))]),
               tolerance = eps)

  expect_true(all(diag(vcov(fit_w, complete = FALSE)) <=
                    diag(vcov(fit_w_hc0, complete = FALSE))))

  # Dropping the collinear column also changes nothing on the M-estimation path,
  # where a wrongly chosen column would give wrong standard errors rather than
  # an error
  expect_no_condition({
    fit_w_red <- coxph_weightit(survival::Surv(time, event) ~ A + X1 + X2,
                                data = test_data, weightit = W)
  })

  expect_equal(coef(fit_w)[!is.na(coef(fit_w))], coef(fit_w_red), tolerance = eps)
  expect_equal(vcov(fit_w, complete = FALSE), vcov(fit_w_red),
               ignore_attr = TRUE, tolerance = eps)
})

test_that("a model with no events returns NA coefficients rather than erroring", {
  skip_on_cran()
  skip_if_not_installed("survival")

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))
  test_data <- .censor_Y_S(test_data)

  # Censor every unit, so no events are observed
  test_data$event <- 0

  expect_warning({
    fit <- coxph_weightit(survival::Surv(time, event) ~ A + X1 + X2,
                          data = test_data)
  }, "no events", ignore.case = TRUE)

  expect_s3_class(fit, "coxph_weightit")
  expect_identical(fit$nevent, 0)

  # `survival::coxph()` also returns `NA` coefficients in this case
  fit_ref <- survival::coxph(survival::Surv(time, event) ~ A + X1 + X2,
                             data = test_data)

  expect_true(all(is.na(coef(fit))))
  expect_identical(names(coef(fit)), names(coef(fit_ref)))
  expect_type(coef(fit), "double")

  # Unlike coxph(), which reports a variance of 0, the variance is NA, matching
  # the coefficients
  V <- vcov(fit)
  expect_identical(dim(V), c(3L, 3L))
  expect_true(all(is.na(V)))
  expect_identical(colnames(V), names(coef(fit)))
  expect_identical(rownames(V), names(coef(fit)))

  # Requesting the variance after the fact is also NA rather than an error
  expect_no_condition({
    V_const <- vcov(fit, vcov = "const")
  })

  expect_true(all(is.na(V_const)))

  # print(), summary(), and confint() work on an unestimable model
  expect_no_condition(invisible(capture.output(print(fit))))
  expect_no_condition(invisible(capture.output(print(summary(fit)))))
  expect_no_condition(ci <- confint(fit))
  expect_true(all(is.na(ci)))

  # The returned object is otherwise complete; the no-events path used to
  # return a stripped-down object missing the model frame, formula, and more
  expect_true(is_not_null(fit$model))
  expect_true(is_not_null(fit$formula))
  expect_equal(nrow(fit$model), nrow(test_data))

  # A null model with no events (both degenerate cases at once)
  expect_warning({
    fit_null <- coxph_weightit(survival::Surv(time, event) ~ 1,
                               data = test_data)
  }, "no events", ignore.case = TRUE)

  expect_length(coef(fit_null), 0L)
  expect_identical(dim(vcov(fit_null)), c(0L, 0L))

  # Weights don't change any of this, for any variance type
  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5, data = test_data,
                  method = "glm", estimand = "ATE")
  })

  expect_warning({
    fit_w <- coxph_weightit(survival::Surv(time, event) ~ A + X1 + X2,
                            data = test_data, weightit = W)
  }, "no events", ignore.case = TRUE)

  expect_identical(fit_w$vcov_type, "asympt")
  expect_true(all(is.na(coef(fit_w))))
  expect_true(all(is.na(vcov(fit_w))))

  for (v in c("HC0", "BS", "FWB")) {
    expect_warning({
      fit_v <- coxph_weightit(survival::Surv(time, event) ~ A + X1 + X2,
                              data = test_data, weightit = W, vcov = v, R = 5)
    }, "no events", ignore.case = TRUE)

    # `vcov_type` carries R/fwb.args attributes for the bootstrap types
    expect_equal(fit_v$vcov_type, v, ignore_attr = TRUE)
    expect_true(all(is.na(vcov(fit_v))))
  }
})

test_that("anova() against a null (intercept-only) model works", {
  skip_on_cran()
  skip_if_not_installed("survival")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))
  test_data <- .censor_Y_S(test_data)

  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5, data = test_data,
                  method = "glm", estimand = "ATE")
  })

  expect_no_condition({
    fit <- coxph_weightit(survival::Surv(time, event) ~ A + X1 + X2,
                          data = test_data, weightit = W)
  })

  expect_no_condition({
    fit_null <- coxph_weightit(survival::Surv(time, event) ~ 1,
                               data = test_data, weightit = W)
  })

  expect_no_condition({
    a <- anova(fit, fit_null)
  })

  expect_s3_class(a, "anova")

  b <- coef(fit)
  V <- vcov(fit)

  # Testing all coefficients against a null model is the joint Wald test
  expect_identical(a[["Df"]][2L], length(b))
  expect_equal(a[["Chisq"]][2L], drop(t(b) %*% solve(V, b)), tolerance = eps)
  expect_equal(a[["Pr(>Chisq)"]][2L],
               pchisq(a[["Chisq"]][2L], length(b), lower.tail = FALSE),
               tolerance = eps)
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

# `.compute_vcov()` and `estfun()` reduce the model matrix to its estimable
# columns, so a stored `gradient` (which has one column per model.matrix column)
# would not line up. `coxph_weightit()` therefore stores none and lets them
# evaluate `psi` themselves.
test_that("nobs() counts units, and AIC()/BIC() still use the number of events", {
  skip_on_cran()
  skip_if_not_installed("survival")

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))
  test_data <- .censor_Y_S(test_data)

  W <- weightit(A ~ X1 + X2, data = test_data, method = "glm")

  fit <- coxph_weightit(survival::Surv(time, event) ~ A + X1,
                        data = test_data, weightit = W)

  # WeightIt's method is dispatched to, not `survival:::nobs.coxph()`, which
  # returns the number of events
  expect_identical(stats::nobs(fit), sum(fit$weights != 0))
  expect_identical(stats::nobs(fit), nrow(test_data))
  expect_false(identical(stats::nobs(fit), fit$nevent))

  # `logLik()` carries its own "nobs" attribute (the number of events), which is
  # what AIC() and BIC() consume, so they are unaffected by the above
  expect_identical(attr(logLik(fit), "nobs"), fit$nevent)

  ref <- survival::coxph(survival::Surv(time, event) ~ A + X1,
                         data = test_data, weights = W$weights)

  expect_equal(AIC(fit), AIC(ref))
  expect_equal(BIC(fit), BIC(ref))
})

test_that("anova.glm_weightit", {
  skip_on_cran()

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                data = test_data, method = "glm", estimand = "ATE")

  #object = larger model, object2 = smaller (nested) model
  fit_full <- glm_weightit(Y_B ~ A * (X1 + X2 + X3 + X4 + X5),
                           data = test_data, family = binomial,
                           weightit = W, vcov = "asympt")
  fit_reduced <- glm_weightit(Y_B ~ A + X1 + X2 + X3 + X4 + X5,
                              data = test_data, family = binomial,
                              weightit = W, vcov = "asympt")

  a <- anova(fit_full, fit_reduced)
  expect_s3_class(a, "anova")
  expect_s3_class(a, "data.frame")
  expect_equal(nrow(a), 2L)
  expect_true(all(c("Res.Df", "Df", "Chisq", "Pr(>Chisq)") %in% colnames(a)))
  expect_true(is.finite(a[["Chisq"]][2L]))
  expect_true(a[["Pr(>Chisq)"]][2L] >= 0 && a[["Pr(>Chisq)"]][2L] <= 1)

  #test = "Chisq" is documented as the only currently-allowed option
  expect_no_condition(anova(fit_full, fit_reduced, test = "Chisq"))
  expect_error(anova(fit_full, fit_reduced, test = "F"))

  #method = "Wald" is documented as the only currently-allowed option
  expect_no_condition(anova(fit_full, fit_reduced, method = "Wald"))
  expect_error(anova(fit_full, fit_reduced, method = "LRT"))

  #vcov override should generally change the test statistic
  a_hc0 <- anova(fit_full, fit_reduced, vcov = "HC0")
  expect_s3_class(a_hc0, "anova")
  expect_not_equal(a[["Chisq"]][2L], a_hc0[["Chisq"]][2L])

  #vcov = "none" has no variance matrix to use and should error
  fit_full_none <- glm_weightit(Y_B ~ A * (X1 + X2 + X3 + X4 + X5),
                                data = test_data, family = binomial,
                                weightit = W, vcov = "none")
  fit_reduced_none <- glm_weightit(Y_B ~ A + X1 + X2 + X3 + X4 + X5,
                                   data = test_data, family = binomial,
                                   weightit = W, vcov = "none")
  expect_error(anova(fit_full_none, fit_reduced_none),
              "no variance matrix", ignore.case = TRUE)

  #models must be fit with the same outcome/units
  fit_diff_y <- glm_weightit(Y_C ~ A + X1 + X2 + X3 + X4 + X5,
                             data = test_data, weightit = W, vcov = "asympt")
  expect_error(anova(fit_full, fit_diff_y))
})

test_that("anova.multinom_weightit", {
  skip_on_cran()

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))
  test_data$Y_M <- factor(test_data$Y_O, ordered = FALSE)

  W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                data = test_data, method = "glm", estimand = "ATE")

  #object = larger model, object2 = smaller (nested) model
  fit_full <- multinom_weightit(Y_M ~ A * (X1 + X2 + X3 + X4 + X5),
                                data = test_data, weightit = W, vcov = "asympt")
  fit_reduced <- multinom_weightit(Y_M ~ A + X1 + X2 + X3 + X4 + X5,
                                   data = test_data, weightit = W, vcov = "asympt")

  a <- anova(fit_full, fit_reduced)
  expect_s3_class(a, "anova")
  expect_s3_class(a, "data.frame")
  expect_equal(nrow(a), 2L)
  expect_true(all(c("Res.Df", "Df", "Chisq", "Pr(>Chisq)") %in% colnames(a)))
  expect_true(is.finite(a[["Chisq"]][2L]))
  expect_true(a[["Pr(>Chisq)"]][2L] >= 0 && a[["Pr(>Chisq)"]][2L] <= 1)

  #test = "Chisq" is documented as the only currently-allowed option
  expect_no_condition(anova(fit_full, fit_reduced, test = "Chisq"))
  expect_error(anova(fit_full, fit_reduced, test = "F"))

  #method = "Wald" is documented as the only currently-allowed option
  expect_no_condition(anova(fit_full, fit_reduced, method = "Wald"))
  expect_error(anova(fit_full, fit_reduced, method = "LRT"))

  #vcov override should generally change the test statistic
  a_hc0 <- anova(fit_full, fit_reduced, vcov = "HC0")
  expect_not_equal(a[["Chisq"]][2L], a_hc0[["Chisq"]][2L])
})

test_that("anova.ordinal_weightit", {
  skip_on_cran()

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                data = test_data, method = "glm", estimand = "ATE")

  #object = larger model, object2 = smaller (nested) model
  fit_full <- ordinal_weightit(Y_O ~ A * (X1 + X2 + X3 + X4 + X5),
                               data = test_data, weightit = W, vcov = "asympt")
  fit_reduced <- ordinal_weightit(Y_O ~ A + X1 + X2 + X3 + X4 + X5,
                                  data = test_data, weightit = W, vcov = "asympt")

  a <- anova(fit_full, fit_reduced)
  expect_s3_class(a, "anova")
  expect_s3_class(a, "data.frame")
  expect_equal(nrow(a), 2L)
  expect_true(all(c("Res.Df", "Df", "Chisq", "Pr(>Chisq)") %in% colnames(a)))
  expect_true(is.finite(a[["Chisq"]][2L]))
  expect_true(a[["Pr(>Chisq)"]][2L] >= 0 && a[["Pr(>Chisq)"]][2L] <= 1)

  #test = "Chisq" is documented as the only currently-allowed option
  expect_no_condition(anova(fit_full, fit_reduced, test = "Chisq"))
  expect_error(anova(fit_full, fit_reduced, test = "F"))

  #method = "Wald" is documented as the only currently-allowed option
  expect_no_condition(anova(fit_full, fit_reduced, method = "Wald"))
  expect_error(anova(fit_full, fit_reduced, method = "LRT"))

  #vcov override should generally change the test statistic
  a_hc0 <- anova(fit_full, fit_reduced, vcov = "HC0")
  expect_not_equal(a[["Chisq"]][2L], a_hc0[["Chisq"]][2L])
})

test_that("anova.coxph_weightit", {
  skip_on_cran()
  skip_if_not_installed("survival")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                data = test_data, method = "glm", estimand = "ATE")

  #object = larger model, object2 = smaller (nested) model
  fit_full <- coxph_weightit(survival::Surv(Y_S) ~ A * (X1 + X2 + X3 + X4 + X5),
                             data = test_data, weightit = W, vcov = "asympt")
  fit_reduced <- coxph_weightit(survival::Surv(Y_S) ~ A + X1 + X2 + X3 + X4 + X5,
                                data = test_data, weightit = W, vcov = "asympt")

  a <- anova(fit_full, fit_reduced)
  expect_s3_class(a, "anova")
  expect_s3_class(a, "data.frame")
  expect_equal(nrow(a), 2L)
  expect_true(all(c("Res.Df", "Df", "Chisq", "Pr(>Chisq)") %in% colnames(a)))
  expect_true(is.finite(a[["Chisq"]][2L]))
  expect_true(a[["Pr(>Chisq)"]][2L] >= 0 && a[["Pr(>Chisq)"]][2L] <= 1)

  #test = "Chisq" is documented as the only currently-allowed option
  expect_no_condition(anova(fit_full, fit_reduced, test = "Chisq"))
  expect_error(anova(fit_full, fit_reduced, test = "F"))

  #method = "Wald" is documented as the only currently-allowed option
  expect_no_condition(anova(fit_full, fit_reduced, method = "Wald"))
  expect_error(anova(fit_full, fit_reduced, method = "LRT"))

  #vcov override changes the test statistic, exactly as for the other model classes
  a_hc0 <- anova(fit_full, fit_reduced, vcov = "HC0")
  expect_not_equal(a[["Chisq"]][2L], a_hc0[["Chisq"]][2L])

  #cluster-robust vcov override
  set.seed(123)
  clus <- sample(1:50, nrow(test_data), replace = TRUE)
  a_clus <- anova(fit_full, fit_reduced, vcov = "HC0", cluster = clus)
  expect_s3_class(a_clus, "anova")
  expect_not_equal(a_hc0[["Chisq"]][2L], a_clus[["Chisq"]][2L])
})

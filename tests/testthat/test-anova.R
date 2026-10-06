test_that("anova() works for each model class", {
  skip_on_cran()
  skip_if_not_installed("patrick")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  patrick::with_parameters_test_that(
    "anova.{cls}_weightit",
    {
      if (cls == "coxph") {
        skip_if_not_installed("survival")
      }

      test_data <- readRDS(test_path("fixtures", "test_data.rds"))
      test_data <- prep(test_data)

      W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                    data = test_data, method = "glm", estimand = "ATE")

      # Each call is built from the fitter's name and unevaluated arguments, and evaluated
      # here, so `$call` reads exactly as if it had been typed
      fit_call <- function(rhs) {
        cl <- rlang::call2(paste0(cls, "_weightit"), call("~", outcome, rhs),
                           data = quote(test_data))

        if (binomial_family) {
          cl$family <- quote(binomial)
        }

        cl$weightit <- quote(W)
        cl$vcov <- "asympt"
        cl
      }

      #object = larger model, object2 = smaller (nested) model
      fit_full <- eval(fit_call(quote(A * (X1 + X2 + X3 + X4 + X5))))
      fit_reduced <- eval(fit_call(quote(A + X1 + X2 + X3 + X4 + X5)))

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

      if (check_hc0_class) {
        expect_s3_class(a_hc0, "anova")
      }

      expect_not_equal(a[["Chisq"]][2L], a_hc0[["Chisq"]][2L])
    },
    # `prep`: builds the outcome that `outcome` names.
    # `binomial_family`: whether the fits pass `family = binomial`.
    # `check_hc0_class`: whether the class of the `vcov = "HC0"` result is checked.
    .cases = patrick::cases(
      list(cls = "glm", outcome = quote(Y_B),
           prep = identity,
           binomial_family = TRUE, check_hc0_class = TRUE),
      list(cls = "multinom", outcome = quote(Y_M),
           prep = function(d) {
             d$Y_M <- factor(d$Y_O, ordered = FALSE)
             d
           },
           binomial_family = FALSE, check_hc0_class = FALSE),
      list(cls = "ordinal", outcome = quote(Y_O),
           prep = identity,
           binomial_family = FALSE, check_hc0_class = FALSE),
      list(cls = "coxph", outcome = quote(survival::Surv(Y_S)),
           prep = identity,
           binomial_family = FALSE, check_hc0_class = FALSE)
    )
  )
})

test_that("anova.glm_weightit needs a variance matrix and a common outcome", {
  skip_on_cran()

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                data = test_data, method = "glm", estimand = "ATE")

  fit_full <- glm_weightit(Y_B ~ A * (X1 + X2 + X3 + X4 + X5),
                           data = test_data, family = binomial,
                           weightit = W, vcov = "asympt")

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

test_that("anova.coxph_weightit accepts a cluster-robust vcov", {
  skip_on_cran()
  skip_if_not_installed("survival")

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                data = test_data, method = "glm", estimand = "ATE")

  fit_full <- coxph_weightit(survival::Surv(Y_S) ~ A * (X1 + X2 + X3 + X4 + X5),
                             data = test_data, weightit = W, vcov = "asympt")
  fit_reduced <- coxph_weightit(survival::Surv(Y_S) ~ A + X1 + X2 + X3 + X4 + X5,
                                data = test_data, weightit = W, vcov = "asympt")

  a_hc0 <- anova(fit_full, fit_reduced, vcov = "HC0")

  #cluster-robust vcov override
  set.seed(123)
  clus <- sample(1:50, nrow(test_data), replace = TRUE)
  a_clus <- anova(fit_full, fit_reduced, vcov = "HC0", cluster = clus)
  expect_s3_class(a_clus, "anova")
  expect_not_equal(a_hc0[["Chisq"]][2L], a_clus[["Chisq"]][2L])
})

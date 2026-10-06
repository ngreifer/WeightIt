# Weights supplied to `weights` are fixed: they enter the fit as `glm()`'s `weights`
# do, the default variance treats them as fixed, and the bootstrap holds them fixed
# rather than re-estimating them. Internally they are stored as the sampling weights
# of a placeholder `weightit` object, so a fit with `weights` should match one with
# such an object supplied to `weightit` (as users previously had to do).

.prep_weights_data <- function() {
  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  cutoff <- quantile(test_data$Y_S, .8)
  test_data$event <- as.numeric(test_data$Y_S < cutoff)
  test_data$time <- pmin(test_data$Y_S, cutoff)
  test_data$Y_M <- factor(test_data$Y_O, ordered = FALSE)

  test_data
}

test_that("weights can be supplied to every model fitting function", {
  skip_on_cran()
  skip_if_not_installed("patrick")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- .prep_weights_data()

  fake <- as.weightit(test_data$SW, treat = test_data$A, estimand = "ATE",
                      s.weights = rep.int(1, nrow(test_data)))

  patrick::with_parameters_test_that(
    "{fitter}",
    {
      if (fitter == "coxph_weightit") {
        skip_if_not_installed("survival")
      }

      fit_call <- function(...) {
        do.call(fitter, c(fit_args, list(data = quote(test_data), ...)))
      }

      expect_no_condition({
        fit <- fit_call(weights = quote(SW))
      })

      # The variable name, its name as a string, and the vector itself are equivalent
      expect_equal(coef(fit_call(weights = "SW")), coef(fit))
      expect_equal(coef(fit_call(weights = test_data$SW)), coef(fit))

      # The weights are treated as fixed
      expect_identical(fit$vcov_type, "HC0")
      expect_error(fit_call(weights = quote(SW), vcov = "asympt"))

      fit_fake <- fit_call(weightit = fake)

      expect_equal(coef(fit), coef(fit_fake), tolerance = eps)
      expect_equal(vcov(fit), vcov(fit_fake), tolerance = eps)

      # They are not sampling weights, so only `(weights)` is added to the model frame
      expect_equal(model.frame(fit)[["(weights)"]], test_data$SW)
      expect_null(model.frame(fit)[["(s.weights)"]])
    },
    .cases = patrick::cases(
      list(fitter = "glm_weightit",
           fit_args = list(quote(Y_B ~ A + X1), family = quote(binomial))),
      list(fitter = "lm_weightit",
           fit_args = list(quote(Y_C ~ A + X1))),
      list(fitter = "multinom_weightit",
           fit_args = list(quote(Y_M ~ A + X1))),
      list(fitter = "ordinal_weightit",
           fit_args = list(quote(Y_O ~ A + X1))),
      list(fitter = "coxph_weightit",
           fit_args = list(quote(survival::Surv(time, event) ~ A + X1)))
    )
  )
})

test_that("weights give the same fit as glm(), lm(), and survival::coxph()", {
  skip_on_cran()
  skip_if_not_installed("sandwich")
  skip_if_not_installed("survival")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- .prep_weights_data()

  fit <- glm_weightit(Y_B ~ A + X1, data = test_data, family = binomial, weights = SW)

  #Warnings are about the non-integer weights
  fit_g <- suppressWarnings({
    glm(Y_B ~ A + X1, data = test_data, family = binomial, weights = SW)
  })

  expect_equal(coef(fit), coef(fit_g), tolerance = eps)
  expect_equal(vcov(fit), sandwich::sandwich(fit_g), tolerance = eps)

  # With fixed weights, the model-based variance is valid, so it gives no warning
  expect_no_warning({
    fit_const <- glm_weightit(Y_B ~ A + X1, data = test_data, family = binomial,
                              weights = SW, vcov = "const")
  })

  expect_equal(vcov(fit_const), vcov(fit_g), tolerance = eps)

  fit_lm <- lm_weightit(Y_C ~ A + X1, data = test_data, weights = SW)
  fit_l <- lm(Y_C ~ A + X1, data = test_data, weights = SW)

  expect_equal(coef(fit_lm), coef(fit_l), tolerance = eps)
  expect_equal(vcov(fit_lm), sandwich::vcovHC(fit_l, type = "HC0"), tolerance = eps)

  fit_cox <- coxph_weightit(survival::Surv(time, event) ~ A + X1, data = test_data,
                            weights = SW)
  fit_c <- survival::coxph(survival::Surv(time, event) ~ A + X1, data = test_data,
                           weights = SW, robust = TRUE)

  expect_equal(unname(coef(fit_cox)), unname(coef(fit_c)), tolerance = eps)
  expect_equal(unname(vcov(fit_cox)), unname(vcov(fit_c)), tolerance = eps)
})

test_that("a weightit object supplied to weights is treated as though supplied to weightit", {
  skip_on_cran()

  test_data <- .prep_weights_data()

  W <- weightit(A ~ X1 + X2 + X3, data = test_data, method = "glm", estimand = "ATE")

  fit_w <- glm_weightit(Y_B ~ A + X1, data = test_data, family = binomial, weights = W)
  fit <- glm_weightit(Y_B ~ A + X1, data = test_data, family = binomial, weightit = W)

  expect_identical(coef(fit_w), coef(fit))
  expect_identical(vcov(fit_w), vcov(fit))
  expect_identical(fit_w$vcov_type, "asympt")
  expect_identical(fit_w$weightit, W)

  # The stored call names the object as `weightit`, so `update()` refits it
  expect_identical(fit_w$call[["weightit"]], quote(W))
  expect_null(fit_w$call[["weights"]])
})

test_that("weights and weightit cannot both be supplied, and weights are checked", {
  skip_on_cran()

  test_data <- .prep_weights_data()

  W <- weightit(A ~ X1 + X2 + X3, data = test_data, method = "glm", estimand = "ATE")

  expect_error(glm_weightit(Y_B ~ A, data = test_data, weightit = W, weights = SW),
               "Only one of")

  test_data$neg <- test_data$X1
  test_data$miss <- test_data$SW
  is.na(test_data$miss[1L]) <- TRUE

  expect_error(glm_weightit(Y_B ~ A, data = test_data, weights = neg), "weights")
  expect_error(glm_weightit(Y_B ~ A, data = test_data, weights = miss), "weights")
  expect_error(glm_weightit(Y_B ~ A, data = test_data, weights = "not_a_var"),
               "name of a variable")
  expect_error(glm_weightit(Y_B ~ A, data = test_data, weights = test_data$SW[-1L]),
               "same length")

  # `NULL` weights give an unweighted fit
  expect_equal(coef(glm_weightit(Y_B ~ A, data = test_data, weights = NULL)),
               coef(glm_weightit(Y_B ~ A, data = test_data)))
})

test_that("the bootstrap holds weights fixed", {
  skip_on_cran()
  skip_if_not_installed("fwb")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- .prep_weights_data()

  set.seed(123)
  fit <- glm_weightit(Y_B ~ A + X1, data = test_data, family = binomial,
                      weights = SW, vcov = "FWB", R = 20)

  # The same draws, applied to the fixed weights by hand
  set.seed(123)
  fwb_out <- suppressWarnings({
    fwb::fwb(data.frame(SW = test_data$SW), R = 20, verbose = FALSE,
             statistic = function(data, w) {
               coef(glm(Y_B ~ A + X1, data = test_data, family = binomial,
                        weights = test_data$SW * w))
             })
  })

  expect_equal(unname(vcov(fit)), unname(cov(fwb_out$t)), tolerance = eps)

  set.seed(123)
  fit_bs <- glm_weightit(Y_B ~ A + X1, data = test_data, family = binomial,
                         weights = SW, vcov = "BS", R = 20)

  expect_true(all(is.finite(sqrt(diag(vcov(fit_bs))))))
})

test_that("update() works with weights", {
  skip_on_cran()

  test_data <- .prep_weights_data()

  W <- weightit(A ~ X1 + X2 + X3, data = test_data, method = "glm", estimand = "ATE")

  fit <- glm_weightit(Y_B ~ A + X1, data = test_data, family = binomial, weights = SW)

  # New weights are evaluated in `data`, as in the original call
  expect_equal(update(fit, weights = X8^2),
               glm_weightit(Y_B ~ A + X1, data = test_data, family = binomial,
                            weights = X8^2),
               ignore_formula_env = TRUE)

  expect_equal(coef(update(fit, data = test_data[1:400, ])),
               coef(glm_weightit(Y_B ~ A + X1, data = test_data[1:400, ],
                                 family = binomial, weights = SW)))

  # `NULL` removes the weights, and a new `weightit` object replaces them
  fit_unw <- update(fit, weights = NULL)
  expect_null(fit_unw$weightit)
  expect_equal(coef(fit_unw),
               coef(glm_weightit(Y_B ~ A + X1, data = test_data, family = binomial)))

  fit_w <- update(fit, weightit = W)
  expect_null(fit_w$call[["weights"]])
  expect_equal(coef(fit_w),
               coef(glm_weightit(Y_B ~ A + X1, data = test_data, family = binomial,
                                 weightit = W)))

  # With a weightit object present, new weights are still its sampling weights
  fit_W <- glm_weightit(Y_B ~ A + X1, data = test_data, family = binomial, weightit = W)
  expect_equal(coef(update(fit_W, weights = "SW")),
               coef(glm_weightit(Y_B ~ A + X1, data = test_data, family = binomial,
                                 weightit = update(W, s.weights = "SW"))))
})

test_that("weights not in data are found in the environment of the formula", {
  skip_on_cran()

  test_data <- .prep_weights_data()

  fit_fun <- function(dat) {
    y <- dat$Y_B
    a <- dat$A
    w <- dat$SW

    glm_weightit(y ~ a, family = binomial, weights = w)
  }

  expect_equal(unname(coef(fit_fun(test_data))),
               unname(coef(glm_weightit(Y_B ~ A, data = test_data, family = binomial,
                                        weights = SW))))
})

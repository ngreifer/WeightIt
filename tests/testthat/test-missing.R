# Tests for the `missing` argument of weightit(), across methods that
# exercise different code paths:
#   - "ind"  (missingness indicators): glm, ebal, gbm
#   - "saem" (misaem-based EM, glm only)
#   - unspecified (NULL) default, and disallowed explicit values
#
# Missingness is injected ad hoc into in-memory copies of test_data via the
# inject_missingness() helper (see helpers.R). The fixture .rds and its
# generator script are never modified.

eps <- if (capabilities("long.double")) 1e-5 else 1e-3

test_data <- readRDS(test_path("fixtures", "test_data.rds"))

# Inject 10% missingness into two numeric covariates, reproducibly
test_data_na <- inject_missingness(test_data, c("X1", "X2"), prop = 0.1, seed = 4321)

test_that("missing = 'ind' (explicit) fits with glm and covs retains original NAs", {
  skip_on_cran()

  expect_no_error({
    W <- weightit(A ~ X1 + X2 + X3 + X4, data = test_data_na,
                  method = "glm", missing = "ind")
  })

  expect_s3_class(W, "weightit")
  expect_true(all(is.finite(W$weights)))
  expect_true(all(W$weights > 0))

  # weightit()'s documented `covs` value: "the covariates used in the fitting
  # ... may have been altered in the fitting process." In practice, `covs` is
  # extracted from the original formula/data *before* any method-internal
  # missing-data handling (median imputation, indicator columns, etc.), so it
  # should still contain the original NAs even though the method used a
  # median-imputed/indicator-augmented version internally to fit the model.
  expect_true(anyNA(W$covs[, c("X1", "X2")]))
  expect_equal(W$covs[, c("X1", "X2")], test_data_na[, c("X1", "X2")],
              ignore_attr = TRUE)
})

test_that("missing = 'ind' (explicit) fits with ebal and covs retains original NAs", {
  skip_on_cran()
  skip_if_not_installed("cobalt")

  expect_no_error({
    W <- weightit(A ~ X1 + X2 + X3 + X4, data = test_data_na,
                  method = "ebal", missing = "ind")
  })

  expect_s3_class(W, "weightit")
  expect_true(all(is.finite(W$weights)))
  expect_true(all(W$weights > 0))

  expect_true(anyNA(W$covs[, c("X1", "X2")]))
  expect_equal(W$covs[, c("X1", "X2")], test_data_na[, c("X1", "X2")],
              ignore_attr = TRUE)
})

test_that("missing = 'ind' (explicit) fits with gbm and covs retains original NAs", {
  skip_on_cran()
  skip_if_not_installed("gbm")
  skip_if_not_installed("cobalt")

  set.seed(1)
  expect_no_error({
    W <- weightit(A ~ X1 + X2 + X3 + X4, data = test_data_na,
                  method = "gbm", missing = "ind", criterion = "smd.mean",
                  n.trees = 500)
  })

  expect_s3_class(W, "weightit")
  expect_true(all(is.finite(W$weights)))
  expect_true(all(W$weights > 0))

  expect_true(anyNA(W$covs[, c("X1", "X2")]))
  expect_equal(W$covs[, c("X1", "X2")], test_data_na[, c("X1", "X2")],
              ignore_attr = TRUE)
})

test_that("missing unspecified (NULL) warns and defaults to the method's first allowable value", {
  skip_on_cran()

  # .process_missing()'s allowable.missings[1L] fallback: for "glm", this is
  # "ind" (missing = c("ind", "saem") in .weightit_methods).
  expect_warning({
    W <- weightit(A ~ X1 + X2 + X3 + X4, data = test_data_na, method = "glm")
  }, "missing values are present", ignore.case = TRUE)

  expect_s3_class(W, "weightit")
  expect_true(all(is.finite(W$weights)))
  expect_true(all(W$weights > 0))
  # `missing` field reflects the *resolved* value chosen by .process_missing()'s
  # allowable.missings[1L] fallback, not the NULL originally supplied.
  expect_equal(W$missing, "ind")
})

test_that("an explicit, disallowed missing value errors informatively", {
  skip_on_cran()

  # ebal's .weightit_methods entry only lists missing = "ind"; requesting
  # "surr" (a gbm-only option) should hit the `missing %nin% allowable.missings`
  # branch of .process_missing() and error before fitting.
  expect_error({
    weightit(A ~ X1 + X2 + X3 + X4, data = test_data_na,
             method = "ebal", missing = "surr")
  }, "only.*allowed.*missing", ignore.case = TRUE)
})

test_that("missing = 'saem' with glm: binary treatment, link = 'logit' fits", {
  skip_on_cran()
  skip_if_not_installed("misaem")

  expect_no_error({
    W <- weightit(A ~ X1 + X2 + X3 + X4, data = test_data_na,
                 method = "glm", missing = "saem", link = "logit")
  })

  expect_s3_class(W, "weightit")
  expect_true(all(is.finite(W$weights)))
  expect_true(all(W$weights > 0))
  expect_true(anyNA(W$covs[, c("X1", "X2")]))
})

test_that("missing = 'saem' with glm: binary treatment, non-logit link errors", {
  skip_on_cran()
  skip_if_not_installed("misaem")

  # Per R/weightit2glm.R: acceptable.links <- "logit" when missing == "saem"
  # for binary treatments; any other link should error.
  expect_error({
    weightit(A ~ X1 + X2 + X3 + X4, data = test_data_na,
            method = "glm", missing = "saem", link = "probit")
  }, "only.*logit.*allowed", ignore.case = TRUE)
})

test_that("missing = 'saem' with glm: continuous treatment fits", {
  skip_on_cran()
  skip_if_not_installed("misaem")

  expect_no_error({
    W <- weightit(Ac ~ X1 + X2 + X3 + X4, data = test_data_na,
                 method = "glm", missing = "saem")
  })

  expect_s3_class(W, "weightit")
  expect_true(all(is.finite(W$weights)))
  expect_true(all(W$weights > 0))
  expect_true(anyNA(W$covs[, c("X1", "X2")]))
})

test_that("missing = 'saem' with glm: s.weights errors", {
  skip_on_cran()
  skip_if_not_installed("misaem")

  # Per R/weightit2glm.R (binary/multi/continuous paths all check
  # !all_the_same(s.weights) and error with this message when missing = "saem").
  expect_error({
    weightit(A ~ X1 + X2 + X3 + X4, data = test_data_na,
            method = "glm", missing = "saem", s.weights = "SW")
  }, "sampling weights cannot be used.*saem", ignore.case = TRUE)
})

test_that("missing = 'saem' with glm produces no M-estimation parts", {
  skip_on_cran()
  skip_if_not_installed("misaem")

  # Documented at R/weightit2glm.R:92 as unsupported for M-estimation with
  # missing = "saem" (binary treatment path).
  W_bin <- weightit(A ~ X1 + X2 + X3 + X4, data = test_data_na,
                    method = "glm", missing = "saem", link = "logit")
  expect_null(attr(W_bin, "Mparts", exact = TRUE))

  # Continuous treatment path
  W_cont <- weightit(Ac ~ X1 + X2 + X3 + X4, data = test_data_na,
                     method = "glm", missing = "saem")
  expect_null(attr(W_cont, "Mparts", exact = TRUE))
})

test_that("missing = 'saem' with glm warns if multi.method is explicitly set (multi-category treatment)", {
  skip_on_cran()
  skip_if_not_installed("misaem")

  # multi.method only matters for multi-category treatments (Am); the warning
  # fires for any explicit multi.method other than "saem"/"glm" (see
  # R/weightit2glm.R:595-601 -- multi.method = "glm" is silently accepted
  # without warning and still forced to "saem").
  expect_warning({
    W <- weightit(Am ~ X1 + X2 + X3 + X4, data = test_data_na,
                 method = "glm", missing = "saem", multi.method = "weightit")
  }, "multi\\.method.*ignored", ignore.case = TRUE)

  expect_s3_class(W, "weightit")
  expect_null(attr(W, "Mparts", exact = TRUE))
})

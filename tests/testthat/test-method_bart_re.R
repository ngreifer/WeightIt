#Small sampler settings with set.seed() before each fit for speed and
#reproducibility. bartisan fits the random intercepts itself, alongside the trees.

test_that("Binary treatment with random intercepts", {
  skip_on_cran()
  skip_if_not_installed("bartisan")
  skip_if_not_installed("cobalt")

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  #`cluster` is a proxy for the omitted confounder X1 in the fixture
  set.seed(123)
  expect_no_condition({
    W <- weightit(A ~ X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + (1 | cluster),
                  data = test_data, method = "bart", estimand = "ATE",
                  num_trees = 20L, num_draws = 100L, num_burn = 100L,
                  include.obj = TRUE)
  })

  expect_true(is.numeric(W$ps))
  expect_true(all(W$ps > 0 & W$ps < 1))
  expect_true(all(is.finite(W$weights) & W$weights > 0))
  expect_s3_class(W$obj, "bartisan_fit")
  expect_false(is_null(W$obj$ranef))
  expect_null(attr(W, "Mparts", exact = TRUE))

  #Distinguishable from a fit that ignores clustering (X1 also omitted)
  set.seed(123)
  W_none <- weightit(A ~ X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                     data = test_data, method = "bart", estimand = "ATE",
                     num_trees = 20L, num_draws = 100L, num_burn = 100L)
  expect_not_equal(W$weights, W_none$weights)

  #Balance on the omitted confounder X1 improves when cluster stands in for it
  smd_re <- abs(cobalt::col_w_smd(test_data$X1, test_data$A, W$weights))
  smd_none <- abs(cobalt::col_w_smd(test_data$X1, test_data$A, W_none$weights))
  expect_lt(smd_re, smd_none)
})

test_that("Continuous treatment with random intercepts", {
  skip_on_cran()
  skip_if_not_installed("bartisan")

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  set.seed(123)
  expect_no_condition({
    W <- weightit(Ac ~ X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + (1 | cluster),
                  data = test_data, method = "bart",
                  num_trees = 20L, num_draws = 100L, num_burn = 100L,
                  include.obj = TRUE)
  })

  expect_true(all(is.finite(W$weights) & W$weights > 0))
  expect_s3_class(W$obj, "bartisan_fit")
  expect_false(is_null(W$obj$ranef))
  expect_null(attr(W, "Mparts", exact = TRUE))
})

test_that("Multi-category treatment with random intercepts", {
  skip_on_cran()
  skip_if_not_installed("bartisan")

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  set.seed(123)
  expect_no_condition({
    W <- weightit(Am ~ X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + (1 | cluster),
                  data = test_data, method = "bart", estimand = "ATE",
                  num_trees = 20L, num_draws = 100L, num_burn = 100L,
                  include.obj = TRUE)
  })

  expect_true(all(is.finite(W$weights) & W$weights > 0))
  expect_s3_class(W$obj, "bartisan_fit")
  expect_false(is_null(W$obj$ranef))
  expect_null(attr(W, "Mparts", exact = TRUE))
})

test_that("Multiple grouping factors are supported, random slopes are not", {
  skip_on_cran()
  skip_if_not_installed("bartisan")

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  set.seed(123)
  expect_no_condition({
    W <- weightit(A ~ X2 + X3 + X7 + (1 | cluster) + (1 | X6),
                  data = test_data, method = "bart", estimand = "ATE",
                  num_trees = 20L, num_draws = 100L, num_burn = 100L,
                  include.obj = TRUE)
  })

  expect_true(all(W$ps > 0 & W$ps < 1))
  expect_true(all(is.finite(W$weights) & W$weights > 0))
  expect_s3_class(W$obj, "bartisan_fit")

  #bartisan refuses a random slope, and its error is relayed
  expect_error({
    weightit(A ~ X2 + X3 + (1 + X2 | cluster), data = test_data,
             method = "bart", num_draws = 20L, num_burn = 20L)
  }, "bartisan")
})

test_that("Random effects method gating", {
  skip_on_cran()

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  #re_ok gating: a method with re_ok = FALSE errors on a random effects formula
  expect_error(weightit(A ~ X2 + (1 | cluster), data = test_data, method = "cbps"),
               "random effects", ignore.case = TRUE)
})

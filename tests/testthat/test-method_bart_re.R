#Small, single-threaded, seeded sampler settings throughout for speed and
#reproducibility. stan4bart's iter/warmup/chains/cores control the sampling;
#BART hyperparameters go in bart_args. stan4bart may emit benign Stan diagnostic
#warnings with these tiny settings, so fits are wrapped in suppressWarnings().

test_that("Binary treatment with random effects (stan4bart)", {
  skip_on_cran()
  skip_if_not_installed("stan4bart")
  skip_if_not_installed("cobalt")

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  #`cluster` is a proxy for the omitted confounder X1 in the fixture
  expect_no_error(suppressWarnings({
    W <- weightit(A ~ X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + (1 | cluster),
                  data = test_data, method = "bart", estimand = "ATE",
                  iter = 40, warmup = 20, chains = 1, cores = 1, seed = 123,
                  bart_args = list(n.trees = 20), include.obj = TRUE)
  }))

  expect_true(is.numeric(W$ps))
  expect_true(all(W$ps > 0 & W$ps < 1))
  expect_true(all(is.finite(W$weights) & W$weights > 0))
  expect_s3_class(W$obj, "stan4bartFit")
  expect_null(attr(W, "Mparts", exact = TRUE))

  #bart2-style control arguments are translated to their stan4bart equivalents
  expect_no_error(suppressWarnings({
    Wt <- weightit(A ~ X2 + X3 + X7 + (1 | cluster),
                   data = test_data, method = "bart", estimand = "ATE",
                   n.chains = 1, n.threads = 1, n.trees = 20, seed = 123)
  }))
  expect_true(all(is.finite(Wt$weights)))

  #Distinguishable from a fit that ignores clustering (X1 also omitted)
  set.seed(123)
  W_none <- weightit(A ~ X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                     data = test_data, method = "bart", estimand = "ATE",
                     n.trees = 20, n.samples = 40, n.burn = 40,
                     n.chains = 1, n.threads = 1)
  expect_not_equal(W$weights, W_none$weights)

  #Balance on the omitted confounder X1 improves when cluster stands in for it
  smd_re <- abs(cobalt::col_w_smd(test_data$X1, test_data$A, W$weights))
  smd_none <- abs(cobalt::col_w_smd(test_data$X1, test_data$A, W_none$weights))
  expect_lt(smd_re, smd_none)
})

test_that("Continuous treatment with random effects (stan4bart)", {
  skip_on_cran()
  skip_if_not_installed("stan4bart")

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  expect_no_error(suppressWarnings({
    W <- weightit(Ac ~ X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + (1 | cluster),
                  data = test_data, method = "bart",
                  iter = 40, warmup = 20, chains = 1, cores = 1, seed = 123,
                  bart_args = list(n.trees = 20), include.obj = TRUE)
  }))

  expect_true(all(is.finite(W$weights) & W$weights > 0))
  expect_s3_class(W$obj, "stan4bartFit")
  expect_null(attr(W, "Mparts", exact = TRUE))
})

test_that("Multi-category treatment with random effects (stan4bart)", {
  skip_on_cran()
  skip_if_not_installed("stan4bart")

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  expect_no_error(suppressWarnings({
    W <- weightit(Am ~ X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + (1 | cluster),
                  data = test_data, method = "bart", estimand = "ATE",
                  iter = 40, warmup = 20, chains = 1, cores = 1, seed = 123,
                  bart_args = list(n.trees = 20), include.obj = TRUE)
  }))

  expect_true(all(is.finite(W$weights) & W$weights > 0))
  expect_type(W$obj, "list")
  expect_s3_class(W$obj[[1L]], "stan4bartFit")
  expect_null(attr(W, "Mparts", exact = TRUE))
})

test_that("Multiple grouping factors are supported (stan4bart)", {
  skip_on_cran()
  skip_if_not_installed("stan4bart")

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  #Two grouping factors, allowed with stan4bart (unlike the old rbart_vi backend)
  expect_no_error(suppressWarnings({
    W <- weightit(A ~ X2 + X3 + X7 + (1 | cluster) + (1 | X6),
                  data = test_data, method = "bart", estimand = "ATE",
                  iter = 40, warmup = 20, chains = 1, cores = 1, seed = 123,
                  bart_args = list(n.trees = 20), include.obj = TRUE)
  }))

  expect_true(all(W$ps > 0 & W$ps < 1))
  expect_true(all(is.finite(W$weights) & W$weights > 0))
  expect_s3_class(W$obj, "stan4bartFit")
})

test_that("Random effects method gating", {
  skip_on_cran()

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  #re_ok gating: a method with re_ok = FALSE errors on a random effects formula
  expect_error(weightit(A ~ X2 + (1 | cluster), data = test_data, method = "cbps"),
               "random effects", ignore.case = TRUE)
})

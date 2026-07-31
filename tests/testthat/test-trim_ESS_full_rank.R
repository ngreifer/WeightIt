# `trim()`, `ESS()`, and `make_full_rank()` are exported and produce numbers, but had
# essentially no assertions: `trim()` was only ever called once with default options,
# `ESS()` only as a test helper, and `make_full_rank()` only indirectly.

skip_on_cran()

eps <- if (capabilities("long.double")) 1e-5 else 1e-3

test_data <- readRDS(test_path("fixtures", "test_data.rds"))

# ---- ESS() -----------------------------------------------------------------

test_that("ESS() is Kish's effective sample size", {
  w <- c(1, 1, 2, 4, 8)

  expect_equal(ESS(w), sum(w)^2 / sum(w^2))

  # Equal weights give the full sample size, and scale invariance holds
  expect_equal(ESS(rep.int(1, 20L)), 20)
  expect_equal(ESS(rep.int(7, 20L)), 20)
  expect_equal(ESS(w), ESS(3 * w))

  # Unequal weights always lose information
  expect_lt(ESS(w), length(w))

  # Zero weights drop out entirely
  expect_equal(ESS(c(w, 0, 0)), ESS(w))

  expect_error(ESS())
  expect_error(ESS("a"))
})

# ---- trim() ----------------------------------------------------------------

test_that("trim.weightit() winsorizes at a quantile and records what it did", {
  d <- test_data

  W <- weightit(A ~ X1 + X2 + X3, data = d, method = "glm", estimand = "ATE")

  expect_message(Wt <- trim(W, at = .9), "Trimming weights to 90%", fixed = TRUE)

  # Upper-only by default: the top decile is pulled down to the 90th percentile and
  # nothing below it moves
  cutoff <- quantile(W$weights, .9, type = 3)
  expect_equal(max(Wt$weights), unname(cutoff), tolerance = eps)
  expect_equal(unname(Wt$weights[W$weights <= cutoff]),
               unname(W$weights[W$weights <= cutoff]), tolerance = eps)
  expect_equal(min(Wt$weights), min(W$weights), tolerance = eps)

  expect_identical(attr(Wt, "trim")$at, .9)
  expect_false(attr(Wt, "trim")$lower)
  expect_false(attr(Wt, "trim")$drop)

  expect_lt(ESS(Wt$weights), nrow(d))
  expect_output(print(Wt), "trimmed")

  # `at = 0` and `at = NULL` are no-ops
  expect_identical(trim(W, at = 0), W)
  expect_identical(trim(W, at = NULL), W)

  expect_error(trim(W, at = -1))
})

test_that("trim(lower = TRUE) trims both tails", {
  d <- test_data

  W <- weightit(A ~ X1 + X2 + X3, data = d, method = "glm", estimand = "ATE")

  expect_message(Wt <- trim(W, at = .9, lower = TRUE),
                 "Trimming weights to 10% and 90%", fixed = TRUE)

  lo <- quantile(W$weights, .1, type = 3)
  hi <- quantile(W$weights, .9, type = 3)

  expect_equal(min(Wt$weights), unname(lo), tolerance = eps)
  expect_equal(max(Wt$weights), unname(hi), tolerance = eps)
  expect_true(attr(Wt, "trim")$lower)

  # Strictly more trimming than the upper-only version
  Wu <- suppressMessages(trim(W, at = .9))
  expect_gt(min(Wt$weights), min(Wu$weights))
})

test_that("trim(drop = TRUE) zeroes the trimmed units instead of winsorizing", {
  d <- test_data

  W <- weightit(A ~ X1 + X2 + X3, data = d, method = "glm", estimand = "ATE")

  expect_message(Wd <- trim(W, at = .9, drop = TRUE),
                 "Setting weights beyond 90% to 0", fixed = TRUE)

  cutoff <- quantile(W$weights, .9, type = 3)

  expect_true(all(Wd$weights[W$weights > cutoff] == 0))
  expect_equal(unname(Wd$weights[W$weights <= cutoff]),
               unname(W$weights[W$weights <= cutoff]), tolerance = eps)
  expect_true(attr(Wd, "trim")$drop)

  # Winsorizing keeps every unit; dropping sets the trimmed ones to 0
  Wu <- suppressMessages(trim(W, at = .9))
  expect_identical(sum(Wu$weights == 0), 0L)
  expect_identical(sum(Wd$weights == 0), sum(W$weights > cutoff))

  # Kish's ESS actually rises when the extreme weights are removed rather than pulled
  # in, since it depends on the variance of the weights that remain
  expect_gt(ESS(Wd$weights), ESS(Wu$weights))
})

test_that("at >= 1 is treated as a count of units rather than a quantile", {
  d <- test_data

  W <- weightit(A ~ X1 + X2 + X3, data = d, method = "glm", estimand = "ATE")

  expect_message(Wt <- trim(W, at = 10), "Trimming the top 10 weights",
                 fixed = TRUE)

  # The 10 largest weights are pulled down to the 11th largest
  o <- order(W$weights, decreasing = TRUE)
  expect_equal(unname(Wt$weights[o[1:10]]),
               rep(unname(W$weights[o[11L]]), 10L), tolerance = eps)

  # A count larger than the sample cannot be honored
  expect_warning(trim(W, at = nrow(d) + 1), "must be less than")
})

test_that("trim.weightit() leaves the focal group alone", {
  d <- test_data

  W <- weightit(A ~ X1 + X2 + X3, data = d, method = "glm", estimand = "ATT")

  expect_message(Wt <- trim(W, at = .9), "where treat is not")

  # ATT weights are 1 in the focal group and must stay there
  expect_true(all(Wt$weights[W$treat == W$focal] == 1))
  expect_ATT_weights_okay(Wt, tolerance = eps)

  # Only the non-focal group is affected
  non_focal <- W$treat != W$focal
  expect_equal(max(Wt$weights[non_focal]),
               unname(quantile(W$weights[non_focal], .9, type = 3)),
               tolerance = eps)
})

test_that("trim.default() works on a bare weight vector", {
  d <- test_data

  W <- weightit(A ~ X1 + X2 + X3, data = d, method = "glm", estimand = "ATE")
  w <- W$weights

  expect_message(wt <- trim(w, at = .9), "Trimming weights to 90%", fixed = TRUE)
  expect_equal(unname(wt), unname(suppressMessages(trim(W, at = .9))$weights),
               tolerance = eps)

  # Supplying `treat` reproduces the focal-group protection of the weightit method
  W_att <- weightit(A ~ X1 + X2 + X3, data = d, method = "glm", estimand = "ATT")
  expect_message(wt_att <- trim(W_att$weights, at = .9, treat = d$A),
                 "where treat is not")
  expect_true(all(wt_att[d$A == 1] == 1))

  # Constant weights cannot be trimmed, and say so rather than failing
  expect_warning(out <- trim(rep.int(1, 100L), at = .9), "all the same")
  expect_identical(out, rep.int(1, 100L))

  expect_error(trim("a", at = .9), "numeric vector")
})

# ---- make_full_rank() ------------------------------------------------------

test_that("make_full_rank() drops redundant columns", {
  X <- cbind(a = c(1, 2, 3, 4), b = c(2, 4, 6, 8), c = c(1, 0, 1, 0))

  # `b` is a multiple of `a`
  out <- make_full_rank(X)
  expect_identical(colnames(out), c("a", "c"))
  expect_identical(qr(out)$rank, ncol(out))

  # A data frame keeps its class
  out_df <- make_full_rank(as.data.frame(X))
  expect_s3_class(out_df, "data.frame")
  expect_named(out_df, c("a", "c"))

  # An already-full-rank matrix is returned unchanged
  expect_identical(make_full_rank(X[, c("a", "c")]), X[, c("a", "c")])
})

test_that("make_full_rank(with.intercept =) controls whether a constant is redundant", {
  # `d` is collinear with the intercept: two dummies that sum to 1
  X <- cbind(a = c(1, 2, 3, 4), d0 = c(1, 1, 0, 0), d1 = c(0, 0, 1, 1))

  # With an implicit intercept, one of the two dummies is redundant
  with_int <- make_full_rank(X, with.intercept = TRUE)
  expect_identical(ncol(with_int), 2L)

  # Without one, both dummies are estimable
  without_int <- make_full_rank(X, with.intercept = FALSE)
  expect_identical(ncol(without_int), 3L)
  expect_identical(colnames(without_int), colnames(X))
})

test_that("make_full_rank() rejects bad input", {
  expect_error(make_full_rank("a"), "numeric matrix or data frame")
  expect_error(make_full_rank(matrix(letters[1:4], 2L, 2L)), "numeric matrix")
  expect_error(make_full_rank(data.frame(a = 1:2, b = letters[1:2])),
               "must be numeric")
  expect_error(make_full_rank(matrix(c(1, NA, 3, 4), 2L, 2L)))
})

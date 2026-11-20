# Small internal helpers that support the M-estimation machinery. None of them is
# exported, and all three were previously exercised only through their default
# arguments, which hid bugs that surface as soon as a caller deviates from those
# defaults.

skip_on_cran()

# ---- .gradient() -----------------------------------------------------------

test_that(".gradient() differentiates an arbitrary subset of the parameters", {
  skip_if_not_installed("patrick")

  #Both components have gradients that are easy to write down
  f <- function(x) c(sum(x^2), prod(x))
  x <- c(1, 2, 3)

  analytic <- rbind(2 * x,
                    c(x[2L] * x[3L], x[1L] * x[3L], x[1L] * x[2L]))

  patrick::with_parameters_test_that(
    ".method = {method}",
    {
      #All parameters, the default
      expect_equal(unname(.gradient(f, x, .method = method)), analytic,
                   tolerance = 1e-5, info = method)

      #A subset that does not start at the first parameter. This used to error
      #("object 'jacob' not found") for "fd" and to return an extra all-NA column
      #for "richardson".
      expect_equal(unname(.gradient(f, x, .parm = 2:3, .method = method)),
                   analytic[, 2:3],
                   tolerance = 1e-5, info = method)

      #An out-of-order subset: the columns follow `.parm`, not `x`
      expect_equal(unname(.gradient(f, x, .parm = c(3L, 1L), .method = method)),
                   analytic[, c(3L, 1L)],
                   tolerance = 1e-5, info = method)

      #A single parameter
      expect_equal(unname(.gradient(f, x, .parm = 2L, .method = method)),
                   analytic[, 2L, drop = FALSE],
                   tolerance = 1e-5, info = method)
    },
    .cases = data.frame(method = c("fd", "richardson"))
  )
})

test_that(".gradient() names its columns after the differentiated parameters", {
  f <- function(x) sum(x^2)
  x <- c(a = 1, b = 2, c = 3)

  expect_identical(colnames(.gradient(f, x, .parm = c(3L, 1L))), c("c", "a"))
  expect_identical(colnames(.gradient(f, x, .parm = c(3L, 1L), .method = "richardson")),
                   c("c", "a"))
})

# ---- .vec2list() -----------------------------------------------------------

test_that(".vec2list() splits a vector into blocks of the given lengths", {
  expect_identical(.vec2list(1:6, c(2L, 3L, 1L)),
                   list(1:2, 3:5, 6L))

  #A zero-length block yields a zero-length element. Building the index as
  #`start:end` instead would count backwards and return two elements.
  expect_identical(.vec2list(1:5, c(2L, 0L, 3L)),
                   list(1:2, integer(), 3:5))

  expect_identical(.vec2list(integer(), 0L), list(integer()))
})

# ---- .block_diag() ---------------------------------------------------------

test_that(".block_diag() stacks matrices block-diagonally", {
  out <- .block_diag(list(matrix(1, 2L, 2L), matrix(3, 1L, 1L)))

  expect_identical(dim(out), c(3L, 3L))
  expect_equal(out, matrix(c(1, 1, 0,
                             1, 1, 0,
                             0, 0, 3), 3L, 3L))

  #A block with no columns contributes rows but no columns, and must not write
  #into its neighbors
  out0 <- .block_diag(list(matrix(1, 2L, 2L),
                           matrix(0, 2L, 0L),
                           matrix(3, 1L, 1L)))

  expect_identical(dim(out0), c(5L, 3L))
  expect_true(all(out0[3:4, ] == 0))
  expect_equal(out0[1:2, 1:2], matrix(1, 2L, 2L))
  expect_equal(out0[5L, 3L], 3)
})

# ---- w.quantile() ----------------------------------------------------------

test_that("w.quantile() is the weighted version of quantile(type = 2)", {
  set.seed(123)
  x <- rnorm(101)
  p <- seq(0, 1, by = .05)

  #Unweighted, it agrees with `quantile(type = 2)`, including at ties and at
  #probabilities that fall in a flat region of the empirical CDF
  expect_equal(w.quantile(x, p), unname(quantile(x, p, type = 2)))
  expect_equal(w.quantile(1:4, .5), 2.5)
  expect_equal(w.quantile(c(3, 1, 2, 2, 4), c(.2, .5, .6, .9)),
               unname(quantile(c(3, 1, 2, 2, 4), c(.2, .5, .6, .9), type = 2)))

  #Integer weights are the same as replicating units. The probabilities compared
  #with `quantile()` are exact in binary: `quantile()` compares `n * p` with its
  #floor with almost no tolerance, so, e.g., `.15 * 180` misses the tie that
  #`w.quantile()` correctly finds
  w <- sample(1:3, length(x), replace = TRUE)
  expect_equal(w.quantile(x, p, w), w.quantile(rep(x, w), p))

  p2 <- c(0, .125, .25, .5, .75, 1)
  expect_equal(w.quantile(x, p2, w), unname(quantile(rep(x, w), p2, type = 2)))

  #A flat region from a zero weight spans to the next value with weight
  expect_equal(w.quantile(1:4, .5, w = c(1, 1, 0, 2)), 3)

  #Missing values are removed before the empirical CDF is computed
  expect_equal(w.quantile(c(1, 2, NA), c(0, .5, 1), na.rm = TRUE), c(1, 1.5, 2))
  expect_equal(w.quantile(c(1, 2, NA), .5), NA_real_)
})

# ---- .isoregw() ------------------------------------------------------------

test_that(".isoregw() is weighted isotonic regression", {
  set.seed(123)
  x <- runif(500)
  y <- rbinom(500, 1, x)

  #Unweighted, it agrees with `stats::isoreg()`
  expect_equal(.isoregw(x, y), isoreg(x, y)$yf[order(order(x))])

  #Integer weights are the same as replicating units
  w <- sample(1:3, length(x), replace = TRUE)
  rep_fit <- isoreg(rep(x, w), rep(y, w))$yf[order(order(rep(x, w)))]
  expect_equal(.isoregw(x, y, w), rep_fit[cumsum(w)])

  #Zero weights, including runs of them, give finite, monotone fits
  w0 <- runif(length(x))
  w0[sample.int(length(x), 200)] <- 0
  f <- .isoregw(x, y, w0)

  expect_true(all(is.finite(f)))
  expect_false(is.unsorted(f[order(x)]))
})

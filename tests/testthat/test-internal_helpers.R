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

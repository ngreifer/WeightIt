test_that("Binary treatment", {
  skip_on_cran()
  skip_if_not_installed("rootSolve")
  skip_if_not_installed("cobalt")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  expect_no_condition({
    W0 <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                   data = test_data, method = "cbps", estimand = "ATE",
                   include.obj = TRUE, solver = "optim")
  })

  sw.opts <- c(FALSE, TRUE)
  over.opts <- c("exact", "twostep", "cont")
  estimand.opts <- c("ATE", "ATT", "ATC", "ATO")
  link.opts <- c("logit", "probit", "loglog", "cauchit", "softplus")

  # The full link.opts grid is only crossed with the cheap "exact" solver,
  # since link doesn't affect the (expensive) over-identified GMM machinery
  # exercised by "twostep"/"cont"; those use a single representative link.
  weight.mat <- matrix(nrow = nrow(test_data),
                       ncol = length(sw.opts) * length(estimand.opts) *
                         (length(link.opts) + 2L))
  colnames(weight.mat) <- rep("", ncol(weight.mat))
  k <- 1

  for (sw in sw.opts) {
    for (over in over.opts) {
      links_to_test <- if (over == "exact") link.opts else "logit"

      for (estimand in estimand.opts) {
        for (link in links_to_test) {
          test_that(sprintf("CBPS: sw = %s, over = %s, estimand = %s, link = %s", sw, over, estimand, link), {
            suppressWarnings({
            W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                          data = test_data, method = "cbps", estimand = estimand,
                          over = over != "exact", twostep = if (over == "exact") NULL else over == "twostep",
                          link = link,
                          s.weights = if (sw) "SW" else NULL,
                          include.obj = TRUE, solver = "multiroot")
            })

            if (over == "exact") {
              expect_M_parts_okay(W, tolerance = eps)
              expect_equal(cobalt::col_w_smd(W$covs, W$treat, W$weights,
                                             s.weights = W$s.weights),
                           0 * cobalt::col_w_smd(W$covs, W$treat,
                                                 s.weights = W$s.weights),
                           expected.label = "all 0s",
                           tolerance = eps)
            }

            expect_true(is.numeric(W$ps))
            expect_false(is_null(W$obj))

            if (estimand %in% c("ATT", "ATC")) {
              expect_ATT_weights_okay(W, tolerance = eps)
            }

            for (i in 0:1) {
              e <- {
                if (estimand == "ATT" && i == 1) expect_equal
                else if (estimand == "ATC" && i == 0) expect_equal
                else expect_not_equal
              }

              e(unname(W$weights[W$treat == i]),
                rep(1, sum(W$treat == i)),
                label = sprintf("%s weights", i),
                expected.label = "all 1s",
                tolerance = eps)
            }

            if (link != "logit" || estimand != "ATO") {
              for (i in seq_len(k)) {
                expect_not_equal(unname(W$weights), weight.mat[,i],
                                 expected.label = sprintf("weights for model %s", colnames(weight.mat)[i]),
                                 tolerance = min(1e-8, eps))
              }
            }

            n <- sprintf("W_%s_%s_%s_%s", sw, over, estimand, link)
            weight.mat[,k] <<- W$weights
            colnames(weight.mat)[k] <<- n
            k <<- k + 1
          })
        }
      }
    }
  }

  # Estimands
  expect_error({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "cbps", estimand = "ATM")
  }, "not an allowable estimand", ignore.case = TRUE)

  #Non-full rank
  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 +
                    I(1 - X5) + I(X9 * 2),
                  data = test_data, method = "cbps", estimand = "ATE",
                  include.obj = TRUE)
  })

  expect_equal(W$weights, W0$weights, tolerance = eps)
})

test_that("Multi-category treatment", {
  skip_on_cran()
  skip_if_not_installed("rootSolve")
  skip_if_not_installed("cobalt")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  expect_no_condition({
    W0 <- weightit(Am ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                   data = test_data, method = "cbps", estimand = "ATE",
                   include.obj = TRUE, solver = "optim")
  })

  sw.opts <- c(FALSE, TRUE)
  over.opts <- c("exact", "twostep", "cont")
  estimand.opts <- c("ATE", "ATT", "ATO")

  weight.mat <- matrix(nrow = nrow(test_data),
                       ncol = length(sw.opts) * length(over.opts) *
                         length(estimand.opts))
  colnames(weight.mat) <- rep("", ncol(weight.mat))
  k <- 1

  for (sw in sw.opts) {
    for (over in over.opts) {
      for (estimand in estimand.opts) {
        test_that(sprintf("CBPS: sw = %s, over = %s, estimand = %s", sw, over, estimand), {
          suppressWarnings({
          W <- weightit(Am ~ X1 + X2 + X3 + X4 + X5,
                        data = test_data, method = "cbps", estimand = estimand,
                        over = over != "exact", twostep = if (over == "exact") NULL else over == "twostep",
                        focal = if (estimand == "ATT") "T" else NULL,
                        s.weights = if (sw) "SW" else NULL,
                        include.obj = TRUE, solver = "multiroot")
          })

          if (over == "exact") {
            expect_M_parts_okay(W, tolerance = eps)
            for (tt in combn(levels(W$treat), 2, simplify = FALSE)) {
              in_tt <- W$treat %in% tt
              expect_equal(cobalt::col_w_smd(W$covs[in_tt,], W$treat[in_tt], W$weights[in_tt],
                                             s.weights = W$s.weights[in_tt]),
                           0 * cobalt::col_w_smd(W$covs[in_tt,], W$treat[in_tt],
                                                 s.weights = W$s.weights[in_tt]),
                           label = sprintf("SMDs for %s", paste(tt, collapse = " vs. ")),
                           expected.label = "all 0s",
                           tolerance = eps)
            }
          }

          expect_true(is_null(W$ps))
          expect_false(is_null(W$obj))

          if (estimand %in% c("ATT", "ATC")) {
            expect_ATT_weights_okay(W, tolerance = eps)
          }

          for (i in levels(W$treat)) {
            e <- {
              if (estimand == "ATT" && i == W$focal) expect_equal
              else expect_not_equal
            }

            e(unname(W$weights[W$treat == i]),
              rep(1, sum(W$treat == i)),
              label = sprintf("%s weights", i),
              expected.label = "all 1s",
              tolerance = eps)
          }

          for (i in seq_len(k - 1)) {
            expect_not_equal(unname(W$weights), weight.mat[,i],
                             expected.label = sprintf("weights for %s", colnames(weight.mat)[i]),
                             tolerance = min(1e-8, eps))
          }

          n <- sprintf("W_%s_%s_%s", sw, over, estimand)
          colnames(weight.mat)[k] <<- n
          weight.mat[,k] <<- W$weights
          k <<- k + 1
        })
      }
    }
  }
})

test_that("Continuous treatment", {
  skip_on_cran()
  skip_if_not_installed("rootSolve")
  skip_if_not_installed("cobalt")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  W0 <- expect_no_unexpected_warning(
    weightit(Ac ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
             data = test_data, method = "cbps",
             include.obj = TRUE, solver = "optim")
  )

  sw.opts <- c(FALSE, TRUE)
  over.opts <- c("exact", "twostep", "cont")

  weight.mat <- matrix(nrow = nrow(test_data),
                       ncol = length(sw.opts) *
                         length(over.opts))
  colnames(weight.mat) <- rep("", ncol(weight.mat))

  k <- 1

  for (sw in sw.opts) {
    for (over in over.opts) {
      test_that(sprintf("CBPS: sw = %s, over = %s", sw, over), {
        suppressWarnings({
        W <- weightit(Ac ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                      data = test_data, method = "cbps",
                      over = over != "exact", twostep = if (over == "exact") NULL else over == "twostep",
                      s.weights = if (sw) "SW" else NULL,
                      include.obj = TRUE, solver = "multiroot")
        })

        if (over == "exact") {
          # For a continuous treatment the balance conditions cannot always be solved:
          # at some samples the estimating equations have no root, the optimizer settles
          # at a stationary point where the objective (the norm of the mean estimating
          # function) is far from 0, and the weighted covariances are not exact. Whether
          # that happens depends on the sample and the specification -- with this
          # generator the full-covariate spec solves at n = 900 and n = 1000 but not at
          # n = 750 or n = 1200 -- so the exact-balance assertion is conditional on the
          # objective saying the conditions were actually met. `weightit()` warns in the
          # other branch; see `.check_cbps_converged()`.
          if (W$obj$value < 1e-8) {
            expect_equal(cobalt::col_w_cov(W$covs, W$treat, W$weights, std = TRUE,
                                           s.weights = W$s.weights),
                         0 * cobalt::col_w_cov(W$covs, W$treat, std = TRUE,
                                               s.weights = W$s.weights),
                         expected.label = "all 0s",
                         tolerance = eps)
          }
          else {
            # Unsolved, but the weights must still be usable and must improve on the
            # unweighted covariances
            expect_true(all(is.finite(W$weights) & W$weights > 0))
            expect_lt(max(abs(cobalt::col_w_cov(W$covs, W$treat, W$weights, std = TRUE,
                                                s.weights = W$s.weights))),
                      max(abs(cobalt::col_w_cov(W$covs, W$treat, std = TRUE,
                                                s.weights = W$s.weights))))
          }
        }

        expect_true(is_null(W$ps))
        expect_false(is_null(W$obj))

        for (i in seq_len(k - 1)) {
          expect_not_equal(unname(W$weights), weight.mat[,i],
                           expected.label = sprintf("weights for %s", colnames(weight.mat)[i]),
                           tolerance = min(1e-8, eps))
        }

        n <- sprintf("W_%s_%s", sw, over)
        colnames(weight.mat)[k] <<- n
        weight.mat[,k] <<- W$weights
        k <<- k + 1
      })
    }
  }

  #Non-full rank
  W <- expect_no_unexpected_warning(
    weightit(Ac ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 +
               I(1 - X5) + I(X9 * 2),
             data = test_data, method = "cbps",
             include.obj = TRUE, solver = "optim")
  )

  expect_equal(W$weights, W0$weights, tolerance = eps)
})

test_that("unsolved balance conditions are reported rather than silently returned", {
  skip_on_cran()
  skip_if_not_installed("cobalt")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  # With `over = FALSE` the system is exactly identified, so the objective -- the norm of
  # the mean estimating function -- has to reach ~0 for the conditions to be solved.
  # `optim()` reports `convergence == 0` whether or not it got there, so the objective is
  # the only signal, and it used to be ignored. Regression test for that: whenever the
  # objective says the conditions were not met, there must be a warning saying so.
  obj_value <- function(W) W$obj$value

  specs <- list(
    full = Ac ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
    small = Ac ~ X1 + X2 + X3,
    binary = A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
    multi = Am ~ X1 + X2 + X3 + X4 + X5
  )

  n_warned <- 0L

  for (nm in names(specs)) {
    ws <- character()
    W <- withCallingHandlers(
      weightit(specs[[nm]], data = test_data, method = "cbps", include.obj = TRUE),
      warning = function(w) {
        ws <<- c(ws, conditionMessage(w))
        invokeRestart("muffleWarning")
      })

    warned <- any(grepl("could not be solved", ws, fixed = TRUE))
    stalled <- obj_value(W) > 1e-5

    expect_identical(warned, stalled,
                     label = sprintf("warning matches the objective for %s (obj = %.2e)",
                                     nm, obj_value(W)))

    if (stalled) n_warned <- n_warned + 1L

    # The objective is not just a solver diagnostic: it tracks the balance actually
    # achieved. Checked on the continuous specs, where `col_w_cov()` is the relevant
    # statistic (binary and multi-category use standardized mean differences instead).
    if (nm %in% c("full", "small")) {
      achieved <- max(abs(cobalt::col_w_cov(W$covs, W$treat, W$weights, std = TRUE,
                                            s.weights = W$s.weights)))

      if (stalled) expect_gt(achieved, eps)
      else expect_lt(achieved, eps)
    }
  }

  # The continuous full-covariate spec is the case that motivated this check, so at least
  # one of the four above must exercise the warning path; otherwise this test is vacuous
  expect_gt(n_warned, 0L)

  # `over = TRUE` is over-identified, so a non-zero objective is expected and must not
  # warn
  expect_no_warning(weightit(specs$full, data = test_data, method = "cbps",
                             over = TRUE, include.obj = TRUE))

  # Hitting the iteration limit is still reported separately, with its own message
  expect_warning(weightit(specs$full, data = test_data, method = "cbps",
                          maxit = 1L, solver = "optim"),
                 "higher value of")
})

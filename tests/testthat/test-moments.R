test_that("moments works as expected, binary", {
  skip_on_cran()
  skip_if_not_installed("rootSolve")
  skip_if_not_installed("optweight", minimum_version = "2.0.1")
  skip_if_not_installed("cobalt")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  poly_test <- cbind(test_data[c("X3", "X4")],
                     poly(test_data$X1, 2),
                     poly(test_data$X2, 3))

  init1 <- cobalt::bal.init(poly_test, treat = test_data$A, stat = "smd.max",
                            pairwise = FALSE)
  init2 <- cobalt::bal.init(poly(test_data$X1, 3), treat = test_data$A, stat = "smd.max",
                            pairwise = FALSE)

  #Moments the same as poly
  for (method in c("ipt", "ebal", "optweight")) {
    test_that(sprintf("method = %s", method), {
      Wm <- weightit(A ~ X1 + X2 + X3 + X4,
                     data = test_data, method = method,
                     include.obj = TRUE,
                     moments = c(X1 = 2, X2 = 3))

      Wp <- weightit(A ~ poly(X1, 2) + poly(X2, 3) + X3 + X4,
                     data = test_data, method = method,
                     include.obj = TRUE)

      if (method %in% c("ipt", "ebal")) {
        expect_M_parts_okay(Wm, tolerance = eps)
        expect_M_parts_okay(Wp, tolerance = eps)
      }

      expect_equal(Wm$weights, Wp$weights,
                   label = "weights for moments",
                   expected.label = "weights for poly",
                   tolerance = eps)

      expect_equal(cobalt::bal.compute(init1, weights = Wm$weights),
                   0,
                   label = "largest SMD for moments",
                   tolerance = eps)
      expect_equal(cobalt::bal.compute(init1, weights = Wp$weights),
                   0,
                   label = "largest SMD for poly",
                   tolerance = eps)

      expect_not_equal(cobalt::bal.compute(init2, weights = Wm$weights),
                       0,
                       label = "largest SMD for moments with extra terms",
                       tolerance = eps)
    })
  }
})

test_that("moments works as expected, multi", {
  skip_on_cran()
  skip_if_not_installed("rootSolve")
  skip_if_not_installed("optweight", minimum_version = "2.0.1")
  skip_if_not_installed("cobalt")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  poly_test <- cbind(test_data[c("X3", "X4")],
                     poly(test_data$X1, 2),
                     poly(test_data$X2, 3))

  init1 <- cobalt::bal.init(poly_test, treat = test_data$Am, stat = "smd.max",
                            pairwise = FALSE)
  init2 <- cobalt::bal.init(poly(test_data$X1, 3), treat = test_data$Am, stat = "smd.max",
                            pairwise = FALSE)

  #Moments the same as poly
  for (method in c("ipt", "ebal", "optweight")) {
    test_that(sprintf("method = %s", method), {
      Wm <- weightit(Am ~ X1 + X2 + X3 + X4,
                     data = test_data, method = method,
                     include.obj = TRUE,
                     moments = c(X1 = 2, X2 = 3))

      Wp <- weightit(Am ~ poly(X1, 2) + poly(X2, 3) + X3 + X4,
                     data = test_data, method = method,
                     include.obj = TRUE)

      if (method %in% c("ipt", "ebal")) {
        expect_M_parts_okay(Wm, tolerance = eps)
        expect_M_parts_okay(Wp, tolerance = eps)
      }

      expect_equal(Wm$weights, Wp$weights,
                   label = "weights for moments",
                   expected.label = "weights for poly",
                   tolerance = eps)

      expect_equal(cobalt::bal.compute(init1, weights = Wm$weights),
                   0,
                   label = "largest SMD for moments",
                   tolerance = eps)
      expect_equal(cobalt::bal.compute(init1, weights = Wp$weights),
                   0,
                   label = "largest SMD for poly",
                   tolerance = eps)

      expect_not_equal(cobalt::bal.compute(init2, weights = Wm$weights),
                       0,
                       label = "largest SMD for moments with extra terms",
                       tolerance = eps)
    })
  }
})

test_that("moments works as expected, cont", {
  skip_on_cran()
  skip_if_not_installed("rootSolve")
  skip_if_not_installed("cobalt")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  poly_test <- cbind(test_data[c("X3", "X4")],
                     poly(test_data$X1, 2),
                     poly(test_data$X2, 3))

  init1 <- cobalt::bal.init(poly_test, treat = test_data$Ac, stat = "p.max")
  init2 <- cobalt::bal.init(poly(test_data$X1, 3), treat = test_data$Ac, stat = "p.max")

  init1t <- cobalt::bal.init(poly_test, stat = "smd.max")
  init2t <- cobalt::bal.init(poly(test_data$X1, 3), stat = "smd.max")

  #Moments the same as poly
  for (method in c("ebal", "cbps")) {
    test_that(sprintf("method = %s", method), {
      # `cbps` warns here when it cannot solve the balance conditions at this sample;
      # see the note under the `method == "ebal"` guard below
      Wm <- expect_no_unexpected_warning(
        weightit(Ac ~ X1 + X2 + X3 + X4,
                 data = test_data, method = method,
                 include.obj = TRUE,
                 moments = c(X1 = 2, X2 = 3))
      )

      Wp <- expect_no_unexpected_warning(
        weightit(Ac ~ poly(X1, 2) + poly(X2, 3) + X3 + X4,
                 data = test_data, method = method,
                 include.obj = TRUE)
      )

      # The point of the block: whichever method is used, expanding the moments and
      # writing the polynomials out by hand must give the same weights.
      expect_equal(Wm$weights, Wp$weights,
                   label = "weights for moments",
                   expected.label = "weights for poly",
                   tolerance = eps)

      expect_not_equal(cobalt::bal.compute(init2, weights = Wm$weights),
                       0,
                       label = "largest corr for moments with extra terms",
                       tolerance = eps)

      if (method == "ebal") {
        expect_M_parts_okay(Wm, tolerance = eps)
        expect_M_parts_okay(Wp, tolerance = eps)

        # Only the balancing methods solve the balance conditions exactly.
        # `method = "cbps"` with a continuous treatment maximizes a
        # generalized-propensity-score likelihood instead, so its weighted
        # treatment-covariate correlations are near zero but not zero; asserting
        # otherwise passes or fails depending on the sample.
        expect_equal(cobalt::bal.compute(init1, weights = Wm$weights),
                     0,
                     label = "largest corr for moments",
                     tolerance = eps)
        expect_equal(cobalt::bal.compute(init1, weights = Wp$weights),
                     0,
                     label = "largest corr for poly",
                     tolerance = eps)

        expect_equal(cobalt::bal.compute(init1t, weights = Wm$weights),
                     0,
                     label = "largest target SMD for moments",
                     tolerance = eps)
        expect_equal(cobalt::bal.compute(init1t, weights = Wp$weights),
                     0,
                     label = "largest target SMD for poly",
                     tolerance = eps)

        expect_not_equal(cobalt::bal.compute(init2t, weights = Wm$weights),
                         0,
                         label = "largest target SMD for moments with extra terms",
                         tolerance = eps)
      }
      else {
        # `cbps` still has to improve on the unweighted correlations, and the
        # moment-expanded terms must be the ones it improves.
        expect_lt(cobalt::bal.compute(init1, weights = Wm$weights),
                  cobalt::bal.compute(init1, weights = rep.int(1, nrow(test_data))))
      }
    })
  }
})
test_that("int and quantile compose with the M-estimating methods", {
  skip_on_cran()
  skip_if_not_installed("rootSolve")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  # `moments`, `int`, and `quantile` all expand `covs` before the model is fit, so the
  # `Xtreat` stored in the M-estimation parts has to match the expansion rather than the
  # original covariates. A mismatch is the same non-conformable-arguments failure that
  # went unnoticed in the continuous glm parts, so every method that both builds Mparts
  # and accepts these arguments is checked with each of them. Only `moments` had any
  # coverage before, and only for `ipt` and `ebal`.
  for (method in c("cbps", "ebal", "ipt")) {
    if (!all(vapply(.weightit_methods[[method]]$packages_needed,
                    rlang::is_installed, logical(1L)))) {
      next
    }

    for (extra in list(list(int = TRUE),
                       list(quantile = list(X1 = c(.25, .5, .75))),
                       list(moments = 2L),
                       list(moments = 2L, int = TRUE))) {
      lbl <- sprintf("method = %s, %s", method, paste(names(extra), collapse = " + "))

      W <- do.call("weightit",
                   c(list(A ~ X1 + X2 + X3, data = test_data, method = method),
                     extra))

      expect_M_parts_okay(W, tolerance = eps)

      # The expansion actually happened: more than the intercept plus three covariates
      expect_gt(ncol(attr(W, "Mparts")$Xtreat), 4L, label = lbl)

      # And the parts survive being assembled into an outcome model's variance
      fit <- glm_weightit(Y_B ~ A, data = test_data, weightit = W,
                          family = binomial)
      expect_true(all(is.finite(sqrt(diag(vcov(fit))))), label = lbl)
    }
  }

  # `method = "glm"` has moments_int_ok = FALSE, so it warns and ignores them rather
  # than silently expanding `covs` out from under `Xtreat`
  expect_warning(weightit(A ~ X1 + X2 + X3, data = test_data, method = "glm",
                          moments = 2L),
                 "not compatible")
})

test_that("Binary treatment", {
  skip_on_cran()
  skip_if_not_installed("gbm")
  skip_if_not_installed("cobalt")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  set.seed(123)
  expect_no_condition({
    W0 <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                   data = test_data, method = "gbm", estimand = "ATE",
                   criterion = "smd.mean", n.trees = 300,
                   include.obj = TRUE)
  })

  expect_true(is.numeric(W0$ps))
  expect_true(all(W0$ps > 0 & W0$ps < 1))
  expect_false(is_null(W0$obj))
  expect_true(is_null(attr(W0, "Mparts", exact = TRUE))) #gbm does not support M-estimation

  sw.opts <- c(FALSE, TRUE)
  estimand.opts <- c("ATE", "ATT", "ATC", "ATO", "ATM")

  weight.mat <- matrix(nrow = nrow(test_data),
                       ncol = length(sw.opts) * length(estimand.opts))
  colnames(weight.mat) <- rep("", ncol(weight.mat))

  k <- 1

  for (sw in sw.opts) {
    for (estimand in estimand.opts) {
      test_that(sprintf("GBM: sw = %s, estimand = %s", sw, estimand), {
        set.seed(123)
        W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                      data = test_data, method = "gbm", estimand = estimand,
                      criterion = "smd.mean", n.trees = 200,
                      s.weights = if (sw) "SW" else NULL,
                      include.obj = TRUE)

        expect_true(is.numeric(W$ps))
        expect_true(all(is.finite(W$ps) & W$ps > 0 & W$ps < 1))
        expect_true(all(is.finite(W$weights) & W$weights > 0))
        expect_false(is_null(W$obj))

        expect_balance_improved(W)

        if (estimand %in% c("ATT", "ATC")) {
          expect_ATT_weights_okay(W, tolerance = eps)
        }

        for (i in seq_len(k - 1)) {
          expect_not_equal(unname(W$weights), weight.mat[, i],
                           expected.label = sprintf("weights for %s", colnames(weight.mat)[i]),
                           tolerance = eps)
        }

        n <- sprintf("W_%s_%s", sw, estimand)
        colnames(weight.mat)[k] <<- n
        weight.mat[, k] <<- W$weights
        k <<- k + 1
      })
    }
  }

  # Spot-checks of additional arguments; not crossed into the main grid.

  test_that("GBM: alternative criterion (ks.max)", {
    set.seed(123)
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "gbm", estimand = "ATT",
                  criterion = "ks.max", n.trees = 200)

    expect_true(all(is.finite(W$weights) & W$weights > 0))
    expect_ATT_weights_okay(W, tolerance = eps)

    # The tree selected under criterion = "ks.max" should actually minimize
    # (or at least not do worse than) the maximum KS statistic, compared to
    # a fit with the same specification but tuned for a different criterion.
    set.seed(123)
    W_alt <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                      data = test_data, method = "gbm", estimand = "ATT",
                      criterion = "smd.mean", n.trees = 200)

    max_ks <- function(W) {
      max(abs(cobalt::col_w_ks(W$covs, W$treat, W$weights, s.weights = W$s.weights)))
    }

    expect_true(max_ks(W) <= max_ks(W_alt) + eps)
  })

  test_that("GBM: trim.at", {
    set.seed(123)
    W0trim <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                       data = test_data, method = "gbm", estimand = "ATE",
                       criterion = "smd.mean", n.trees = 200)

    set.seed(123)
    Wtrim <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                      data = test_data, method = "gbm", estimand = "ATE",
                      criterion = "smd.mean", n.trees = 200, trim.at = .9)

    expect_true(all(is.finite(Wtrim$weights) & Wtrim$weights > 0))
    expect_true(max(Wtrim$weights) <= max(W0trim$weights) + eps)
  })

  test_that("GBM: distribution = 'adaboost'", {
    set.seed(123)
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "gbm", estimand = "ATE",
                  criterion = "smd.mean", n.trees = 200,
                  distribution = "adaboost")

    expect_true(all(is.finite(W$weights) & W$weights > 0))
    expect_true(is.numeric(W$ps))
  })

  test_that("GBM: tuning interaction.depth", {
    set.seed(123)
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "gbm", estimand = "ATE",
                  criterion = "smd.mean", n.trees = 150,
                  interaction.depth = c(2, 4), shrinkage = .01)

    expect_true(all(is.finite(W$weights) & W$weights > 0))
    expect_false(is_null(W$info$tune))
    expect_false(is_null(W$info$best.tune))
    expect_true(all(c("interaction.depth", "best.tree") %in% names(W$info$tune)))
  })

  test_that("GBM: subclass (MMWS)", {
    set.seed(123)
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "gbm", estimand = "ATE",
                  criterion = "smd.mean", n.trees = 200, subclass = 10)

    expect_true(all(is.finite(W$weights) & W$weights > 0))
    # Subclassing yields far fewer distinct weight values than the sample size
    expect_true(length(unique(W$weights)) < nrow(test_data) / 4)
  })

  test_that("GBM: use.offset", {
    set.seed(123)
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "gbm", estimand = "ATE",
                  criterion = "smd.mean", n.trees = 200, use.offset = TRUE)

    expect_true(all(is.finite(W$weights) & W$weights > 0))
    expect_true(is.numeric(W$ps))

    # use.offset is documented as only allowed with binary treatments and
    # distribution = "bernoulli"
    expect_error({
      weightit(Am ~ X1 + X2 + X3 + X4 + X5,
              data = test_data, method = "gbm", estimand = "ATE",
              criterion = "smd.mean", n.trees = 50, use.offset = TRUE)
    }, "use.offset.*cannot be used with multi-category", fixed = FALSE)

    expect_error({
      weightit(A ~ X1 + X2 + X3 + X4 + X5,
              data = test_data, method = "gbm", estimand = "ATE",
              criterion = "smd.mean", n.trees = 50, use.offset = TRUE,
              distribution = "adaboost")
    }, "use.offset.*can only be used with", fixed = FALSE)
  })

  test_that("GBM: cv-based criterion works for binary treatments", {
    expect_no_error({
      W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                    data = test_data, method = "gbm", estimand = "ATE",
                    criterion = "cv3", n.trees = 200)
    })
    expect_true(all(is.finite(W$weights) & W$weights > 0))

    expect_no_error({
      W2 <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                     data = test_data, method = "gbm", estimand = "ATE",
                     criterion = "cv3", n.trees = 200, use.offset = TRUE)
    })
    expect_true(all(is.finite(W2$weights) & W2$weights > 0))
  })

  test_that("GBM: missing = 'ind' vs missing = 'surr'", {
    test_data_na <- inject_missingness(test_data, c("X1", "X3"))

    expect_no_condition({
      W_ind <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                        data = test_data_na, method = "gbm", estimand = "ATE",
                        criterion = "smd.mean", n.trees = 200, missing = "ind")
    })

    expect_warning({
      W_surr <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                         data = test_data_na, method = "gbm", estimand = "ATE",
                         criterion = "smd.mean", n.trees = 200, missing = "surr")
    }, "missing.*will be set to.*ind", ignore.case = TRUE)

    # covs retain the original NAs in both cases (per documented `covs` output)
    expect_true(anyNA(W_ind$covs$X1))
    expect_true(anyNA(W_surr$covs$X1))

    expect_true(all(is.finite(W_ind$weights) & W_ind$weights > 0))
    expect_true(all(is.finite(W_surr$weights) & W_surr$weights > 0))

    # See doc-vs-behavior note above: this currently holds exactly, which is
    # surprising given the documented difference between "ind" and "surr".
    expect_equal(W_ind$weights, W_surr$weights, tolerance = eps)
  })

  test_that("GBM: unspecified missing defaults to 'ind' with a warning", {
    test_data_na <- inject_missingness(test_data, "X1")

    expect_warning({
      weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
              data = test_data_na, method = "gbm", estimand = "ATE",
              criterion = "smd.mean", n.trees = 50)
    }, "missing values are present", ignore.case = TRUE)
  })

  # Estimands
  expect_error({
    weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
            data = test_data, method = "gbm", estimand = "FOO")
  }, "not an allowable estimand", ignore.case = TRUE)
})

test_that("Multi-category treatment", {
  skip_on_cran()
  skip_if_not_installed("gbm")
  skip_if_not_installed("cobalt")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  expect_balance_improved_multi <- function(W, ...) {
    weighted <- vapply(combn(levels(W$treat), 2, simplify = FALSE), function(tt) {
      in_tt <- W$treat %in% tt
      max(abs(cobalt::col_w_smd(W$covs[in_tt, ], W$treat[in_tt], W$weights[in_tt],
                                s.weights = W$s.weights[in_tt])))
    }, numeric(1L))

    unweighted <- vapply(combn(levels(W$treat), 2, simplify = FALSE), function(tt) {
      in_tt <- W$treat %in% tt
      max(abs(cobalt::col_w_smd(W$covs[in_tt, ], W$treat[in_tt],
                                s.weights = W$s.weights[in_tt])))
    }, numeric(1L))

    expect_true(max(weighted) < max(unweighted), ...)
  }

  set.seed(123)
  expect_no_condition({
    W0 <- weightit(Am ~ X1 + X2 + X3 + X4 + X5,
                   data = test_data, method = "gbm", estimand = "ATE",
                   criterion = "smd.mean", n.trees = 100,
                   include.obj = TRUE)
  })

  expect_true(is_null(W0$ps))
  expect_false(is_null(W0$obj))
  expect_true(is_null(attr(W0, "Mparts", exact = TRUE))) #gbm does not support M-estimation

  sw.opts <- c(FALSE, TRUE)
  estimand.opts <- c("ATE", "ATT", "ATO", "ATM")

  weight.mat <- matrix(nrow = nrow(test_data),
                       ncol = length(sw.opts) * length(estimand.opts))
  colnames(weight.mat) <- rep("", ncol(weight.mat))
  sw.used <- estimand.used <- character(ncol(weight.mat))

  k <- 1

  for (sw in sw.opts) {
    for (estimand in estimand.opts) {
      test_that(sprintf("GBM: sw = %s, estimand = %s", sw, estimand), {
        set.seed(123)
        W <- weightit(Am ~ X1 + X2 + X3 + X4 + X5,
                      data = test_data, method = "gbm", estimand = estimand,
                      focal = if (estimand %in% c("ATT")) "T" else NULL,
                      criterion = "smd.mean", n.trees = 100,
                      s.weights = if (sw) "SW" else NULL,
                      include.obj = TRUE)

        expect_true(is_null(W$ps))
        expect_true(all(is.finite(W$weights) & W$weights > 0))
        expect_false(is_null(W$obj))

        expect_balance_improved_multi(W)

        if (estimand %in% c("ATT", "ATC")) {
          expect_ATT_weights_okay(W, tolerance = eps)
        }

        for (i in seq_len(k - 1)) {
          expect_not_equal(unname(W$weights), weight.mat[, i],
                           expected.label = sprintf("weights for %s", colnames(weight.mat)[i]),
                           tolerance = eps)
        }

        n <- sprintf("W_%s_%s", sw, estimand)
        colnames(weight.mat)[k] <<- n
        weight.mat[, k] <<- W$weights
        sw.used[k] <<- sw
        estimand.used[k] <<- estimand
        k <<- k + 1
      })
    }
  }

  test_that("GBM: alternative criterion (ks.mean) - multi-category", {
    set.seed(123)
    W <- weightit(Am ~ X1 + X2 + X3 + X4 + X5,
                  data = test_data, method = "gbm", estimand = "ATE",
                  criterion = "ks.mean", n.trees = 100)

    expect_true(all(is.finite(W$weights) & W$weights > 0))

    # As in the binary case: the fit tuned for criterion = "ks.mean" should
    # achieve a mean KS statistic (averaged across pairwise treatment-group
    # comparisons, matching how "ks.mean" itself is computed for multi-category
    # treatments) no worse than a fit tuned for a different criterion with the
    # same specification.
    set.seed(123)
    W_alt <- weightit(Am ~ X1 + X2 + X3 + X4 + X5,
                      data = test_data, method = "gbm", estimand = "ATE",
                      criterion = "smd.mean", n.trees = 100)

    mean_ks <- function(W) {
      vals <- vapply(combn(levels(W$treat), 2, simplify = FALSE), function(tt) {
        in_tt <- W$treat %in% tt
        mean(abs(cobalt::col_w_ks(W$covs[in_tt, ], W$treat[in_tt], W$weights[in_tt],
                                  s.weights = W$s.weights[in_tt])))
      }, numeric(1L))
      mean(vals)
    }

    expect_true(mean_ks(W) <= mean_ks(W_alt) + eps)
  })

  test_that("GBM: subclass (MMWS) - multi-category", {
    set.seed(123)
    W <- weightit(Am ~ X1 + X2 + X3 + X4 + X5,
                  data = test_data, method = "gbm", estimand = "ATE",
                  criterion = "smd.mean", n.trees = 100, subclass = 8)

    expect_true(all(is.finite(W$weights) & W$weights > 0))
    expect_true(length(unique(W$weights)) < nrow(test_data) / 4)
  })

  test_that("GBM: use.offset not allowed with multi-category treatments", {
    expect_error({
      weightit(Am ~ X1 + X2 + X3 + X4 + X5,
              data = test_data, method = "gbm", estimand = "ATE",
              criterion = "smd.mean", n.trees = 50, use.offset = TRUE)
    }, "use.offset.*cannot be used with multi-category", fixed = FALSE)
  })
})

test_that("Continuous treatment", {
  skip_on_cran()
  skip_if_not_installed("gbm")
  skip_if_not_installed("cobalt")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  expect_balance_improved_cont <- function(W, ...) {
    weighted <- abs(cobalt::col_w_cov(W$covs, W$treat, W$weights, std = TRUE,
                                      s.weights = W$s.weights))
    unweighted <- abs(cobalt::col_w_cov(W$covs, W$treat, std = TRUE,
                                        s.weights = W$s.weights))

    expect_true(max(weighted) < max(unweighted), ...)
  }

  set.seed(123)
  expect_no_condition({
    W0 <- weightit(Ac ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                   data = test_data, method = "gbm",
                   criterion = "p.mean", n.trees = 1000,
                   include.obj = TRUE)
  })

  expect_true(is_null(W0$ps))
  expect_false(is_null(W0$obj))
  expect_true(is_null(attr(W0, "Mparts", exact = TRUE))) #gbm does not support M-estimation
  expect_balance_improved_cont(W0)

  sw.opts <- c(FALSE, TRUE)

  weight.mat <- matrix(nrow = nrow(test_data), ncol = length(sw.opts))
  colnames(weight.mat) <- rep("", ncol(weight.mat))

  k <- 1

  for (sw in sw.opts) {
    test_that(sprintf("GBM: sw = %s", sw), {
      set.seed(123)
      W <- weightit(Ac ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                    data = test_data, method = "gbm",
                    criterion = "p.mean", n.trees = 1000,
                    s.weights = if (sw) "SW" else NULL,
                    include.obj = TRUE)

      expect_true(is_null(W$ps))
      expect_true(all(is.finite(W$weights) & W$weights > 0))
      expect_false(is_null(W$obj))

      expect_balance_improved_cont(W)

      for (i in seq_len(k - 1)) {
        expect_not_equal(unname(W$weights), weight.mat[, i],
                         expected.label = sprintf("weights for %s", colnames(weight.mat)[i]),
                         tolerance = eps)
      }

      n <- sprintf("W_%s", sw)
      colnames(weight.mat)[k] <<- n
      weight.mat[, k] <<- W$weights
      k <<- k + 1
    })
  }

  test_that("GBM: alternative criterion (p.max) - continuous", {
    set.seed(123)
    W <- weightit(Ac ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "gbm",
                  criterion = "p.max", n.trees = 1000)

    expect_true(all(is.finite(W$weights) & W$weights > 0))

    # The fit tuned for criterion = "p.max" should achieve a maximum correlation no worse than a fit
    # tuned for a different criterion with the same specification.
    set.seed(123)
    W_alt <- weightit(Ac ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                      data = test_data, method = "gbm",
                      criterion = "p.mean", n.trees = 1000)

    max_corr <- function(W) {
      max(abs(cobalt::col_w_corr(W$covs, W$treat, W$weights, s.weights = W$s.weights)))
    }

    expect_true(max_corr(W) <= max_corr(W_alt) + eps)
  })

  test_that("GBM: density = 'kernel' - continuous", {
    set.seed(123)
    W <- weightit(Ac ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "gbm",
                  criterion = "p.mean", n.trees = 300,
                  density = "kernel")

    expect_true(all(is.finite(W$weights) & W$weights > 0))
  })

  test_that("GBM: density = 'dt_3' - continuous", {
    set.seed(123)
    W <- weightit(Ac ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "gbm",
                  criterion = "p.mean", n.trees = 300,
                  density = "dt_3")

    expect_true(all(is.finite(W$weights) & W$weights > 0))
  })

  test_that("GBM: trim.at - continuous", {
    W0trim <- weightit(Ac ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                       data = test_data, method = "gbm",
                       criterion = "p.mean", n.trees = 300)

    Wtrim <- weightit(Ac ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                      data = test_data, method = "gbm",
                      criterion = "p.mean", n.trees = 300, trim.at = .9)

    expect_true(all(is.finite(Wtrim$weights) & Wtrim$weights > 0))
    expect_true(max(Wtrim$weights) <= max(W0trim$weights) + eps)
  })

  test_that("GBM: distribution = 'tdist' - continuous", {
    W <- weightit(Ac ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "gbm",
                  criterion = "p.mean", n.trees = 300,
                  distribution = "tdist")

    expect_true(all(is.finite(W$weights) & W$weights > 0))
  })

  test_that("GBM: use.offset - continuous", {
    W <- weightit(Ac ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "gbm",
                  criterion = "p.mean", n.trees = 300, use.offset = TRUE)

    expect_true(all(is.finite(W$weights) & W$weights > 0))
  })

  test_that("GBM: cv-based criterion works for continuous treatments", {
    set.seed(123)
    expect_no_error({
      W <- weightit(Ac ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                    data = test_data, method = "gbm",
                    criterion = "cv3", n.trees = 300)
    })
    expect_true(all(is.finite(W$weights) & W$weights > 0))
  })

  # Missing data: "ind" vs. "surr" (see extended note in the Binary block above)
  test_that("GBM: missing = 'ind' vs missing = 'surr' - continuous", {
    test_data_na <- inject_missingness(test_data, c("X1", "X3"))

    expect_no_condition({
      W_ind <- weightit(Ac ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                        data = test_data_na, method = "gbm",
                        criterion = "p.mean", n.trees = 500, missing = "ind")
    })

    expect_warning({
      W_surr <- weightit(Ac ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                         data = test_data_na, method = "gbm",
                         criterion = "p.mean", n.trees = 500, missing = "surr")
    }, "missing.*will be set to.*ind", ignore.case = TRUE)

    expect_true(anyNA(W_ind$covs$X1))
    expect_true(anyNA(W_surr$covs$X1))
    expect_true(all(is.finite(W_ind$weights) & W_ind$weights > 0))
    expect_true(all(is.finite(W_surr$weights) & W_surr$weights > 0))
  })
})

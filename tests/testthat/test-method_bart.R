# Small sampler settings throughout to keep fits fast (`bartisan_control()`'s
# defaults are num_trees = 50, num_draws = 800, num_burn = 200). bartisan uses R's
# random number generator, so `set.seed()` before every fit makes it reproducible
# (see the *Reproducibility* section of `?method_bart`).

# For a continuous treatment fit with `include.obj = TRUE`: the weights times the
# model's conditional densities (the generalized propensity scores) are the
# marginal density of the treatment, which integrates to nearly 1 over the
# observed range of the treatment
expect_marginal_density_okay <- function(W) {
  gps <- predict(W$obj, type = "density", weights = rep(1, length(W$treat)))
  ord <- order(W$treat)
  a <- as.numeric(W$treat)[ord]
  f_a <- unname(W$weights * gps)[ord]
  area <- sum(diff(a) * (f_a[-1L] + f_a[-length(f_a)]) / 2)

  expect_gt(area, .9)
  expect_lt(area, 1.01)
}

test_that("Binary treatment", {
  skip_on_cran()
  skip_if_not_installed("bartisan")
  skip_if_not_installed("cobalt")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  num_trees <- 20L
  num_draws <- 100L
  num_burn <- 100L

  set.seed(123)
  expect_no_condition({
    W0 <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                   data = test_data, method = "bart", estimand = "ATE",
                   num_trees = num_trees, num_draws = num_draws, num_burn = num_burn,
                   include.obj = TRUE)
  })

  expect_true(is.numeric(W0$ps))
  expect_true(all(W0$ps > 0 & W0$ps < 1))
  expect_true(all(is.finite(W0$weights) & W0$weights >= 0))

  # The model is bartisan's default for a binary response: binomial with a logit link
  expect_s3_class(W0$obj, "bartisan_fit")
  expect_identical(W0$obj$family$family, "binomial")
  expect_identical(W0$obj$family$link, "logit")

  # BART-based weighting does not support M-estimation, unlike
  # cbps/ebal/ipt/glm -- confirm no Mparts are attached.
  expect_null(attr(W0, "Mparts", exact = TRUE))

  expect_balance_improved(W0)

  # `set.seed()` reproduces exactly
  set.seed(123)
  W0b <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "bart", estimand = "ATE",
                  num_trees = num_trees, num_draws = num_draws, num_burn = num_burn)

  expect_equal(W0$weights, W0b$weights)

  # A different seed gives different weights
  set.seed(124)
  W0c <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "bart", estimand = "ATE",
                  num_trees = num_trees, num_draws = num_draws, num_burn = num_burn)

  expect_not_equal(W0$weights, W0c$weights)

  # Chains run one after another without a future plan and stay reproducible
  set.seed(123)
  expect_no_condition({
    W.ch <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                     data = test_data, method = "bart", estimand = "ATE",
                     num_trees = num_trees, num_draws = num_draws, num_burn = num_burn,
                     chains = 2L)
  })

  set.seed(123)
  W.ch2 <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                    data = test_data, method = "bart", estimand = "ATE",
                    num_trees = num_trees, num_draws = num_draws, num_burn = num_burn,
                    chains = 2L)

  expect_equal(W.ch$weights, W.ch2$weights)
  expect_not_equal(W.ch$weights, W0$weights)

  # Estimands
  expect_error({
    weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
             data = test_data, method = "bart", estimand = "XYZ",
             num_trees = num_trees, num_draws = num_draws, num_burn = num_burn)
  }, "not an allowable estimand", ignore.case = TRUE)

  # Non-full rank
  set.seed(123)
  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 +
                    I(1 - X5) + I(X9 * 2),
                  data = test_data, method = "bart", estimand = "ATE",
                  num_trees = num_trees, num_draws = num_draws, num_burn = num_burn)
  })

  expect_true(is.numeric(W$ps))
  expect_true(all(W$ps > 0 & W$ps < 1))

  # `subclass` computes MMWS weights instead of standard IPW; restricted to
  # ATE/ATT/ATC. Spot-check via the sharp drop in the number of unique
  # weights.
  set.seed(123)
  expect_no_condition({
    W.sub <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                      data = test_data, method = "bart", estimand = "ATE",
                      num_trees = num_trees, num_draws = num_draws, num_burn = num_burn,
                      subclass = 10)
  })

  expect_true(length(unique(W.sub$weights)) < nrow(test_data) / 10)
  expect_not_equal(W.sub$weights, W0$weights)

  # `use.offset` supplies the linear predictor of a logistic regression as an offset
  set.seed(123)
  expect_no_condition({
    W.off <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                      data = test_data, method = "bart", estimand = "ATE",
                      num_trees = num_trees, num_draws = num_draws, num_burn = num_burn,
                      use.offset = TRUE, include.obj = TRUE)
  })

  expect_true(isTRUE(W.off$obj$has_offset))
  expect_true(all(W.off$ps > 0 & W.off$ps < 1))
  expect_not_equal(W.off$weights, W0$weights)

  # Missing data with missingness indicators; the default, `missing = "mia"`, is
  # tested below
  data_na <- test_data
  set.seed(4321)
  is.na(data_na$X1) <- sample(nrow(data_na), round(.05 * nrow(data_na)))

  set.seed(123)
  expect_no_condition({
    W.na <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                     data = data_na, method = "bart", estimand = "ATE",
                     missing = "ind",
                     num_trees = num_trees, num_draws = num_draws, num_burn = num_burn)
  })

  expect_true(anyNA(W.na$covs$X1))
  expect_true(is.numeric(W.na$ps))
  expect_true(all(W.na$ps > 0 & W.na$ps < 1))

  expect_error({
    weightit(A ~ X1 + X2 + X3, data = data_na, method = "bart",
             missing = "surr",
             num_trees = num_trees, num_draws = num_draws, num_burn = num_burn)
  }, "only.*allowed for.*missing", ignore.case = TRUE)
})

test_that("Binary treatment: sampling weights", {
  skip_on_cran()
  skip_if_not_installed("bartisan")
  skip_if_not_installed("cobalt")

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  num_trees <- 20L
  num_draws <- 100L
  num_burn <- 100L

  set.seed(123)
  expect_no_condition({
    W.sw <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                     data = test_data, method = "bart", estimand = "ATE",
                     s.weights = "SW",
                     num_trees = num_trees, num_draws = num_draws, num_burn = num_burn,
                     include.obj = TRUE)
  })

  expect_true(all(is.finite(W.sw$weights) & W.sw$weights >= 0))

  # The sampling weights are passed to bartisan as prior weights with a mean of 1
  expect_equal(unname(W.sw$obj$prior_weights), test_data$SW / mean(test_data$SW))

  # so multiplying them by a constant does not change the fit
  set.seed(123)
  W.sw10 <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                     data = test_data, method = "bart", estimand = "ATE",
                     s.weights = 10 * test_data$SW,
                     num_trees = num_trees, num_draws = num_draws, num_burn = num_burn)

  expect_equal(W.sw$weights, W.sw10$weights)

  # Varying sampling weights change the weights relative to none
  set.seed(123)
  W.none <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                     data = test_data, method = "bart", estimand = "ATE",
                     num_trees = num_trees, num_draws = num_draws, num_burn = num_burn)

  expect_not_equal(W.sw$weights, W.none$weights)

  # Balance assessed with the sampling weights improves
  expect_balance_improved(W.sw)
})

test_that("Binary treatment: configurations", {
  skip_on_cran()
  skip_if_not_installed("bartisan")
  skip_if_not_installed("cobalt")
  skip_if_not_installed("patrick")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  num_trees <- 20L
  num_draws <- 100L
  num_burn <- 100L

  # Weights from the configurations already run, so each new one can be checked
  # against all of them
  seen <- new.env()

  patrick::with_parameters_test_that(
    "BART: estimand = {estimand}",
    {
      set.seed(123)
      W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                    data = test_data, method = "bart", estimand = estimand,
                    num_trees = num_trees, num_draws = num_draws, num_burn = num_burn,
                    include.obj = TRUE)

      expect_true(is.numeric(W$ps))
      expect_true(all(W$ps > 0 & W$ps < 1))
      expect_false(is_null(W$obj))
      expect_true(all(is.finite(W$weights) & W$weights >= 0))

      # BART is a machine-learning PS method; it approximates but does not
      # solve exactly for balance, so we check improvement over the
      # unweighted sample rather than exact-zero SMDs.
      expect_balance_improved(W)

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

      for (other in ls(seen)) {
        expect_not_equal(unname(W$weights), seen[[other]],
                         expected.label = sprintf("weights for %s", other),
                         tolerance = eps)
      }

      seen[[sprintf("estimand = %s", estimand)]] <- unname(W$weights)
    },
    estimand = c("ATE", "ATT", "ATC", "ATO", "ATM", "ATOS")
  )
})

test_that("Multi-category treatment", {
  skip_on_cran()
  skip_if_not_installed("bartisan")
  skip_if_not_installed("cobalt")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  num_trees <- 20L
  num_draws <- 100L
  num_burn <- 100L

  set.seed(123)
  expect_no_condition({
    W0 <- weightit(Am ~ X1 + X2 + X3 + X4 + X5,
                   data = test_data, method = "bart", estimand = "ATE",
                   num_trees = num_trees, num_draws = num_draws, num_burn = num_burn,
                   include.obj = TRUE)
  })

  expect_true(is_null(W0$ps)) #ps not directly returned for multi-category bart
  expect_true(all(is.finite(W0$weights) & W0$weights > 0))

  # A single fit of bartisan's default model for a multi-category response,
  # whose probabilities sum to 1 for each unit
  expect_s3_class(W0$obj, "bartisan_fit")
  expect_identical(W0$obj$family$family, "multinomial")
  expect_identical(W0$obj$family$link, "logit")

  p <- predict(W0$obj, type = "prob")
  expect_identical(dim(p), c(nrow(test_data), 3L))
  expect_equal(unname(rowSums(p)), rep(1, nrow(test_data)), tolerance = eps)

  # Documented behavior: ATC == ATT for multi-category treatments with the same
  # focal.
  set.seed(123)
  W.att <- weightit(Am ~ X1 + X2 + X3 + X4 + X5,
                    data = test_data, method = "bart", estimand = "ATT",
                    focal = "T",
                    num_trees = num_trees, num_draws = num_draws, num_burn = num_burn)
  set.seed(123)
  W.atc <- weightit(Am ~ X1 + X2 + X3 + X4 + X5,
                    data = test_data, method = "bart", estimand = "ATC",
                    focal = "T",
                    num_trees = num_trees, num_draws = num_draws, num_burn = num_burn)

  expect_equal(W.att$weights, W.atc$weights)

  # `ATOS` is not an allowable estimand for multi-category treatments
  expect_error({
    weightit(Am ~ X1 + X2 + X3 + X4 + X5,
             data = test_data, method = "bart", estimand = "ATOS",
             num_trees = num_trees, num_draws = num_draws, num_burn = num_burn)
  }, "not an allowable estimand", ignore.case = TRUE)

  expect_error({
    weightit(Am ~ X1 + X2 + X3 + X4 + X5,
             data = test_data, method = "bart", use.offset = TRUE,
             num_trees = num_trees, num_draws = num_draws, num_burn = num_burn)
  }, "use.offset")

  # Sampling weights are passed to bartisan with a mean of 1
  set.seed(123)
  expect_no_condition({
    W.sw <- weightit(Am ~ X1 + X2 + X3 + X4 + X5,
                     data = test_data, method = "bart", estimand = "ATE",
                     s.weights = "SW",
                     num_trees = num_trees, num_draws = num_draws, num_burn = num_burn,
                     include.obj = TRUE)
  })

  expect_true(all(is.finite(W.sw$weights) & W.sw$weights > 0))
  expect_equal(unname(W.sw$obj$prior_weights), test_data$SW / mean(test_data$SW))
  expect_not_equal(W.sw$weights, W0$weights)
})

test_that("Multi-category treatment: configurations", {
  skip_on_cran()
  skip_if_not_installed("bartisan")
  skip_if_not_installed("cobalt")
  skip_if_not_installed("patrick")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  num_trees <- 20L
  num_draws <- 100L
  num_burn <- 100L

  # Weights from the configurations already run, so each new one can be checked
  # against all of them
  seen <- new.env()

  patrick::with_parameters_test_that(
    "BART: estimand = {estimand}",
    {
      set.seed(123)
      W <- weightit(Am ~ X1 + X2 + X3 + X4 + X5,
                    data = test_data, method = "bart", estimand = estimand,
                    focal = if (estimand == "ATT") "T" else NULL,
                    num_trees = num_trees, num_draws = num_draws, num_burn = num_burn,
                    include.obj = TRUE)

      expect_true(is_null(W$ps))
      expect_false(is_null(W$obj))
      expect_true(all(is.finite(W$weights) & W$weights >= 0))

      # Only check pairs involving "T": by construction (see
      # fixtures/make_test_data.R), "C1" and "C2" are randomly assigned
      # among untreated units and share no true relationship with the
      # covariates, so a C1-vs-C2 "balance improved" check would compare
      # pure noise and can fail by chance.
      for (tt in combn(levels(W$treat), 2, simplify = FALSE)) {
        if ("T" %nin% tt) next

        in_tt <- W$treat %in% tt
        W_sub <- list(covs = W$covs[in_tt, , drop = FALSE],
                     treat = factor(W$treat[in_tt]),
                     weights = W$weights[in_tt],
                     s.weights = W$s.weights[in_tt])
        expect_balance_improved(W_sub,
                                label = sprintf("SMDs for %s", paste(tt, collapse = " vs. ")))
      }

      if (estimand == "ATT") {
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

      for (other in ls(seen)) {
        expect_not_equal(unname(W$weights), seen[[other]],
                         expected.label = sprintf("weights for %s", other),
                         tolerance = eps)
      }

      seen[[sprintf("estimand = %s", estimand)]] <- unname(W$weights)
    },
    estimand = c("ATE", "ATT", "ATO", "ATM")
  )
})

test_that("Continuous treatment", {
  skip_on_cran()
  skip_if_not_installed("bartisan")
  skip_if_not_installed("cobalt")

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  num_trees <- 20L
  num_draws <- 100L
  num_burn <- 100L

  set.seed(123)
  expect_no_condition({
    W0 <- weightit(Ac ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                   data = test_data, method = "bart",
                   num_trees = num_trees, num_draws = num_draws, num_burn = num_burn,
                   include.obj = TRUE)
  })

  expect_true(all(is.finite(W0$weights) & W0$weights > 0))

  # By default, the conditional density comes from bartisan's default model for a
  # continuous response, a Dirichlet process mixture for the errors
  expect_s3_class(W0$obj, "bartisan_fit")
  expect_identical(W0$obj$family$family, "dpm")

  expect_marginal_density_okay(W0)

  # ML-based GPS is approximate; check improvement in weighted
  # treatment-covariate correlation rather than exact-zero.
  weighted <- abs(cobalt::col_w_cov(W0$covs, W0$treat, W0$weights, std = TRUE))
  unweighted <- abs(cobalt::col_w_cov(W0$covs, W0$treat, std = TRUE))
  expect_lt(max(weighted), max(unweighted))

  # Stabilized GPS weights have a mean near 1
  expect_lt(abs(mean(W0$weights) - 1), .25)

  # Non-full rank
  set.seed(123)
  expect_no_condition({
    W <- weightit(Ac ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 +
                    I(1 - X5) + I(X9 * 2),
                  data = test_data, method = "bart",
                  num_trees = num_trees, num_draws = num_draws, num_burn = num_burn)
  })

  expect_true(all(is.finite(W$weights) & W$weights > 0))

  # The density arguments of the other methods are ignored with a warning
  set.seed(123)
  expect_warning({
    W.dens <- weightit(Ac ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                       data = test_data, method = "bart", density = "kernel",
                       num_trees = num_trees, num_draws = num_draws, num_burn = num_burn)
  }, "density.*ignored")

  expect_equal(W.dens$weights, W0$weights)

  # `use.offset` supplies the linear predictor of a linear regression as an offset
  set.seed(123)
  expect_no_condition({
    W.off <- weightit(Ac ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                      data = test_data, method = "bart", use.offset = TRUE,
                      num_trees = num_trees, num_draws = num_draws, num_burn = num_burn,
                      include.obj = TRUE)
  })

  expect_true(isTRUE(W.off$obj$has_offset))
  expect_true(all(is.finite(W.off$weights) & W.off$weights > 0))
  expect_not_equal(W.off$weights, W0$weights)
})

test_that("Continuous treatment: family", {
  skip_on_cran()
  skip_if_not_installed("bartisan")

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  num_trees <- 20L
  num_draws <- 100L
  num_burn <- 100L

  set.seed(123)
  W0 <- weightit(Ac ~ X1 + X2 + X3 + X4 + X5,
                 data = test_data, method = "bart",
                 num_trees = num_trees, num_draws = num_draws, num_burn = num_burn)

  # A normal model whose mean and standard deviation both depend on the covariates
  set.seed(123)
  expect_no_condition({
    W.ls <- weightit(Ac ~ X1 + X2 + X3 + X4 + X5,
                     data = test_data, method = "bart",
                     family = bartisan::gaussian_ls(),
                     num_trees = num_trees, num_draws = num_draws, num_burn = num_burn,
                     include.obj = TRUE)
  })

  expect_identical(W.ls$obj$family$family, "gaussian_ls")
  expect_true(all(is.finite(W.ls$weights) & W.ls$weights > 0))
  expect_marginal_density_okay(W.ls)
  expect_not_equal(W.ls$weights, W0$weights)

  # A family given by name
  set.seed(123)
  expect_no_condition({
    W.g <- weightit(Ac ~ X1 + X2 + X3 + X4 + X5,
                    data = test_data, method = "bart", family = "gaussian",
                    num_trees = num_trees, num_draws = num_draws, num_burn = num_burn,
                    include.obj = TRUE)
  })

  expect_identical(W.g$obj$family$family, "gaussian")
  expect_true(all(is.finite(W.g$weights) & W.g$weights > 0))
  expect_marginal_density_okay(W.g)

  # Arguments to the family are supplied by constructing it
  set.seed(123)
  expect_no_condition({
    W.nu <- weightit(Ac ~ X1 + X2 + X3 + X4 + X5,
                     data = test_data, method = "bart", family = bartisan::dpm(nu = 5),
                     num_trees = num_trees, num_draws = num_draws, num_burn = num_burn)
  })

  expect_not_equal(W.nu$weights, W0$weights)

  expect_error({
    weightit(Ac ~ X1 + X2 + X3, data = test_data, method = "bart", family = 3)
  }, "family")

  # With an ordinal model, each value of the treatment is a category, and the
  # weights are P(A = a) / P(A = a | X), as for a discrete treatment, with
  # P(A = a) the observed proportion of units with that value
  d <- test_data
  d$Ad <- round(d$Ac)
  d$SW2 <- d$SW / mean(d$SW)

  set.seed(123)
  expect_no_condition({
    W.ord <- weightit(Ad ~ X1 + X2 + X3 + X4 + X5, data = d, method = "bart",
                      family = "ordinal", s.weights = "SW2",
                      num_trees = num_trees, num_draws = num_draws, num_burn = num_burn,
                      include.obj = TRUE)
  })

  expect_identical(W.ord$obj$family$family, "ordinal")

  p_cond <- predict(W.ord$obj, type = "density", weights = rep(1, nrow(d)))
  p_marg <- ave(d$SW2, d$Ad, FUN = sum) / sum(d$SW2)

  expect_equal(unname(W.ord$weights), unname(p_marg / p_cond))
})

test_that("Continuous treatment: sampling weights", {
  skip_on_cran()
  skip_if_not_installed("bartisan")

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  num_trees <- 20L
  num_draws <- 100L
  num_burn <- 100L

  set.seed(123)
  W0 <- weightit(Ac ~ X1 + X2 + X3 + X4 + X5,
                 data = test_data, method = "bart",
                 num_trees = num_trees, num_draws = num_draws, num_burn = num_burn)

  # dpm(), the default family, does not take prior weights, so sampling weights
  # that vary are an error, whether dpm() is the default or supplied
  expect_error({
    weightit(Ac ~ X1 + X2 + X3 + X4 + X5,
             data = test_data, method = "bart", s.weights = "SW",
             num_trees = num_trees, num_draws = num_draws, num_burn = num_burn)
  }, "does not accept sampling weights")

  expect_error({
    weightit(Ac ~ X1 + X2 + X3 + X4 + X5,
             data = test_data, method = "bart", s.weights = "SW",
             family = bartisan::dpm(nu = 5),
             num_trees = num_trees, num_draws = num_draws, num_burn = num_burn)
  }, "does not accept sampling weights")

  # Constant sampling weights change nothing
  set.seed(123)
  expect_no_condition({
    W.const <- weightit(Ac ~ X1 + X2 + X3 + X4 + X5,
                        data = test_data, method = "bart",
                        s.weights = rep(2, nrow(test_data)),
                        num_trees = num_trees, num_draws = num_draws, num_burn = num_burn)
  })

  expect_equal(W.const$weights, W0$weights)

  # A family that takes prior weights gets the sampling weights with a mean of 1
  set.seed(123)
  expect_no_condition({
    W.ls <- weightit(Ac ~ X1 + X2 + X3 + X4 + X5,
                     data = test_data, method = "bart", s.weights = "SW",
                     family = bartisan::gaussian_ls(),
                     num_trees = num_trees, num_draws = num_draws, num_burn = num_burn,
                     include.obj = TRUE)
  })

  expect_true(all(is.finite(W.ls$weights) & W.ls$weights > 0))
  expect_equal(unname(W.ls$obj$prior_weights), test_data$SW / mean(test_data$SW))

  # The generalized propensity scores are the conditional densities themselves,
  # not raised to the power of the prior weights
  expect_marginal_density_okay(W.ls)

  # Multiplying the sampling weights by a constant does not change the fit
  set.seed(123)
  W.ls10 <- weightit(Ac ~ X1 + X2 + X3 + X4 + X5,
                     data = test_data, method = "bart", s.weights = 10 * test_data$SW,
                     family = bartisan::gaussian_ls(),
                     num_trees = num_trees, num_draws = num_draws, num_burn = num_burn)

  expect_equal(W.ls$weights, W.ls10$weights)
})

test_that("Binary and multi-category treatments: link", {
  skip_on_cran()
  skip_if_not_installed("bartisan")

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  f <- A ~ X1 + X2 + X3

  set.seed(1)
  W0 <- weightit(f, data = test_data, method = "bart",
                 num_trees = 10L, num_draws = 60L, num_burn = 40L,
                 include.obj = TRUE)

  # `link` sets the link of the binomial model, and `family` is ignored
  set.seed(1)
  expect_no_condition({
    W.probit <- weightit(f, data = test_data, method = "bart", link = "probit",
                         num_trees = 10L, num_draws = 60L, num_burn = 40L,
                         include.obj = TRUE)
  })

  expect_identical(W.probit$obj$family$family, "binomial")
  expect_identical(W.probit$obj$family$link, "probit")
  expect_true(all(W.probit$ps > 0 & W.probit$ps < 1))
  expect_not_equal(W.probit$weights, W0$weights)

  set.seed(1)
  expect_no_condition({
    W.cll <- weightit(f, data = test_data, method = "bart", link = "cloglog",
                      family = gaussian(),
                      num_trees = 10L, num_draws = 60L, num_burn = 40L,
                      include.obj = TRUE)
  })

  expect_identical(W.cll$obj$family$family, "binomial")
  expect_identical(W.cll$obj$family$link, "cloglog")

  # `use.offset` fits the GLM with the same link
  set.seed(1)
  expect_no_condition({
    W.off <- weightit(f, data = test_data, method = "bart", link = "probit",
                      use.offset = TRUE,
                      num_trees = 10L, num_draws = 60L, num_burn = 40L)
  })

  expect_true(all(W.off$ps > 0 & W.off$ps < 1))

  expect_error({
    weightit(f, data = test_data, method = "bart", link = "cauchit")
  }, "link")

  # Censoring weights use the binary model, link included
  d <- test_data
  set.seed(55)
  d$C <- rbinom(nrow(d), 1, plogis(-1.5 + .5 * d$X1))

  set.seed(1)
  expect_no_condition({
    W.cens <- weightit(.cens(C) ~ X1 + X2 + X3, data = d, method = "bart",
                       link = "probit",
                       num_trees = 10L, num_draws = 60L, num_burn = 40L,
                       include.obj = TRUE)
  })

  expect_identical(W.cens$obj$family$link, "probit")

  # The multinomial model takes the logit and probit links
  set.seed(1)
  expect_no_condition({
    W.mprobit <- weightit(Am ~ X1 + X2 + X3, data = test_data, method = "bart",
                          link = "probit",
                          num_trees = 10L, num_draws = 60L, num_burn = 40L,
                          include.obj = TRUE)
  })

  expect_identical(W.mprobit$obj$family$family, "mnp")
  expect_true(all(is.finite(W.mprobit$weights) & W.mprobit$weights > 0))

  expect_error({
    weightit(Am ~ X1 + X2 + X3, data = test_data, method = "bart", link = "cloglog")
  }, "link")
})

test_that("Covariates reach bartisan as they appear in the formula", {
  skip_on_cran()
  skip_if_not_installed("bartisan")

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  d <- test_data
  d$Xf <- factor(cut(d$X2, c(-Inf, -1, 0, 1, Inf)), labels = c("a", "b", "c", "d"))
  d$Xc <- as.character(d$Xf)

  num_trees <- 20L
  num_draws <- 100L
  num_burn <- 100L

  # A factor or character variable is one column holding a factor, not a set of
  # dummies; a transformed term is kept; an interaction is represented by its
  # components, which trees can interact on their own
  set.seed(123)
  expect_no_condition({
    W <- weightit(A ~ X1 + Xf + Xc + I(X3^2) + X4:X5, data = d, method = "bart",
                  include.obj = TRUE,
                  num_trees = num_trees, num_draws = num_draws, num_burn = num_burn)
  })

  fit_data <- W$obj$call$data
  X <- fit_data[startsWith(names(fit_data), ".x")]

  expect_identical(ncol(X), 6L)
  expect_identical(sum(vapply(X, is.factor, logical(1L))), 2L)
  expect_identical(levels(X[[2L]]), c("a", "b", "c", "d"))
  expect_true(all(W$ps > 0 & W$ps < 1))

  # With `missing = "ind"`, missing values of a factor become a category of
  # their own
  d_na <- d
  set.seed(4321)
  is.na(d_na$Xf) <- sample(nrow(d_na), 40L)

  set.seed(123)
  expect_no_condition({
    W.ind <- weightit(A ~ X1 + Xf, data = d_na, method = "bart", missing = "ind",
                      include.obj = TRUE,
                      num_trees = num_trees, num_draws = num_draws, num_burn = num_burn)
  })

  fit_data <- W.ind$obj$call$data

  expect_false(anyNA(fit_data))
  expect_identical(levels(fit_data[[".x2"]]), c("a", "b", "c", "d", "<NA>"))
  expect_identical(sum(startsWith(names(fit_data), ".x")), 2L)
})

test_that("Missing covariates: missingness incorporated in attributes", {
  skip_on_cran()
  skip_if_not_installed("bartisan")

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  data_na <- test_data
  set.seed(4321)
  is.na(data_na$X1) <- sample(nrow(data_na), round(.1 * nrow(data_na)))
  is.na(data_na$X3) <- sample(nrow(data_na), round(.05 * nrow(data_na)))

  num_trees <- 20L
  num_draws <- 100L
  num_burn <- 100L

  # "mia" is the default, with the warning every method gives when `missing` is
  # not supplied
  set.seed(123)
  expect_warning({
    W.def <- weightit(A ~ X1 + X2 + X3 + X4 + X5, data = data_na, method = "bart",
                      num_trees = num_trees, num_draws = num_draws, num_burn = num_burn)
  }, "missing values are present", ignore.case = TRUE)

  set.seed(123)
  expect_no_condition({
    W.mia <- weightit(A ~ X1 + X2 + X3 + X4 + X5, data = data_na, method = "bart",
                      missing = "mia", include.obj = TRUE,
                      num_trees = num_trees, num_draws = num_draws, num_burn = num_burn)
  })

  expect_identical(W.def$missing, "mia")
  expect_identical(W.mia$missing, "mia")
  expect_equal(W.def$weights, W.mia$weights)
  expect_true(all(W.mia$ps > 0 & W.mia$ps < 1))
  expect_true(anyNA(W.mia$covs$X1))
  expect_output(print(W.mia), "missingness incorporated in attributes")

  # The trees get the missing values themselves, with no indicators added
  fit_data <- W.mia$obj$call$data
  expect_true(anyNA(fit_data))
  expect_identical(sum(startsWith(names(fit_data), ".x")), 5L)

  set.seed(123)
  W.ind <- weightit(A ~ X1 + X2 + X3 + X4 + X5, data = data_na, method = "bart",
                    missing = "ind", include.obj = TRUE,
                    num_trees = num_trees, num_draws = num_draws, num_burn = num_burn)

  expect_false(anyNA(W.ind$obj$call$data))
  expect_identical(sum(startsWith(names(W.ind$obj$call$data), ".x")), 7L)
  expect_not_equal(W.mia$weights, W.ind$weights)

  # `use.offset` fits its GLM with missingness indicators
  set.seed(123)
  expect_no_condition({
    W.off <- weightit(A ~ X1 + X2 + X3 + X4 + X5, data = data_na, method = "bart",
                      missing = "mia", use.offset = TRUE,
                      num_trees = num_trees, num_draws = num_draws, num_burn = num_burn)
  })

  expect_true(all(W.off$ps > 0 & W.off$ps < 1))

  # Multi-category and continuous treatments
  set.seed(123)
  expect_no_condition({
    W.multi <- weightit(Am ~ X1 + X2 + X3 + X4 + X5, data = data_na, method = "bart",
                        missing = "mia",
                        num_trees = num_trees, num_draws = num_draws, num_burn = num_burn)
  })

  expect_true(all(is.finite(W.multi$weights) & W.multi$weights > 0))

  set.seed(123)
  expect_no_condition({
    W.cont <- weightit(Ac ~ X1 + X2 + X3 + X4 + X5, data = data_na, method = "bart",
                       missing = "mia", use.offset = TRUE,
                       num_trees = num_trees, num_draws = num_draws, num_burn = num_burn)
  })

  expect_true(all(is.finite(W.cont$weights) & W.cont$weights > 0))

  # `weightit.fit()` uses the same default
  set.seed(123)
  W.fit <- weightit.fit(as.matrix(data_na[c("X1", "X2", "X3", "X4", "X5")]),
                        treat = data_na$A, method = "bart",
                        num_trees = num_trees, num_draws = num_draws, num_burn = num_burn)

  expect_identical(W.fit$missing, "mia")
})

test_that("Arguments from the dbarts-based method are translated", {
  skip_on_cran()
  skip_if_not_installed("bartisan")

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  f <- A ~ X1 + X2 + X3

  set.seed(1)
  W.native <- weightit(f, data = test_data, method = "bart",
                       num_trees = 10L, num_draws = 60L, num_burn = 40L,
                       include.obj = TRUE)

  # `n.threads` is dropped silently
  set.seed(1)
  expect_no_condition({
    W.old <- weightit(f, data = test_data, method = "bart",
                      n.trees = 10L, n.samples = 60L, n.burn = 40L, n.chains = 1L,
                      n.threads = 2L, include.obj = TRUE)
  })

  expect_equal(W.old$weights, W.native$weights)
  expect_identical(W.old$obj$control$num_trees, 10L)
  expect_identical(nrow(W.old$obj$eta[[1L]]), 60L)

  # A native name wins over its translation
  set.seed(1)
  W.both <- weightit(f, data = test_data, method = "bart",
                     n.trees = 30L, num_trees = 10L, num_draws = 60L, num_burn = 40L,
                     include.obj = TRUE)

  expect_identical(W.both$obj$control$num_trees, 10L)

  # `seed` is ignored with a warning, since `set.seed()` controls the sampler
  set.seed(1)
  expect_warning({
    W.seed <- weightit(f, data = test_data, method = "bart", seed = 5L,
                       num_trees = 10L, num_draws = 60L, num_burn = 40L)
  }, "seed.*ignored")

  expect_equal(W.seed$weights, W.native$weights)
})

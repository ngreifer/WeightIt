test_that(".find_re_bars() and .no_re_bars() parse lme4-style formulas", {
  #No bars
  expect_null(.find_re_bars(a ~ x1 + x2))

  #One bar
  b1 <- .find_re_bars(a ~ x1 + x2 + (1 | g))
  expect_length(b1, 1L)
  expect_identical(deparse1(b1[[1L]]), "1 | g")
  expect_identical(deparse1(.no_re_bars(a ~ x1 + x2 + (1 | g))),
                   deparse1(a ~ x1 + x2))

  #Random slope
  b2 <- .find_re_bars(a ~ x1 + (x1 | g))
  expect_length(b2, 1L)
  expect_setequal(all.vars(b2[[1L]]), c("x1", "g"))

  #Multiple grouping factors
  b3 <- .find_re_bars(a ~ x1 + (1 | g1) + (1 | g2))
  expect_length(b3, 2L)

  #Only random effects -> fixed part becomes intercept
  expect_identical(deparse1(.no_re_bars(a ~ (1 | g))),
                   deparse1(a ~ 1))
})

test_that("Binary treatment with random effects (lme4::glmer)", {
  skip_on_cran()
  skip_if_not_installed("lme4")
  skip_if_not_installed("cobalt")

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  #ATE
  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + (1 | cluster),
                  data = test_data, method = "glm", estimand = "ATE",
                  include.obj = TRUE)
  })

  expect_true(is.numeric(W$ps))
  expect_true(all(W$ps > 0 & W$ps < 1))
  expect_true(all(is.finite(W$weights) & W$weights > 0))
  expect_s4_class(W$obj, "glmerMod")
  expect_balance_improved(W)

  #No M-estimation parts for mixed models
  expect_null(attr(W, "Mparts"))
  expect_null(attr(W, "Mparts.list"))

  #Weights differ from the fixed-effects-only fit
  W_fe <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                   data = test_data, method = "glm", estimand = "ATE")
  expect_not_equal(W$weights, W_fe$weights)

  #ATT
  expect_no_condition({
    W2 <- weightit(A ~ X1 + X2 + X3 + X7 + (1 | cluster),
                   data = test_data, method = "glm", estimand = "ATT")
  })
  expect_ATT_weights_okay(W2)

  #Link other than logit
  expect_no_condition({
    Wp <- weightit(A ~ X1 + X2 + X7 + (1 | cluster),
                   data = test_data, method = "glm", estimand = "ATE",
                   link = "probit")
  })
  expect_true(all(Wp$ps > 0 & Wp$ps < 1))

  #Unsupported link errors
  expect_error(weightit(A ~ X1 + (1 | cluster), data = test_data,
                        method = "glm", link = "identity"))
})

test_that("Continuous treatment with random effects", {
  skip_on_cran()
  skip_if_not_installed("lme4")
  skip_if_not_installed("cobalt")

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  expect_no_condition({
    W <- weightit(Ac ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + (1 | cluster),
                  data = test_data, method = "glm", include.obj = TRUE)
  })

  expect_true(all(is.finite(W$weights) & W$weights > 0))
  expect_s4_class(W$obj, "lmerMod")

  #Balance (treatment-covariate correlations) improves for continuous treatment
  weighted <- abs(cobalt::col_w_corr(W$covs, W$treat, W$weights,
                                     s.weights = W$s.weights))
  unweighted <- abs(cobalt::col_w_corr(W$covs, W$treat,
                                       s.weights = W$s.weights))
  expect_true(max(weighted) < max(unweighted))

  expect_null(attr(W, "Mparts"))

  #Unsupported link errors
  expect_error(weightit(Ac ~ X1 + (1 | cluster), data = test_data,
                        method = "glm", link = "sqrt"))
})

test_that("Multi-category treatment with random effects (mclogit::mblogit)", {
  skip_on_cran()
  skip_if_not_installed("mclogit")
  skip_if_not_installed("cobalt")

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  #mclogit's IWLS can emit benign convergence warnings on this data
  expect_no_error(suppressWarnings({
    W <- weightit(Am ~ X1 + X2 + X3 + X7 + (1 | cluster),
                  data = test_data, method = "glm", estimand = "ATE",
                  include.obj = TRUE)
  }))

  expect_true(all(is.finite(W$weights) & W$weights > 0))
  expect_s3_class(W$obj, "mmblogit")
  expect_null(attr(W, "Mparts"))

  #Legacy `random` argument (mclogit-style formula) still works
  expect_no_error(suppressWarnings({
    W2 <- weightit(Am ~ X1 + X2 + X7, data = test_data, method = "glm",
                   estimand = "ATE", multi.method = "mclogit",
                   random = ~ 1 | cluster)
  }))
  expect_true(all(is.finite(W2$weights)))
})

test_that("Random effects standard-error fallback and guardrails", {
  skip_on_cran()
  skip_if_not_installed("lme4")

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  W <- weightit(A ~ X1 + X2 + X7 + (1 | cluster), data = test_data,
                method = "glm", estimand = "ATE")

  #SEs fall back to HC0; asymptotic (M-estimation) SEs are unavailable
  expect_no_error({
    m <- glm_weightit(Y_B ~ A, data = test_data, weightit = W,
                      family = binomial)
  })
  expect_error(glm_weightit(Y_B ~ A, data = test_data, weightit = W,
                            family = binomial, vcov = "asympt"))

  #Random effects only allowed with method = "glm"
  expect_error(weightit(A ~ X1 + (1 | cluster), data = test_data,
                        method = "cbps"))
})

test_that("Random effects substitute for a correlated omitted predictor", {
  skip_on_cran()
  skip_if_not_installed("lme4")
  skip_if_not_installed("cobalt")

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  #In the fixture, `cluster` is a fine quantile-binning of X1 (a confounder), so
  #including (1 | cluster) with X1 omitted should adjust for the X1 confounding.

  #Cluster-aware fit (X1 omitted, cluster in its place)
  expect_no_condition({
    W_cluster <- weightit(A ~ X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + (1 | cluster),
                          data = test_data, method = "glm", estimand = "ATE",
                          include.obj = TRUE)
  })

  #Fit ignoring clustering (X1 also omitted)
  W_none <- weightit(A ~ X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                     data = test_data, method = "glm", estimand = "ATE")

  #Random-effect variance is non-negligible (cluster carries the X1 signal)
  re.var <- as.numeric(lme4::VarCorr(W_cluster$obj)[["cluster"]])
  expect_gt(re.var, 0.1)

  #The two fits are distinguishable by the weights they produce
  expect_not_equal(W_cluster$weights, W_none$weights)

  #Balance on the omitted confounder X1 improves when cluster stands in for it
  smd_cluster <- abs(cobalt::col_w_smd(test_data$X1, test_data$A, W_cluster$weights))
  smd_none <- abs(cobalt::col_w_smd(test_data$X1, test_data$A, W_none$weights))
  expect_lt(smd_cluster, smd_none)
})

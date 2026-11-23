# Families for continuous treatments with `method = "glm"` other than the default
# Gaussian one (`family = "gaussian_ls"`, `"poisson"`, `"negbin"`, and `"Gamma"`).
# See R/glm_cont_families.R and _dev/glm-continuous-families-notes.md.

skip_on_cran()
skip_if_not_installed("patrick")

eps <- if (capabilities("long.double")) 1e-5 else 1e-3

test_data <- local({
  d <- readRDS(test_path("fixtures", "test_data.rds"))

  set.seed(20261120)
  lp <- with(d, .4 * X1 - .3 * X5 + .3 * X3)

  d$A_pois <- rpois(nrow(d), exp(1 + lp))
  d$A_nb <- rnbinom(nrow(d), size = 1.5, mu = exp(1 + lp))
  d$A_gam <- rgamma(nrow(d), shape = 2, rate = 2 / exp(.5 + lp))
  d$A_ls <- lp + exp(-.2 + .4 * d$X1) * rnorm(nrow(d))

  d
})

patrick::with_parameters_test_that(
  "family = {family}, link = {link}, s.weights = {sw}",
  {
    if (family == "negbin") skip_if_not_installed("MASS")

    W <- expect_no_condition({
      weightit(formula, data = test_data, method = "glm", family = family,
               link = link, s.weights = if (sw) "SW")
    })

    expect_true(all(is.finite(W$weights)))
    expect_true(all(W$weights > 0))

    expect_M_parts_okay(W, tolerance = eps)

    fit <- lm_weightit(Y_C ~ ., data = cbind(test_data["Y_C"], A = W$treat),
                       weightit = W)
    expect_true(all(is.finite(sqrt(diag(vcov(fit))))))
  },
  .cases = {
    specs <- list(
      list(family = "gaussian_ls", link = "identity", formula = A_ls ~ X1 + X3 + X5),
      list(family = "gaussian_ls", link = "log", formula = I(A_ls + 10) ~ X1 + X3 + X5),
      list(family = "poisson", link = "log", formula = A_pois ~ X1 + X3 + X5),
      list(family = "poisson", link = "sqrt", formula = A_pois ~ X1 + X3 + X5),
      list(family = "negbin", link = "log", formula = A_nb ~ X1 + X3 + X5),
      list(family = "Gamma", link = "log", formula = A_gam ~ X1 + X3 + X5)
    )

    do.call(patrick::cases, unlist(lapply(specs, function(s) {
      lapply(c(FALSE, TRUE), function(sw) c(s, sw = sw))
    }), recursive = FALSE))
  }
)

test_that("the Gaussian family is the default and is unchanged", {
  #`Ac` gives some extreme weights, which is not what this test is about
  suppressWarnings({
    W0 <- weightit(Ac ~ X1 + X2 + X3, data = test_data, method = "glm")
    W1 <- weightit(Ac ~ X1 + X2 + X3, data = test_data, method = "glm",
                   family = "gaussian")
    W2 <- weightit(Ac ~ X1 + X2 + X3, data = test_data, method = "glm",
                   family = gaussian())
  })

  expect_identical(W1$weights, W0$weights)
  expect_identical(W2$weights, W0$weights)

  #A family object supplies its link unless `link` is given (the square root link
  #is misspecified here and gives some extreme weights)
  suppressWarnings({
    W3 <- weightit(A_pois ~ X1 + X3, data = test_data, method = "glm",
                   family = poisson(link = "sqrt"))
    W4 <- weightit(A_pois ~ X1 + X3, data = test_data, method = "glm",
                   family = "poisson", link = "sqrt")
  })
  expect_identical(W3$weights, W4$weights)

  W5 <- weightit(A_pois ~ X1 + X3, data = test_data, method = "glm",
                 family = poisson(link = "sqrt"), link = "log")
  W6 <- weightit(A_pois ~ X1 + X3, data = test_data, method = "glm",
                 family = "poisson")
  expect_identical(W5$weights, W6$weights)

  #Case-insensitive names
  W7 <- weightit(A_gam ~ X1 + X3, data = test_data, method = "glm",
                 family = "gamma")
  W8 <- weightit(A_gam ~ X1 + X3, data = test_data, method = "glm",
                 family = Gamma(link = "log"))
  expect_identical(W7$weights, W8$weights)
})

test_that("weights use the exact marginal density for discrete families and the location-scale model", {
  #Poisson: f(a) is the sampling-weighted average of the conditional pmfs
  W <- weightit(A_pois ~ X1 + X3, data = test_data, method = "glm",
                family = "poisson", s.weights = "SW", include.obj = TRUE)
  mu <- W$obj$fitted.values
  A <- test_data$A_pois
  SW <- test_data$SW

  num <- vapply(A, function(a) sum(SW * dpois(a, mu)) / sum(SW), numeric(1L))
  expect_equal(W$weights, num / dpois(A, mu), tolerance = 1e-10, ignore_attr = TRUE)

  #Location-scale: the same, with each unit's own SD
  W <- weightit(A_ls ~ X1 + X3, data = test_data, method = "glm",
                family = "gaussian_ls", include.obj = TRUE)
  mu <- W$obj$fitted.values
  s <- W$obj$sd
  A <- test_data$A_ls

  num <- vapply(A, function(a) mean(dnorm(a, mu, s)), numeric(1L))
  expect_equal(W$weights, num / dnorm(A, mu, s), tolerance = 1e-10, ignore_attr = TRUE)

  #Gamma: interpolated on a log-scale grid, close to the exact average
  W <- weightit(A_gam ~ X1 + X3, data = test_data, method = "glm",
                family = "Gamma")
  M <- attr(W, "Mparts")
  shape <- exp(M$btreat[["log(shape)"]])
  mu <- exp(drop(M$Xtreat %*% M$btreat[-1L]))
  A <- test_data$A_gam

  num <- vapply(A, function(a) mean(dgamma(a, shape, rate = shape / mu)), numeric(1L))
  expect_equal(W$weights, num / dgamma(A, shape, rate = shape / mu), tolerance = 1e-3,
               ignore_attr = TRUE)
})

test_that("the Gamma shape and negative binomial theta are maximum likelihood estimates", {
  skip_if_not_installed("MASS")

  W <- weightit(A_gam ~ X1 + X3, data = test_data, method = "glm",
                family = "Gamma", include.obj = TRUE)
  shape <- exp(attr(W, "Mparts")$btreat[["log(shape)"]])

  expect_equal(shape, MASS::gamma.shape(W$obj, it.lim = 100, eps.max = 1e-12)$alpha,
               tolerance = 1e-6, ignore_attr = TRUE)

  W <- weightit(A_nb ~ X1 + X3, data = test_data, method = "glm",
                family = "negbin", include.obj = TRUE)
  theta <- exp(attr(W, "Mparts")$btreat[["log(theta)"]])

  expect_equal(theta, W$obj$theta, tolerance = 1e-4, ignore_attr = TRUE)
})

test_that("families compose with stabilize, weightitMSM(), and random effects", {
  skip_if_not_installed("rootSolve")

  W <- expect_no_condition({
    weightit(A_pois ~ X1 + X3 + X5, data = test_data, method = "glm",
             family = "poisson", stabilize = ~X5)
  })
  expect_M_parts_okay(W, tolerance = eps)

  test_data$A_pois2 <- rpois(nrow(test_data), exp(.5 + .2 * test_data$A_pois / 5))

  W <- expect_no_condition({
    weightitMSM(list(A_pois ~ X1 + X3, A_pois2 ~ X1 + X3 + A_pois),
                data = test_data, method = "glm", family = "poisson")
  })
  expect_M_parts_okay(W, tolerance = eps)

  skip_if_not_installed("lme4")

  test_data$clus <- factor(rep_len(1:15, nrow(test_data)))

  for (f in c("poisson", "negbin", "Gamma")) {
    A <- switch(f, poisson = "A_pois", negbin = "A_nb", Gamma = "A_gam")

    W <- expect_no_error({
      suppressWarnings(weightit(reformulate(c("X1", "X3", "(1 | clus)"), A),
                                data = test_data, method = "glm", family = f))
    })
    expect_true(all(is.finite(W$weights)))
    expect_null(attr(W, "Mparts"))
  }

  #For Gamma, the shape is the maximum likelihood estimate given the fitted
  #conditional means
  W <- suppressWarnings(weightit(A_gam ~ X1 + X3 + (1 | clus), data = test_data,
                                 method = "glm", family = "Gamma",
                                 include.obj = TRUE))
  mu <- as.numeric(predict(W$obj, type = "response"))
  shape <- .gamma_shape_ml(test_data$A_gam, mu, rep(1, nrow(test_data)))
  A <- test_data$A_gam
  num <- vapply(A, function(a) mean(dgamma(a, shape, rate = shape / mu)), numeric(1L))
  expect_equal(W$weights, num / dgamma(A, shape, rate = shape / mu),
               tolerance = 1e-3, ignore_attr = TRUE)
})

test_that("a model with no covariates gives weights of 1", {
  W <- weightit(A_pois ~ 1, data = test_data, method = "glm", family = "poisson")
  expect_equal(unname(W$weights), rep(1, nrow(test_data)))
})

test_that("families give informative errors and warnings", {
  expect_error(weightit(Ac ~ X1, data = test_data, method = "glm", family = "poisson"),
               "nonnegative integer")
  expect_error(weightit(Ac ~ X1, data = test_data, method = "glm", family = "Gamma"),
               "positive values")
  expect_error(weightit(Ac ~ X1, data = test_data, method = "glm", family = "beta"),
               "must be one of")
  expect_error(weightit(Ac ~ X1, data = test_data, method = "glm", family = binomial()),
               "must be one of")
  #Bias reduction is not available for continuous treatments
  expect_error(weightit(Ac ~ X1, data = test_data, method = "glm", link = "br.identity"),
               "must be one of")
  expect_error(weightit(A_pois ~ X1, data = test_data, method = "glm",
                        family = "poisson", link = "br.log"),
               "must be one of")
  expect_error(weightit(A_pois ~ X1, data = test_data, method = "glm",
                        family = "poisson", link = "inverse"),
               "must be one of")
  expect_warning(weightit(A_pois ~ X1, data = test_data, method = "glm",
                          family = "poisson", density = "dt_4"),
                 "ignored")

  skip_if_not_installed("lme4")
  test_data$clus <- factor(rep_len(1:15, nrow(test_data)))

  expect_error(weightit(A_ls ~ X1 + (1 | clus), data = test_data, method = "glm",
                        family = "gaussian_ls"),
               "random effects", ignore.case = TRUE)
})

test_that("another method's family is not passed to the stabilization numerator", {
  skip_if_not_installed("bartisan")

  #`family = "dpm"` is a bartisan family, not one `method = "glm"` has; the
  #numerator is fit by `method = "glm"`, which must use its own default rather
  #than reject it
  expect_no_error({
    W <- weightit(Ac ~ X1 + X3, data = test_data[1:200, ], method = "bart",
                  family = "dpm", stabilize = ~X5,
                  n.trees = 20, n.burn = 50, n.samples = 50, n.chains = 1)
  })
})

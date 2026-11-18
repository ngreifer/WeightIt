#`hatvalues()` for `multinom_weightit` and `ordinal_weightit` objects returns the
#generalized leverage tr(F^{-1} F_i), where F_i is unit i's contribution to the
#expected information F. The brute-force version builds each F_i from finite
#differences of the fitted category probabilities, F_i = w_i J_i' diag(1 / p_i) J_i,
#so it shares nothing with the closed form but the fitted model.
brute_hatvalues <- function(fit, X, offset) {
  theta <- na.rem(coef(fit))
  w <- fit$weights

  get_p <- function(B) {
    unname(fit$get_p(B, X, offset))
  }

  P <- get_p(theta)

  eps <- 1e-6

  J <- lapply(seq_along(theta), function(j) {
    d <- replace(numeric(length(theta)), j, eps)
    (get_p(theta + d) - get_p(theta - d)) / (2 * eps)
  })

  F_i <- lapply(seq_len(nrow(P)), function(i) {
    Ji <- do.call("cbind", lapply(J, function(Jj) Jj[i, ]))
    w[i] * crossprod(Ji / sqrt(P[i, ]))
  })

  Finv <- solve(Reduce("+", F_i))

  vapply(F_i, function(Fi) sum(Finv * t(Fi)), numeric(1L))
}

test_that("hatvalues() gives the generalized leverage", {
  skip_if_not_installed("patrick")

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  set.seed(1234)
  test_data <- test_data[sample.int(nrow(test_data), 100L), ]
  test_data$Y_M <- factor(test_data$Y_O, ordered = FALSE)
  test_data$off <- runif(nrow(test_data))

  test_data$w <- test_data$SW
  test_data$w[1:2] <- 0

  patrick::with_parameters_test_that(
    "{model}, link = {link}",
    {
      fit <- {
        if (model == "multinom") {
          multinom_weightit(Y_M ~ A + X1 + X5 + offset(off), data = test_data,
                            weights = w)
        }
        else {
          ordinal_weightit(Y_O ~ A + X1 + X5 + offset(off), data = test_data,
                           weights = w, link = link)
        }
      }

      h <- hatvalues(fit)

      expect_length(h, nrow(sandwich::estfun(fit)))
      expect_identical(names(h), rownames(test_data))
      expect_identical(unname(h[1:2]), c(0, 0))
      expect_true(all(h >= 0))

      #The leverages sum to the number of coefficients
      expect_equal(sum(h), length(coef(fit)), tolerance = 1e-10)

      expect_equal(unname(h),
                   brute_hatvalues(fit, X = model.matrix(fit), offset = test_data$off),
                   tolerance = 1e-6)
    },
    .cases = data.frame(model = c("multinom", "ordinal", "ordinal", "ordinal"),
                        link = c("logit", "logit", "probit", "cloglog"))
  )
})

#With two categories, both models are binomial GLMs, whose hat values `glm()` gives.
#The fits differ at the level of their convergence criteria.
test_that("hatvalues() agrees with glm() for a binary outcome", {
  skip_if_not_installed("patrick")

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))
  test_data$Y_O2 <- factor(test_data$Y_B, ordered = TRUE)

  patrick::with_parameters_test_that(
    "{model}, link = {link}",
    {
      fit <- {
        if (model == "multinom") {
          multinom_weightit(Y_O2 ~ A + X1 + X2 + X5, data = test_data,
                            weights = SW)
        }
        else {
          ordinal_weightit(Y_O2 ~ A + X1 + X2 + X5, data = test_data,
                           weights = SW, link = link)
        }
      }

      fit_g <- glm(Y_O2 ~ A + X1 + X2 + X5, data = test_data, weights = SW,
                   family = quasibinomial(link),
                   control = list(epsilon = 1e-12, maxit = 100L))

      expect_equal(hatvalues(fit), hatvalues(fit_g), tolerance = 1e-5)
    },
    .cases = data.frame(model = c("multinom", "ordinal", "ordinal"),
                        link = c("logit", "logit", "probit"))
  )
})

test_that("hatvalues() ignores aliased coefficients", {
  test_data <- readRDS(test_path("fixtures", "test_data.rds"))
  test_data$Y_M <- factor(test_data$Y_O, ordered = FALSE)
  test_data$A2 <- 2 * test_data$A

  fit <- ordinal_weightit(Y_O ~ A + A2 + X1, data = test_data, weights = SW)
  fit0 <- ordinal_weightit(Y_O ~ A + X1, data = test_data, weights = SW)

  expect_equal(hatvalues(fit), hatvalues(fit0), tolerance = 1e-8)

  fit <- multinom_weightit(Y_M ~ A + A2 + X1, data = test_data, weights = SW)
  fit0 <- multinom_weightit(Y_M ~ A + X1, data = test_data, weights = SW)

  expect_equal(hatvalues(fit), hatvalues(fit0), tolerance = 1e-8)

  #With no covariates, only the thresholds are estimated
  fit <- ordinal_weightit(Y_O ~ 1, data = test_data, weights = SW)

  expect_equal(sum(hatvalues(fit)), nlevels(test_data$Y_O) - 1L,
               tolerance = 1e-10)
})

# Bias-reduced fits scale the weights to have a mean of 1 among the units with a
# nonzero weight before fitting (Mukhopadhyay, 2020). The score grows with the
# weights but the bias-reducing adjustment does not, so without this, multiplying
# the weights by a constant would change the estimates. With it, bias-reduced
# estimates and their standard errors are invariant to the scale of the weights,
# as those of unpenalized fits are.

test_that("bias-reduced outcome models are invariant to the scale of the weights", {
  skip_on_cran()
  skip_if_not_installed("patrick")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  # A smaller sample, so that the adjustment moves the estimates noticeably
  test_data <- readRDS(test_path("fixtures", "test_data.rds"))[1:200, ]
  test_data$Y_M <- factor(test_data$Y_O, ordered = FALSE)

  cutoff <- quantile(test_data$Y_S, .8)
  test_data$event <- as.numeric(test_data$Y_S < cutoff)
  test_data$time <- pmin(test_data$Y_S, cutoff)

  W <- weightit(A ~ X1 + X2 + X3, data = test_data, method = "glm", estimand = "ATE")

  W10 <- W
  W10$weights <- 10 * W$weights

  patrick::with_parameters_test_that(
    "{fitter}",
    {
      if (fitter == "glm_weightit") {
        skip_if_not_installed("brglm2")
      }

      if (fitter == "coxph_weightit") {
        skip_if_not_installed("survival")
      }

      fit <- do.call(fitter, c(fit_args, list(data = quote(test_data), weightit = quote(W),
                                              br = TRUE)))

      fit10 <- do.call(fitter, c(fit_args, list(data = quote(test_data), weightit = quote(W10),
                                                br = TRUE)))

      # The adjustment changes the estimates, so the checks below can fail
      fit_ml <- do.call(fitter, c(fit_args, list(data = quote(test_data), weightit = quote(W))))
      expect_not_equal(coef(fit), coef(fit_ml), tolerance = eps)

      expect_equal(coef(fit10), coef(fit), tolerance = eps)
      expect_equal(vcov(fit10), vcov(fit), tolerance = eps)
      expect_equal(vcov(fit10, vcov = "HC0"), vcov(fit, vcov = "HC0"), tolerance = eps)
    },
    .cases = patrick::cases(
      list(fitter = "glm_weightit",
           fit_args = list(quote(Y_B ~ A + X1), family = quote(binomial))),
      list(fitter = "multinom_weightit",
           fit_args = list(quote(Y_M ~ A + X1))),
      list(fitter = "ordinal_weightit",
           fit_args = list(quote(Y_O ~ A + X1))),
      list(fitter = "coxph_weightit",
           fit_args = list(quote(survival::Surv(time, event) ~ A + X1)))
    )
  )
})

test_that("bias-reduced propensity score models are invariant to the scale of the sampling weights", {
  skip_on_cran()
  skip_if_not_installed("patrick")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))[1:200, ]
  test_data$Ao <- factor(test_data$Am, ordered = TRUE)
  test_data$SW10 <- 10 * test_data$SW

  patrick::with_parameters_test_that(
    "{case}",
    {
      skip_if_not_installed(pkg)

      W <- do.call("weightit", c(wargs, list(data = quote(test_data), method = "glm",
                                             s.weights = "SW")))

      W10 <- do.call("weightit", c(wargs, list(data = quote(test_data), method = "glm",
                                               s.weights = "SW10")))

      expect_equal(W10$weights, W$weights, tolerance = eps)
    },
    .cases = patrick::cases(
      list(case = "binary, br.logit",
           wargs = list(quote(A ~ X1 + X2 + X3), link = "br.logit"),
           pkg = "brglm2"),
      list(case = "binary, flic",
           wargs = list(quote(A ~ X1 + X2 + X3), link = "flic"),
           pkg = "logistf"),
      list(case = "binary, flac",
           wargs = list(quote(A ~ X1 + X2 + X3), link = "flac"),
           pkg = "logistf"),
      list(case = "multinomial, br.logit",
           wargs = list(quote(Am ~ X1 + X2 + X3), link = "br.logit"),
           pkg = "stats"),
      list(case = "ordinal, br.logit",
           wargs = list(quote(Ao ~ X1 + X2 + X3), link = "br.logit"),
           pkg = "stats")
    )
  )
})

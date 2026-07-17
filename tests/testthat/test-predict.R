test_that("predict.glm_weightit", {
  skip_on_cran()

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                data = test_data, method = "glm", estimand = "ATE")

  fit <- glm_weightit(Y_B ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                      data = test_data, family = binomial, weightit = W,
                      vcov = "none")

  #Default type is "response", which differs from stats::predict.glm()
  p_response <- predict(fit)
  expect_length(p_response, nrow(test_data))
  expect_true(all(p_response >= 0 & p_response <= 1))
  expect_equal(unname(p_response), unname(fit$fitted.values), tolerance = eps)

  #"probs" is an alias for "response"
  expect_equal(predict(fit, type = "probs"), p_response, tolerance = eps)

  #"link"/"lp" is the linear predictor scale
  p_link <- predict(fit, type = "link")
  expect_equal(unname(p_link), unname(fit$linear.predictors), tolerance = eps)
  expect_equal(predict(fit, type = "lp"), p_link, tolerance = eps)

  #response = linkinv(link)
  expect_equal(unname(p_response), unname(plogis(p_link)), tolerance = eps)

  #newdata with the same data should reproduce the fitted values
  p_newdata <- predict(fit, newdata = test_data)
  expect_equal(unname(p_newdata), unname(p_response), tolerance = eps)

  #newdata with a subset
  sub <- test_data[1:50, ]
  p_sub <- predict(fit, newdata = sub)
  expect_length(p_sub, nrow(sub))
  expect_true(all(p_sub >= 0 & p_sub <= 1))
  expect_equal(unname(p_sub), unname(p_response[1:50]), tolerance = eps)

  #G-computation-style newdata (setting treatment to a fixed value for all units)
  p0 <- predict(fit, newdata = transform(test_data, A = 0))
  p1 <- predict(fit, newdata = transform(test_data, A = 1))
  expect_length(p0, nrow(test_data))
  expect_length(p1, nrow(test_data))
  expect_true(all(p0 >= 0 & p0 <= 1))
  expect_true(all(p1 >= 0 & p1 <= 1))

  #Only "response"/"probs" and "link"/"lp" are allowed for glm_weightit
  expect_error(predict(fit, type = "class"))
  expect_error(predict(fit, type = "mean"))
})

test_that("predict.multinom_weightit", {
  skip_on_cran()

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))
  test_data$Y_M <- factor(test_data$Y_O, ordered = FALSE)

  W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                data = test_data, method = "glm", estimand = "ATE")

  fit <- multinom_weightit(Y_M ~ A * (X1 + X2 + X3 + X4 + X5),
                           data = test_data, weightit = W, vcov = "none")

  outcome_levels <- levels(test_data$Y_M)

  #Default type is "response": a matrix with one column per outcome level
  p_response <- predict(fit)
  expect_true(is.matrix(p_response))
  expect_equal(dim(p_response), c(nrow(test_data), nlevels(test_data$Y_M)))
  expect_equal(colnames(p_response), outcome_levels)
  expect_true(all(p_response >= 0 & p_response <= 1))
  expect_equal(unname(rowSums(p_response)), rep(1, nrow(test_data)), tolerance = eps)
  expect_equal(p_response, fit$fitted.values, tolerance = eps)

  #"probs" is an alias for "response"
  expect_equal(predict(fit, type = "probs"), p_response, tolerance = eps)

  #"level" restricts the response matrix to a single column
  p_level <- predict(fit, level = outcome_levels[1L])
  expect_equal(unname(p_level), unname(p_response[, 1L]), tolerance = eps)

  p_level_idx <- predict(fit, level = 1L)
  expect_equal(p_level_idx, p_level, tolerance = eps)

  #"class": modal predicted category, an unordered factor with the outcome's levels
  p_class <- predict(fit, type = "class")
  expect_s3_class(p_class, "factor")
  expect_false(is.ordered(p_class))
  expect_equal(levels(p_class), outcome_levels)
  expect_length(p_class, nrow(test_data))
  expect_equal(as.integer(p_class), max.col(p_response, ties.method = "first"))

  #"mean": expected value of the outcome; here the levels are numeric-like
  #("1"-"4") so `values` can be inferred automatically
  p_mean <- predict(fit, type = "mean")
  expect_true(is.numeric(p_mean))
  expect_length(p_mean, nrow(test_data))
  auto_values <- setNames(as.numeric(outcome_levels), outcome_levels)
  expect_equal(unname(p_mean),
               unname(drop(p_response %*% auto_values[outcome_levels])),
               tolerance = eps)

  #"mean" with explicit `values`
  values <- setNames(c(10, 20, 30, 40), outcome_levels)
  p_mean2 <- predict(fit, type = "mean", values = values)
  expect_equal(unname(p_mean2),
               unname(drop(p_response %*% values[outcome_levels])),
               tolerance = eps)

  #"link"/"lp" is documented as not usable with multinomial models
  expect_error(predict(fit, type = "link"))
  expect_error(predict(fit, type = "lp"))

  #"stdlv" is documented as ordinal-only
  expect_error(predict(fit, type = "stdlv"))

  #newdata with the same data should reproduce the fitted values
  p_newdata <- predict(fit, newdata = test_data)
  expect_equal(unname(p_newdata), unname(p_response), tolerance = eps)

  #newdata with a subset
  sub <- test_data[1:100, ]
  p_sub <- predict(fit, newdata = sub)
  expect_equal(dim(p_sub), c(nrow(sub), nlevels(test_data$Y_M)))
  expect_equal(unname(p_sub), unname(p_response[1:100, ]), tolerance = eps)

  p_sub_class <- predict(fit, newdata = sub, type = "class")
  expect_equal(p_sub_class, p_class[1:100], ignore_attr = TRUE)

  p_sub_mean <- predict(fit, newdata = sub, type = "mean")
  expect_equal(unname(p_sub_mean), unname(p_mean[1:100]), tolerance = eps)
})

test_that("predict.ordinal_weightit", {
  skip_on_cran()

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                data = test_data, method = "glm", estimand = "ATE")

  fit <- ordinal_weightit(Y_O ~ A * (X1 + X2 + X3 + X4 + X5),
                          data = test_data, weightit = W, vcov = "none")

  outcome_levels <- levels(test_data$Y_O)

  #Default type is "response": a matrix with one column per outcome level
  p_response <- predict(fit)
  expect_true(is.matrix(p_response))
  expect_equal(dim(p_response), c(nrow(test_data), nlevels(test_data$Y_O)))
  expect_true(all(p_response >= 0 & p_response <= 1))
  expect_equal(unname(rowSums(p_response)), rep(1, nrow(test_data)), tolerance = eps)
  expect_equal(p_response, fit$fitted.values, tolerance = eps)

  #"link"/"lp": linear predictor without thresholds -- allowed for ordinal models
  p_link <- predict(fit, type = "link")
  expect_equal(unname(p_link), unname(fit$linear.predictors), tolerance = eps)
  expect_equal(predict(fit, type = "lp"), p_link, tolerance = eps)
  expect_length(p_link, nrow(test_data))

  #object$linear.predictors previously used a mismatched scale (the
  #standardized design matrix multiplied by already-destandardized
  #coefficients), so predict(fit, type = "link") with no `newdata` disagreed
  #with predict(fit, newdata = <the same data>, type = "link"). Now fixed in
  #R/ordinal_weightit.R -- both should agree, and both should match a
  #from-scratch computation on the original (unstandardized) design matrix.
  tt <- terms(fit)
  Terms <- delete.response(tt)
  mf <- model.frame(Terms, test_data, na.action = na.pass, xlev = fit$xlevels)
  Xmat <- model.matrix(Terms, mf, contrasts.arg = fit$contrasts)
  Xmat <- Xmat[, colnames(Xmat) != "(Intercept)", drop = FALSE]
  beta <- fit$coefficients[seq_len(ncol(Xmat))]

  p_link_correct <- predict(fit, newdata = test_data, type = "link")
  expect_equal(unname(p_link_correct), unname(drop(Xmat %*% beta)), tolerance = eps)
  expect_equal(unname(p_link), unname(p_link_correct), tolerance = eps)

  #"class": modal predicted category, an *ordered* factor for ordinal models
  p_class <- predict(fit, type = "class")
  expect_s3_class(p_class, "factor")
  expect_true(is.ordered(p_class))
  expect_equal(levels(p_class), outcome_levels)
  expect_equal(as.integer(p_class), max.col(p_response, ties.method = "first"))

  #"mean": expected value of the outcome
  p_mean <- predict(fit, type = "mean")
  expect_true(is.numeric(p_mean))
  auto_values <- setNames(as.numeric(outcome_levels), outcome_levels)
  expect_equal(unname(p_mean),
               unname(drop(p_response %*% auto_values[outcome_levels])),
               tolerance = eps)

  #"stdlv": standardized latent variable; defined for the default "logit" link
  p_stdlv <- predict(fit, type = "stdlv")
  expect_true(is.numeric(p_stdlv))
  expect_length(p_stdlv, nrow(test_data))

  #"stdlv" derives from object$linear.predictors -- now fixed, the no-newdata
  #and newdata-based (from-scratch) computations should agree.
  p_stdlv_correct <- predict(fit, newdata = test_data, type = "stdlv")
  expect_equal(unname(p_stdlv), unname(p_stdlv_correct), tolerance = eps)

  #"stdlv" also works with the "probit" link
  fit_probit <- ordinal_weightit(Y_O ~ A * (X1 + X2 + X3 + X4 + X5),
                                 data = test_data, weightit = W, vcov = "none",
                                 link = "probit")
  expect_no_condition(predict(fit_probit, type = "stdlv"))

  #"stdlv" errors for links outside {probit, logit, cloglog, loglog} (per docs
  #and the explicit check in predict.multinom_weightit()/predict.ordinal_weightit())
  fit_cauchit <- ordinal_weightit(Y_O ~ A * (X1 + X2 + X3 + X4 + X5),
                                  data = test_data, weightit = W, vcov = "none",
                                  link = "cauchit")
  expect_error(predict(fit_cauchit, type = "stdlv"), "cannot be used")

  #newdata with the same data should reproduce the fitted values
  p_newdata <- predict(fit, newdata = test_data)
  expect_equal(unname(p_newdata), unname(p_response), tolerance = eps)

  #newdata with a subset, across all types
  sub <- test_data[1:100, ]

  p_sub <- predict(fit, newdata = sub)
  expect_equal(dim(p_sub), c(nrow(sub), nlevels(test_data$Y_O)))
  expect_equal(unname(p_sub), unname(p_response[1:100, ]), tolerance = eps)

  #Compare against the newdata-based full-data predictions computed above --
  #the newdata code path should be internally self-consistent regardless of subset.
  p_sub_link <- predict(fit, newdata = sub, type = "link")
  expect_equal(unname(p_sub_link), unname(p_link_correct[1:100]), tolerance = eps)

  p_sub_stdlv <- predict(fit, newdata = sub, type = "stdlv")
  expect_equal(unname(p_sub_stdlv), unname(p_stdlv_correct[1:100]), tolerance = eps)

  p_sub_mean <- predict(fit, newdata = sub, type = "mean")
  expect_equal(unname(p_sub_mean), unname(p_mean[1:100]), tolerance = eps)
})

test_that("predict.coxph_weightit", {
  skip_on_cran()
  skip_if_not_installed("survival")

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                data = test_data, method = "glm", estimand = "ATE")

  fit <- coxph_weightit(survival::Surv(Y_S) ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                        data = test_data, weightit = W, vcov = "none")

  #WeightIt does NOT define its own predict.coxph_weightit method (confirmed:
  #no S3method(predict, coxph_weightit) is registered in the package NAMESPACE).
  #coxph_weightit objects inherit class "coxph" in addition to "coxph_weightit"
  #(confirmed below), so predict() dispatches via ordinary S3 method resolution
  #to survival's (unexported but registered) predict.coxph() method -- this is
  #genuine S3 dispatch on the inherited class, not a WeightIt-specific method.
  expect_true(inherits(fit, "coxph"))
  expect_null(getS3method("predict", "coxph_weightit", optional = TRUE))
  expect_false(is.null(getS3method("predict", "coxph", optional = TRUE)))

  #Predictions from predict() should be identical to calling survival's method directly
  p_default <- predict(fit)
  p_default_direct <- survival:::predict.coxph(fit)
  expect_equal(p_default, p_default_direct)

  for (tt in c("lp", "risk", "expected", "terms")) {
    p_tt <- predict(fit, type = tt)
    p_tt_direct <- survival:::predict.coxph(fit, type = tt)
    expect_equal(p_tt, p_tt_direct)
  }

  #newdata
  sub <- test_data[1:50, ]
  p_sub <- predict(fit, newdata = sub, type = "lp")
  p_sub_direct <- survival:::predict.coxph(fit, newdata = sub, type = "lp")
  expect_equal(p_sub, p_sub_direct)
  expect_length(p_sub, nrow(sub))
})

test_that("Binary treatment", {
  skip_on_cran()
  skip_if_not_installed("rootSolve")
  skip_if_not_installed("cobalt")
  skip_if_not_installed("brglm2")
  skip_if_not_installed("logistf")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  expect_no_condition({
    W0 <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                   data = test_data, method = "glm", estimand = "ATE",
                   include.obj = TRUE)
  })

  expect_M_parts_okay(W0, tolerance = eps)

  expect_true(is.numeric(W0$ps))

  # quick

  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "glm", estimand = "ATE",
                  quick = TRUE, include.obj = TRUE)
  })

  expect_M_parts_okay(W, tolerance = eps)

  expect_equal(W$weights, W0$weights, tolerance = eps)

  expect_false(is_null(W$obj))
  expect_false(is_null(W0$obj))
  expect_not_equal(W$obj, W0$obj)

  # Estimands

  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "glm", estimand = "ATT")
  })

  expect_M_parts_okay(W, tolerance = eps)

  expect_equal(unname(W$weights[W$treat == 1]),
               rep(1, sum(W$treat == 1)),
               tolerance = eps)

  expect_ATT_weights_okay(W, tolerance = eps)

  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "glm", estimand = "ATC")
  })

  expect_M_parts_okay(W, tolerance = eps)

  expect_equal(unname(W$weights[W$treat == 0]),
               rep(1, sum(W$treat == 0)),
               tolerance = eps)

  expect_ATT_weights_okay(W, tolerance = eps)

  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "glm", estimand = "ATO")
  })

  expect_M_parts_okay(W, tolerance = eps)

  expect_equal(unname(cobalt::col_w_smd(W$covs, W$treat, W$weights)),
               rep(0, 12),
               tolerance = eps)

  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "glm", estimand = "ATM")
  })

  expect_M_parts_okay(W, tolerance = eps)

  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "glm", estimand = "ATOS")
  })

  expect_M_parts_okay(W, tolerance = eps)

  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "glm", estimand = "ATE",
                  link = "probit")
  })

  expect_M_parts_okay(W, tolerance = eps)

  # brglm2
  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "glm", estimand = "ATE",
                  link = "br.logit", epsilon = 1e-10)
  })

  expect_M_parts_okay(W, tolerance = eps)

  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "glm", estimand = "ATE",
                  link = "br.probit", type = "AS_median", epsilon = 1e-10)
  })

  expect_M_parts_okay(W, tolerance = eps)

  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "glm", estimand = "ATE",
                  link = "br.logit", type = "correction", epsilon = 1e-10)
  })

  expect_null(attr(W, "Mparts"))

  # logistf
  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "glm", estimand = "ATE",
                  link = "flic")
  })

  expect_null(attr(W, "Mparts"))

  # s.weights

  expect_no_condition({
    WS<- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "glm", estimand = "ATE",
                  s.weights = "SW")
  })

  expect_equal(test_data$SW, WS$s.weights)

  expect_M_parts_okay(WS, tolerance = eps)

  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X5 + X6,
                  data = test_data, method = "glm", estimand = "ATE",
                  s.weights = "SW", link = "log")
  })

  expect_M_parts_okay(W, tolerance = eps)

  # No warning for non-integer #successes

  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "glm", estimand = "ATE",
                  link = "br.logit", s.weights = "SW", epsilon = 1e-10)
  })
  expect_M_parts_okay(W, tolerance = eps)

  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "glm", estimand = "ATE",
                  link = "flac", s.weights = "SW")
  })

  expect_equal(test_data$SW, W$s.weights)

  #Stabilization

  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "glm", estimand = "ATE",
                  include.obj = TRUE, stabilize = TRUE)
  })

  expect_M_parts_okay(W, tolerance = eps)

  expect_null(attr(W, "Mparts", exact = TRUE))
  expect_false(is_null(attr(W, "Mparts.list", exact = TRUE)))

  expect_equal(cobalt::col_w_smd(W$covs, W$treat, W$weights),
               cobalt::col_w_smd(W0$covs, W0$treat, W0$weights),
               tolerance = eps)

  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "glm", estimand = "ATE",
                  include.obj = TRUE, stabilize = ~X1)
  })

  expect_M_parts_okay(W, tolerance = eps)

  expect_not_equal(cobalt::col_w_smd(W$covs, W$treat, W$weights),
                   cobalt::col_w_smd(W0$covs, W0$treat, W0$weights),
                   tolerance = eps)

  #Stab + s.weights
  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "glm", estimand = "ATE",
                  s.weights = "SW", stabilize = TRUE)
  })

  expect_M_parts_okay(W, tolerance = eps)

  expect_equal(cobalt::col_w_smd(W$covs, W$treat, W$weights, s.weights = W$s.weights),
               cobalt::col_w_smd(WS$covs, WS$treat, WS$weights, s.weights = WS$s.weights),
               tolerance = eps)

  expect_no_condition({
    WSs <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                    data = test_data, method = "glm", estimand = "ATE",
                    s.weights = "SW", stabilize = ~X1)
  })

  expect_M_parts_okay(WSs, tolerance = eps)

  expect_not_equal(cobalt::col_w_smd(WSs$covs, WSs$treat, WSs$weights, s.weights = WSs$s.weights),
                   cobalt::col_w_smd(W$covs, W$treat, W$weights, s.weights = W$s.weights),
                   tolerance = eps)

  #Non-full rank
  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 +
                    I(1 - X5) + I(X9 * 2),
                  data = test_data, method = "glm", estimand = "ATE",
                  include.obj = TRUE)
  })

  expect_equal(W$weights, W0$weights, tolerance = eps)

  # Separation
  set.seed(123)
  test_data$Xx <- rbinom(nrow(test_data), 1, .01 + .99 * test_data$A)

  expect_warning({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 +
                    Xx,
                  data = test_data, method = "glm", estimand = "ATE",
                  include.obj = TRUE)
  }, "Propensity scores numerically equal to 0 or 1 were estimated",
  ignore.case = TRUE)

  # expect_failure(expect_M_parts_okay(W))

  test_data$Xx <- NULL
})

test_that("Treatment guessing works for non-0/1 treatment", {
  skip_on_cran()
  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  # Non-0/1 treatment
  expect_message({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = transform(test_data, A = factor(A, levels = 0:1, labels = c("A", "B"))),
                  method = "glm", estimand = "ATT")
  }, '"B" is the treated')

  expect_ATT_weights_okay(W, focal = "B", tolerance = eps)

  expect_message({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = transform(test_data, A = factor(A, levels = 0:1, labels = c("B", "A"))),
                  method = "glm", estimand = "ATT")
  }, '"A" is the treated')

  expect_ATT_weights_okay(W, focal = "A", tolerance = eps)

  #When character, Z should always be guessed as treatment
  expect_message({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = transform(test_data, A = as.character(factor(A, levels = 0:1, labels = c("Z", "O")))),
                  method = "glm", estimand = "ATT")
  }, '"Z" is the treated')

  expect_ATT_weights_okay(W, focal = "Z", tolerance = eps)

  expect_message({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = transform(test_data, A = as.character(factor(A, levels = 0:1, labels = c("O", "Z")))),
                  method = "glm", estimand = "ATT")
  }, '"Z" is the treated')

  expect_ATT_weights_okay(W, focal = "Z", tolerance = eps)

  #ATC
  expect_message({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = transform(test_data, A = factor(A, levels = 0:1, labels = c("A", "B"))),
                  method = "glm", estimand = "ATC")
  }, '"A" is the control')

  expect_ATT_weights_okay(W, focal = "A", tolerance = eps)

  expect_message({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = transform(test_data, A = factor(A, levels = 0:1, labels = c("B", "A"))),
                  method = "glm", estimand = "ATC")
  }, '"B" is the control')

  expect_ATT_weights_okay(W, focal = "B", tolerance = eps)

  #When character, Z should always be guessed as treatment
  expect_message({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = transform(test_data, A = as.character(factor(A, levels = 0:1, labels = c("Z", "O")))),
                  method = "glm", estimand = "ATC")
  }, '"O" is the control')

  expect_ATT_weights_okay(W, focal = "O", tolerance = eps)

  expect_message({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = transform(test_data, A = as.character(factor(A, levels = 0:1, labels = c("O", "Z")))),
                  method = "glm", estimand = "ATC")
  }, '"O" is the control')

  expect_ATT_weights_okay(W, focal = "O", tolerance = eps)

  #Using "treat" and "control" synonyms should override other rules
  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = transform(test_data, A = as.character(factor(A, levels = 0:1, labels = c("unexposed", "exposed")))),
                  method = "glm", estimand = "ATT")
  })

  expect_ATT_weights_okay(W, focal = "exposed", tolerance = eps)

  expect_message({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = transform(test_data, A = as.character(factor(A, levels = 0:1, labels = c("unexposed", "control")))),
                  method = "glm", estimand = "ATT")
  }, '"unexposed" is the treated')

  expect_ATT_weights_okay(W, focal = "unexposed", tolerance = eps)

  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = transform(test_data, A = as.character(factor(A, levels = 0:1, labels = c("unexposed", "exposed")))),
                  method = "glm", estimand = "ATC")
  })

  expect_ATT_weights_okay(W, focal = "unexposed", tolerance = eps)

  expect_message({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = transform(test_data, A = as.character(factor(A, levels = 0:1, labels = c("unexposed", "control")))),
                  method = "glm", estimand = "ATC")
  }, '"control" is the control')

  expect_ATT_weights_okay(W, focal = "control", tolerance = eps)
})

test_that("Ordinal treatment", {
  skip_on_cran()
  skip_if_not_installed("rootSolve")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))
  test_data$Ao <- ordered(findInterval(test_data$Ac, quantile(test_data$Ac, seq(0, 1, length.out = 5)),
                                       all.inside = TRUE))

  expect_no_condition({
    W0 <- weightit(Ao ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                   data = test_data, method = "glm", estimand = "ATE",
                   include.obj = TRUE)
  })

  expect_M_parts_okay(W0, tolerance = eps)

  expect_no_condition({
    W <- weightit(Ao ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "glm", estimand = "ATE",
                  multi.method = "weightit",
                  include.obj = TRUE)
  })

  expect_M_parts_okay(W, tolerance = eps)

  #A "br." link requests bias reduction from the in-house cumulative link fitter,
  #which (unlike `brglm2::bracl()`, formerly used here) supports M-estimation
  for (link in c("br.logit", "br.probit", "br.cloglog", "br.loglog", "br.cauchit")) {
    expect_no_condition({
      W_br <- weightit(Ao ~ X1 + X2 + X3, data = test_data, method = "glm",
                       estimand = "ATE", link = link, include.obj = TRUE)
    })

    expect_true(W_br$obj$br)
    expect_M_parts_okay(W_br, tolerance = 1e-4)

    W_ml <- weightit(Ao ~ X1 + X2 + X3, data = test_data, method = "glm",
                     estimand = "ATE", link = sub("^br\\.", "", link),
                     include.obj = TRUE)

    expect_false(isTRUE(W_ml$obj$br))
    expect_true(max(abs(coef(W_br$obj) - coef(W_ml$obj))) > 1e-5)
  }

  #`multi.method = "polr"` cannot do bias reduction, so a "br." link routes to the
  #in-house fitter and gives exactly what `multi.method = "weightit"` gives
  expect_no_condition({
    W_polr_br <- weightit(Ao ~ X1 + X2 + X3, data = test_data, method = "glm",
                          estimand = "ATE", multi.method = "polr",
                          link = "br.logit", include.obj = TRUE)
  })

  W_wi_br <- weightit(Ao ~ X1 + X2 + X3, data = test_data, method = "glm",
                      estimand = "ATE", multi.method = "weightit",
                      link = "br.logit", include.obj = TRUE)

  expect_equal(W_polr_br$weights, W_wi_br$weights, tolerance = eps)

  #Without a "br." link, `multi.method = "polr"` still uses MASS::polr()
  skip_if_not_installed("MASS")

  W_polr <- weightit(Ao ~ X1 + X2 + X3, data = test_data, method = "glm",
                     estimand = "ATE", multi.method = "polr", include.obj = TRUE)

  expect_s3_class(W_polr$obj, "polr")

  #`multi.method = "bracl"` is deprecated and reroutes to the in-house fitter
  expect_warning({
    W_dep <- weightit(Ao ~ X1 + X2 + X3, data = test_data, method = "glm",
                      estimand = "ATE", multi.method = "bracl",
                      include.obj = TRUE)
  }, "deprecated")

  expect_true(W_dep$obj$br)
  expect_equal(W_dep$weights, W_wi_br$weights, tolerance = eps)
})

test_that("Multi-category treatment with a bias-reduced model", {
  skip_on_cran()
  skip_if_not_installed("rootSolve")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  expect_no_condition({
    W_br <- weightit(Am ~ X1 + X2 + X3, data = test_data, method = "glm",
                     estimand = "ATE", link = "br.logit", include.obj = TRUE)
  })

  expect_true(W_br$obj$br)
  expect_M_parts_okay(W_br, tolerance = 1e-4)

  W_ml <- weightit(Am ~ X1 + X2 + X3, data = test_data, method = "glm",
                   estimand = "ATE", include.obj = TRUE)

  expect_false(isTRUE(W_ml$obj$br))
  expect_true(max(abs(coef(W_br$obj) - coef(W_ml$obj))) > 1e-5)

  #Only the logit link is available for unordered multi-category treatments
  expect_error({
    weightit(Am ~ X1 + X2 + X3, data = test_data, method = "glm",
             estimand = "ATE", link = "br.probit")
  }, "logit")

  #`multi.method = "brmultinom"` is deprecated and reroutes to the in-house fitter
  expect_warning({
    W_dep <- weightit(Am ~ X1 + X2 + X3, data = test_data, method = "glm",
                      estimand = "ATE", multi.method = "brmultinom",
                      include.obj = TRUE)
  }, "deprecated")

  expect_true(W_dep$obj$br)
  expect_equal(W_dep$weights, W_br$weights, tolerance = eps)
})
